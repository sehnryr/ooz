#include "kraken.h"

const uint8_t *Kraken_ParseQuantumHeader(OozQuantumHeader *hdr,
                                         const uint8_t *p, bool use_checksum) {
  uint32_t v = (p[0] << 16) | (p[1] << 8) | p[2];
  uint32_t size = v & 0x3FFFF;
  if (size != 0x3ffff) {
    hdr->compressed_size = size + 1;
    hdr->flag1 = (v >> 18) & 1;
    hdr->flag2 = (v >> 19) & 1;
    if (use_checksum) {
      hdr->checksum = (p[3] << 16) | (p[4] << 8) | p[5];
      return p + 6;
    } else {
      return p + 3;
    }
  }
  v >>= 18;
  if (v == 1) {
    // memset
    hdr->checksum = p[3];
    hdr->compressed_size = 0;
    hdr->whole_match_distance = 0;
    return p + 4;
  }
  return NULL;
}

bool Kraken_ReadLzTable(int mode, const uint8_t *src, const uint8_t *src_end,
                        uint8_t *dst, int dst_size, int offset,
                        uint8_t *scratch, uint8_t *scratch_end,
                        KrakenLzTable *lztable) {
  uint8_t *out;
  int decode_count, n;
  uint8_t *packed_offs_stream, *packed_len_stream;

  if (mode > 1)
    return false;

  if (src_end - src < 13)
    return false;

  if (offset == 0) {
    COPY_64(dst, src);
    dst += 8;
    src += 8;
  }

  if (*src & 0x80) {
    uint8_t flag = *src++;
    if ((flag & 0xc0) != 0x80)
      return false; // reserved flag set

    return false; // excess bytes not supported
  }

  // Disable no copy optimization if source and dest overlap
  bool force_copy = dst <= src_end && src <= dst + dst_size;

  // Decode lit stream, bounded by dst_size
  out = scratch;
  n = Ooz_Decodebytes(&out, src, src_end, &decode_count,
                      MIN(scratch_end - scratch, dst_size), force_copy, scratch,
                      scratch_end);
  if (n < 0)
    return false;
  src += n;
  lztable->lit_stream = out;
  lztable->lit_stream_size = decode_count;
  scratch += decode_count;

  // Decode command stream, bounded by dst_size
  out = scratch;
  n = Ooz_Decodebytes(&out, src, src_end, &decode_count,
                      MIN(scratch_end - scratch, dst_size), force_copy, scratch,
                      scratch_end);
  if (n < 0)
    return false;
  src += n;
  lztable->cmd_stream = out;
  lztable->cmd_stream_size = decode_count;
  scratch += decode_count;

  // Check if to decode the multistuff crap
  if (src_end - src < 3)
    return false;

  int offs_scaling = 0;
  uint8_t *packed_offs_stream_extra = NULL;

  if (src[0] & 0x80) {
    // uses the mode where distances are coded with 2 tables
    offs_scaling = src[0] - 127;
    src++;

    packed_offs_stream = scratch;
    n = Ooz_Decodebytes(&packed_offs_stream, src, src_end,
                        &lztable->offs_stream_size,
                        MIN(scratch_end - scratch, lztable->cmd_stream_size),
                        false, scratch, scratch_end);
    if (n < 0)
      return false;
    src += n;
    scratch += lztable->offs_stream_size;

    if (offs_scaling != 1) {
      packed_offs_stream_extra = scratch;
      n = Ooz_Decodebytes(&packed_offs_stream_extra, src, src_end,
                          &decode_count,
                          MIN(scratch_end - scratch, lztable->offs_stream_size),
                          false, scratch, scratch_end);
      if (n < 0 || decode_count != lztable->offs_stream_size)
        return false;
      src += n;
      scratch += decode_count;
    }
  } else {
    // Decode packed offset stream, it's bounded by the command length.
    packed_offs_stream = scratch;
    n = Ooz_Decodebytes(&packed_offs_stream, src, src_end,
                        &lztable->offs_stream_size,
                        MIN(scratch_end - scratch, lztable->cmd_stream_size),
                        false, scratch, scratch_end);
    if (n < 0)
      return false;
    src += n;
    scratch += lztable->offs_stream_size;
  }

  // Decode packed litlen stream. It's bounded by 1/4 of dst_size.
  packed_len_stream = scratch;
  n = Ooz_Decodebytes(
      &packed_len_stream, src, src_end, &lztable->len_stream_size,
      MIN(scratch_end - scratch, dst_size >> 2), false, scratch, scratch_end);
  if (n < 0)
    return false;
  src += n;
  scratch += lztable->len_stream_size;

  // Reserve memory for final dist stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->offs_stream = (int *)scratch;
  scratch += lztable->offs_stream_size * 4;

  // Reserve memory for final len stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->len_stream = (int *)scratch;
  scratch += lztable->len_stream_size * 4;

  if (scratch + 64 > scratch_end)
    return false;

  return Ooz_UnpackOffsets(src, src_end, packed_offs_stream,
                           packed_offs_stream_extra, lztable->offs_stream_size,
                           offs_scaling, packed_len_stream,
                           lztable->len_stream_size, lztable->offs_stream,
                           lztable->len_stream, 0, 0);
}

// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type0(KrakenLzTable *lzt, uint8_t *dst,
                                uint8_t *dst_end, uint8_t *dst_start) {
  const uint8_t *cmd_stream = lzt->cmd_stream,
                *cmd_stream_end = cmd_stream + lzt->cmd_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *len_stream_end = lzt->len_stream + lzt->len_stream_size;
  const uint8_t *lit_stream = lzt->lit_stream;
  const uint8_t *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *offs_stream = lzt->offs_stream;
  const int *offs_stream_end = lzt->offs_stream + lzt->offs_stream_size;
  const uint8_t *copyfrom;
  uint32_t final_len;
  int32_t offset;
  int32_t recent_offs[7];
  int32_t last_offset;

  recent_offs[3] = -8;
  recent_offs[4] = -8;
  recent_offs[5] = -8;
  last_offset = -8;

  while (cmd_stream < cmd_stream_end) {
    uint32_t f = *cmd_stream++;
    uint32_t litlen = f & 3;
    uint32_t offs_index = f >> 6;
    uint32_t matchlen = (f >> 2) & 0xF;

    // use cmov
    uint32_t next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = (litlen == 3) ? next_len_stream : len_stream;
    litlen = (litlen == 3) ? next_long_length : litlen;
    recent_offs[6] = *offs_stream;

    COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
    if (litlen > 8) {
      COPY_64_ADD(dst + 8, lit_stream + 8, &dst[last_offset + 8]);
      if (litlen > 16) {
        COPY_64_ADD(dst + 16, lit_stream + 16, &dst[last_offset + 16]);
        if (litlen > 24) {
          do {
            COPY_64_ADD(dst + 24, lit_stream + 24, &dst[last_offset + 24]);
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    offset = recent_offs[offs_index + 3];
    recent_offs[offs_index + 3] = recent_offs[offs_index + 2];
    recent_offs[offs_index + 2] = recent_offs[offs_index + 1];
    recent_offs[offs_index + 1] = recent_offs[offs_index + 0];
    recent_offs[3] = offset;
    last_offset = offset;

    offs_stream = (int *)((intptr_t)offs_stream + ((offs_index + 1) & 4));

    if ((uintptr_t)offset < (uintptr_t)(dst_start - dst))
      return false; // offset out of bounds

    copyfrom = dst + offset;
    if (matchlen != 15) {
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      dst += matchlen + 2;
    } else {
      matchlen = 14 + *len_stream++; // why is the value not 16 here, the above
                                     // case copies up to 16 bytes.
      if ((uintptr_t)matchlen > (uintptr_t)(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      COPY_64(dst + 16, copyfrom + 16);
      do {
        COPY_64(dst + 24, copyfrom + 24);
        matchlen -= 8;
        dst += 8;
        copyfrom += 8;
      } while (matchlen > 24);
      dst += matchlen;
    }
  }

  // check for incorrect input
  if (offs_stream != offs_stream_end || len_stream != len_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 8) {
    do {
      COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
      dst += 8, lit_stream += 8, final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst = *lit_stream++ + dst[last_offset];
    } while (dst++, --final_len);
  }
  return true;
}

// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type1(KrakenLzTable *lzt, uint8_t *dst,
                                uint8_t *dst_end, uint8_t *dst_start) {
  const uint8_t *cmd_stream = lzt->cmd_stream,
                *cmd_stream_end = cmd_stream + lzt->cmd_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *len_stream_end = lzt->len_stream + lzt->len_stream_size;
  const uint8_t *lit_stream = lzt->lit_stream;
  const uint8_t *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *offs_stream = lzt->offs_stream;
  const int *offs_stream_end = lzt->offs_stream + lzt->offs_stream_size;
  const uint8_t *copyfrom;
  uint32_t final_len;
  int32_t offset;
  int32_t recent_offs[7];

  recent_offs[3] = -8;
  recent_offs[4] = -8;
  recent_offs[5] = -8;

  while (cmd_stream < cmd_stream_end) {
    uint32_t f = *cmd_stream++;
    uint32_t litlen = f & 3;
    uint32_t offs_index = f >> 6;
    uint32_t matchlen = (f >> 2) & 0xF;

    // use cmov
    uint32_t next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = (litlen == 3) ? next_len_stream : len_stream;
    litlen = (litlen == 3) ? next_long_length : litlen;
    recent_offs[6] = *offs_stream;

    COPY_64(dst, lit_stream);
    if (litlen > 8) {
      COPY_64(dst + 8, lit_stream + 8);
      if (litlen > 16) {
        COPY_64(dst + 16, lit_stream + 16);
        if (litlen > 24) {
          do {
            COPY_64(dst + 24, lit_stream + 24);
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    offset = recent_offs[offs_index + 3];
    recent_offs[offs_index + 3] = recent_offs[offs_index + 2];
    recent_offs[offs_index + 2] = recent_offs[offs_index + 1];
    recent_offs[offs_index + 1] = recent_offs[offs_index + 0];
    recent_offs[3] = offset;

    offs_stream = (int *)((intptr_t)offs_stream + ((offs_index + 1) & 4));

    if ((uintptr_t)offset < (uintptr_t)(dst_start - dst))
      return false; // offset out of bounds

    copyfrom = dst + offset;
    if (matchlen != 15) {
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      dst += matchlen + 2;
    } else {
      matchlen = 14 + *len_stream++; // why is the value not 16 here, the above
                                     // case copies up to 16 bytes.
      if ((uintptr_t)matchlen > (uintptr_t)(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      COPY_64(dst + 16, copyfrom + 16);
      do {
        COPY_64(dst + 24, copyfrom + 24);
        matchlen -= 8;
        dst += 8;
        copyfrom += 8;
      } while (matchlen > 24);
      dst += matchlen;
    }
  }

  // check for incorrect input
  if (offs_stream != offs_stream_end || len_stream != len_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 64) {
    do {
      COPY_64_bytes(dst, lit_stream);
      dst += 64, lit_stream += 64, final_len -= 64;
    } while (final_len >= 64);
  }
  if (final_len >= 8) {
    do {
      COPY_64(dst, lit_stream);
      dst += 8, lit_stream += 8, final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst++ = *lit_stream++;
    } while (--final_len);
  }
  return true;
}

bool Kraken_ProcessLzRuns(int mode, uint8_t *dst, int dst_size, int offset,
                          KrakenLzTable *lztable) {
  uint8_t *dst_end = dst + dst_size;

  if (mode == 1)
    return Kraken_ProcessLzRuns_Type1(lztable, dst + (offset == 0 ? 8 : 0),
                                      dst_end, dst - offset);

  if (mode == 0)
    return Kraken_ProcessLzRuns_Type0(lztable, dst + (offset == 0 ? 8 : 0),
                                      dst_end, dst - offset);

  return false;
}

// Decode one 256kb big quantum block. It's divided into two 128k blocks
// internally that are compressed separately but with a shared history.
int Kraken_DecodeQuantum(uint8_t *dst, uint8_t *dst_end, uint8_t *dst_start,
                         const uint8_t *src, const uint8_t *src_end,
                         uint8_t *scratch, uint8_t *scratch_end) {
  const uint8_t *src_in = src;
  int mode, chunkhdr, dst_count, src_used, written_bytes;

  while (dst_end - dst != 0) {
    dst_count = dst_end - dst;
    if (dst_count > 0x20000)
      dst_count = 0x20000;
    if (src_end - src < 4)
      return -1;
    chunkhdr = src[2] | src[1] << 8 | src[0] << 16;
    if (!(chunkhdr & 0x800000)) {
      // Stored as entropy without any match copying.
      uint8_t *out = dst;
      src_used = Ooz_Decodebytes(&out, src, src_end, &written_bytes, dst_count,
                                 false, scratch, scratch_end);
      if (src_used < 0 || written_bytes != dst_count)
        return -1;
    } else {
      src += 3;
      src_used = chunkhdr & 0x7FFFF;
      mode = (chunkhdr >> 19) & 0xF;
      if (src_end - src < src_used)
        return -1;
      if (src_used < dst_count) {
        size_t scratch_usage = MIN(MIN(3 * dst_count + 32 + 0xd000, 0x6C000),
                                   scratch_end - scratch);
        if (scratch_usage < sizeof(KrakenLzTable))
          return -1;
        if (!Kraken_ReadLzTable(
                mode, src, src + src_used, dst, dst_count, dst - dst_start,
                scratch + sizeof(KrakenLzTable), scratch + scratch_usage,
                (KrakenLzTable *)scratch))
          return -1;
        if (!Kraken_ProcessLzRuns(mode, dst, dst_count, dst - dst_start,
                                  (KrakenLzTable *)scratch))
          return -1;
      } else if (src_used > dst_count || mode != 0) {
        return -1;
      } else {
        memmove(dst, src, dst_count);
      }
    }
    src += src_used;
    dst += dst_count;
  }
  return src - src_in;
}
