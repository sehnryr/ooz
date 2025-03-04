#pragma once

#include "bitreader.h"
#include "huffman.h"
#include "stdafx.h"

#define ALIGN_POINTER(p, align)                                                \
  ((uint8_t *)(((uintptr_t)(p) + (align - 1)) & ~(align - 1)))
#define ALIGN_16(x) (((x) + 15) & ~15)
#define COPY_64(d, s)                                                          \
  { *(uint64_t *)(d) = *(uint64_t *)(s); }
#define COPY_64_bytes(d, s)                                                    \
  {                                                                            \
    _mm_storeu_si128((__m128i *)d + 0, _mm_loadu_si128((__m128i *)s + 0));     \
    _mm_storeu_si128((__m128i *)d + 1, _mm_loadu_si128((__m128i *)s + 1));     \
    _mm_storeu_si128((__m128i *)d + 2, _mm_loadu_si128((__m128i *)s + 2));     \
    _mm_storeu_si128((__m128i *)d + 3, _mm_loadu_si128((__m128i *)s + 3));     \
  }
#define COPY_64_ADD(d, s, t)                                                   \
  _mm_storel_epi64((__m128i *)(d),                                             \
                   _mm_add_epi8(_mm_loadl_epi64((__m128i *)(s)),               \
                                _mm_loadl_epi64((__m128i *)(t))))

static const uint32_t kRiceCodeBits2Value[256] = {
    0x80000000, 0x00000007, 0x10000006, 0x00000006, 0x20000005, 0x00000105,
    0x10000005, 0x00000005, 0x30000004, 0x00000204, 0x10000104, 0x00000104,
    0x20000004, 0x00010004, 0x10000004, 0x00000004, 0x40000003, 0x00000303,
    0x10000203, 0x00000203, 0x20000103, 0x00010103, 0x10000103, 0x00000103,
    0x30000003, 0x00020003, 0x10010003, 0x00010003, 0x20000003, 0x01000003,
    0x10000003, 0x00000003, 0x50000002, 0x00000402, 0x10000302, 0x00000302,
    0x20000202, 0x00010202, 0x10000202, 0x00000202, 0x30000102, 0x00020102,
    0x10010102, 0x00010102, 0x20000102, 0x01000102, 0x10000102, 0x00000102,
    0x40000002, 0x00030002, 0x10020002, 0x00020002, 0x20010002, 0x01010002,
    0x10010002, 0x00010002, 0x30000002, 0x02000002, 0x11000002, 0x01000002,
    0x20000002, 0x00000012, 0x10000002, 0x00000002, 0x60000001, 0x00000501,
    0x10000401, 0x00000401, 0x20000301, 0x00010301, 0x10000301, 0x00000301,
    0x30000201, 0x00020201, 0x10010201, 0x00010201, 0x20000201, 0x01000201,
    0x10000201, 0x00000201, 0x40000101, 0x00030101, 0x10020101, 0x00020101,
    0x20010101, 0x01010101, 0x10010101, 0x00010101, 0x30000101, 0x02000101,
    0x11000101, 0x01000101, 0x20000101, 0x00000111, 0x10000101, 0x00000101,
    0x50000001, 0x00040001, 0x10030001, 0x00030001, 0x20020001, 0x01020001,
    0x10020001, 0x00020001, 0x30010001, 0x02010001, 0x11010001, 0x01010001,
    0x20010001, 0x00010011, 0x10010001, 0x00010001, 0x40000001, 0x03000001,
    0x12000001, 0x02000001, 0x21000001, 0x01000011, 0x11000001, 0x01000001,
    0x30000001, 0x00000021, 0x10000011, 0x00000011, 0x20000001, 0x00001001,
    0x10000001, 0x00000001, 0x70000000, 0x00000600, 0x10000500, 0x00000500,
    0x20000400, 0x00010400, 0x10000400, 0x00000400, 0x30000300, 0x00020300,
    0x10010300, 0x00010300, 0x20000300, 0x01000300, 0x10000300, 0x00000300,
    0x40000200, 0x00030200, 0x10020200, 0x00020200, 0x20010200, 0x01010200,
    0x10010200, 0x00010200, 0x30000200, 0x02000200, 0x11000200, 0x01000200,
    0x20000200, 0x00000210, 0x10000200, 0x00000200, 0x50000100, 0x00040100,
    0x10030100, 0x00030100, 0x20020100, 0x01020100, 0x10020100, 0x00020100,
    0x30010100, 0x02010100, 0x11010100, 0x01010100, 0x20010100, 0x00010110,
    0x10010100, 0x00010100, 0x40000100, 0x03000100, 0x12000100, 0x02000100,
    0x21000100, 0x01000110, 0x11000100, 0x01000100, 0x30000100, 0x00000120,
    0x10000110, 0x00000110, 0x20000100, 0x00001100, 0x10000100, 0x00000100,
    0x60000000, 0x00050000, 0x10040000, 0x00040000, 0x20030000, 0x01030000,
    0x10030000, 0x00030000, 0x30020000, 0x02020000, 0x11020000, 0x01020000,
    0x20020000, 0x00020010, 0x10020000, 0x00020000, 0x40010000, 0x03010000,
    0x12010000, 0x02010000, 0x21010000, 0x01010010, 0x11010000, 0x01010000,
    0x30010000, 0x00010020, 0x10010010, 0x00010010, 0x20010000, 0x00011000,
    0x10010000, 0x00010000, 0x50000000, 0x04000000, 0x13000000, 0x03000000,
    0x22000000, 0x02000010, 0x12000000, 0x02000000, 0x31000000, 0x01000020,
    0x11000010, 0x01000010, 0x21000000, 0x01001000, 0x11000000, 0x01000000,
    0x40000000, 0x00000030, 0x10000020, 0x00000020, 0x20000010, 0x00001010,
    0x10000010, 0x00000010, 0x30000000, 0x00002000, 0x10001000, 0x00001000,
    0x20000000, 0x00100000, 0x10000000, 0x00000000,
};

static const uint8_t kRiceCodeBits2Len[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
};

static uint32_t bitmasks[32] = {
    0x1,        0x3,       0x7,       0xf,       0x1f,       0x3f,
    0x7f,       0xff,      0x1ff,     0x3ff,     0x7ff,      0xfff,
    0x1fff,     0x3fff,    0x7fff,    0xffff,    0x1ffff,    0x3ffff,
    0x7ffff,    0xfffff,   0x1fffff,  0x3fffff,  0x7fffff,   0xffffff,
    0x1ffffff,  0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff,
    0x7fffffff, 0xffffffff};

// Header in front of each 256k block
typedef struct OozHeader {
  // Type of decoder used, 6 means kraken
  int decoder_type;

  // Whether to restart the decoder
  bool restart_decoder;

  // Whether this block is uncompressed
  bool uncompressed;

  // Whether this block uses checksums.
  bool use_checksums;
} OozHeader;

// Additional header in front of each 256k block ("quantum").
typedef struct OozQuantumHeader {
  // The compressed size of this quantum. If this value is 0 it means
  // the quantum is a special quantum such as memset.
  uint32_t compressed_size;
  // If checksums are enabled, holds the checksum.
  uint32_t checksum;
  // Two flags
  uint8_t flag1;
  uint8_t flag2;
  // Whether the whole block matched a previous block
  uint32_t whole_match_distance;
} OozQuantumHeader;

typedef struct OozDecoder {
  // Updated after the |*_DecodeStep| function completes to hold
  // the number of bytes read and written.
  int src_used, dst_used;

  // Pointer to a 256k buffer that holds the intermediate state
  // in between decode phase 1 and 2.
  uint8_t *scratch;
  size_t scratch_size;

  OozHeader hdr;
} OozDecoder;

void *MallocAligned(size_t size, size_t alignment);
void FreeAligned(void *p);

uint32_t BSR(uint32_t x);
uint32_t BSF(uint32_t x);

int CountLeadingZeros(uint32_t bits);
int Log2RoundUp(uint32_t v);

bool DecodeGolombRiceLengths(uint8_t *dst, size_t size, BitReader2 *br);
bool DecodeGolombRiceBits(uint8_t *dst, uint32_t size, uint32_t bitcount,
                          BitReader2 *br);

void FillByteOverflow16(uint8_t *dst, uint8_t v, size_t n);

void CombineScaledOffsetArrays(int *offs_stream, size_t offs_stream_size,
                               int scale, const uint8_t *low_bits);

OozDecoder *Ooz_Create();
void Ooz_Destroy(OozDecoder *decoder);
const uint8_t *Ooz_ParseHeader(OozHeader *hdr, const uint8_t *p);
int Ooz_Decodebytes(uint8_t **output, const uint8_t *src,
                    const uint8_t *src_end, int *decoded_size,
                    size_t output_size, bool force_memmove, uint8_t *scratch,
                    uint8_t *scratch_end);
int Ooz_Decodebytes_Type12(const uint8_t *src, size_t src_size, uint8_t *output,
                           int output_size, int type);
bool Ooz_DecodebytesCore(HuffReader *hr, HuffRevLut *lut);
int Ooz_DecodeRecursive(const uint8_t *src, size_t src_size, uint8_t *output,
                        int output_size, uint8_t *scratch,
                        uint8_t *scratch_end);
int Ooz_DecodeMultiArray(const uint8_t *src, const uint8_t *src_end,
                         uint8_t *dst, uint8_t *dst_end, uint8_t **array_data,
                         int *array_lens, int array_count, int *total_size_out,
                         bool force_memmove, uint8_t *scratch,
                         uint8_t *scratch_end);
int Ooz_GetBlockSize(const uint8_t *src, const uint8_t *src_end, int *dest_size,
                     int dest_capacity);
int Ooz_DecodeRLE(const uint8_t *src, size_t src_size, uint8_t *dst,
                  int dst_size, uint8_t *scratch, uint8_t *scratch_end);
int Ooz_DecodeTans(const uint8_t *src, size_t src_size, uint8_t *dst,
                   int dst_size, uint8_t *scratch, uint8_t *scratch_end);
bool Ooz_UnpackOffsets(const uint8_t *src, const uint8_t *src_end,
                       const uint8_t *packed_offs_stream,
                       const uint8_t *packed_offs_stream_extra,
                       int packed_offs_stream_size, int multi_dist_scale,
                       const uint8_t *packed_litlen_stream,
                       int packed_litlen_stream_size, int *offs_stream,
                       int *len_stream, bool excess_flag, int excess_bytes);
uint32_t Ooz_GetCrc(const uint8_t *p, size_t p_size);
void Ooz_CopyWholeMatch(uint8_t *dst, uint32_t offset, size_t length);
bool Ooz_DecodeStep(struct OozDecoder *dec, uint8_t *dst_start, int offset,
                    size_t dst_bytes_left_in, const uint8_t *src,
                    size_t src_bytes_left);
int Ooz_Decompress(const uint8_t *src, size_t src_len, uint8_t *dst,
                   size_t dst_len);
