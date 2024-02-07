#pragma once

#include "decompress.h"
#include "stdafx.h"

// Kraken decompression happens in two phases, first one decodes
// all the literals and copy lengths using huffman and second
// phase runs the copy loop. This holds the tables needed by stage 2.
typedef struct KrakenLzTable {
  // Stream of (literal, match) pairs. The flag byte contains
  // the length of the match, the length of the literal and whether
  // to use a recent offset.
  uint8_t *cmd_stream;
  int cmd_stream_size;

  // Holds the actual distances in case we're not using a recent
  // offset.
  int *offs_stream;
  int offs_stream_size;

  // Holds the sequence of literals. All literal copying happens from
  // here.
  uint8_t *lit_stream;
  int lit_stream_size;

  // Holds the lengths that do not fit in the flag stream. Both literal
  // lengths and match length are stored in the same array.
  int *len_stream;
  int len_stream_size;
} KrakenLzTable;

const uint8_t *Kraken_ParseQuantumHeader(OozQuantumHeader *hdr,
                                         const uint8_t *p, bool use_checksum);
bool Kraken_ReadLzTable(int mode, const uint8_t *src, const uint8_t *src_end,
                        uint8_t *dst, int dst_size, int offset,
                        uint8_t *scratch, uint8_t *scratch_end,
                        KrakenLzTable *lztable);
bool Kraken_ProcessLzRuns_Type0(KrakenLzTable *lzt, uint8_t *dst,
                                uint8_t *dst_end, uint8_t *dst_start);
bool Kraken_ProcessLzRuns_Type1(KrakenLzTable *lzt, uint8_t *dst,
                                uint8_t *dst_end, uint8_t *dst_start);
bool Kraken_ProcessLzRuns(int mode, uint8_t *dst, int dst_size, int offset,
                          KrakenLzTable *lztable);
int Kraken_DecodeQuantum(uint8_t *dst, uint8_t *dst_end, uint8_t *dst_start,
                         const uint8_t *src, const uint8_t *src_end,
                         uint8_t *scratch, uint8_t *scratch_end);
