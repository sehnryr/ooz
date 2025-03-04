#pragma once

#include "bitreader.h"
#include "stdafx.h"

typedef struct {
  uint32_t A_used;
  uint32_t B_used;
  uint8_t A[256];
  uint32_t B[256];
} TansData;

typedef struct {
  uint32_t x;
  uint8_t bits_x;
  uint8_t symbol;
  uint16_t w;
} TansLutEnt;

typedef struct {
  TansLutEnt *lut;
  uint8_t *dst, *dst_end;
  const uint8_t *ptr_f, *ptr_b;
  uint32_t bits_f, bits_b;
  int bitpos_f, bitpos_b;
  uint32_t state_0, state_1, state_2, state_3, state_4;
} TansDecoderParams;

void SimpleSort(uint8_t *p, uint8_t *pend);
bool Tans_DecodeTable(BitReader *bits, int L_bits, TansData *tans_data);
void Tans_InitLut(TansData *tans_data, int L_bits, TansLutEnt *lut);
bool Tans_Decode(TansDecoderParams *params);
