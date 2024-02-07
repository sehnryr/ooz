#pragma once

#include "stdafx.h"

typedef struct {
  // |p| holds the current byte and |p_end| the end of the buffer.
  const uint8_t *p, *p_end;
  // Bits accumulated so far
  uint32_t bits;
  // Next byte will end up in the |bitpos| position in |bits|.
  int bitpos;
} BitReader;

typedef struct {
  const uint8_t *p, *p_end;
  uint32_t bitpos;
} BitReader2;

void BitReader_Refill(BitReader *bits);
void BitReader_RefillBackwards(BitReader *bits);
int BitReader_ReadBit(BitReader *bits);
int BitReader_ReadBitNoRefill(BitReader *bits);
int BitReader_ReadBitsNoRefill(BitReader *bits, int n);
int BitReader_ReadBitsNoRefillZero(BitReader *bits, int n);
uint32_t BitReader_ReadMoreThan24Bits(BitReader *bits, int n);
uint32_t BitReader_ReadMoreThan24BitsB(BitReader *bits, int n);
int BitReader_ReadGamma(BitReader *bits);
int BitReader_ReadGammaX(BitReader *bits, int forced);
uint32_t BitReader_ReadDistance(BitReader *bits, uint32_t v);
uint32_t BitReader_ReadDistanceB(BitReader *bits, uint32_t v);
bool BitReader_ReadLength(BitReader *bits, uint32_t *v);
bool BitReader_ReadLengthB(BitReader *bits, uint32_t *v);
int BitReader_ReadFluff(BitReader *bits, int num_symbols);
