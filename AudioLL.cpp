#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <chrono>
#include <cstdint>
#include <iostream>

#include "../Untitled1.h"
#include "../defs.h"
#include "../wave.h"
#include "kissfft/kiss_fft.h"
#include "structs.h"
using namespace std;

extern uint32_t g_bitMask32[];
extern huffdata g_huffDecBandQs[];
extern huffdata g_huffDecBandQsDiff[];
extern huffdata g_huffDecBandQsLayerQ2[];
extern huffdata g_huffDecBandQsLayerQ3[];
extern int32_t g_mdctHraWindow480With10MsQ30[];
extern int32_t g_mdctHraWindow960With10MsQ30[];

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef INT32_MAX_F
#define INT32_MAX_F 2147483647.0
#endif
uint32_t readbits_times = 0;

typedef struct {
  uint32_t* stream_buffer;        // 0
  int32_t index_in_stream;        // 20
  int32_t index_in_DWORD;         // 24
  int32_t total_bits;  // 28
  int32_t read_bits;              // 32
} read_index2;
FILE* bit_record = NULL;
uint64_t ReadBitsInDWORD(read_index2* i32, int nbits) {
  uint64_t result = 0;
  if (!nbits) {
    return result;
  }
  if ((i32->total_bits) < (i32->read_bits + nbits)) {
    // puts("bits is not enough");
    return 0;
  }

  if (32 - i32->index_in_DWORD < nbits) {
    int Remaining_bits = nbits - (32 - i32->index_in_DWORD);
    result = i32->stream_buffer[i32->index_in_stream] & g_bitMask32[32 - i32->index_in_DWORD];
    result <<= Remaining_bits;
    i32->index_in_stream += 1;

    result |= i32->stream_buffer[i32->index_in_stream] >> (64 - i32->index_in_DWORD - nbits);
    i32->index_in_DWORD = Remaining_bits;
  } else {
    result = (i32->stream_buffer[i32->index_in_stream] & g_bitMask32[32 - i32->index_in_DWORD]) >>
             (32 - i32->index_in_DWORD - nbits);
    i32->index_in_DWORD += nbits;
  }
  i32->read_bits += nbits;
  readbits_times += 1;
  // fwrite(&i32->read_bits,1,4,bit_record);
  return result;
}
void AudioAlignReadHeadToByte(read_index2* i32) {
  auto tmp = i32->read_bits - (i32->index_in_DWORD & 7) + 8;
  i32->read_bits = tmp;
  i32->index_in_stream = tmp >> 5;
  i32->index_in_DWORD = tmp & 0x1F;
}
void subi32(read_index2* i32, int nbits) {
  i32->read_bits -= nbits;
  i32->index_in_DWORD -= nbits;
  if (i32->index_in_DWORD < 0) {
    i32->index_in_DWORD += 32;
    i32->index_in_stream -= 1;
  }
}

uint8_t stream_buffer[4000];

void AudioBitstreamRearrangeByte(uint32_t* stream_buffer, size_t s) {
  int32_t DWORD_NUM;  // x8
  DWORD_NUM = s / 4;
  if (s % 4) {
    DWORD_NUM += 1;
  }
  for (size_t i = 0; i < DWORD_NUM; i++) {
    stream_buffer[i] = __builtin_bswap32(stream_buffer[i]);
  }
}

int32_t smpRate_table[4] = {44100, 48000, 88200, 96000};
int32_t frLength_table[4] = {120, 240, 480, 960};

void AudioGetBandId(int32* a1, int Nband, int frLength, double a4) {
  double v8;   // d0
  double tmp;  // d9
  int tmp2;    // w8

  *a1 = 0;
  if (Nband >= 1) {
    tmp = (double)Nband / pow(frLength, 1.0 / a4);
    for (int Nbandi = 1; Nbandi < Nband; Nbandi++) {
      tmp2 = (int)pow((double)Nbandi / tmp, a4);
      if (Nbandi > tmp2) {
        tmp2 = Nbandi;
      }
      a1[Nbandi] = tmp2;
    }
  }
  a1[Nband] = frLength;
}

int32_t* g_mdctHraWindow = NULL;
int32_t* RotationTable_sin = NULL;
int32_t* RotationTable_tan = NULL;
int32_t* FoldingParamA = NULL;
int32_t* FoldingParamB = NULL;
int32_t* FoldingParamC = NULL;
int32_t* _TwParamA = NULL;
int32_t* _TwParamB = NULL;
int32_t* _TwParamC = NULL;
int32_t* _TwParamD = NULL;

int32_t* ImdctPrevious[2];

int32_t* is_set_signed[1];  // 符号位是否已赋值
int32_t* had_signed[1];     // 符号位
int32_t* unpack_data[1];    // 数据

int32_t* ImdctOutput[2];
int32_t* DataFromStream[2];

int32_t* fft_block1 = NULL;
int32_t* fft_block2 = NULL;
int32_t* mdct_block1 = NULL;
int32_t* mdct_block2 = NULL;
// int32_t (*quantScale)[0];

FILE* wave_file_out = NULL;
int32_t wave_pcm_buf[2 * 960];

int32_t BandId[64];

int32_t old_frLength = 0;
int32_t old_band_split_scale = 0;
kiss_fft_cfg fft_cfg = NULL;
int32_t print_times = 0;

void LLinit(int32_t frLength, int32_t intrinsic_delay) {
  int32_t half_intrinsic_delay = intrinsic_delay >> 1;
  int32_t half_frLength = frLength >> 1;
  int32_t half_active_length = (frLength - intrinsic_delay) >> 1;
  size_t i;

  if (frLength == 480 && intrinsic_delay == 120) {
    g_mdctHraWindow = (int32_t*)malloc(960 * sizeof(int32_t));
    for (i = 0; i < 960; i++) {
      g_mdctHraWindow[i] = 2 * g_mdctHraWindow480With10MsQ30[i];
    }

  } else if (frLength == 960 || intrinsic_delay == 240) {
    g_mdctHraWindow = (int32_t*)malloc(1920 * sizeof(int32_t));
    for (i = 0; i < 1920; i++) {
      g_mdctHraWindow[i] = 2 * g_mdctHraWindow960With10MsQ30[i];
    }
  } else {
    return;
  }
  if (g_mdctHraWindow == NULL) {
    return;
  }
  fft_cfg = kiss_fft_alloc(frLength / 4, false, 0, 0);
  ImdctPrevious[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  ImdctPrevious[1] = (int32_t*)malloc(frLength * sizeof(int32_t));
  memset(ImdctPrevious[0], 0, frLength * sizeof(int32_t));
  memset(ImdctPrevious[1], 0, frLength * sizeof(int32_t));

  RotationTable_sin = (int32_t*)malloc(half_frLength * sizeof(int32_t));
  RotationTable_tan = (int32_t*)malloc(half_frLength * sizeof(int32_t));
  for (i = 0; i < half_frLength; i++) {
    auto tmp_sin = sin((double)(i * 2 + 1) * M_PI / (double)(4 * frLength));
    RotationTable_sin[i] = (tmp_sin * INT32_MAX_F);
    auto tmp_cos = cos((double)(i * 2 + 1) * M_PI / (double)(4 * frLength));
    RotationTable_tan[i] = ((tmp_cos + -1.0) / tmp_sin * INT32_MAX_F);
  }
  FoldingParamA = (int32_t*)malloc(half_frLength * sizeof(int32_t));
  FoldingParamB = (int32_t*)malloc(half_frLength * sizeof(int32_t));
  FoldingParamC = (int32_t*)malloc(half_frLength * sizeof(int32_t));
  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamB[i] = -g_mdctHraWindow[half_frLength - half_intrinsic_delay + i];
  }

  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamA[i] = g_mdctHraWindow[frLength + half_frLength - half_intrinsic_delay + i] - 0x7FFFFFFF;
  }

  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamA[i] = (int32_t)((double)FoldingParamA[i] / (double)FoldingParamB[i] * INT32_MAX_F);
  }
  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamC[i] = g_mdctHraWindow[half_frLength + half_intrinsic_delay - i - 1];
  }
  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamC[i] -= 0x7FFFFFFF;
  }
  for (i = 0; i < half_intrinsic_delay; i++) {
    FoldingParamC[i] = (int32_t)((double)FoldingParamC[i] / (double)FoldingParamB[i] * INT32_MAX_F);
  }

  for (i = 0; i < half_active_length; i++) {
    FoldingParamB[half_intrinsic_delay + i] = -g_mdctHraWindow[i];
  }

  for (i = 0; i < half_active_length; i++) {
    FoldingParamA[half_intrinsic_delay + i] = g_mdctHraWindow[frLength + i] - 0x7FFFFFFF;
  }

  for (i = 0; i < half_active_length; i++) {
    FoldingParamA[half_intrinsic_delay + i] = (int32_t)((double)FoldingParamA[half_intrinsic_delay + i] /
                                                        (double)FoldingParamB[half_intrinsic_delay + i] * INT32_MAX_F);
    if (FoldingParamA[half_intrinsic_delay + i] == 0x80000000 && half_intrinsic_delay + i >= half_frLength / 2) {
      FoldingParamA[half_intrinsic_delay + i] = 0x7FFFFFFF;
    }
  }

  for (i = 0; i < half_active_length; i++) {
    FoldingParamC[half_intrinsic_delay + i] = g_mdctHraWindow[frLength - i - 1];
  }
  for (i = 0; i < half_active_length; i++) {
    FoldingParamC[half_intrinsic_delay + i] -= 0x7FFFFFFF;
  }
  for (i = 0; i < half_active_length; i++) {
    FoldingParamC[half_intrinsic_delay + i] = (int32_t)((double)FoldingParamC[half_intrinsic_delay + i] /
                                                        (double)FoldingParamB[half_intrinsic_delay + i] * INT32_MAX_F);
    if (FoldingParamC[half_intrinsic_delay + i] == 0x80000000 && half_intrinsic_delay + i < half_frLength / 2) {
      FoldingParamC[half_intrinsic_delay + i] = 0x7FFFFFFF;
    }
  }
  _TwParamA = (int32_t*)malloc(frLength / 4 * sizeof(int32_t));
  _TwParamB = (int32_t*)malloc(frLength / 4 * sizeof(int32_t));
  _TwParamC = (int32_t*)malloc(frLength / 2 * sizeof(int32_t));
  _TwParamD = (int32_t*)malloc(frLength / 2 * sizeof(int32_t));
  for (i = 0; i < frLength / 4; i++) {
    _TwParamA[i] = cos((double)(i * 8 + 1) * -M_PI / (double)(frLength * 4)) * INT32_MAX_F;
  }

  for (i = 0; i < frLength / 4; i++) {
    _TwParamB[i] = sin((double)(i * 8 + 1) * -M_PI / (double)(frLength * 4)) * INT32_MAX_F;
  }
  for (i = 0; i < frLength / 2; i++) {
    _TwParamC[i] = cos((double)(i * 8 + 1) * -M_PI / (double)(frLength * 4)) * INT32_MAX_F;
  }

  for (i = 0; i < frLength / 2; i++) {
    _TwParamD[i] = sin((double)(i * 8 + 1) * -M_PI / (double)(frLength * 4)) * INT32_MAX_F;
  }
  /**/
  int test_mem = 4;
  is_set_signed[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  had_signed[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  unpack_data[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  ImdctOutput[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  ImdctOutput[1] = (int32_t*)malloc(frLength * sizeof(int32_t));
  DataFromStream[0] = (int32_t*)malloc(frLength * sizeof(int32_t));
  DataFromStream[1] = (int32_t*)malloc(frLength * sizeof(int32_t));

  fft_block1 = (int32_t*)malloc(frLength * sizeof(int32_t));
  fft_block2 = (int32_t*)malloc(frLength * sizeof(int32_t));
  mdct_block1 = (int32_t*)malloc(2 * frLength * sizeof(int32_t));
  mdct_block2 = (int32_t*)malloc(2 * frLength * sizeof(int32_t));
}
void LLdeinit() {
  if (g_mdctHraWindow == NULL) {
    return;
  }
  kiss_fft_free(fft_cfg);
  free(g_mdctHraWindow);
  free(RotationTable_sin);
  free(RotationTable_tan);
  free(FoldingParamA);
  free(FoldingParamB);
  free(FoldingParamC);
  free(_TwParamA);
  free(_TwParamB);
  free(_TwParamC);
  free(_TwParamD);

  free(ImdctPrevious[0]);
  free(ImdctPrevious[1]);

  free(is_set_signed[0]);
  free(had_signed[0]);
  free(unpack_data[0]);

  free(ImdctOutput[0]);
  free(ImdctOutput[1]);
  free(DataFromStream[0]);
  free(DataFromStream[1]);
  free(fft_block1);
  free(fft_block2);
  free(mdct_block1);
  free(mdct_block2);
  fft_cfg = NULL;
  g_mdctHraWindow = NULL;

  RotationTable_sin = NULL;
  RotationTable_tan = NULL;

  FoldingParamA = NULL;
  FoldingParamB = NULL;
  FoldingParamC = NULL;

  _TwParamA = NULL;
  _TwParamB = NULL;
  _TwParamC = NULL;
  _TwParamD = NULL;

  ImdctPrevious[0] = NULL;
  ImdctPrevious[1] = NULL;
  
  is_set_signed[0] = NULL;
  had_signed[0] = NULL;
  unpack_data[0] = NULL;

  ImdctOutput[0] = NULL;
  ImdctOutput[1] = NULL;
  DataFromStream[0] = NULL;
  DataFromStream[1] = NULL;
  
  fft_block1 = NULL;
  fft_block2 = NULL;
  mdct_block1 = NULL;
  mdct_block2 = NULL;
}

void AudioWaveOutputFromInt24(int32_t** in, int32_t frLength, int32_t channels, int bitPerSample, void* pcm_out);
void LLunpack(int one_pack_size, int pack_index) {
  uint32_t* stream_buffer_s32 = (uint32_t*)stream_buffer;
  AudioBitstreamRearrangeByte(stream_buffer_s32, one_pack_size);

  read_index2 i32 = {
      .stream_buffer = stream_buffer_s32,
      .index_in_stream = 0,
      .index_in_DWORD = 0,
      .total_bits = one_pack_size * 8,
      .read_bits = 0};

  int32_t codectype = 0;                   // 0 0
  int32_t smpRate2 = 96000;                // 4 1
  int32_t bitPerSample = 32;               // 8 2
  int32_t frLength = 960;                  // 12 3
  int32_t channels = 2;                    // 16 4
  int32_t stream_bitrate = 960;            // 20 5  15b
  int32_t UsedCBR = 1;                     // 24 6  2b
  int32_t intrinsic_delay = frLength / 4;  // 28 7  10b
  int32_t a8 = 2;                          // 32 8  2b flag
  int32_t a9 = 1;                          // 36 9  1b flag
  int32_t band_split_scale = 5;            // 40 10 4b
  int32_t a11 = 0;                         // 44 11 1b flag
  int32_t a12 = 0;                         // 48 12 1b
  int32_t a13 = 1;                         // 52 13 1b
  int32_t a14 = 1;                         // 56 14
  int32_t ext_cfg = 1;                     // 60 15 1b

  int32_t isInterleaveStream;  // 56
  int32_t UsedInt16;           // 60
  int32_t Nband = 0;
  codectype = ReadBitsInDWORD(&i32, 2);
  smpRate2 = smpRate_table[ReadBitsInDWORD(&i32, 2)];
  Nband = sqrt((double)smpRate2 / 50.0);
  channels = ReadBitsInDWORD(&i32, 1) + 1;
  frLength = frLength_table[ReadBitsInDWORD(&i32, 2)];
  stream_bitrate = ReadBitsInDWORD(&i32, 15);
  ext_cfg = ReadBitsInDWORD(&i32, 1);
  UsedInt16 = ReadBitsInDWORD(&i32, 1);
  isInterleaveStream = ReadBitsInDWORD(&i32, 1);
  if (codectype == 3) {
    UsedCBR = 2;
  } else if (codectype == 2) {
    UsedCBR = 1;
  } else if (codectype == 1) {
    puts("is l2hcv2");
  } else {
    UsedCBR = 0;
  }

  if (ext_cfg) {
    UsedCBR = ReadBitsInDWORD(&i32, 2);
    intrinsic_delay = ReadBitsInDWORD(&i32, 10);  // AudioGetL2hcWindow
    a8 = ReadBitsInDWORD(&i32, 2);                // sub_206DC
    band_split_scale = ReadBitsInDWORD(&i32, 4);
    a9 = ReadBitsInDWORD(&i32, 1);   // AudioLLDecode
    a11 = ReadBitsInDWORD(&i32, 1);  // AudioLLDecode
    a12 = ReadBitsInDWORD(&i32, 1);
    a13 = ReadBitsInDWORD(&i32, 1);
  }

  if (old_frLength != frLength) {
    LLdeinit();
    LLinit(frLength, intrinsic_delay);

    old_frLength = frLength;
    wave_write_header(wave_file_out, bitPerSample, (bitPerSample >> 3), smpRate2, channels, frLength);
  }
  int32_t half_active_length = (frLength - intrinsic_delay) / 2;
  if (old_band_split_scale != band_split_scale) {
    AudioGetBandId(BandId, Nband, frLength, (double)band_split_scale);
    band_split_scale = old_band_split_scale;
  }

  int32_t sel_FoldingParamA = 0;
  int32_t sel_FoldingParamB = 1;
  int32_t sel_FoldingParamC = 0;
  uint32_t out_channels = 2;

  int32_t(*quantScale)[Nband];
  *(int32_t**)&quantScale = new int32_t[channels * Nband];
  memset(quantScale, 0, sizeof(int32_t) * channels * Nband);
  AudioAlignReadHeadToByte(&i32);
  memset(ImdctOutput[0], 0, sizeof(int32_t) * frLength);
  memset(ImdctOutput[1], 0, sizeof(int32_t) * frLength);
  memset(DataFromStream[0], 0, sizeof(int32_t) * frLength);
  memset(DataFromStream[1], 0, sizeof(int32_t) * frLength);
  if (a8 != 2) {  // sub_206DC
    puts("a8 must equal 2");
  }
  int32_t BitTarget = 0;
  if (UsedCBR == 2) {
    BitTarget = one_pack_size * 8;
  } else if (UsedCBR == 1) {
    BitTarget = (double)(stream_bitrate * frLength) / ((double)smpRate2 / 1000.0);
    BitTarget &= 0xFFFFFFF8;
  } else {
    BitTarget = 0x10000;
  }
  if (a11 == 1) {
    puts("a11 cannot equal 1");
  }
  int32_t Remaining_Bits_Per_Channel = (BitTarget - i32.read_bits) / channels;

  int32_t read_times = 0;
  int32_t const1 = 1;  // about channel
  size_t const1_i;
  size_t Nband_i;

  if (UsedCBR == 1 && a12 == 1 && out_channels < channels && !isInterleaveStream) {
  } else if (channels >= 1) {
    for (size_t channels_i = 0; channels_i < channels; channels_i++) {
      memset(is_set_signed[0], 0, sizeof(int32_t) * frLength * const1);
      memset(had_signed[0], 0, sizeof(int32_t) * frLength * const1);
      memset(unpack_data[0], 0, sizeof(int32_t) * frLength * const1);
      if (UsedCBR != 1)  // if not CBR
      {
        // VBR
        // int32_t diffFlag[const1];
        int32_t diffFlag[1];

        for (const1_i = 0; const1_i < const1; const1_i++) {
          diffFlag[const1_i] = ReadBitsInDWORD(&i32, 1);
        }

        int32_t const0 = 0;
        for (const1_i = 0; const1_i < const1; const1_i++) {
          for (Nband_i = const0; Nband_i < Nband; Nband_i++) {
            if (diffFlag[const1_i]) {
              if (Nband_i == 0) {
                auto tmp_p = g_huffDecBandQs[ReadBitsInDWORD(&i32, 8)];
                quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p.data;
                subi32(&i32, 8 - tmp_p.nbits);
              } else {
                auto tmp_p2 = g_huffDecBandQsDiff[ReadBitsInDWORD(&i32, 7)];
                quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p2.data;
                quantScale[channels_i + const1_i][Nband_i + const0] +=
                    quantScale[channels_i + const1_i][Nband_i + const0 - 1] - 10;
                subi32(&i32, 7 - tmp_p2.nbits);
              }
            } else {
              auto tmp_p = g_huffDecBandQs[ReadBitsInDWORD(&i32, 8)];
              quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p.data;
              subi32(&i32, 8 - tmp_p.nbits);
            }
            printf("%d ,", quantScale[channels_i][Nband_i]);
          }
        }
        printf(" \n");
        if (const1 >= 1) {
          for (const1_i = 0; const1_i < const1; const1_i++) {
            for (Nband_i = 0; Nband_i < Nband; Nband_i++) {
              if (quantScale[const1_i + channels_i][Nband_i] > 32) {
                printf("error \n");
              }
            }
          }
        }
        for (Nband_i = 0; Nband_i < Nband; Nband_i++) {
          for (const1_i = 0; const1_i < const1; const1_i++) {
            if (quantScale[channels_i + const1_i][Nband_i] <= 0) {
              continue;
            }
            if (quantScale[channels_i + const1_i][Nband_i] == 2) {
              int32_t tmpBandId = BandId[Nband_i];
              while (tmpBandId < BandId[Nband_i + 1]) {
                auto tmp_p3 = g_huffDecBandQsLayerQ2[ReadBitsInDWORD(&i32, 3)];
                subi32(&i32, 3 - tmp_p3.nbits);
                unpack_data[const1_i][tmpBandId] += tmp_p3.data;

                if (!is_set_signed[const1_i][tmpBandId] && tmp_p3.data >= 1) {
                  is_set_signed[const1_i][tmpBandId] = 1;
                  had_signed[const1_i][tmpBandId] = ReadBitsInDWORD(&i32, 1);
                }
                tmpBandId += 1;
              }
            } else if (quantScale[channels_i + const1_i][Nband_i] == 1) {
              int32_t tmpBandsize = BandId[Nband_i + 1] - BandId[Nband_i];
              int32_t band_size_div3 = tmpBandsize / 3;
              int32_t band_size_round3 = 3 * band_size_div3;
              int32_t band_size_remain3 = tmpBandsize % 3;
              if (tmpBandsize >= 3) {
                int32_t BandIdoffset = 2;
                for (size_t i = 0; i < band_size_div3; BandIdoffset += 3, i++) {
                  auto tmp_p4 = g_huffDecBandQsLayerQ3[ReadBitsInDWORD(&i32, 5)];
                  int32_t huff_data_2bit = (tmp_p4.data >> 2) & 1;
                  int32_t huff_data_1bit = (tmp_p4.data >> 1) & 1;
                  int32_t huff_data_0bit = (tmp_p4.data >> 0) & 1;

                  int32_t huff_data_0bit_index = BandId[Nband_i] + BandIdoffset - 0;
                  int32_t huff_data_1bit_index = BandId[Nband_i] + BandIdoffset - 1;
                  int32_t huff_data_2bit_index = BandId[Nband_i] + BandIdoffset - 2;
                  subi32(&i32, 5 - tmp_p4.nbits);

                  unpack_data[const1_i][huff_data_2bit_index] += huff_data_2bit;
                  unpack_data[const1_i][huff_data_1bit_index] += huff_data_1bit;
                  unpack_data[const1_i][huff_data_0bit_index] += huff_data_0bit;
                  if (!is_set_signed[const1_i][huff_data_2bit_index] && huff_data_2bit != 0) {
                    is_set_signed[const1_i][huff_data_2bit_index] = 1;
                    had_signed[const1_i][huff_data_2bit_index] = ReadBitsInDWORD(&i32, 1);
                  }
                  if (!is_set_signed[const1_i][huff_data_1bit_index] && huff_data_1bit != 0) {
                    is_set_signed[const1_i][huff_data_1bit_index] = 1;
                    had_signed[const1_i][huff_data_1bit_index] = ReadBitsInDWORD(&i32, 1);
                  }
                  if (!is_set_signed[const1_i][huff_data_0bit_index] && huff_data_0bit != 0) {
                    is_set_signed[const1_i][huff_data_0bit_index] = 1;
                    had_signed[const1_i][huff_data_0bit_index] = ReadBitsInDWORD(&i32, 1);
                  }
                }
              }
              if (band_size_remain3 >= 1) {
                for (size_t band_size_remain3_index = 0; band_size_remain3_index < band_size_remain3;
                     band_size_remain3_index++) {
                  auto tmp_p5 = ReadBitsInDWORD(&i32, 1);
                  auto tmp_index = band_size_round3 + band_size_remain3_index + BandId[Nband_i];
                  unpack_data[const1_i][tmp_index] += tmp_p5;
                  if (!is_set_signed[const1_i][tmp_index] && tmp_p5 >= 1) {
                    had_signed[const1_i][tmp_index] = ReadBitsInDWORD(&i32, 1);
                    is_set_signed[const1_i][tmp_index] = 1;
                  }
                }
              }
            } else if (quantScale[channels_i + const1_i][Nband_i] >= 3) {
              int32_t need_read = quantScale[channels_i + const1_i][Nband_i] - 3;
              int32_t tmpBandId = BandId[Nband_i];

              while (tmpBandId < BandId[Nband_i + 1]) {
                auto tmp_p6 = g_huffDecBandQsLayerQ3[ReadBitsInDWORD(&i32, 5)];
                subi32(&i32, 5 - tmp_p6.nbits);
                auto tmp_data = tmp_p6.data << need_read;
                read_times += 1;
                if (quantScale[channels_i + const1_i][Nband_i] >= 4) {
                  tmp_data += ReadBitsInDWORD(&i32, need_read);
                }
                unpack_data[const1_i][tmpBandId] += tmp_data;
                if (!is_set_signed[const1_i][tmpBandId] && tmp_data >= 1) {
                  had_signed[const1_i][tmpBandId] = ReadBitsInDWORD(&i32, 1);
                  is_set_signed[const1_i][tmpBandId] = 1;
                }
                tmpBandId += 1;
              }
            }
          }
        }
        for (const1_i = 0; const1_i < const1; const1_i++) {
          for (size_t frLength_index = 0; frLength_index < frLength; frLength_index++) {
            if (had_signed[const1_i][frLength_index] > 0) {
              unpack_data[const1_i][frLength_index] = -unpack_data[const1_i][frLength_index];
            }
          }
          memcpy(DataFromStream[const1_i + channels_i], unpack_data[const1_i], frLength * sizeof(int32_t));
        }
      } else {
        // CBR
        // int32_t diffFlag[const1];
        int32_t already_Readbits = i32.read_bits;
        int32_t diffFlag[1];
        for (const1_i = 0; const1_i < const1; const1_i++) {
          diffFlag[const1_i] = ReadBitsInDWORD(&i32, 1);
        }
        auto tmp_f1 = fmin(fmax(sqrt((double)Remaining_Bits_Per_Channel / (double)(frLength * const1)), 0.2), 1.0);
        int32_t newNband = 0;
        int32_t pre_val = BandId[0], next_val = 0;
        for (newNband = 0; newNband < Nband; newNband++) {
          next_val = BandId[newNband + 1];
          if (pre_val < next_val) {
            int iVar6 = (int)((double)frLength * tmp_f1);
            if (iVar6 <= pre_val) {
              iVar6 = pre_val;
            }
            iVar6 = iVar6 - pre_val;
            pre_val = next_val - pre_val;
            if (iVar6 < pre_val || iVar6 == 0) {
              break;
            }
          }
          pre_val = next_val;
        }

        int32_t const0 = 0;
        for (const1_i = 0; const1_i < const1; const1_i++) {
          for (Nband_i = const0; Nband_i < newNband; Nband_i++) {
            if (diffFlag[const1_i]) {
              if (Nband_i == 0) {
                auto tmp_p = g_huffDecBandQs[ReadBitsInDWORD(&i32, 8)];
                quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p.data;
                subi32(&i32, 8 - tmp_p.nbits);
              } else {
                auto tmp_p2 = g_huffDecBandQsDiff[ReadBitsInDWORD(&i32, 7)];
                quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p2.data;
                quantScale[channels_i + const1_i][Nband_i + const0] +=
                    quantScale[channels_i + const1_i][Nband_i + const0 - 1] - 10;
                subi32(&i32, 7 - tmp_p2.nbits);
              }
            } else {
              auto tmp_p = g_huffDecBandQs[ReadBitsInDWORD(&i32, 8)];
              quantScale[channels_i + const1_i][Nband_i + const0] = tmp_p.data;
              subi32(&i32, 8 - tmp_p.nbits);
            }
          }
        }
        for (const1_i = 0; const1_i < const1; const1_i++) {
          for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
            if (quantScale[const1_i + channels_i][Nband_i] > 32) {
              printf("error \n");
            }
          }
        }

        int32_t(*psyScalefactor)[Nband];
        *(int32_t**)&psyScalefactor = new int32_t[const1 * Nband];
        int32_t(*psyScalefactor_bak)[Nband];
        *(int32_t**)&psyScalefactor_bak = new int32_t[const1 * Nband];

        memcpy(psyScalefactor_bak, quantScale[channels_i], sizeof(int32_t) * Nband * const1);
        memcpy(psyScalefactor, quantScale[channels_i], sizeof(int32_t) * Nband * const1);
        memset(quantScale[channels_i], 0, sizeof(int32_t) * Nband * const1);
        auto AudioBandPsyAcoustic = [](int* p_psyScalefactor, int nband, int x, int* out) {
          int max = p_psyScalefactor[0];
          for (size_t i = 0; i < nband; i++) {
            if (max < p_psyScalefactor[i]) {
              max = p_psyScalefactor[i];
            }
          }
          if (max < 1) {
            memset(out, 0, sizeof(int) * nband);
            return 0;
          }
          if (max == 1) {
            memcpy(out, p_psyScalefactor, sizeof(int) * nband);
          } else {
            int min = max - x;
            if (min <= 0) {
              min = 0;
            }
            for (size_t i = 0; i < nband; i++) {
              out[i] = (p_psyScalefactor[i] > min);
            }
          }
          return 1;
        };
        int32_t(*inner_mem_pool3)[Nband];
        *(int32_t**)&inner_mem_pool3 = new int32_t[const1 * Nband];
        int continue_flag = 0;
        int loop_index = 0;

        do {
          int x1 = 3 * ((++loop_index) % 3);
          if (x1 < 1) {
            x1 = 1;
          }
          continue_flag = 0;
          for (const1_i = 0; const1_i < const1; const1_i++) {
            continue_flag |= AudioBandPsyAcoustic(psyScalefactor[const1_i], Nband, x1, inner_mem_pool3[const1_i]);

            for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
              psyScalefactor[const1_i][Nband_i] -= inner_mem_pool3[const1_i][Nband_i];
              quantScale[const1_i + channels_i][Nband_i] =
                  psyScalefactor_bak[const1_i][Nband_i] - psyScalefactor[const1_i][Nband_i];
            }
          }
          // AudioEstimateBitCount
          int EstimateBitCount = 0;
          for (const1_i = 0; const1_i < const1; const1_i++) {
            for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
              EstimateBitCount +=
                  (BandId[Nband_i + 1] - BandId[Nband_i]) * (quantScale[const1_i + channels_i][Nband_i] + 2);
            }
          }

          if ((EstimateBitCount >= (Remaining_Bits_Per_Channel + already_Readbits - i32.read_bits)) || !continue_flag) {
            for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
              print_times += 1;
              for (const1_i = 0; const1_i < const1; const1_i++) {
                auto NoiseFloorScale = psyScalefactor[const1_i][Nband_i];
                if (quantScale[channels_i + const1_i][Nband_i] <= 0) {
                  continue;
                }
                if (quantScale[channels_i + const1_i][Nband_i] == 2) {
                  int32_t tmpBandId = BandId[Nband_i];
                  while (tmpBandId < BandId[Nband_i + 1]) {
                    if (i32.read_bits >= (Remaining_Bits_Per_Channel + already_Readbits - 12)) {
                      // printf("out of read\n");
                      continue_flag = 0;
                      break;
                    }
                    auto tmp_p3 = g_huffDecBandQsLayerQ2[ReadBitsInDWORD(&i32, 3)];
                    subi32(&i32, 3 - tmp_p3.nbits);
                    unpack_data[const1_i][tmpBandId] += (tmp_p3.data << NoiseFloorScale);

                    if (!is_set_signed[const1_i][tmpBandId] && tmp_p3.data >= 1) {
                      is_set_signed[const1_i][tmpBandId] = 1;
                      had_signed[const1_i][tmpBandId] = ReadBitsInDWORD(&i32, 1);
                    }
                    tmpBandId += 1;
                  }
                } else if (quantScale[channels_i + const1_i][Nband_i] == 1) {
                  int32_t tmpBandsize = BandId[Nband_i + 1] - BandId[Nband_i];
                  int32_t band_size_div3 = tmpBandsize / 3;
                  int32_t band_size_round3 = 3 * band_size_div3;
                  int32_t band_size_remain3 = tmpBandsize % 3;
                  if (tmpBandsize >= 3) {
                    int32_t BandIdoffset = 2;
                    for (size_t i = 0; i < band_size_div3; BandIdoffset += 3, i++) {
                      if (i32.read_bits >= (Remaining_Bits_Per_Channel + already_Readbits - 12)) {
                        // printf("out of read\n");
                        continue_flag = 0;
                        break;
                      }
                      auto tmp_p4 = g_huffDecBandQsLayerQ3[ReadBitsInDWORD(&i32, 5)];
                      int32_t huff_data_2bit = ((tmp_p4.data >> 2) & 1) << NoiseFloorScale;
                      int32_t huff_data_1bit = ((tmp_p4.data >> 1) & 1) << NoiseFloorScale;
                      int32_t huff_data_0bit = ((tmp_p4.data >> 0) & 1) << NoiseFloorScale;

                      int32_t huff_data_0bit_index = BandId[Nband_i] + BandIdoffset - 0;
                      int32_t huff_data_1bit_index = BandId[Nband_i] + BandIdoffset - 1;
                      int32_t huff_data_2bit_index = BandId[Nband_i] + BandIdoffset - 2;
                      subi32(&i32, 5 - tmp_p4.nbits);

                      unpack_data[const1_i][huff_data_2bit_index] += huff_data_2bit;
                      unpack_data[const1_i][huff_data_1bit_index] += huff_data_1bit;
                      unpack_data[const1_i][huff_data_0bit_index] += huff_data_0bit;
                      if (!is_set_signed[const1_i][huff_data_2bit_index] && huff_data_2bit != 0) {
                        is_set_signed[const1_i][huff_data_2bit_index] = 1;
                        had_signed[const1_i][huff_data_2bit_index] = ReadBitsInDWORD(&i32, 1);
                      }
                      if (!is_set_signed[const1_i][huff_data_1bit_index] && huff_data_1bit != 0) {
                        is_set_signed[const1_i][huff_data_1bit_index] = 1;
                        had_signed[const1_i][huff_data_1bit_index] = ReadBitsInDWORD(&i32, 1);
                      }
                      if (!is_set_signed[const1_i][huff_data_0bit_index] && huff_data_0bit != 0) {
                        is_set_signed[const1_i][huff_data_0bit_index] = 1;
                        had_signed[const1_i][huff_data_0bit_index] = ReadBitsInDWORD(&i32, 1);
                      }
                    }
                  }
                  if (band_size_remain3 >= 1) {
                    for (size_t band_size_remain3_index = 0; band_size_remain3_index < band_size_remain3;
                         band_size_remain3_index++) {
                      if (i32.read_bits >= (Remaining_Bits_Per_Channel + already_Readbits - 12)) {
                        // printf("out of read\n");
                        continue_flag = 0;
                        break;
                      }
                      auto tmp_p5 = ReadBitsInDWORD(&i32, 1) << NoiseFloorScale;
                      auto tmp_index = band_size_round3 + band_size_remain3_index + BandId[Nband_i];
                      unpack_data[const1_i][tmp_index] += tmp_p5;
                      if (!is_set_signed[const1_i][tmp_index] && tmp_p5 >= 1) {
                        had_signed[const1_i][tmp_index] = ReadBitsInDWORD(&i32, 1);
                        is_set_signed[const1_i][tmp_index] = 1;
                      }
                    }
                  }
                } else if (quantScale[channels_i + const1_i][Nband_i] >= 3) {
                  int32_t need_read = quantScale[channels_i + const1_i][Nband_i] - 3;
                  int32_t tmpBandId = BandId[Nband_i];

                  while (tmpBandId < BandId[Nband_i + 1]) {
                    if (i32.read_bits >= (Remaining_Bits_Per_Channel + already_Readbits - 12)) {
                      // printf("out of read\n");
                      continue_flag = 0;
                      break;
                    }
                    auto tmp_p6 = g_huffDecBandQsLayerQ3[ReadBitsInDWORD(&i32, 5)];
                    subi32(&i32, 5 - tmp_p6.nbits);
                    auto tmp_data = tmp_p6.data << need_read;
                    read_times += 1;
                    if (quantScale[channels_i + const1_i][Nband_i] >= 4) {
                      tmp_data += ReadBitsInDWORD(&i32, need_read);
                    }
                    unpack_data[const1_i][tmpBandId] += tmp_data << NoiseFloorScale;
                    if (!is_set_signed[const1_i][tmpBandId] && tmp_data >= 1) {
                      had_signed[const1_i][tmpBandId] = ReadBitsInDWORD(&i32, 1);
                      is_set_signed[const1_i][tmpBandId] = 1;
                    }
                    tmpBandId += 1;
                  }
                }
              }
            }
            memcpy(psyScalefactor_bak, psyScalefactor, sizeof(int32_t) * Nband * const1);
          } else {
            if ((EstimateBitCount >= (Remaining_Bits_Per_Channel + already_Readbits - i32.read_bits))) {
              continue_flag = 0;
            }
          }
        } while (continue_flag);
        /*
        for (size_t j = 0; j < 480; j++)
        {
          printf("%d ,",unpack_data[0][j]);
        }
        printf("\n");
        */
        auto AudioGetBandQsTotal = [](int* psyScalefactor, int Nband, int* pBandId) {
          uint16_t res = 0;
          for (size_t Nband_i = 0; Nband_i < Nband; Nband_i++) {
            res += (pBandId[Nband_i + 1] - pBandId[Nband_i]) * psyScalefactor[Nband_i];
            // printf(" %d,",psyScalefactor[Nband_i]);
          }
          // printf("\n");
          return res;
        };
        for (const1_i = 0; const1_i < const1; const1_i++) {
          int32_t BandQsTotal = AudioGetBandQsTotal(psyScalefactor[const1_i], newNband, BandId);
          int32_t BandQsTotal_bak = BandQsTotal;
          for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
            auto pSf = psyScalefactor[const1_i][Nband_i];
            if (pSf <= 1) pSf = 0;
            if (pSf != 0 && Nband_i < (newNband >> 1)) {
              pSf -= 1;
            }
            for (size_t j = BandId[Nband_i]; j < BandId[Nband_i + 1]; j++) {
              auto t114 = unpack_data[const1_i][j];
              if (t114 >= 1) {
                BandQsTotal = (int16_t)(0x7C4D * BandQsTotal + 0x3619);
                auto tmp_mask = (~(0xffffffff << pSf));
                auto gain = (0x7FFF * BandQsTotal) & tmp_mask;
                unpack_data[const1_i][j] += gain;
              }
            }
          }
          if (a13 == 1) {
            for (Nband_i = 0; Nband_i < newNband; Nband_i++) {
              auto pSf = psyScalefactor[const1_i][Nband_i] - 3;
              if (pSf <= 0) {
                continue;
              }
              for (size_t j = BandId[Nband_i]; j < BandId[Nband_i + 1]; j++) {
                BandQsTotal_bak = (int16_t)(0x7C4D * BandQsTotal_bak + 0x3619);
                if (!unpack_data[const1_i][j]) {
                  unpack_data[const1_i][j] = ~(-1 << pSf);
                  had_signed[const1_i][j] = BandQsTotal_bak < 0x8000;
                }
              }
            }
          }
          for (size_t frLength_index = 0; frLength_index < frLength; frLength_index++) {
            if (had_signed[const1_i][frLength_index] > 0) {
              unpack_data[const1_i][frLength_index] = -unpack_data[const1_i][frLength_index];
            }
          }
          memcpy(DataFromStream[const1_i + channels_i], unpack_data[const1_i], frLength * sizeof(int32_t));
        }
        operator delete[](psyScalefactor);
        operator delete[](psyScalefactor_bak);
        operator delete[](inner_mem_pool3);
      }
    }
  }
  operator delete[](quantScale);
  /*
  for (size_t i = 0; i < 2; i++)
  {
    for (size_t j = 0; j < 480; j++)
    {
      printf("%d ,", DataFromStream[i][j]);
    }
    printf("\n\n");
  }
  */
  if (abs(i32.total_bits - i32.read_bits) > 8) {
    printf("error %d\n", pack_index);
  }
  AudioAlignReadHeadToByte(&i32);

  auto gen_fixpt = [](double in) { return (uint64_t)(in * (double)0x7FFFFFFF) & 0x7FFFFFFFFFFFFFFF; };

  auto AudioIntInvMs = [](int32_t* ch1, int32_t* ch2, int32_t infrLength) {
    for (size_t i = 0; i < infrLength; i++) {
      ch1[i] -= (int64_t)(0x7FFFFFFFCAFB0CCDll * ((int64_t)ch2[i])) >> 31;  // *= 1/-(1+sqrt(2))
      ch2[i] -= (int64_t)(0x000000005A827999ll * ((int64_t)ch1[i])) >> 31;  // *= sqrt(2)/2
      ch1[i] -= (int64_t)(0x7FFFFFFFCAFB0CCDll * ((int64_t)ch2[i])) >> 31;  // *= 1/-(1+sqrt(2))
      // ch1=(sqrt(2)/2)*(ch2+ch1)
      // ch2=(sqrt(2)/2)*(ch2-ch1)
    }
  };

  if (isInterleaveStream == 1) {
    AudioIntInvMs(DataFromStream[0], DataFromStream[1], frLength);
  }

  if (a9 != 1) {
    puts("a9 must equal 1");
  }
  auto AudioDct4 =
      [](int32_t* in, int32_t len, int32_t* out, int32_t* inTwParamA, int32_t* inTwParamB, kiss_fft_cfg cfg) {
        int32_t _max_ = abs(in[0]);
        memset(fft_block1, 0, 2 * len * sizeof(int32_t));
        memset(fft_block2, 0, 2 * len * sizeof(int32_t));
        size_t i;
        for (i = 0; i < len; i++) {
          if (abs(in[i]) > _max_) {
            _max_ = abs(in[i]);
          }
        }

        int32_t _max_bits = 0;
        while ((_max_ >> _max_bits) != 0) {
          _max_bits += 1;
        }
        int32_t len_bits = 0;
        while ((len >> len_bits) != 0) {
          len_bits += 1;
        }

        int32_t over_flowing_bits = len_bits + _max_bits;

        over_flowing_bits -= 31;
        if (over_flowing_bits > 2) {
          over_flowing_bits = 2;
        } else if (over_flowing_bits < 0) {
          over_flowing_bits = 0;
        }

        for (i = 0; i < len; i++) {
          fft_block1[i] = (in[i] >> over_flowing_bits);  // * (len >> 1);//* (len >> 1)是因为kissfft的算法问题
        }
        for (size_t i = 1, j = len - 1; i < len / 2; i += 2, j -= 2) {
          auto tmp = fft_block1[i];
          fft_block1[i] = fft_block1[j];
          fft_block1[j] = tmp;
        }
        for (i = 0; i < len; i += 2) {
          auto tmp1 = fft_block1[i];
          auto tmp2 = fft_block1[i + 1];
          fft_block1[i] =
              (((uint64_t)((int64_t)tmp1 * (int64_t)inTwParamA[i >> 1]) >> 31) -
               ((uint64_t)((int64_t)tmp2 * (int64_t)inTwParamB[i >> 1]) >> 31));
          fft_block1[i + 1] =
              (((uint64_t)((int64_t)tmp1 * (int64_t)inTwParamB[i >> 1]) >> 31) +
               ((uint64_t)((int64_t)tmp2 * (int64_t)inTwParamA[i >> 1]) >> 31));
        }

        kiss_fft(cfg, (kiss_fft_cpx*)fft_block1, (kiss_fft_cpx*)fft_block2);

        for (i = 0; i < len; i += 2) {
          auto tmp1 = fft_block2[i];
          auto tmp2 = fft_block2[i + 1];
          out[i] = (((uint64_t)((int64_t)tmp1 * (int64_t)inTwParamA[i >> 1]) >> 31) -
                    ((uint64_t)((int64_t)tmp2 * (int64_t)inTwParamB[i >> 1]) >> 31)) *
                   2;
          out[i + 1] = (((uint64_t)((int64_t)tmp1 * (int64_t)inTwParamB[i >> 1]) >> 31) +
                        ((uint64_t)((int64_t)tmp2 * (int64_t)inTwParamA[i >> 1]) >> 31)) *
                       2;
        }
        size_t j;
        for (i = 1, j = len - 1; i < len / 2; i += 2, j -= 2) {
          auto tmp = -out[i];
          out[i] = -out[j];
          out[j] = tmp;
        }
        int64_t normal = sqrt(0.5 / (double)len) * INT32_MAX_F;
        for (i = 0; i < len; i++) {
          out[i] = ((normal * out[i]) >> 31) << over_flowing_bits;
        }
        return;
      };
  // AudioMonoIntMdctSyn start

  memcpy(ImdctOutput[0], DataFromStream[0], sizeof(int32_t) * frLength);
  memcpy(ImdctOutput[1], DataFromStream[1], sizeof(int32_t) * frLength);
  // decode per channel
  // sub_21814 start
  for (size_t channels_i = 0; channels_i < channels; channels_i++) {
    memset(mdct_block1, 0, sizeof(int32_t) * 2 * frLength);
    auto half_frLength = frLength >> 1;
    auto last_index = frLength - 1;
    auto ImdctOutput_chn = ImdctOutput[channels_i];
    size_t i;
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[half_frLength + i] = -ImdctOutput_chn[half_frLength + i];
    }
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[i] -= (uint64_t)((int64_t)RotationTable_tan[i] * (int64_t)ImdctOutput_chn[last_index - i]) >> 31;
    }
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[last_index - i] -= (uint64_t)((int64_t)RotationTable_sin[i] * (int64_t)ImdctOutput_chn[i]) >> 31;
    }
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[i] -= (uint64_t)((int64_t)RotationTable_tan[i] * (int64_t)ImdctOutput_chn[last_index - i]) >> 31;
    }

    AudioDct4(&ImdctOutput_chn[half_frLength], half_frLength, mdct_block1, _TwParamA, _TwParamB, fft_cfg);

    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[i] -= (uint64_t)(0x5A827999LL * (int64_t)mdct_block1[i]) >> 31;
    }
    AudioDct4(&ImdctOutput_chn[0], half_frLength, mdct_block1, _TwParamA, _TwParamB, fft_cfg);
    for (i = 0; i < half_frLength; i++) { /*
       ImdctOutput_chn[half_frLength + i] -= ((uint64_t)(0x7FFFFFFFA57D8668LL * (int64_t)mdct_block1[i]) >> 31)
                                             << 1;*/
      ImdctOutput_chn[half_frLength + i] -= ((uint64_t)(0x7FFFFFFFA57D8668LL * (int64_t)mdct_block1[i]) >> 30);
    }
    AudioDct4(&ImdctOutput_chn[half_frLength], half_frLength, mdct_block1, _TwParamA, _TwParamB, fft_cfg);
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[i] -= (uint64_t)(0x5A827999LL * (int64_t)mdct_block1[i]) >> 31;
      ImdctOutput_chn[i] -= (uint64_t)(0x7FFFFFFFC0000001LL * (int64_t)ImdctOutput_chn[half_frLength + i]) >> 31;
      ImdctOutput_chn[half_frLength + i] -= ImdctOutput_chn[i];
    }
    memcpy(mdct_block1, ImdctOutput_chn, sizeof(int32_t) * frLength);
    for (i = 0; i < half_frLength; i++) {
      ImdctOutput_chn[i * 2] = mdct_block1[i];
      ImdctOutput_chn[i * 2 + 1] = mdct_block1[half_frLength + i];
    }

    // sub_21D20 start
    for (i = 0; i < frLength; i += 4) {
      auto backup1 = ImdctOutput_chn[i + 2];
      ImdctOutput_chn[i + 2] = ImdctOutput_chn[i + 3];
      ImdctOutput_chn[i + 3] = backup1;
    }
    // sub_21D20 end
    // sub_21814 end

    for (i = 0; i < half_frLength; i++) {
      auto backup1 = -ImdctOutput_chn[i];
      ImdctOutput_chn[i] = -ImdctOutput_chn[last_index - i];
      ImdctOutput_chn[last_index - i] = backup1;
    }

    //  AudioIntInvWinTdac start
    memset(mdct_block1, 0, sizeof(int32_t) * 2 * frLength);

    // memset(mdct_block2, 0, sizeof(int32_t) * frLength);
    int32_t half_intrinsic_delay = intrinsic_delay >> 1;
    memcpy(mdct_block2, ImdctOutput_chn, sizeof(int32_t) * frLength);
    memcpy(mdct_block1, &ImdctPrevious[channels_i][half_intrinsic_delay], sizeof(int32_t) * half_active_length);
    memcpy(&mdct_block1[half_active_length], &mdct_block2[half_intrinsic_delay], sizeof(int32_t) * half_active_length);
    memcpy(
        &ImdctPrevious[channels_i][half_intrinsic_delay],
        &mdct_block2[half_intrinsic_delay + half_active_length],
        sizeof(int32_t) * half_active_length);

    // sub_214B0
    {
      for (i = 0; i < half_active_length / 2; i++) {
        auto backup1 = mdct_block1[half_active_length + i];
        mdct_block1[half_active_length + i] = mdct_block1[2 * half_active_length - i - 1];
        mdct_block1[2 * half_active_length - i - 1] = backup1;
      }
      for (i = 0; i < half_active_length * sel_FoldingParamA; i++) {
        mdct_block1[i] -= (uint64_t)((int64_t)mdct_block1[half_active_length + i] *
                                     (int64_t)FoldingParamA[half_intrinsic_delay + i]) >>
                          31;
      }

      for (i = 0; i < half_active_length * sel_FoldingParamB; i++) {
        mdct_block1[half_active_length + i] -=
            (uint64_t)((int64_t)mdct_block1[i] * (int64_t)FoldingParamB[half_intrinsic_delay + i]) >> 31;
      }

      for (i = 0; i < half_active_length * sel_FoldingParamC; i++) {
        mdct_block1[i] -= (uint64_t)((int64_t)mdct_block1[half_active_length + i] *
                                     (int64_t)FoldingParamC[half_intrinsic_delay + i]) >>
                          31;
      }

      for (i = 0; i < half_active_length / 2; i++) {
        auto backup1 = mdct_block1[half_active_length + i];
        mdct_block1[half_active_length + i] = mdct_block1[2 * half_active_length - i - 1];
        mdct_block1[2 * half_active_length - i - 1] = backup1;
      }
    }
    for (i = 0; i < half_active_length; i++) {  //?
      ImdctOutput_chn[intrinsic_delay + i] = mdct_block1[half_active_length + i];
    }

    memcpy(
        &ImdctOutput_chn[intrinsic_delay + half_active_length],
        &mdct_block2[half_intrinsic_delay + half_active_length],
        sizeof(int32_t) * half_active_length);

    for (i = 0; i < half_intrinsic_delay; i++) {
      ImdctOutput_chn[i] = ImdctPrevious[channels_i][i];
    }
    memcpy(
        ImdctPrevious[channels_i],
        &mdct_block2[frLength - half_intrinsic_delay],
        half_intrinsic_delay * sizeof(int32_t));
    memcpy(&ImdctOutput_chn[half_intrinsic_delay], mdct_block2, half_intrinsic_delay * sizeof(int32_t));
    // sub_214B0
    {
      for (i = 0; i < half_intrinsic_delay / 2; i++) {
        auto backup1 = ImdctOutput_chn[half_intrinsic_delay + i];
        ImdctOutput_chn[half_intrinsic_delay + i] = ImdctOutput_chn[2 * half_intrinsic_delay - 1 - i];
        ImdctOutput_chn[2 * half_intrinsic_delay - 1 - i] = backup1;
      }
      for (i = 0; i < half_intrinsic_delay; i++) {
        ImdctOutput_chn[i] -=
            (uint64_t)((int64_t)ImdctOutput_chn[half_intrinsic_delay + i] * (int64_t)FoldingParamA[i]) >> 31;
      }

      for (i = 0; i < half_intrinsic_delay; i++) {
        ImdctOutput_chn[half_intrinsic_delay + i] -=
            (uint64_t)((int64_t)ImdctOutput_chn[i] * (int64_t)FoldingParamB[i]) >> 31;
      }

      for (i = 0; i < half_intrinsic_delay; i++) {
        ImdctOutput_chn[i] -=
            (uint64_t)((int64_t)ImdctOutput_chn[half_intrinsic_delay + i] * (int64_t)FoldingParamC[i]) >> 31;
      }

      for (i = 0; i < half_intrinsic_delay / 2; i++) {
        auto backup1 = ImdctOutput_chn[half_intrinsic_delay + i];
        ImdctOutput_chn[half_intrinsic_delay + i] = ImdctOutput_chn[2 * half_intrinsic_delay - 1 - i];
        ImdctOutput_chn[2 * half_intrinsic_delay - 1 - i] = backup1;
      }
    }
  }
  memset(wave_pcm_buf, 0, sizeof(wave_pcm_buf));

  AudioWaveOutputFromInt24((int32_t**)ImdctOutput, frLength, channels, bitPerSample, wave_pcm_buf);

  fwrite(wave_pcm_buf, 1, frLength * channels * bitPerSample >> 3, wave_file_out);

  return;
}
void AudioWaveOutputFromInt24(int32_t** in, int32_t frLength, int32_t channels, int bitPerSample, void* pcm_out) {
  switch (bitPerSample) {
    case -32: {
      for (size_t i = 0; i < channels; i++) {
        for (size_t j = 0; j < frLength; j++) {
          auto sample = in[i][j];
          if (sample < -0x800000) {
            sample = -0x800000;
          } else if (sample > 0x800000) {
            sample = 0x800000;
          }
          *((float*)pcm_out + i + j * channels) = (double)sample / ((double)0x800000);
        }
      }
    } break;
    case 16: {
      for (size_t i = 0; i < channels; i++) {
        for (size_t j = 0; j < frLength; j++) {
          auto sample = in[i][j];
          if (sample < -0x800000) {
            sample = -0x800000;
          } else if (sample > 0x800000) {
            sample = 0x800000;
          }
          sample >>= 8;
          *((int16_t*)pcm_out + i + j * channels) = sample;
        }
      }
    } break;
    case 24: {
      for (size_t i = 0; i < channels; i++) {
        for (size_t j = 0; j < frLength; j++) {
          auto sample = in[i][j];
          if (sample < -0x800000) {
            sample = -0x800000;
          } else if (sample > 0x800000) {
            sample = 0x800000;
          }
          ((uint8_t*)pcm_out + (i + j * channels) * 3)[0] = ((uint8_t*)(&sample))[0];
          ((uint8_t*)pcm_out + (i + j * channels) * 3)[1] = ((uint8_t*)(&sample))[1];
          ((uint8_t*)pcm_out + (i + j * channels) * 3)[2] = ((uint8_t*)(&sample))[2];
        }
      }
    } break;
    case 32: {
      for (size_t i = 0; i < channels; i++) {
        for (size_t j = 0; j < frLength; j++) {
          auto sample = in[i][j];
          if (sample < -0x800000) {
            sample = -0x800000;
          } else if (sample > 0x800000) {
            sample = 0x800000;
          }
          sample <<= 8;
          *(((int32_t*)pcm_out) + i + j * channels) = sample;
        }
      }
    } break;
    default:
      break;
  }
}

char x[100];
int main() {
  // bit_record = fopen("E:/codec/L2HC/bit_record3.bin", "wb");
  FILE* fp = NULL;
  fp = fopen("./8dw_enc.bin", "rb");

  int read_count = 0;
  int one_pack_size = 0;

  fseek(fp, 0, SEEK_END);

  int file_size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
#if 1
  size_t i = 0;
  wave_file_out = fopen("./8dw_tests.wav", "wb");
  memset(x, ' ', sizeof(x));
  x[5 + 52] = 0;
  int old_rate = -1;
  for (i = 0; read_count < file_size; i++) {
    fread(&one_pack_size, 1, 4, fp);
    fread(stream_buffer, 1, one_pack_size, fp);
    read_count += one_pack_size + 4;

    int rate = (int)(((double)read_count) * 100. / (double)file_size);
    if (rate != old_rate) {
      sprintf(&x[0], "%3d", rate);
      x[3] = '%';
      x[4] = '[';
      x[5 + (rate >> 1)] = '=';
      x[5 + 51] = ']';
      printf("\r%s", x);
      fflush(stdout);
      old_rate = rate;
    }
    LLunpack(one_pack_size, i);
    /*
    if(i==2){
      break;
    }*/
  }
  LLdeinit();
  printf("\n");
  fseek(wave_file_out, 0, SEEK_SET);

  wave_write_header(wave_file_out, 32, 4, 48000, 2, i * 480);

  fclose(wave_file_out);
#else
  fread(&one_pack_size, 1, 4, fp);
  fread(stream_buffer, 1, one_pack_size, fp);
  read_count += one_pack_size + 4;
  LLunpack(one_pack_size, 0);
#endif
  fclose(fp);

  // Unpack(&test);
  return 0;
}