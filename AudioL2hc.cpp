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
extern uint64_t BIT_MASK[];
__int16 HUF_DEC_DIFF_SF[32][3] = {{4, 0, 4},  {4, 1, 4},  {5, 2, 5},  {6, 3, 5},  {3, 4, 3},  {3, 5, 3},  {3, 6, 3},
                                  {3, 7, 3},  {2, 8, 2},  {2, 9, 2},  {2, 10, 2}, {2, 11, 2}, {2, 12, 2}, {2, 13, 2},
                                  {2, 14, 2}, {2, 15, 2}, {1, 16, 2}, {1, 17, 2}, {1, 18, 2}, {1, 19, 2}, {1, 20, 2},
                                  {1, 21, 2}, {1, 22, 2}, {1, 23, 2}, {0, 24, 2}, {0, 25, 2}, {0, 26, 2}, {0, 27, 2},
                                  {0, 28, 2}, {0, 29, 2}, {0, 30, 2}, {0, 31, 2}};

uint64_t stream_buffer_in_QWORD[500];
AudioL2hc test = {0};
typedef struct {
  int i;  // QWORD 剩余可读个数
  int j;
  /* data */
} read_index;

uint64_t ReadBitsInQWORD(read_index* i64, uint64_t* stream_buffer, int nbits) {
  uint64_t result = 0;
  if (!nbits) {
    return result;
  }

  if (i64->i < nbits) {
    int Remaining_bits = nbits - i64->i;
    result = stream_buffer[i64->j] & BIT_MASK[i64->i];
    result <<= Remaining_bits;
    if ((i64->j + 1) * sizeof(uint64_t) <= 1376) {
      i64->j += 1;
      i64->i = 64 - Remaining_bits;
      result |= (stream_buffer[i64->j] >> i64->i) & BIT_MASK[i64->i];
    } else {
      i64->i = 0;
    }
  } else {
    i64->i -= nbits;
    result = (stream_buffer[i64->j] >> i64->i) & BIT_MASK[nbits];
  }
  return result;
}

__int64 BytesToWords(uint8_t* stream_buffer, int stream_buffer_size, uint64_t* p_inner_stream_buffer, read_index* a4) {
  int stream_buffer_size_div_8;  // w8
  int index;                     // w9
  unsigned __int64 index_x8;     // x10
  unsigned __int8* v12;          // x11
  __int64 result;                // x0

  memset(p_inner_stream_buffer, 0, 0x560uLL);
  a4->i = 64;
  a4->j = 0;
  if (((stream_buffer_size & 0x8000u) != 0) | (((__int16)stream_buffer_size & 0x80000007) == 0))
    stream_buffer_size_div_8 =
        (int)((stream_buffer_size + (((unsigned int)(__int16)stream_buffer_size >> 28) & 7)) << 16) >> 19;
  else
    stream_buffer_size_div_8 =
        ((int)((stream_buffer_size + (((unsigned int)(__int16)stream_buffer_size >> 28) & 7)) << 16) >> 19) + 1;
  if (stream_buffer_size_div_8 >= 1) {
    index = 0;
    do {
      index_x8 = 8LL * index;  // < pcm_buffer_size
      index = (__int16)(index + 1);
      v12 = (unsigned __int8*)&stream_buffer[(int)index_x8];
      p_inner_stream_buffer[index_x8 / 8] = ((unsigned __int64)v12[0] << 56) | ((unsigned __int64)v12[1] << 48) |
                                            ((unsigned __int64)v12[2] << 40) | ((unsigned __int64)v12[3] << 32) |
                                            ((unsigned __int64)v12[4] << 24) | ((unsigned __int64)v12[5] << 16) |
                                            ((unsigned __int64)v12[6] << 8) | v12[7];
    } while (stream_buffer_size_div_8 > index);
  }
  result = 0LL;

  return result;
}

#if 0
void MdctSNS1(void *result, void *pAudioL2hc_add32)
{
  int channels_0; // w24
  __int64 bandNum; // x26
  _WORD *pAudioL2hc_add4440; // x20
  __int64 channels_0_1; // x27
  __int16 bandNum_1; // w25
  __int16 *p_sf; // x22
  __int16 (*p_sf_inner)[32]; // x23
  __int64 channels_0_index_1; // x28
  float v11; // s0
  __int64 bandIndex; // x8
  __int16 *v13; // x10
  __int64 i; // x11
  __int16 *v15; // x14
  int sf_1ch; // w13
  int sf_2ch; // w15
  int v18; // w12
  __int64 channels_0_index; // x8
  __int64 v20; // x11
  __int16 (*v21)[32]; // x13
  __int16 sf_sum; // w14
  __int16 (*v23)[32]; // x15
  __int64 bandIndex_1; // x16
  __int16 v25; // t1
  __int16 sf_avg_sub1; // w15
  __int16 *v27; // x14
  int v28; // w16
  int sf_avg; // w15
  __int16 v30; // w0
  float maskScale; // s4
  int v32; // w0
  int v33; // w17
  int v34; // w16
  __int64 v35; // x17
  __int64 v36; // x1
  __int16 v37; // w2
  int v38; // w2
  int v39; // w16
  int v40; // w0
  float maskScale_1; // s4
  __int64 v42; // x15
  __int16 v43; // w14
  __int16 v44; // w16
  __int64 v45; // x15
  int v46; // w14
  double v47; // d4
  int v48; // w16
  double v49; // d5
  __int64 channels_0_index_2; // x8
  _WORD *p_dr_adjust; // x9
  __int16 (*v52)[32]; // x10
  __int16 (*v53)[32]; // x11
  _WORD *v54; // x12
  __int64 v55; // x13
  __int64 bandIndex_2; // x16
  __int16 v57; // t1
  __int16 sf_new[2][32]; // [xsp+0h] [xbp-80h] BYREF


  channels_0 = *((__int16 *)pAudioL2hc_add32 + 2);
  bandNum = *((unsigned __int16 *)pAudioL2hc_add32 + 103);
  pAudioL2hc_add4440 = result;
  channels_0_1 = *((unsigned __int16 *)pAudioL2hc_add32 + 2);
  bandNum_1 = *((_WORD *)pAudioL2hc_add32 + 103);
  memset(sf_new, 0, sizeof(sf_new));
  if ( channels_0 >= 1 )
  {
    p_sf = (__int16 *)((char *)result + 7680);
    p_sf_inner = sf_new;
    channels_0_index_1 = channels_0_1;
    do
    {
      if ( (__int16)bandNum >= 1 )
        result = memcpy(p_sf_inner, p_sf, 2 * bandNum);
      --channels_0_index_1;
      ++p_sf_inner;
      p_sf += 32;
    }
    while ( channels_0_index_1 );
  }
  v11 = 0.75;
  if ( !pAudioL2hc_add4440[5894] )              // msFlag
    v11 = 0.5;
  if ( (__int16)bandNum >= 1 )
  {
    bandIndex = 0LL;
    v13 = sf_new[1];
    for ( i = bandNum; i; --i )
    {
      sf_1ch = *(v13 - 32);
      sf_2ch = *v13;
      v18 = sf_1ch - sf_2ch;
      if ( sf_1ch > sf_2ch )
      {
        v15 = v13;
      }
      else
      {
        v18 = sf_2ch - sf_1ch;
        if ( sf_2ch <= sf_1ch )
          goto LABEL_12;
        v15 = &sf_new[0][bandIndex];
        LOWORD(sf_1ch) = *v13;
      }
      *v15 = sf_1ch - (int)(float)(v11 * (float)v18);
LABEL_12:
      ++bandIndex;
      ++v13;
    }
  }
  if ( channels_0 >= 1 )
  {
    channels_0_index = 0LL;
    v20 = (unsigned int)((__int16)bandNum - 1);
    v21 = sf_new;
    do                                          // get sf_new
    {
      if ( (__int16)bandNum >= 1 )
      {
        sf_sum = 0;
        v23 = v21;
        bandIndex_1 = bandNum;
        do
        {
          v25 = *(_WORD *)v23;
          v23 = (__int16 (*)[32])((char *)v23 + 2);
          --bandIndex_1;
          sf_sum += v25;
        }
        while ( bandIndex_1 );
        sf_avg_sub1 = sf_sum / (__int16)bandNum;
        v27 = sf_new[channels_0_index];
        v28 = *v27;
        sf_avg = (__int16)(sf_avg_sub1 + 1);
        v30 = v28 - 2;
        if ( sf_avg >= v28 )
          maskScale = 0.25;
        else
          maskScale = 0.375;
        if ( (__int16)bandNum != 1 )
          v30 = v27[1];
        v32 = v30 - v28;
        v33 = (__int16)(v28 - 2) - v28;
        if ( v32 <= 0 )
          v32 = 0;
        if ( v33 <= 0 )
          v33 = 0;
        v34 = (int)(float)((float)(maskScale * (float)v33) + (float)((float)(maskScale * (float)v32) + (float)v28));
        *v27 = v34;
        if ( (__int16)bandNum >= 2 )
        {
          v35 = 1LL;
          do
          {
            v40 = v27[v35];
            if ( sf_avg >= v40 )
              maskScale_1 = 0.25;
            else
              maskScale_1 = 0.375;
            if ( v20 == v35 )
            {
              v37 = v40 - 2;
              v36 = v20 + 1;
            }
            else
            {
              v36 = v35 + 1;
              v37 = v27[v35 + 1];
            }
            v38 = v37 - v40;
            v39 = (__int16)v34 - v40;
            if ( v38 <= 0 )
              v38 = 0;
            if ( v39 <= 0 )
              v39 = 0;
            v34 = (int)(float)((float)(maskScale_1 * (float)v39)
                             + (float)((float)(maskScale_1 * (float)v38) + (float)v40));
            v27[v35] = v34;
            v35 = v36;
          }
          while ( v36 < (__int16)bandNum );
        }
        v42 = 0LL;
        v43 = 0;
        do
        {
          v44 = (*v21)[v42++];
          v43 += v44;
        }
        while ( bandNum != v42 );
        v45 = 0LL;
        v46 = (__int16)(v43 / (__int16)bandNum + 1);
        v47 = (double)v46;
        do
        {
          v48 = (*v21)[v45];
          result = (void *)(unsigned int)(v48 - v46);
          v49 = v47 + (double)(v46 - v48) * -0.5;
          if ( v46 < v48 )
            v49 = (double)(int)result * 0.375 + v47;
          (*v21)[v45++] = (int)v49;
        }
        while ( (__int16)bandNum != v45 );
      }
      ++channels_0_index;
      ++v21;
    }
    while ( channels_0_index != channels_0_1 );
    channels_0_index_2 = 0LL;
    p_dr_adjust = pAudioL2hc_add4440 + 3904;
    v52 = sf_new;
    do
    {
      if ( bandNum_1 >= 1 )
      {
        v53 = v52;
        v54 = p_dr_adjust;
        v55 = 3840LL;
        do
        {
          bandIndex_2 = v55 - 3839;
          v57 = *(_WORD *)v53;
          v53 = (__int16 (*)[32])((char *)v53 + 2);
          ++v55;
          *v54 = *(v54 - 64) - v57;             // dr_adjust[b]=sf[b]-sf_new[b]
          ++v54;
        }
        while ( bandIndex_2 < *((__int16 *)pAudioL2hc_add32 + 103) );
        bandNum_1 = *((_WORD *)pAudioL2hc_add32 + 103);
        LOWORD(channels_0) = *((_WORD *)pAudioL2hc_add32 + 2);
      }
      ++channels_0_index_2;
      p_dr_adjust += 32;
      ++v52;
    }
    while ( channels_0_index_2 < (__int16)channels_0 );
  }
}

void MdctSNS2() {}
#endif

void Unpack(AudioL2hc* a) {
  FILE* fp = NULL;
  fp = fopen("E:/codec/L2HC/48kS32w_enc.bin", "rb");
  int one_pack_size = 0;
  fread(&one_pack_size, 1, 4, fp);
  fread(stream_buffer, 1, one_pack_size, fp);
  read_index i64 = {.i = 2, .j = 0};
  uint64_t stream[] = {
      0b0000000000000000000000000000000000000000000000000000000000000011,
      0b1110000000000000000000000000000000000000000000000000000000000111,
      0b1100000000000000000000000000000000000000000000000000000000000000};
  printf("test %d\n", ReadBitsInQWORD(&i64, stream, 5));
  printf("test %d\n", ReadBitsInQWORD(&i64, stream, 58));
  printf("test %d\n", ReadBitsInQWORD(&i64, stream, 5));
  BytesToWords(stream_buffer, one_pack_size, stream_buffer_in_QWORD, &i64);
  printf("AudioL2hc size:%ld\n", sizeof(AudioL2hc));
  a->codectype_index_from_stream = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
  printf("codecType   %ld\n", a->codectype_index_from_stream);
  a->smpRate_index_from_stream = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
  printf("sampleRate  %ld\n", smpRate_table[a->smpRate_index_from_stream]);
  a->channels_index_from_stream = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 1);
  a->channels_0 = a->channels_index_from_stream + 1;
  printf("chNum       %ld\n", a->channels_0);
  a->frameLength_index_from_stream = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
  a->frLength_0 = frLength_table[a->frameLength_index_from_stream];
  printf("frameLength %ld\n", a->frLength_0);
  a->lowBrFlag = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 1);
  printf("lowBrFlag   %ld\n", a->lowBrFlag);
  a->msFlag = 0;
  if (a->channels_0 == 2) {
    a->msFlag = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 1);
  }
  printf("msFlag      %ld\n", a->msFlag);
  a->dr = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 5);
  printf("dr          %ld\n", a->dr);
  a->drQuater = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 3);
  printf("drQuater    %ld\n", a->drQuater);
  a->sfId = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 3);
  printf("sfId        %ld\n", a->sfId);
  a->bandNum = 32;
  if (a->lowBrFlag == 1) {
    a->bandNum = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 4);
  }
  printf("bandNum     %ld\n", a->bandNum);
  a->diffFlag = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
  printf("diffFlag    %ld\n", a->diffFlag);
  uint32_t channels_0_index = 0;
  uint32_t bandIndex = 0;
  if (a->diffFlag != 1) {
    for (channels_0_index = 0; channels_0_index < a->channels_0; channels_0_index++) {
      for (bandIndex = 0; bandIndex < a->bandNum; bandIndex++) {
        a->sf[channels_0_index][bandIndex] = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 5);
        a->sf[channels_0_index][bandIndex] -= 8;
      }
    }
  } else {
    for (channels_0_index = 0; channels_0_index < a->channels_0; channels_0_index++) {
      a->sf[channels_0_index][0] = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 5);
      a->sf[channels_0_index][0] -= 8;
      for (bandIndex = 1; bandIndex < a->bandNum; bandIndex++) {
        int64 HUF_DEC_DIFF_SF_index = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 5);
        int sfValueDiff = HUF_DEC_DIFF_SF[HUF_DEC_DIFF_SF_index][0];
        ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, HUF_DEC_DIFF_SF[HUF_DEC_DIFF_SF_index][2]);
        if (sfValueDiff) {
          int64 sfSign = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 1);
          if (sfSign == 0) {
            a->sf[channels_0_index][bandIndex] = a->sf[channels_0_index][bandIndex - 1] - sfValueDiff;
          } else {
            a->sf[channels_0_index][bandIndex] = a->sf[channels_0_index][bandIndex - 1] + sfValueDiff;
          }
        } else {
          a->sf[channels_0_index][bandIndex] = a->sf[channels_0_index][bandIndex - 1];
        }
      }
    }
  }
  if (a->channels_0 >= 1) {
    a->hufTupId_4tuple_1ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
    a->hufTupId_2tuple_1ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
    a->hufTupId_1tuple_1ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
    if (a->channels_0 == 2) {
      a->hufTupId_4tuple_2ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
      a->hufTupId_2tuple_2ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
      a->hufTupId_1tuple_2ch = ReadBitsInQWORD(&i64, stream_buffer_in_QWORD, 2);
    }
  }
  if (a->lowBrFlag < 1) {
    memset(a->dr_adjust, 0, sizeof(a->dr_adjust));
  } else if (a->channels_0 == 2) {
    // MdctSNS1(&a->sth1,&a->smpRate_0);
  } else {
    // MdctSNS2();
  }
}