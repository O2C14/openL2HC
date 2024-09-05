#include <stdint.h>
typedef struct {
  int32_t index;
  int32_t data;
  int32_t nbits;
} huffdata;
#if 1
typedef struct  // aull dec
{
  int32_t Nband;               // 32
  int32_t isInterleaveStream;  // 56 AudioUnpackSideInfo
  int32_t UsedInt16;           // 60 AudioUnpackSideInfo
  int32_t out_channels;        // 64
  int32_t init_flag;           // 68
  int32_t codectype;           // 0  0
  int32_t smpRate2;            // 4  1
  int32_t bitPerSample;        // 8  2
  int32_t frLength;            // 12 3
  int32_t channels;            // 16 4
  int32_t a5;                  // 20 5  15b about CBR bitrate
  int32_t UsedCBR;             // 24 6  2b
  int32_t intrinsic_delay;     // 28 7  10b
  int32_t a8;                  // 32 8  2b
  int32_t a9;                  // 36 9  1b
  int32_t band_split_scale;    // 40 10 4b
  int32_t a11;                 // 44 11 1b
  int32_t a12;                 // 48 12 1b
  int32_t a13;                 // 52 13 1b
  int32_t a14;                 // 56 14
  int32_t ext_cfg;             // 60 15 1b

  int32_t* ImdctOutput;
  int32_t* DataFromStream;
  int32_t* ImdctPrevious;
  // end_about_pcm_buffer
  int32_t* BandId;
  int32_t* huffDecBand;
  // L2hcWindow AudioGetL2hcWindow
  int32_t* mdctHraWindow;  // 32
  // RotationParam
  int32_t* RotationTable_sin;
  int32_t* RotationTable_tan;
  // FoldingParam
  int32_t* FoldingParamA;
  int32_t* FoldingParamB;
  int32_t* FoldingParamC;
  // Dct4GetTwParam
  int32_t* TwiddleTableA;
  int32_t* TwiddleTableB;
  int32_t* TwiddleTableC;
  int32_t* TwiddleTableD;
} AudioLL_dec;

#endif
