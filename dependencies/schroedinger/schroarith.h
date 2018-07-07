
#ifndef _SCHRO_ARITH_H_
#define _SCHRO_ARITH_H_

#define SCHRO_INTERNAL

#include <stdint.h>

#if __cplusplus
extern "C" {
#endif

/* This definition of SchroBuffer is simplified version of the
 * real version for the purpose of this integration test */
typedef struct _SchroBuffer SchroBuffer;
struct _SchroBuffer {
	unsigned char *data;
	unsigned int length;
};


typedef struct _SchroArith SchroArith;

struct _SchroArith {
  SchroBuffer *buffer;
  uint8_t *dataptr;
  uintptr_t offset;

  uint32_t range[2];
  uint32_t code;
  uint32_t range_size;
  int cntr;
  int carry;

  uint16_t lut[512];
};

void schro_arith_decode_init (SchroArith *arith, SchroBuffer *buffer);
void schro_arith_encode_init (SchroArith *arith, SchroBuffer *buffer);
void schro_arith_flush (SchroArith *arith);
void schro_arith_decode_flush (SchroArith *arith);

void schro_arith_encode_bit (SchroArith *arith, uint16_t *probability,  int value);

int schro_arith_decode_bit (SchroArith *arith, uint16_t *probability);

#ifdef SCHRO_ARITH_DEFINE_INLINE
static inline int
_schro_arith_decode_bit (SchroArith *arith, uint16_t *probability)
{
  unsigned int range_x_prob;
  unsigned int value;
  unsigned int lut_index;
  register unsigned int range = arith->range[1];
  register unsigned int code_minus_low = arith->code;

  while (range <= 0x40000000) {

    range <<= 1;
    code_minus_low <<= 1;

    if (!--arith->cntr) {
      arith->offset++;
      if (arith->offset < arith->buffer->length) {
        code_minus_low |= arith->dataptr[arith->offset] << 8;
      } else {
        code_minus_low |= 0xff00;
      }

      arith->offset++;
      if (arith->offset < arith->buffer->length) {
        code_minus_low |= arith->dataptr[arith->offset];
      } else {
        code_minus_low |= 0xff;
      }

      arith->cntr = 16;
    }
  }

  range_x_prob = ((range >> 16) * (*probability)) & 0xFFFF0000;
  lut_index = (*probability)>>7 & ~1;

  value = (code_minus_low  >= range_x_prob);
  (*probability) += arith->lut[lut_index | value];

  if (value) {
    code_minus_low -= range_x_prob;
    range -= range_x_prob;
  } else {
    range = range_x_prob;
  }

  arith->range[1] = range;
  arith->code = code_minus_low;

  return value;
}

static inline void
_schro_arith_encode_bit (SchroArith *arith, uint16_t *probability, int value)
{
  unsigned int range;
  unsigned int probability0;
  unsigned int range_x_prob;

  probability0 = (*probability);
  range = arith->range[1];
  range_x_prob = (range * probability0) >> 16;

  if (value) {
    arith->range[0] = arith->range[0] + range_x_prob;
    arith->range[1] -= range_x_prob;
    (*probability) -= arith->lut[(*probability)>>8];
  } else {
    arith->range[1] = range_x_prob;
    (*probability) += arith->lut[255-((*probability)>>8)];
  }

  while (arith->range[1] <= 0x4000) {
    arith->range[0] <<= 1;
    arith->range[1] <<= 1;
    arith->cntr++;

    if (arith->cntr == 8) {
      if (arith->range[0] < (1<<24) &&
          (arith->range[0] + arith->range[1]) >= (1<<24)) {
        arith->carry++;
      } else {
        if (arith->range[0] >= (1<<24)) {
          arith->dataptr[arith->offset-1]++;
          while (arith->carry) {
            arith->dataptr[arith->offset] = 0x00;
            arith->carry--;
            arith->offset++;
          }
        } else {
          while (arith->carry) {
            arith->dataptr[arith->offset] = 0xff;
            arith->carry--;
            arith->offset++;
          }
        }
        arith->dataptr[arith->offset] = arith->range[0] >> 16;
        arith->offset++;
      }

      arith->range[0] &= 0xffff;
      arith->cntr = 0;
    }
  }
}

static inline int
maxbit (unsigned int x)
{
#if 0
  int i;
  for(i=0;x;i++){
    x >>= 1;
  }
  return i;
#else
  int i = 0;
  if (x == 0) return 0;
  if (x > 0x00ff) {
    i += 8;
    x >>= 8;
  }
  if (x > 0x000f) {
    i += 4;
    x >>= 4;
  }
  if (x > 0x0003) {
    i += 2;
    x >>= 2;
  }
  if (x > 0x0001) {
    i += 1;
    x >>= 1;
  }
  if (x > 0x0000) {
    i += 1;
  }
  return i;
#endif
}

#else /* SCHRO_ARITH_DEFINE_INLINE */
int _schro_arith_decode_bit (SchroArith *arith, uint16_t *probability);
void _schro_arith_encode_bit (SchroArith *arith, uint16_t *probability, int
    value) SCHRO_INTERNAL;

#endif /* SCHRO_ARITH_DEFINE_INLINE */

#if __cplusplus
} // extern "C"
#endif

#endif


