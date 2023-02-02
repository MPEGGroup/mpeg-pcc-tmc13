
#ifndef _SCHRO_ARITH_H_
#define _SCHRO_ARITH_H_

#define SCHRO_INTERNAL

#include <stdint.h>

#if __cplusplus
extern "C" {
#endif

typedef uint8_t (*SchroRdFn)(void *);
typedef void (*SchroWrFn)(uint8_t byte, void *);

typedef struct _SchroArith SchroArith;

struct _SchroArith {
  SchroRdFn read;
  SchroWrFn write;
  void * io_priv;

  uint32_t range[2];
  uint32_t code;
  uint32_t range_size;
  int cntr;
  int carry;

  uint8_t first_byte;
  uint8_t output_byte;

  uint16_t lut[512];
};

void schro_arith_decode_init (SchroArith *arith, SchroRdFn read_fn, void * priv);
void schro_arith_encode_init (SchroArith *arith, SchroWrFn write_fn, void * priv);
void schro_arith_flush (SchroArith *arith);
void schro_arith_decode_flush (SchroArith *arith);

void schro_arith_encode_bit (SchroArith *arith, uint16_t *probability,  int value);

int schro_arith_decode_bit (SchroArith *arith, uint16_t *probability);

void schro_arith_encode_bypass_bit(SchroArith *arith, int value);

int schro_arith_decode_bypass_bit(SchroArith *arith);

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
    if (!--arith->cntr) {
      code_minus_low |= arith->read(arith->io_priv) << 8;
      arith->cntr = 8;
    }

    range <<= 1;
    code_minus_low <<= 1;
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
        // NB: carry cannot occur on the first byte since:
        //  - requires initial low=ff80, range=8000 for seq of False
        //  - requires initial low=8000, range=8000 for seq of True
        //  - requires initial low=0001, range=ffff for seq of True
        arith->carry++;
      } else {
        if (arith->range[0] >= (1<<24)) {
          // NB: output_byte is always valid here since:
          //  - given the initial value of low = 0, range = ffff
          //  - the largest value of low after the first bit=1 is ff00 (p=ff01)
          //  - coding subsequent symbols cannot cause low > 7fff00
          // => minimum initial value of low that can trigger this = 2
          //    (this is not a valid initial condition)
          //  - with range = 8000 (the maximum renormalised range)
          //  - minimum initial value of low to trigger = 8080
          arith->output_byte++;
          while (arith->carry) {
            arith->write(arith->output_byte, arith->io_priv);
            arith->output_byte = 0x00;
            arith->carry--;
          }
        } else {
          while (arith->carry) {
            arith->write(arith->output_byte, arith->io_priv);
            arith->output_byte = 0xff;
            arith->carry--;
          }
        }
        if (!arith->first_byte)
          arith->write(arith->output_byte, arith->io_priv);
        else
          arith->first_byte = 0;
        arith->output_byte = arith->range[0] >> 16;
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

static inline int
_schro_arith_decode_bypass_bit(SchroArith *arith)
{
  unsigned int value;
  register unsigned int range = arith->range[1];
  register unsigned int code_minus_low = arith->code;

  if (!--arith->cntr) {
    code_minus_low |= arith->read(arith->io_priv) << 8;
    arith->cntr = 8;
  }
  code_minus_low <<= 1;

  value = (code_minus_low >= range);

  if (value) {
    code_minus_low -= range;
  }
  arith->range[1] = range;
  arith->code = code_minus_low;

  return value;
}

static inline void
_schro_arith_encode_bypass_bit(SchroArith *arith, int value)
{
  uint32_t r0old = arith->range[0];
  uint32_t r1old = arith->range[1];
  int ctrold = arith->cntr;

  arith->cntr++;
  arith->range[0] <<= 1;
  if (value)
    arith->range[0] += arith->range[1];

  if (arith->cntr == 8) {
    if (arith->range[0] < (1 << 24) &&
      (arith->range[0] + arith->range[1]) >= (1 << 24)) {
      // NB: carry cannot occur on the first byte since:
      //  - requires initial low=ff80, range=8000 for seq of False
      //  - requires initial low=8000, range=8000 for seq of True
      //  - requires initial low=0001, range=ffff for seq of True
      arith->carry++;
    } else {
      if (arith->range[0] >= (1 << 24)) {
        // NB: output_byte is always valid here since:
        //  - given the initial value of low = 0, range = ffff
        //  - the largest value of low after the first bit=1 is ff00 (p=ff01)
        //  - coding subsequent symbols cannot cause low > 7fff00
        // => minimum initial value of low that can trigger this = 2
        //    (this is not a valid initial condition)
        //  - with range = 8000 (the maximum renormalised range)
        //  - minimum initial value of low to trigger = 8080
        arith->output_byte++;
        while (arith->carry) {
          arith->write(arith->output_byte, arith->io_priv);
          arith->output_byte = 0x00;
          arith->carry--;
        }
      }
      else {
        while (arith->carry) {
          arith->write(arith->output_byte, arith->io_priv);
          arith->output_byte = 0xff;
          arith->carry--;
        }
      }
      if (!arith->first_byte)
        arith->write(arith->output_byte, arith->io_priv);
      else
        arith->first_byte = 0;
      arith->output_byte = arith->range[0] >> 16;
    }

    arith->range[0] &= 0xffff;
    arith->cntr = 0;
  }
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


