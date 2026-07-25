#ifndef MICROSIM_ENDIAN_COMPAT_H_
#define MICROSIM_ENDIAN_COMPAT_H_

#include <stdint.h>
#include <string.h>

/*
 * Byte swapping is needed when writing big-endian legacy VTK output.
 * Do not depend on <endian.h>: it is a glibc/Linux header and is not
 * available on macOS.  Keeping the operation local also avoids undefined
 * behaviour from pointer-based type punning between double and int64_t.
 */
static inline uint64_t microsim_bswap64(uint64_t value) {
#if defined(__clang__) || defined(__GNUC__)
  return __builtin_bswap64(value);
#else
  return ((value & UINT64_C(0x00000000000000ff)) << 56) |
         ((value & UINT64_C(0x000000000000ff00)) << 40) |
         ((value & UINT64_C(0x0000000000ff0000)) << 24) |
         ((value & UINT64_C(0x00000000ff000000)) << 8)  |
         ((value & UINT64_C(0x000000ff00000000)) >> 8)  |
         ((value & UINT64_C(0x0000ff0000000000)) >> 24) |
         ((value & UINT64_C(0x00ff000000000000)) >> 40) |
         ((value & UINT64_C(0xff00000000000000)) >> 56);
#endif
}

static inline int microsim_is_big_endian(void) {
  const uint16_t value = UINT16_C(0x0102);
  unsigned char bytes[sizeof(value)];

  memcpy(bytes, &value, sizeof(value));
  return bytes[0] == UINT8_C(0x01);
}

static inline double microsim_swap_double_bytes(double value) {
  uint64_t bits;

  memcpy(&bits, &value, sizeof(bits));
  bits = microsim_bswap64(bits);
  memcpy(&value, &bits, sizeof(value));
  return value;
}

#endif
