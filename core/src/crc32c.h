#ifndef _G3_CRC32C_H
#define _G3_CRC32C_H

/*
 * CRC32C of G3 frames, with no pre- or post-inversion.  Defined by whichever
 * of crc32.c (portable), crc32c_sse.c (x86-64) or crc32c_neon.c (AArch64)
 * the build compiles into the core library.
 */

#ifdef __cplusplus
extern "C" {
#endif

unsigned long crc32c(unsigned long crc, const unsigned char *buf,
    unsigned int len);

/* Name of the implementation used by crc32c(), for diagnostics */
const char *crc32c_impl(void);

#ifdef __cplusplus
}
#endif

#endif
