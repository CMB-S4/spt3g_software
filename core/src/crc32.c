/* NOTE: THIS IS MODIFIED FROM THE ZLIB VERSION IN ORDER TO COMPUTE A CRC32-C
 * INSTEAD OF A CRC32.
 *
 * In addition, an SSE4.2 implementation has been provided to compute the CRC
 * in hardware.
 */

#include <stddef.h>
#include <stdint.h>

#define DYNAMIC_CRC_TABLE
#define STDC
#define FAR
#define OF(x) x
#define Z_NULL NULL
typedef unsigned int uInt;
typedef unsigned long uLong;

/* crc32.c -- compute the CRC-32 of a data stream
 * Copyright (C) 1995-2006, 2010 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
 *
 * Thanks to Rodney Brown <rbrown64@csc.com.au> for his contribution of faster
 * CRC methods: exclusive-oring 32 bits of data at a time, and pre-computing
 * tables for updating the shift register in one step with three exclusive-ors
 * instead of four steps with four exclusive-ors.  This results in about a
 * factor of two increase in speed on a Power PC G4 (PPC7455) using gcc -O3.
 */

/* Copyright (C) 1995-2010 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu


  The data format used by the zlib library is described by RFCs (Request for
  Comments) 1950 to 1952 in the files http://www.ietf.org/rfc/rfc1950.txt
  (zlib format), rfc1951.txt (deflate format) and rfc1952.txt (gzip format).
*/

/* @(#) $Id: crc32.c 1437 2012-11-06 19:57:38Z nega $ */

/*
  Note on the use of DYNAMIC_CRC_TABLE: there is no mutex or semaphore
  protection on the static variables used to control the first-use generation
  of the crc tables.  Therefore, if you #define DYNAMIC_CRC_TABLE, you should
  first call get_crc_table() to initialize the tables before allowing more than
  one thread to use crc32().
 */

#ifdef MAKECRCH
#  include <stdio.h>
#  ifndef DYNAMIC_CRC_TABLE
#    define DYNAMIC_CRC_TABLE
#  endif /* !DYNAMIC_CRC_TABLE */
#endif /* MAKECRCH */

#define local static

/* Find a four-byte integer type for crc32_little() and crc32_big(). */
#ifndef NOBYFOUR
#  ifdef STDC           /* need ANSI C limits.h to determine sizes */
#    include <limits.h>
#    define BYFOUR
#    if (UINT_MAX == 0xffffffffUL)
       typedef unsigned int u4;
#    else
#      if (ULONG_MAX == 0xffffffffUL)
         typedef unsigned long u4;
#      else
#        if (USHRT_MAX == 0xffffffffUL)
           typedef unsigned short u4;
#        else
#          undef BYFOUR     /* can't find a four-byte integer type! */
#        endif
#      endif
#    endif
#  endif /* STDC */
#endif /* !NOBYFOUR */

/* Definitions for doing the crc four data bytes at a time. */
#ifdef BYFOUR
#  define REV(w) ((((w)>>24)&0xff)+(((w)>>8)&0xff00)+ \
                (((w)&0xff00)<<8)+(((w)&0xff)<<24))
   local unsigned long crc32_little OF((unsigned long,
                        const unsigned char FAR *, unsigned));
   local unsigned long crc32_big OF((unsigned long,
                        const unsigned char FAR *, unsigned));
#  define TBLS 8
#else
#  define TBLS 1
#endif /* BYFOUR */

#ifdef DYNAMIC_CRC_TABLE

local volatile int crc_table_empty = 1;
local unsigned long FAR crc_table[TBLS][256];
local void make_crc_table OF((void));
#ifdef MAKECRCH
   local void write_table OF((FILE *, const unsigned long FAR *));
#endif /* MAKECRCH */
/*
  Generate tables for a byte-wise 32-bit CRC calculation on the polynomial:
  x^32+x^28+x^27+x^26+x^25+x^23+x^22+x^20+x^19+x^18+x^14+x^13+x^11+x^10+x^9+x^8+x^6+1.

  Polynomials over GF(2) are represented in binary, one bit per coefficient,
  with the lowest powers in the most significant bit.  Then adding polynomials
  is just exclusive-or, and multiplying a polynomial by x is a right shift by
  one.  If we call the above polynomial p, and represent a byte as the
  polynomial q, also with the lowest power in the most significant bit (so the
  byte 0xb1 is the polynomial x^7+x^3+x+1), then the CRC is (q*x^32) mod p,
  where a mod b means the remainder after dividing a by b.

  This calculation is done using the shift-register method of multiplying and
  taking the remainder.  The register is initialized to zero, and for each
  incoming bit, x^32 is added mod p to the register if the bit is a one (where
  x^32 mod p is p+x^32 = x^26+...+1), and the register is multiplied mod p by
  x (which is shifting right by one and adding x^32 mod p if the bit shifted
  out is a one).  We start with the highest power (least significant bit) of
  q and repeat for all eight bits of q.

  The first table is simply the CRC of all possible eight bit values.  This is
  all the information needed to generate CRCs on data a byte at a time for all
  combinations of CRC register values and incoming bytes.  The remaining tables
  allow for word-at-a-time CRC calculation for both big-endian and little-
  endian machines, where a word is four bytes.
*/

local void make_crc_table()
{
    unsigned long c;
    unsigned int n, k;
    unsigned long poly;                 /* polynomial exclusive-or pattern */
    /* terms of polynomial defining this crc (except x^32): */
    static volatile int first = 1;      /* flag to limit concurrent making */
    static const unsigned char p[] = {0,6,8,9,10,11,13,14,18,19,20,22,23,25,26,27,28};

    /* See if another task is already doing this (not thread-safe, but better
       than nothing -- significantly reduces duration of vulnerability in
       case the advice about DYNAMIC_CRC_TABLE is ignored) */
    if (first) {
        first = 0;

        /* make exclusive-or pattern from polynomial (0x82F63B78UL) */
        poly = 0UL;
        for (n = 0; n < sizeof(p)/sizeof(unsigned char); n++)
            poly |= 1UL << (31 - p[n]);

        /* generate a crc for every 8-bit value */
        for (n = 0; n < 256; n++) {
            c = (unsigned long)n;
            for (k = 0; k < 8; k++)
                c = c & 1 ? poly ^ (c >> 1) : c >> 1;
            crc_table[0][n] = c;
        }

#ifdef BYFOUR
        /* generate crc for each value followed by one, two, and three zeros,
           and then the byte reversal of those as well as the first table */
        for (n = 0; n < 256; n++) {
            c = crc_table[0][n];
            crc_table[4][n] = REV(c);
            for (k = 1; k < 4; k++) {
                c = crc_table[0][c & 0xff] ^ (c >> 8);
                crc_table[k][n] = c;
                crc_table[k + 4][n] = REV(c);
            }
        }
#endif /* BYFOUR */

        crc_table_empty = 0;
    }
    else {      /* not first */
        /* wait for the other guy to finish (not efficient, but rare) */
        while (crc_table_empty)
            ;
    }

#ifdef MAKECRCH
    /* write out CRC tables to crc32.h */
    {
        FILE *out;

        out = fopen("crc32.h", "w");
        if (out == NULL) return;
        fprintf(out, "/* crc32.h -- tables for rapid CRC calculation\n");
        fprintf(out, " * Generated automatically by crc32.c\n */\n\n");
        fprintf(out, "local const unsigned long FAR ");
        fprintf(out, "crc_table[TBLS][256] =\n{\n  {\n");
        write_table(out, crc_table[0]);
#  ifdef BYFOUR
        fprintf(out, "#ifdef BYFOUR\n");
        for (k = 1; k < 8; k++) {
            fprintf(out, "  },\n  {\n");
            write_table(out, crc_table[k]);
        }
        fprintf(out, "#endif\n");
#  endif /* BYFOUR */
        fprintf(out, "  }\n};\n");
        fclose(out);
    }
#endif /* MAKECRCH */
}

#ifdef MAKECRCH
local void write_table(out, table)
    FILE *out;
    const unsigned long FAR *table;
{
    int n;

    for (n = 0; n < 256; n++)
        fprintf(out, "%s0x%08lxUL%s", n % 5 ? "" : "    ", table[n],
                n == 255 ? "\n" : (n % 5 == 4 ? ",\n" : ", "));
}
#endif /* MAKECRCH */

#else /* !DYNAMIC_CRC_TABLE */
/* ========================================================================
 * Tables of CRC-32s of all single-byte values, made by make_crc_table().
 */
#include "crc32.h"
#endif /* DYNAMIC_CRC_TABLE */

/* =========================================================================
 * This function can be used by asm versions of crc32()
 */
const unsigned long FAR * get_crc32c_table()
{
#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
        make_crc_table();
#endif /* DYNAMIC_CRC_TABLE */
    return (const unsigned long FAR *)crc_table;
}

/* ========================================================================= */
#define DO1 crc = crc_table[0][((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8)
#define DO8 DO1; DO1; DO1; DO1; DO1; DO1; DO1; DO1

/* ========================================================================= */
#if defined(__i386__) || defined (__x86_64__)
local unsigned long crc32c_tabular(crc, buf, len)
#else
unsigned long crc32c(crc, buf, len)
#endif
    unsigned long crc;
    const unsigned char FAR *buf;
    uInt len;
{
    if (buf == Z_NULL) return 0UL;

#ifdef DYNAMIC_CRC_TABLE
    if (crc_table_empty)
        make_crc_table();
#endif /* DYNAMIC_CRC_TABLE */

#ifdef BYFOUR
    if (sizeof(void *) == sizeof(ptrdiff_t)) {
        u4 endian;

        endian = 1;
        if (*((unsigned char *)(&endian)))
            return crc32_little(crc, buf, len);
        else
            return crc32_big(crc, buf, len);
    }
#endif /* BYFOUR */

    while (len >= 8) {
        DO8;
        len -= 8;
    }
    if (len) do {
        DO1;
    } while (--len);
    return crc;
}

#ifdef BYFOUR

/* ========================================================================= */
#define DOLIT4 c ^= *buf4++; \
        c = crc_table[3][c & 0xff] ^ crc_table[2][(c >> 8) & 0xff] ^ \
            crc_table[1][(c >> 16) & 0xff] ^ crc_table[0][c >> 24]
#define DOLIT32 DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4; DOLIT4

/* ========================================================================= */
local unsigned long crc32_little(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    unsigned len;
{
    register u4 c;
    register const u4 FAR *buf4;

    c = (u4)crc;
    while (len && ((ptrdiff_t)buf & 3)) {
        c = crc_table[0][(c ^ *buf++) & 0xff] ^ (c >> 8);
        len--;
    }

    buf4 = (const u4 FAR *)(const void FAR *)buf;
    while (len >= 32) {
        DOLIT32;
        len -= 32;
    }
    while (len >= 4) {
        DOLIT4;
        len -= 4;
    }
    buf = (const unsigned char FAR *)buf4;

    if (len) do {
        c = crc_table[0][(c ^ *buf++) & 0xff] ^ (c >> 8);
    } while (--len);
    return (unsigned long)c;
}

/* ========================================================================= */
#define DOBIG4 c ^= *++buf4; \
        c = crc_table[4][c & 0xff] ^ crc_table[5][(c >> 8) & 0xff] ^ \
            crc_table[6][(c >> 16) & 0xff] ^ crc_table[7][c >> 24]
#define DOBIG32 DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4; DOBIG4

/* ========================================================================= */
local unsigned long crc32_big(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    unsigned len;
{
    register u4 c;
    register const u4 FAR *buf4;

    c = REV((u4)crc);
    while (len && ((ptrdiff_t)buf & 3)) {
        c = crc_table[4][(c >> 24) ^ *buf++] ^ (c << 8);
        len--;
    }

    buf4 = (const u4 FAR *)(const void FAR *)buf;
    buf4--;
    while (len >= 32) {
        DOBIG32;
        len -= 32;
    }
    while (len >= 4) {
        DOBIG4;
        len -= 4;
    }
    buf4++;
    buf = (const unsigned char FAR *)buf4;

    if (len) do {
        c = crc_table[4][(c >> 24) ^ *buf++] ^ (c << 8);
    } while (--len);
    return (unsigned long)(REV(c));
}

#endif /* BYFOUR */

#if defined(__i386__) || defined (__x86_64__)
local unsigned long crc32c_sse42(crc, buf, len)
    unsigned int crc;
    const unsigned char FAR *buf;
    unsigned len;
{
    /* From: http://byteworm.com/2010/10/13/crc32 */
    uint32_t q = len/sizeof(long);
    uint32_t r = len % sizeof(long);
    long *p = (long *)buf;

    while (q--) {
#ifdef __x86_64__
        __asm__ __volatile__(".byte 0xf2,0x48,0xf,0x38,0xf1,0xf1;" : "=S"(crc)
            : "0"(crc), "c"(*p));
#else
        __asm__ __volatile__(".byte 0xf2, 0xf, 0x38, 0xf1, 0xf1;" : "=S"(crc)
            : "0"(crc), "c"(*p));
#endif
        p++;
    }

    buf = (unsigned char *)p;
    while (r--) {
        __asm__ __volatile__(".byte 0xf2, 0xf, 0x38, 0xf0, 0xf1;" : "=S"(crc)
            : "0"(crc), "c"(*buf));
        buf++;
    }

    return crc;
}

unsigned long crc32c(crc, buf, len)
    unsigned long crc;
    const unsigned char FAR *buf;
    unsigned len;
{
    static int has_sse42 = -1;

    if (has_sse42 == -1) {
	unsigned int a,b,c,d;
#ifdef __i386__
	/* Workaround use of EBX as the PIC base register */
        __asm__ __volatile__("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1" :
	    "=a"(a), "=r"(b), "=c"(c), "=d"(d) : "0"(1));
#else
        __asm__ __volatile__("cpuid" : "=a"(a), "=b"(b), "=c"(c), "=d"(d) :
            "0"(1));
#endif
        has_sse42 = ((c & 0x00100000) != 0);
    }

    if (has_sse42)
        return crc32c_sse42(crc, buf, len);
    else
        return crc32c_tabular(crc, buf, len);
}
#endif

