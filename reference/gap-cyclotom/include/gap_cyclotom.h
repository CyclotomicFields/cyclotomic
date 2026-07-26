/*
 * Standalone extraction of GAP's cyclotomic arithmetic.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef GAP_CYCLOTOM_REFERENCE_H
#define GAP_CYCLOTOM_REFERENCE_H

#include <stddef.h>
#include <stdint.h>

typedef struct gap_cyc_ctx gap_cyc_ctx;
typedef struct gap_cyc gap_cyc;

gap_cyc_ctx * gap_cyc_ctx_new(void);
void gap_cyc_ctx_free(gap_cyc_ctx * ctx);
const char * gap_cyc_ctx_error(const gap_cyc_ctx * ctx);

gap_cyc * gap_cyc_from_terms(
    gap_cyc_ctx * ctx,
    uint32_t order,
    size_t len,
    const uint32_t * exponents,
    const int64_t * coefficients);

gap_cyc * gap_cyc_root(gap_cyc_ctx * ctx, uint32_t order, uint32_t exponent);
gap_cyc * gap_cyc_add(gap_cyc_ctx * ctx, const gap_cyc * lhs, const gap_cyc * rhs);
gap_cyc * gap_cyc_mul(gap_cyc_ctx * ctx, const gap_cyc * lhs, const gap_cyc * rhs);

void gap_cyc_free(gap_cyc * cyc);
uint32_t gap_cyc_order(const gap_cyc * cyc);
size_t gap_cyc_len(const gap_cyc * cyc);
uint32_t gap_cyc_exponent(const gap_cyc * cyc, size_t index);
int64_t gap_cyc_coefficient(const gap_cyc * cyc, size_t index);
int gap_cyc_equal(const gap_cyc * lhs, const gap_cyc * rhs);

#endif
