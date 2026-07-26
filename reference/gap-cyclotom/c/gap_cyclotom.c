/*
 * Standalone extraction of the algorithms in GAP's src/cyclotom.c.
 *
 * GAP is licensed under GPL-2.0-or-later. This derived implementation keeps
 * the original function decomposition and Zumbroich-basis algorithms while
 * replacing GAP bags and Obj arithmetic with checked int64_t coefficients.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gap_cyclotom.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

struct gap_cyc_ctx {
    int64_t * result;
    size_t capacity;
    const char * error;
};

struct gap_cyc {
    uint32_t order;
    size_t len;
    int64_t * coefficients;
    uint32_t * exponents;
};

static void SetError(gap_cyc_ctx * ctx, const char * error)
{
    ctx->error = error;
}

static int AddChecked(gap_cyc_ctx * ctx, int64_t lhs, int64_t rhs, int64_t * out)
{
    if (__builtin_add_overflow(lhs, rhs, out)) {
        SetError(ctx, "coefficient addition overflow");
        return 0;
    }
    return 1;
}

static int SubChecked(gap_cyc_ctx * ctx, int64_t lhs, int64_t rhs, int64_t * out)
{
    if (__builtin_sub_overflow(lhs, rhs, out)) {
        SetError(ctx, "coefficient subtraction overflow");
        return 0;
    }
    return 1;
}

static int MulChecked(gap_cyc_ctx * ctx, int64_t lhs, int64_t rhs, int64_t * out)
{
    if (__builtin_mul_overflow(lhs, rhs, out)) {
        SetError(ctx, "coefficient multiplication overflow");
        return 0;
    }
    return 1;
}

static uint32_t Gcd(uint32_t lhs, uint32_t rhs)
{
    while (rhs != 0) {
        uint32_t next = lhs % rhs;
        lhs = rhs;
        rhs = next;
    }
    return lhs;
}

static int GrowResultCyc(gap_cyc_ctx * ctx, uint32_t order)
{
    if (ctx->capacity >= order)
        return 1;

    int64_t * grown = realloc(ctx->result, (size_t)order * sizeof(int64_t));
    if (grown == NULL) {
        SetError(ctx, "failed to grow cyclotomic scratch buffer");
        return 0;
    }

    memset(grown + ctx->capacity, 0, ((size_t)order - ctx->capacity) * sizeof(int64_t));
    ctx->result = grown;
    ctx->capacity = order;
    return 1;
}

static int ResetResultCyc(gap_cyc_ctx * ctx, uint32_t order)
{
    if (!GrowResultCyc(ctx, order))
        return 0;
    memset(ctx->result, 0, (size_t)order * sizeof(int64_t));
    return 1;
}

/*
 * ConvertToBase is the same rewriting used by GAP: eliminate every root not
 * in the Zumbroich basis using 1 + e_p + ... + e_p^(p-1) = 0.
 */
static int ConvertToBase(gap_cyc_ctx * ctx, uint32_t n)
{
    uint32_t nn = n;

    for (uint32_t p = 2; p <= nn; p += (p == 2 ? 1 : 2)) {
        if (nn % p != 0)
            continue;

        uint32_t q = p;
        while (nn / q % p == 0)
            q *= p;
        nn /= q;

        int64_t start_bad = p == 2 ? (int64_t)q / 2
                                   : -((int64_t)(q / p) - 1) / 2;
        int64_t end_bad = p == 2 ? (int64_t)q - 1
                                 : ((int64_t)(q / p) - 1) / 2;

        for (int64_t raw = start_bad; raw <= end_bad; raw++) {
            int64_t residue = ((int64_t)(n / q) * raw) % (int64_t)q;
            if (residue < 0)
                residue += q;

            for (uint32_t i = (uint32_t)residue; i < n; i += q) {
                int64_t coefficient = ctx->result[i];
                if (coefficient == 0)
                    continue;

                ctx->result[i] = 0;
                for (uint32_t k = 1; k < p; k++) {
                    uint32_t exponent = (i + (uint64_t)k * n / p) % n;
                    if (!SubChecked(
                            ctx,
                            ctx->result[exponent],
                            coefficient,
                            &ctx->result[exponent]))
                        return 0;
                }
            }
        }
    }
    return 1;
}

static void FieldProperties(
    uint32_t n,
    uint32_t * phi,
    int * squarefree,
    uint32_t * number_of_primes)
{
    *phi = n;
    *squarefree = 1;
    *number_of_primes = 0;
    uint32_t remaining = n;

    for (uint32_t p = 2; p <= remaining; p++) {
        if (remaining % p != 0)
            continue;
        *phi = *phi / p * (p - 1);
        (*number_of_primes)++;
        remaining /= p;
        if (remaining % p == 0)
            *squarefree = 0;
        while (remaining % p == 0)
            remaining /= p;
    }
}

/*
 * Cyclotomic packs ResultCyc and reduces it into the smallest cyclotomic
 * field. The `hint` has GAP's meaning: its prime divisors cannot be removed.
 */
static gap_cyc * Cyclotomic(gap_cyc_ctx * ctx, uint32_t n, uint32_t hint)
{
    size_t len = 0;
    uint32_t exponent_gcd = n;
    int coefficients_equal = 1;
    int have_coefficient = 0;
    int64_t common_coefficient = 0;

    for (uint32_t i = 0; i < n; i++) {
        if (ctx->result[i] == 0)
            continue;
        len++;
        exponent_gcd = Gcd(exponent_gcd, i);
        if (!have_coefficient) {
            common_coefficient = ctx->result[i];
            have_coefficient = 1;
        }
        else if (ctx->result[i] != common_coefficient) {
            coefficients_equal = 0;
        }
    }

    if (exponent_gcd > 1) {
        uint32_t reduced_n = n / exponent_gcd;
        for (uint32_t i = 1; i < reduced_n; i++) {
            ctx->result[i] = ctx->result[i * exponent_gcd];
            ctx->result[i * exponent_gcd] = 0;
        }
        n = reduced_n;
    }

    uint32_t phi;
    uint32_t number_of_primes;
    int squarefree;
    FieldProperties(n, &phi, &squarefree, &number_of_primes);

    if (len == phi && coefficients_equal && squarefree) {
        memset(ctx->result, 0, (size_t)n * sizeof(int64_t));
        if (number_of_primes % 2 != 0) {
            if (!SubChecked(ctx, 0, common_coefficient, &common_coefficient))
                return NULL;
        }
        ctx->result[0] = common_coefficient;
        n = 1;
        len = common_coefficient == 0 ? 0 : 1;
    }

    uint32_t divisor_gcd = Gcd(phi, (uint32_t)len);
    uint32_t nn = n;
    for (uint32_t p = 3; p <= nn && p - 1 <= divisor_gcd; p += 2) {
        if (nn % p != 0)
            continue;
        nn /= p;
        while (nn % p == 0)
            nn /= p;

        if (n % ((uint64_t)p * p) == 0 || len % (p - 1) != 0 || hint % p == 0)
            continue;

        int equal_classes = 1;
        for (uint32_t i = 0; i < n && equal_classes; i += p) {
            int64_t coefficient = ctx->result[(i + n / p) % n];
            for (uint32_t k = i + 2 * n / p; k < i + n; k += n / p) {
                if (ctx->result[k % n] != coefficient) {
                    equal_classes = 0;
                    break;
                }
            }
        }

        if (!equal_classes)
            continue;

        for (uint32_t i = 0; i < n; i += p) {
            int64_t coefficient = ctx->result[(i + n / p) % n];
            if (!SubChecked(ctx, 0, coefficient, &ctx->result[i]))
                return NULL;
            for (uint32_t k = i + n / p; k < i + n; k += n / p)
                ctx->result[k % n] = 0;
        }
        len /= p - 1;

        for (uint32_t i = 1; i < n / p; i++) {
            ctx->result[i] = ctx->result[i * p];
            ctx->result[i * p] = 0;
        }
        n /= p;
    }

    gap_cyc * cyc = calloc(1, sizeof(gap_cyc));
    if (cyc == NULL) {
        SetError(ctx, "failed to allocate cyclotomic");
        return NULL;
    }

    cyc->order = n;
    cyc->len = len;
    if (len != 0) {
        cyc->coefficients = malloc(len * sizeof(int64_t));
        cyc->exponents = malloc(len * sizeof(uint32_t));
        if (cyc->coefficients == NULL || cyc->exponents == NULL) {
            gap_cyc_free(cyc);
            SetError(ctx, "failed to allocate cyclotomic terms");
            return NULL;
        }
    }

    size_t index = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (ctx->result[i] == 0)
            continue;
        cyc->coefficients[index] = ctx->result[i];
        cyc->exponents[index] = i;
        ctx->result[i] = 0;
        index++;
    }
    cyc->len = index;
    return cyc;
}

static int FindCommonField(
    gap_cyc_ctx * ctx,
    uint32_t nl,
    uint32_t nr,
    uint32_t * ml,
    uint32_t * mr,
    uint32_t * n)
{
    uint32_t gcd = Gcd(nl, nr);
    uint64_t common = (uint64_t)nl * (nr / gcd);
    if (common > UINT32_MAX) {
        SetError(ctx, "common cyclotomic field exceeds uint32_t");
        return 0;
    }
    *n = (uint32_t)common;
    *ml = *n / nl;
    *mr = *n / nr;
    return GrowResultCyc(ctx, *n);
}

gap_cyc_ctx * gap_cyc_ctx_new(void)
{
    return calloc(1, sizeof(gap_cyc_ctx));
}

void gap_cyc_ctx_free(gap_cyc_ctx * ctx)
{
    if (ctx == NULL)
        return;
    free(ctx->result);
    free(ctx);
}

const char * gap_cyc_ctx_error(const gap_cyc_ctx * ctx)
{
    return ctx == NULL || ctx->error == NULL ? "" : ctx->error;
}

gap_cyc * gap_cyc_from_terms(
    gap_cyc_ctx * ctx,
    uint32_t order,
    size_t len,
    const uint32_t * exponents,
    const int64_t * coefficients)
{
    if (ctx == NULL || order == 0) {
        if (ctx != NULL)
            SetError(ctx, "cyclotomic order must be positive");
        return NULL;
    }
    ctx->error = NULL;
    if (!ResetResultCyc(ctx, order))
        return NULL;

    for (size_t i = 0; i < len; i++) {
        uint32_t exponent = exponents[i] % order;
        if (!AddChecked(
                ctx,
                ctx->result[exponent],
                coefficients[i],
                &ctx->result[exponent]))
            return NULL;
    }
    if (!ConvertToBase(ctx, order))
        return NULL;
    return Cyclotomic(ctx, order, 1);
}

gap_cyc * gap_cyc_root(gap_cyc_ctx * ctx, uint32_t order, uint32_t exponent)
{
    int64_t coefficient = 1;
    return gap_cyc_from_terms(ctx, order, 1, &exponent, &coefficient);
}

gap_cyc * gap_cyc_add(gap_cyc_ctx * ctx, const gap_cyc * lhs, const gap_cyc * rhs)
{
    if (ctx == NULL || lhs == NULL || rhs == NULL)
        return NULL;
    ctx->error = NULL;

    uint32_t ml, mr, n;
    if (!FindCommonField(ctx, lhs->order, rhs->order, &ml, &mr, &n))
        return NULL;
    memset(ctx->result, 0, (size_t)n * sizeof(int64_t));

    for (size_t i = 0; i < lhs->len; i++)
        ctx->result[lhs->exponents[i] * ml] = lhs->coefficients[i];

    for (size_t i = 0; i < rhs->len; i++) {
        uint32_t exponent = rhs->exponents[i] * mr;
        if (!AddChecked(
                ctx,
                ctx->result[exponent],
                rhs->coefficients[i],
                &ctx->result[exponent]))
            return NULL;
    }

    if (lhs->order % ml != 0 || rhs->order % mr != 0) {
        if (!ConvertToBase(ctx, n))
            return NULL;
    }
    return Cyclotomic(ctx, n, ml * mr);
}

gap_cyc * gap_cyc_mul(gap_cyc_ctx * ctx, const gap_cyc * lhs, const gap_cyc * rhs)
{
    if (ctx == NULL || lhs == NULL || rhs == NULL)
        return NULL;
    ctx->error = NULL;

    uint32_t ml, mr, n;
    if (!FindCommonField(ctx, lhs->order, rhs->order, &ml, &mr, &n))
        return NULL;
    memset(ctx->result, 0, (size_t)n * sizeof(int64_t));

    for (size_t i = 0; i < lhs->len; i++) {
        for (size_t j = 0; j < rhs->len; j++) {
            uint32_t exponent =
                ((uint64_t)lhs->exponents[i] * ml + (uint64_t)rhs->exponents[j] * mr) % n;
            int64_t product;
            if (!MulChecked(ctx, lhs->coefficients[i], rhs->coefficients[j], &product))
                return NULL;
            if (!AddChecked(ctx, ctx->result[exponent], product, &ctx->result[exponent]))
                return NULL;
        }
    }

    if (!ConvertToBase(ctx, n))
        return NULL;
    return Cyclotomic(ctx, n, ml * mr);
}

void gap_cyc_free(gap_cyc * cyc)
{
    if (cyc == NULL)
        return;
    free(cyc->coefficients);
    free(cyc->exponents);
    free(cyc);
}

uint32_t gap_cyc_order(const gap_cyc * cyc)
{
    return cyc->order;
}

size_t gap_cyc_len(const gap_cyc * cyc)
{
    return cyc->len;
}

uint32_t gap_cyc_exponent(const gap_cyc * cyc, size_t index)
{
    return cyc->exponents[index];
}

int64_t gap_cyc_coefficient(const gap_cyc * cyc, size_t index)
{
    return cyc->coefficients[index];
}

int gap_cyc_equal(const gap_cyc * lhs, const gap_cyc * rhs)
{
    if (lhs->order != rhs->order || lhs->len != rhs->len)
        return 0;
    for (size_t i = 0; i < lhs->len; i++) {
        if (lhs->exponents[i] != rhs->exponents[i]
            || lhs->coefficients[i] != rhs->coefficients[i])
            return 0;
    }
    return 1;
}
