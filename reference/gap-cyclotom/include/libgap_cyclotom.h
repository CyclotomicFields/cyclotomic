#ifndef LIBGAP_CYCLOTOM_H
#define LIBGAP_CYCLOTOM_H

#include <stddef.h>
#include <stdint.h>

int libgap_cyc_init(const char * gap_root);
const char * libgap_cyc_error(void);

int libgap_cyc_from_terms(
    size_t slot,
    uint32_t order,
    size_t len,
    const uint32_t * exponents,
    const int64_t * numerators,
    const int64_t * denominators);
int libgap_cyc_add(size_t output, size_t lhs, size_t rhs);
int libgap_cyc_mul(size_t output, size_t lhs, size_t rhs);
int libgap_cyc_quo(size_t output, size_t lhs, size_t rhs);
int libgap_cyc_release(size_t slot);
int libgap_cyc_equal(size_t lhs, size_t rhs, int * equal);

int libgap_cyc_order(size_t slot, uint32_t * order);
int libgap_cyc_coefficient(
    size_t slot,
    uint32_t exponent,
    int64_t * numerator,
    int64_t * denominator);

#endif
