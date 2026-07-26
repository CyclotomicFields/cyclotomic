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

int libgap_character_table_dimensions(
    uint32_t table,
    size_t * rows,
    size_t * columns);
int libgap_character_table(
    uint32_t table,
    const size_t * slots,
    size_t slots_len,
    int64_t * class_sizes,
    int64_t * group_order);
int libgap_character_tensor_decomposition(
    uint32_t table,
    size_t lhs,
    size_t rhs,
    int64_t * multiplicities,
    size_t multiplicities_len);

#endif
