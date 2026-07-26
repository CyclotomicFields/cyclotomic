#include "libgap_cyclotom.h"

#include <libgap-api.h>

#include <limits.h>

static Obj Roots;
static Obj FuncE;
static Obj FuncConductor;
static Obj FuncCoeffs;
static Obj FuncNumerator;
static Obj FuncDenominator;
static Obj FuncAlternatingGroup;
static Obj FuncSL;
static Obj FuncPSL;
static Obj FuncCharacterTable;
static Obj FuncIrr;
static Obj FuncValuesOfClassFunction;
static Obj FuncSizesConjugacyClasses;
static Obj FuncSize;
static Obj FuncScalarProduct;
enum { CharacterTableCount = 3 };
static Obj CharacterTables[CharacterTableCount];
static Obj Irreducibles[CharacterTableCount];
static const char * LastError;

static void MarkRootsAndCharacterTables(void)
{
    if (Roots)
        GAP_MarkBag(Roots);
    for (size_t i = 0; i < CharacterTableCount; i++) {
        if (CharacterTables[i])
            GAP_MarkBag(CharacterTables[i]);
        if (Irreducibles[i])
            GAP_MarkBag(Irreducibles[i]);
    }
}

static void RecordError(void)
{
    LastError = "libgap raised an error";
}

static Obj At(size_t slot)
{
    if (Roots == 0 || slot == 0 || slot > GAP_LenList(Roots))
        return 0;
    return GAP_ElmList(Roots, slot);
}

static int Store(size_t slot, Obj value)
{
    if (slot == 0 || value == 0) {
        LastError = "invalid libgap cyclotomic slot or value";
        return 0;
    }
    GAP_AssList(Roots, slot, value);
    return 1;
}

int libgap_cyc_init(const char * gap_root)
{
    if (Roots)
        return 1;
    if (gap_root == 0 || gap_root[0] == '\0') {
        LastError = "GAP root is empty";
        return 0;
    }

    char * argv[] = {
        (char *)"gap",
        (char *)"-l",
        (char *)gap_root,
        (char *)"-q",
        (char *)"-T",
        (char *)"-A",
        0,
    };
    GAP_Initialize(6, argv, MarkRootsAndCharacterTables, RecordError, 0);

    int ok = GAP_Enter();
    if (ok) {
        Roots = GAP_NewPlist(64);
        GAP_AssignGlobalVariable("LIBGAP_CYCLOTOM_RUST_ROOTS", Roots);
        FuncE = GAP_ValueGlobalVariable("E");
        FuncConductor = GAP_ValueGlobalVariable("Conductor");
        FuncCoeffs = GAP_ValueGlobalVariable("COEFFS_CYC");
        FuncNumerator = GAP_ValueGlobalVariable("NumeratorRat");
        FuncDenominator = GAP_ValueGlobalVariable("DenominatorRat");
        FuncAlternatingGroup = GAP_ValueGlobalVariable("AlternatingGroup");
        FuncSL = GAP_ValueGlobalVariable("SL");
        FuncPSL = GAP_ValueGlobalVariable("PSL");
        FuncCharacterTable = GAP_ValueGlobalVariable("CharacterTable");
        FuncIrr = GAP_ValueGlobalVariable("Irr");
        FuncValuesOfClassFunction =
            GAP_ValueGlobalVariable("ValuesOfClassFunction");
        FuncSizesConjugacyClasses =
            GAP_ValueGlobalVariable("SizesConjugacyClasses");
        FuncSize = GAP_ValueGlobalVariable("Size");
        FuncScalarProduct = GAP_ValueGlobalVariable("ScalarProduct");
        ok = Roots && FuncE && FuncConductor && FuncCoeffs
            && FuncNumerator && FuncDenominator && FuncAlternatingGroup
            && FuncSL && FuncPSL && FuncCharacterTable && FuncIrr
            && FuncValuesOfClassFunction && FuncSizesConjugacyClasses
            && FuncSize && FuncScalarProduct;
        if (!ok)
            LastError = "required GAP cyclotomic globals are unavailable";
    }
    GAP_Leave();
    return ok;
}

const char * libgap_cyc_error(void)
{
    return LastError ? LastError : "unknown libgap error";
}

static int EnsureCharacterTable(uint32_t table)
{
    if (table >= CharacterTableCount) {
        LastError = "unknown character table";
        return 0;
    }
    if (CharacterTables[table])
        return 1;

    Obj two = GAP_NewObjIntFromInt(2);
    Obj parameter = GAP_NewObjIntFromInt(table == 2 ? 11 : 5);
    Obj group;
    if (table == 0)
        group = GAP_CallFunc1Args(FuncAlternatingGroup, parameter);
    else if (table == 1)
        group = GAP_CallFunc2Args(FuncSL, two, parameter);
    else
        group = GAP_CallFunc2Args(FuncPSL, two, parameter);
    if (!group) {
        LastError = "failed to construct character table group";
        return 0;
    }
    CharacterTables[table] = GAP_CallFunc1Args(FuncCharacterTable, group);
    Irreducibles[table] = GAP_CallFunc1Args(FuncIrr, CharacterTables[table]);
    if (!CharacterTables[table] || !Irreducibles[table]) {
        LastError = "failed to construct character table";
        return 0;
    }
    return 1;
}

int libgap_cyc_from_terms(
    size_t slot,
    uint32_t order,
    size_t len,
    const uint32_t * exponents,
    const int64_t * numerators,
    const int64_t * denominators)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok) {
        Obj order_obj = GAP_NewObjIntFromInt(order);
        Obj root = GAP_CallFunc1Args(FuncE, order_obj);
        Obj value = GAP_NewObjIntFromInt(0);
        for (size_t i = 0; i < len && ok; i++) {
            if (denominators[i] <= 0) {
                LastError = "rational denominator must be positive";
                ok = 0;
                break;
            }
            Obj numerator = GAP_NewObjIntFromInt(numerators[i]);
            Obj denominator = GAP_NewObjIntFromInt(denominators[i]);
            Obj coefficient = GAP_QUO(numerator, denominator);
            Obj exponent = GAP_NewObjIntFromInt(exponents[i] % order);
            Obj power = GAP_POW(root, exponent);
            value = GAP_SUM(value, GAP_PROD(coefficient, power));
        }
        if (ok)
            ok = Store(slot, value);
    }
    GAP_Leave();
    return ok;
}

static int Binary(size_t output, size_t lhs, size_t rhs, char operation)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok) {
        Obj left = At(lhs);
        Obj right = At(rhs);
        if (!left || !right) {
            LastError = "invalid libgap cyclotomic operand slot";
            ok = 0;
        }
        else {
            Obj result = operation == '+' ? GAP_SUM(left, right)
                : operation == '*' ? GAP_PROD(left, right)
                                   : GAP_QUO(left, right);
            ok = Store(output, result);
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_cyc_add(size_t output, size_t lhs, size_t rhs)
{
    return Binary(output, lhs, rhs, '+');
}

int libgap_cyc_mul(size_t output, size_t lhs, size_t rhs)
{
    return Binary(output, lhs, rhs, '*');
}

int libgap_cyc_quo(size_t output, size_t lhs, size_t rhs)
{
    return Binary(output, lhs, rhs, '/');
}

int libgap_cyc_release(size_t slot)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok && Roots && slot != 0)
        GAP_AssList(Roots, slot, 0);
    GAP_Leave();
    return ok;
}

int libgap_cyc_equal(size_t lhs, size_t rhs, int * equal)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok) {
        Obj left = At(lhs);
        Obj right = At(rhs);
        if (!left || !right || !equal) {
            LastError = "invalid libgap equality operand";
            ok = 0;
        }
        else {
            *equal = GAP_EQ(left, right) != 0;
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_cyc_order(size_t slot, uint32_t * order)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok) {
        Obj value = At(slot);
        Obj conductor = value ? GAP_CallFunc1Args(FuncConductor, value) : 0;
        if (!conductor || !order) {
            LastError = "invalid libgap conductor request";
            ok = 0;
        }
        else {
            Int result = GAP_ValueInt(conductor);
            if (result < 0 || (UInt)result > UINT_MAX) {
                LastError = "cyclotomic conductor does not fit in uint32_t";
                ok = 0;
            }
            else {
                *order = (uint32_t)result;
            }
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_cyc_coefficient(
    size_t slot,
    uint32_t exponent,
    int64_t * numerator,
    int64_t * denominator)
{
    LastError = 0;
    int ok = GAP_Enter();
    if (ok) {
        Obj value = At(slot);
        Obj coefficients = value ? GAP_CallFunc1Args(FuncCoeffs, value) : 0;
        if (!coefficients || !numerator || !denominator) {
            LastError = "invalid libgap coefficient request";
            ok = 0;
        }
        else {
            UInt len = GAP_LenList(coefficients);
            Obj coefficient = exponent < len
                ? GAP_ElmList(coefficients, (UInt)exponent + 1)
                : GAP_NewObjIntFromInt(0);
            Obj num = GAP_CallFunc1Args(FuncNumerator, coefficient);
            Obj den = GAP_CallFunc1Args(FuncDenominator, coefficient);
            Int num_value = GAP_ValueInt(num);
            Int den_value = GAP_ValueInt(den);
            *numerator = (int64_t)num_value;
            *denominator = (int64_t)den_value;
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_character_table_dimensions(
    uint32_t table,
    size_t * rows,
    size_t * columns)
{
    LastError = 0;
    if (!rows || !columns) {
        LastError = "character table dimensions need output pointers";
        return 0;
    }
    int ok = GAP_Enter();
    if (ok) {
        ok = EnsureCharacterTable(table);
        if (ok) {
            *rows = GAP_LenList(Irreducibles[table]);
            Obj first = GAP_ElmList(Irreducibles[table], 1);
            Obj values = GAP_CallFunc1Args(FuncValuesOfClassFunction, first);
            *columns = GAP_LenList(values);
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_character_table(
    uint32_t table,
    const size_t * slots,
    size_t slots_len,
    int64_t * class_sizes,
    int64_t * group_order)
{
    LastError = 0;
    if (!slots || !class_sizes || !group_order) {
        LastError = "character table needs output buffers";
        return 0;
    }
    int ok = GAP_Enter();
    if (ok)
        ok = EnsureCharacterTable(table);
    if (ok) {
        size_t rows = GAP_LenList(Irreducibles[table]);
        Obj first = GAP_ElmList(Irreducibles[table], 1);
        Obj first_values =
            GAP_CallFunc1Args(FuncValuesOfClassFunction, first);
        size_t columns = GAP_LenList(first_values);
        if (slots_len != rows * columns) {
            LastError = "wrong number of character table slots";
            ok = 0;
        }
        Obj sizes =
            GAP_CallFunc1Args(FuncSizesConjugacyClasses, CharacterTables[table]);
        Obj size = GAP_CallFunc1Args(FuncSize, CharacterTables[table]);
        if (ok && !GAP_IsSmallInt(size)) {
            LastError = "character table group order does not fit an integer";
            ok = 0;
        }
        if (ok)
            *group_order = GAP_ValueInt(size);
        for (size_t column = 0; column < columns && ok; column++) {
            Obj class_size = GAP_ElmList(sizes, column + 1);
            if (!GAP_IsSmallInt(class_size)) {
                LastError = "conjugacy class size does not fit an integer";
                ok = 0;
            }
            else {
                class_sizes[column] = GAP_ValueInt(class_size);
            }
        }
        for (size_t row = 0; row < rows && ok; row++) {
            Obj character = GAP_ElmList(Irreducibles[table], row + 1);
            Obj values =
                GAP_CallFunc1Args(FuncValuesOfClassFunction, character);
            for (size_t column = 0; column < columns && ok; column++) {
                ok = Store(
                    slots[row * columns + column],
                    GAP_ElmList(values, column + 1));
            }
        }
    }
    GAP_Leave();
    return ok;
}

int libgap_character_tensor_decomposition(
    uint32_t table,
    size_t lhs,
    size_t rhs,
    int64_t * multiplicities,
    size_t multiplicities_len)
{
    LastError = 0;
    if (!multiplicities) {
        LastError = "tensor decomposition needs an output buffer";
        return 0;
    }
    int ok = GAP_Enter();
    if (ok)
        ok = EnsureCharacterTable(table);
    if (ok) {
        size_t rows = GAP_LenList(Irreducibles[table]);
        if (lhs >= rows || rhs >= rows || multiplicities_len != rows) {
            LastError = "invalid character row or multiplicity length";
            ok = 0;
        }
        Obj left = ok ? GAP_ElmList(Irreducibles[table], lhs + 1) : 0;
        Obj right = ok ? GAP_ElmList(Irreducibles[table], rhs + 1) : 0;
        Obj product = ok ? GAP_PROD(left, right) : 0;
        for (size_t row = 0; row < rows && ok; row++) {
            Obj character = GAP_ElmList(Irreducibles[table], row + 1);
            Obj scalar =
                GAP_CallFunc2Args(FuncScalarProduct, product, character);
            if (!GAP_IsSmallInt(scalar)) {
                LastError = "character multiplicity does not fit an integer";
                ok = 0;
            }
            else {
                multiplicities[row] = GAP_ValueInt(scalar);
            }
        }
    }
    GAP_Leave();
    return ok;
}
