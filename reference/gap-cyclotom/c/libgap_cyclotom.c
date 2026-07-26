#include "libgap_cyclotom.h"

#include <libgap-api.h>

#include <limits.h>

static Obj Roots;
static Obj FuncE;
static Obj FuncConductor;
static Obj FuncCoeffs;
static Obj FuncNumerator;
static Obj FuncDenominator;
static const char * LastError;

static void MarkRoots(void)
{
    if (Roots)
        GAP_MarkBag(Roots);
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
    GAP_Initialize(6, argv, MarkRoots, RecordError, 0);

    int ok = GAP_Enter();
    if (ok) {
        Roots = GAP_NewPlist(64);
        GAP_AssignGlobalVariable("LIBGAP_CYCLOTOM_RUST_ROOTS", Roots);
        FuncE = GAP_ValueGlobalVariable("E");
        FuncConductor = GAP_ValueGlobalVariable("Conductor");
        FuncCoeffs = GAP_ValueGlobalVariable("COEFFS_CYC");
        FuncNumerator = GAP_ValueGlobalVariable("NumeratorRat");
        FuncDenominator = GAP_ValueGlobalVariable("DenominatorRat");
        ok = Roots && FuncE && FuncConductor && FuncCoeffs
            && FuncNumerator && FuncDenominator;
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
