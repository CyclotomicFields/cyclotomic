// These are C libraries I'm using from C++ because C is such trash
#include <flint/arith.h>
#include <antic/nf.h>
#include <antic/nf_elem.h>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <string>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "need n and num_tests!\n";
        return -1;
    }

    size_t n = std::stoi(std::string(argv[1]));
    size_t num_tests = std::stoi(std::string(argv[2]));

    std::cerr << "n = " << n  << "\nnum_tests = " << num_tests << std::endl;

    const size_t num_per_test = 6;
    const size_t max_num_terms = 5;

    fmpz_poly_t cyc_z;
    fmpz_poly_init(cyc_z);
    fmpz_poly_cyclotomic(cyc_z, n);

    fmpq_poly_t cyclotomic_polynomial_n;
    fmpq_poly_init(cyclotomic_polynomial_n);
    fmpq_poly_set_fmpz_poly(cyclotomic_polynomial_n, cyc_z);

    nf_t cyclotomic_field_n;
    nf_init(cyclotomic_field_n, cyclotomic_polynomial_n);

    // some random generator for the test numbers
    flint_rand_t state;
    flint_randinit(state);

    /*
    nf_elem_t** nums = new nf_elem_t*[num_tests];
    for (size_t i = 0; i < num_tests; ++i) {
        nums[i] = new nf_elem_t[num_per_test];
    }
    */
    auto nums = new nf_elem_t[num_tests][num_per_test];

    std::default_random_engine generator;
    std::uniform_int_distribution<int> num_terms_dist(1,max_num_terms-1);
    std::uniform_int_distribution<int> rational_dist(1,10);
    std::uniform_int_distribution<int> exp_dist(1,n-1);

    std::cerr << "generating test data" << std::endl;
    for (size_t i = 0; i < num_tests; ++i) {
        for (size_t j = 0; j < num_per_test; ++j) {
            nf_elem_struct* num = nums[i][j];
            nf_elem_init(num, cyclotomic_field_n);
            int num_terms = num_terms_dist(generator);

            for (int k = 0; k < num_terms; ++k) {
                nf_elem_t term;
                nf_elem_init(term, cyclotomic_field_n);

                fmpq_poly_t pol;
                fmpq_poly_init(pol);

                int exp = exp_dist(generator);
                int numerator = rational_dist(generator);
                int denominator = rational_dist(generator);

                fmpq_t coeff;
                fmpq_init(coeff);
                fmpq_set_si(coeff, numerator, denominator);

                fmpq_poly_set_coeff_fmpq(pol, exp, coeff);

                nf_elem_set_fmpq_poly(term, pol, cyclotomic_field_n);

                nf_elem_t sum;
                nf_elem_init(sum, cyclotomic_field_n);
                nf_elem_add(sum, num, term, cyclotomic_field_n);

                nf_elem_set(num, sum, cyclotomic_field_n);
            }
        }
    }

    /*
    std::cout << "start print\n";
    for (size_t i = 0; i < num_tests; ++i) {
        for (size_t j = 0; j < num_per_test; ++j) {
            std::cout << i << ": " << j << ": ";
            // to check the numbers randomly generated are "reasonable"
            nf_elem_print_pretty(nums[i][j], cyclotomic_field_n, "e");
            std::cout << std::endl;
        }
    }
    std::cout << "end print" << std::endl;
    */

    std::cout << "starting benchmark" << std::endl;

    const auto start = std::chrono::steady_clock::now();

    for (size_t i = 0; i < num_tests; ++i) {
        nf_elem_t* chunk = nums[i];
        nf_elem_t prod1;
        nf_elem_init(prod1, cyclotomic_field_n);
        fmpq_poly_fit_length(prod1->elem, 2*n);
        nf_elem_mul(prod1, chunk[0], chunk[1], cyclotomic_field_n);

        nf_elem_t prod2;
        nf_elem_init(prod2, cyclotomic_field_n);
        fmpq_poly_fit_length(prod2->elem, 2*n);
        nf_elem_mul(prod2, chunk[2], chunk[3], cyclotomic_field_n);

        nf_elem_t prod3;
        nf_elem_init(prod3, cyclotomic_field_n);
        fmpq_poly_fit_length(prod3->elem, 2*n);
        nf_elem_mul(prod3, chunk[4], chunk[5], cyclotomic_field_n);

        nf_elem_t sum1;
        nf_elem_init(sum1, cyclotomic_field_n);
        fmpq_poly_fit_length(sum1->elem, 2*n);
        nf_elem_add(sum1, prod1, prod2, cyclotomic_field_n);

        nf_elem_t sum2;
        nf_elem_init(sum2, cyclotomic_field_n);
        fmpq_poly_fit_length(sum2->elem, 2*n);
        nf_elem_add(sum2, sum1, prod3, cyclotomic_field_n);
    }

    const auto end = std::chrono::steady_clock::now();

    std::cerr << "time elapsed (ms):" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "\n";

    delete[] nums;

    return 0;
}
