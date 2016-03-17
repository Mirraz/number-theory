#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// NULL
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "factorize.h"

void test_round_sqrt() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	fzr_type::test_round_sqrt();
}

struct MyPow {
	uint_fast64_t prime;
	uint_fast8_t exp;
};

struct MyFactors {
	uint_fast8_t pow_count;
	MyPow pows[16];
};

MyFactors my_factorize(uint_fast64_t n) {
	assert(n > 0);
	MyFactors factors;
	factors.pow_count = 0;
	if (n == 1) return factors;
	uint_fast64_t n_sqrt = round(sqrt(n));
	for (uint_fast64_t d=2; d<=n_sqrt; d+=2) {
		if (n % d == 0) {
			uint_fast8_t exp = 0;
			do {
				n /= d;
				++exp;
			} while (n % d == 0);
			assert(factors.pow_count < 16);
			factors.pows[factors.pow_count].prime = d;
			factors.pows[factors.pow_count].exp = exp;
			++factors.pow_count;
			n_sqrt = round(sqrt(n));
		}
		if (d == 2) --d;
	}
	if (n > 1) {
		assert(factors.pow_count == 0 || factors.pows[factors.pow_count-1].prime != n);
		factors.pows[factors.pow_count].prime = n;
		factors.pows[factors.pow_count].exp = 1;
		++factors.pow_count;
	}
	return factors;
}

void test_factorize() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	
	MyFactors factors;
	fzr_type::factorize_cb_type cb = [&factors] (fzr_type::num_type prime, fzr_type::exp_type exp) -> bool {
		factors.pows[factors.pow_count].prime = prime;
		factors.pows[factors.pow_count].exp = exp;
		++factors.pow_count;
		return false;
	};
	fzr_type factorizer(NULL, 0, cb);
	
	for (fzr_type::num_type i=1; i<=1024*256+1; ++i) {
		factors.pow_count = 0;
		factorizer.factorize(i);
		MyFactors my_factors = my_factorize(i);
		assert(factors.pow_count == my_factors.pow_count);
		for (uint_fast8_t i=0; i<factors.pow_count; ++i) {
			assert(factors.pows[i].prime == my_factors.pows[i].prime);
			assert(factors.pows[i].exp == my_factors.pows[i].exp);
		}
	}
}

// sum of two squares theorem
// one of summands can be 0 i.e. if n is square itself then function returns true
typedef uint_fast64_t square_summands_num_type;
bool check_2_square_summands(square_summands_num_type n) {
	if (n <= 1) return true;
	typedef Factorizer<square_summands_num_type> fzr_type;
	bool result = true;
	fzr_type::factorize_cb_type cb = [&result] (fzr_type::num_type prime, fzr_type::exp_type exp) -> bool {
		if ((prime & 3) == 3 && (exp & 1)) {
			result = false;
			return true;
		} else {
			return false;
		}
	};
	fzr_type factorizer(NULL, 0, cb);
	factorizer.factorize(n);
	return result;
}

void test_check_2_square_summands() {
	for (square_summands_num_type i=0; i<1024*4+1; ++i) {
		printf("%" PRIuFAST64 "\t%s\n", i, check_2_square_summands(i) ? "true" : "false");
	}
}

void tests_suite() {
	test_round_sqrt();
	test_factorize();
	//test_check_2_square_summands();
}

int main() {
	tests_suite();
	return 0;
}

