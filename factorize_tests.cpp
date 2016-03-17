#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// NULL
#include <stdint.h>
#include <inttypes.h>
#include "factorize.h"

//void test_round_sqrt() {
//	typedef Factorization<uint_fast64_t> fzn_type;
//	fzn_type::test_round_sqrt();
//}

void test_factorize() {
	typedef Factorization<uint_fast64_t> fzn_type;
	
	bool cb_was_called;
	bool is_prime;
	fzn_type::factorize_cb_type cb = [&cb_was_called, &is_prime] (fzn_type::num_type prime, fzn_type::exp_type exp) -> bool {
		if (cb_was_called) is_prime = false;
		else if (exp > 1) is_prime = false;
		cb_was_called = true;
		printf("%" PRIuFAST64 "^%" PRIuFAST8 " ", prime, exp);
		return false;
	};
	fzn_type::Factorizer factorizer(NULL, 0, cb);
	
	for (fzn_type::num_type i=1; i<=1024*4+1; ++i) {
	//for (fzn_type::num_type i=UINT64_MAX; i>=UINT64_MAX-64; --i) {
		printf("%" PRIuFAST64 " = ", i);
		cb_was_called = false;
		is_prime = true;
		factorizer.factorize(i);
		if (is_prime) printf("\t\tPRIME");
		printf("\n");
	}
}

// sum of two squares theorem
// one of summands can be 0 i.e. if n is square itself then function returns true
typedef uint_fast64_t square_summands_num_type;
bool check_2_square_summands(square_summands_num_type n) {
	if (n <= 1) return true;
	typedef Factorization<square_summands_num_type> fzn_type;
	bool result = true;
	fzn_type::factorize_cb_type cb = [&result] (fzn_type::num_type prime, fzn_type::exp_type exp) -> bool {
		if ((prime & 3) == 3 && (exp & 1)) {
			result = false;
			return true;
		} else {
			return false;
		}
	};
	fzn_type::Factorizer factorizer(NULL, 0, cb);
	factorizer.factorize(n);
	return result;
}

void test_check_2_square_summands() {
	for (square_summands_num_type i=0; i<1024*4+1; ++i) {
		printf("%" PRIuFAST64 "\t%s\n", i, check_2_square_summands(i) ? "true" : "false");
	}
}

void tests_suite() {
	//test_round_sqrt();
	//test_factorize();
	test_check_2_square_summands();
}

int main() {
	tests_suite();
	return 0;
}

