#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// NULL
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

typedef uint_fast64_t num_type;
typedef uint_fast8_t exp_type;
#define round_sqrt(n) ( (num_type)round( sqrt( (double)n ) ) )
#define NUM_MAX UINT64_MAX
#define SQRT_MAX ((num_type)UINT32_MAX)

void test_round_sqrt() {
	const num_type max_sqrt = SQRT_MAX + (num_type)1;
	const num_type max_delta = 1024;
	for (num_type i=max_sqrt; i>=max_sqrt-1024; --i) {
		//printf("%llu\n", i);
		for (num_type j=0; j<=max_delta; ++j) {
			num_type s1, s2;
			if (i == max_sqrt) {
				if (j == 0) s1 = NUM_MAX;
				else s1 = NUM_MAX - j + 1;
				s2 = NUM_MAX;
			} else {
				num_type s = i * i;
				s1 = s - j;
				s2 = s + j;
			}
			//printf("%llu %llu\n", s1, s2);
			num_type r1 = round_sqrt(s1), r2 = round_sqrt(s2);
			assert(r1 == i && r2 == i);
		}
		//printf("\n");
	}
}

typedef bool (*factorize_cb_type)(num_type prime, exp_type exp, bool last_call, void *cb_data);

// if cb returns true, factorize interrupts
// last_call is true if at the end n != 1
void factorize(num_type n, factorize_cb_type cb, void *cb_data) {
	assert(n > 0);
	num_type n_sqrt = round_sqrt(n);
	for (num_type p=2; p<=n_sqrt; p+=2) {
		if (!(n % p)) {
			exp_type exp = 0;
			do {
				n /= p;
				++exp;
			} while (!(n % p));
			if (cb(p, exp, false, cb_data)) break;
			n_sqrt = round_sqrt(n);
		}
		if (p == 2) --p;
	}
	if (n != 1) cb(n, 1, true, cb_data);
}

bool test_factorize_cb(num_type prime, exp_type exp, bool last_call, void *cb_data) {
	if (!last_call) *((bool *)cb_data) = false;
	printf("%"PRIuFAST64"^%"PRIuFAST8" ", prime, exp);
	return false;
}

void test_factorize() {
	/*
	//num_type i = 18446744030759878681llu;
	num_type i = 18446744065119616769llu;
	printf("%"PRIuFAST64" = ", i);
	bool is_prime = true;
	factorize(i, test_factorize_cb, &is_prime);
	if (is_prime) printf("\t\tPRIME");
	printf("\n");
	*/
	for (num_type i=1; i<=1024*4+1; ++i) {
		printf("%"PRIuFAST64" = ", i);
		bool is_prime = true;
		factorize(i, test_factorize_cb, &is_prime);
		if (is_prime) printf("\t\tPRIME");
		printf("\n");
	}
	/*
	for (num_type i=UINT64_MAX; i>=UINT64_MAX-64; --i) {
		printf("%"PRIuFAST64" = ", i);
		bool is_prime = true;
		factorize(i, test_factorize_cb, &is_prime);
		if (is_prime) printf("\t\tPRIME");
		printf("\n");
	}
	*/
}

bool check_2_square_summands_cb(num_type prime, exp_type exp, bool last_call, void *cb_data) {
	(void)last_call;
	if ((prime & 3) == 3 && (exp & 1)) {
		*((bool *)cb_data) = false;
		return true;
	} else {
		return false;
	}
}

// sum of two squares theorem
// one of summands can be 0 i.e. if n is square itself then function returns true
bool check_2_square_summands(num_type n) {
	if (n <= 1) return true;
	bool result = true;
	factorize(n, check_2_square_summands_cb, &result);
	return result;
}

void test_check_2_square_summands() {
	for (num_type i=0; i<1024*4+1; ++i) {
		printf("%"PRIuFAST64"\t%s\n", i, check_2_square_summands(i) ? "true" : "false");
	}
}

void tests_suite() {
	//test_round_sqrt();
	//test_factorize();
	//test_check_2_square_summands();
}

int main() {
	tests_suite();
	return 0;
}

