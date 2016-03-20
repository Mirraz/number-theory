#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "primitive_roots.h"

#define main HIDE_main
#define tests_suite HIDE_tests_suite
#define test_max_primitive_root HIDE_test_max_primitive_root
#include "mul_group_mod_tests.cpp"
#undef main
#undef tests_suite
#undef test_max_primitive_root

void test_is_primitive_root() {
	// primorial(9) <= 2^32-1 < primorial(10)
	typedef PrimitiveRoots<uint_fast32_t, 9, ((uint_fast32_t)1)<<31, uint_fast64_t> prrs_type;
	
	// pi(2^16) = 6542
	prrs_type::num_type primes[6542];
	size_t primes_count = prrs_type::primes_array_type::fill_primes(
		primes,
		sizeof(primes) / sizeof(primes[0]),
		(prrs_type::num_type)UINT16_MAX + 1
	);
	assert(primes_count == sizeof(primes) / sizeof(primes[0]));
	prrs_type::primes_array_type primes_array(primes, primes_count);
	
	for (size_t idx=0; idx<sizeof(primes) / sizeof(primes[0]); ++idx) {
		prrs_type::num_type modulo = primes[idx];
		prrs_type primitive_roots(primes_array, modulo);
		prrs_type::num_type min_root = 0;
		for (prrs_type::num_type i=1; i<modulo; ++i) {
			bool is_primitive_root = primitive_roots.is_primitive_root(i);
			if (is_primitive_root) {min_root = i; break;}
		}
		prrs_type::num_type max_root = 0;
		for (prrs_type::num_type i=modulo-1; i>=1; --i) {
			bool is_primitive_root = primitive_roots.is_primitive_root(i);
			if (is_primitive_root) {max_root = i; break;}
		}
		assert(min_root != 0);
		assert(max_root != 0);
		assert((modulo-1) % 4 != 0 || modulo-max_root == min_root);
		printf("%" PRIuFAST64 "\t %" PRIuFAST64 "\n", modulo, min_root);
		printf("%" PRIuFAST64 "\t-%" PRIuFAST64 "\n", modulo, modulo-max_root);
		printf("\n");
	}
}

void test_max_primitive_root() {
	fill_myprimes();
	fill_max_roots();

	// primorial(9) <= 2^32-1 < primorial(10)
	typedef PrimitiveRoots<uint_fast32_t, 9, ((uint_fast32_t)1)<<31, uint_fast64_t> prrs_type;
	
	// pi(2^16) = 6542
	prrs_type::num_type primes[6542];
	size_t primes_count = prrs_type::primes_array_type::fill_primes(
		primes,
		sizeof(primes) / sizeof(primes[0]),
		(prrs_type::num_type)UINT16_MAX + 1
	);
	assert(primes_count == sizeof(primes) / sizeof(primes[0]));
	prrs_type::primes_array_type primes_array(primes, primes_count);
	
	for (size_t idx=0; idx<sizeof(primes) / sizeof(primes[0]); ++idx) {
		prrs_type::num_type modulo = primes[idx];
		prrs_type primitive_roots(primes_array, modulo);
		prrs_type::num_type max_root = 0;
		for (prrs_type::num_type i=modulo-1; i>=1; --i) {
			bool is_primitive_root = primitive_roots.is_primitive_root(i);
			if (is_primitive_root) {max_root = i; break;}
		}
		assert(max_root != 0);
		
		assert(max_root == myprimes[idx] - max_roots[idx]);
	}
}

void tests_suite() {
	//test_is_primitive_root();
	test_max_primitive_root();
}

int main() {
	tests_suite();
	return 0;
}
