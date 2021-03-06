#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include "square_root_mod.h"
#include "factorize.h"
#include "mul_mod.h"

uint_fast32_t my_min_nonresidue(uint_fast32_t primes[], size_t primes_count, uint_fast32_t p) {
	assert(p > 2);
	typedef MulMod<uint_fast32_t,((uint_fast32_t)1 << 31), uint_fast64_t> mul_mod_type;
	for (size_t i=0; i < primes_count; ++i) {
		uint_fast32_t nr = primes[i];
		if (mul_mod_type::pow_mod(p, nr, (p-1)/2) == p-1) return nr;
	}
	assert(0);
	return 0;
}

void test_least_nonresidue() {
	typedef uint_fast32_t num_type;
	typedef SquareRootMod<num_type, 32, uint_fast64_t> srm_type;
	typedef PrimesArray<num_type> primes_array_type;
	num_type primes[65536];
	size_t primes_count = primes_array_type::fill_primes(primes, 65536, UINT32_MAX);
	assert(primes_count == 65536);
	
	for (size_t idx=1; idx<primes_count; ++idx) {
		num_type p = primes[idx];
		num_type nr = srm_type::least_nonresidue(primes, primes_count, p);
		assert(nr != 0);
		num_type my_nr = my_min_nonresidue(primes, primes_count, p);
		assert(nr == my_nr);
	}
}

void test_square_root_mod_algo_01() {
	typedef uint_fast32_t num_type;
	typedef SquareRootMod<num_type, 32, uint_fast64_t> srm_type;
	typedef PrimesArray<num_type> primes_array_type;
	num_type primes[1024];
	size_t primes_count = primes_array_type::fill_primes(primes, 1024, UINT32_MAX);
	assert(primes_count == 1024);
	
	for (size_t idx=1; idx<primes_count; ++idx) {
		num_type p = primes[idx];
		num_type nr = srm_type::least_nonresidue(primes, primes_count, p);
		assert(nr != 0);
		num_type r_count = 0, nr_count = 0;
		for (num_type a=1; a<p; ++a) {
			num_type r = srm_type::square_root_mod_algo_01(p, nr, a);
			int legendre_symbol = srm_type::legendre_symbol(p, a);
			assert((legendre_symbol != 1) == (r == 0));
			if (r == 0) {
				++nr_count;
			} else {
				assert(srm_type::mul_mod_type::square_mod(p, r) == a);
				++r_count;
			}
		}
		assert(r_count == nr_count);
		assert(r_count == (p-1)/2);
	}
}

void test_tonelli_shanks_algo() {
	typedef uint_fast32_t num_type;
	typedef SquareRootMod<num_type, 32, uint_fast64_t> srm_type;
	typedef PrimesArray<num_type> primes_array_type;
	num_type primes[1024];
	size_t primes_count = primes_array_type::fill_primes(primes, 1024, UINT32_MAX);
	assert(primes_count == 1024);
	
	for (size_t idx=1; idx<primes_count; ++idx) {
		num_type p = primes[idx];
		srm_type square_root_mod(p, primes, primes_count);
		num_type r_count = 0, nr_count = 0;
		for (num_type a=1; a<p; ++a) {
			if (square_root_mod.legendre_symbol(a) != 1) {++nr_count; continue;}
			num_type r = square_root_mod.tonelli_shanks_algo(a);
			assert(srm_type::mul_mod_type::square_mod(p, r) == a);
			++r_count;
		}
		assert(r_count == nr_count);
		assert(r_count == (p-1)/2);
	}
}

void test_tonelli_shanks_algo_02_rand() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned int seed = (unsigned int)tv.tv_sec * 1000000 + tv.tv_usec;
	fprintf(stderr, "seed = %u\n", seed);
	srand(seed);

	typedef uint_fast32_t num_type;
	typedef SquareRootMod<num_type, 32, uint_fast64_t> srm_type;
	typedef PrimesArray<num_type> primes_array_type;
	num_type primes[1024];
	size_t primes_count = primes_array_type::fill_primes(primes, 1024, UINT32_MAX);
	assert(primes_count == 1024);
	PrimeChecker<num_type> prime_checker(primes_array_type(primes, primes_count));
	
	for (num_type p = UINT32_MAX-1024*4;; p+=2) {
		if (prime_checker.is_prime(p)) {
			srm_type square_root_mod(p, primes, primes_count);
			for (uint_fast16_t i=0; i<1024*4; ++i) {
				num_type a = rand() % p;
				if (square_root_mod.legendre_symbol(a) != 1) continue;
				num_type r = square_root_mod.tonelli_shanks_algo(a);
				assert(srm_type::mul_mod_type::square_mod(p, r) == a);
			}
		}
		if (p == UINT32_MAX) break;
	}
}

void tests_suite() {
	//test_least_nonresidue();
	//test_square_root_mod_algo_01();
	test_tonelli_shanks_algo();
	test_tonelli_shanks_algo_02_rand();
}

int main() {
	tests_suite();
	return 0;
}

