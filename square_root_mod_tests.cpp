#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include "square_root_mod.h"
#include "factorize.h"
#include "pow_mod.h"

uint_fast32_t my_min_nonresidue(uint_fast32_t primes[], size_t primes_count, uint_fast32_t p) {
	assert(p > 2);
	typedef PowMod<uint_fast32_t,((uint_fast32_t)1 << 31), uint_fast64_t> pow_mod_type;
	for (size_t i=0; i < primes_count; ++i) {
		uint_fast32_t nr = primes[i];
		if (pow_mod_type::pow_mod(p, nr, (p-1)/2) == p-1) return nr;
	}
	assert(0);
	return 0;
}

void test_least_nonresidue() {
	typedef uint_fast32_t num_type;
	typedef SquareRootMod<num_type, ((num_type)1)<<31, uint_fast64_t> srm_type;
	typedef PrimesArray<num_type> primes_array_type;
	num_type primes[65536];
	size_t primes_count = primes_array_type::fill_primes(primes, 65536, UINT32_MAX);
	assert(primes_count == 65536);
	
	num_type p;
	for (size_t idx=1; idx<primes_count; ++idx) {
		p = primes[idx];
		num_type nr = srm_type::least_nonresidue(primes, primes_count, p);
		assert(nr != 0);
		num_type my_nr = my_min_nonresidue(primes, primes_count, p);
		assert(nr == my_nr);
	}
}

void tests_suite() {
	test_least_nonresidue();
}

int main() {
	tests_suite();
	return 0;
}

