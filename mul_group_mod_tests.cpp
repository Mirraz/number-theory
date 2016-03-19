#include <assert.h>
#include <stdio.h>
#include "mul_group_mod.h"

uint_fast64_t gcd(uint_fast64_t a, uint_fast64_t b) {
	if (a == 0) return b;
	return gcd(b % a, a);
}

void test_order() {
	typedef MulGroupMod<uint_fast32_t, 9, ((uint_fast32_t)1)<<31, uint_fast64_t> mgm_type;
	mgm_type::primes_array_type primes_array;
	mgm_type a(primes_array);
	mgm_type::num_type i, j, k;
	for (i=2; i<1024; ++i) {
		a.assign(i);
		for (j=0; j<i; ++j) {
			if (gcd(i,j) != 1) continue;
			printf("mod = %u element = %u ", (unsigned int)i, (unsigned int)j);
			k = a.element_order(j);
			printf("order = %u\n", (unsigned int)k);
		}
	}
}

void tests_suite() {
	test_order();
}

int main() {
	tests_suite();
	return 0;
}

