#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "canonic_factors.h"

void test_constructor_and_value() {
	typedef CanonicFactors<uint_fast64_t, 15> canonic_factors_type;
	canonic_factors_type::primes_array_type primes_array;
	canonic_factors_type a(primes_array);
	canonic_factors_type::num_type i, j;
	for (i=1; i<=1024+1; ++i) {
		a.assign(i);
		//printf("%u = ", i); a.dump(); printf("\n");
		j = a.value();
		assert(i == j);
	}
	for (i=UINT32_MAX; i>=UINT32_MAX-1024; --i) {
		a.assign(i);
		j = a.value();
		assert(i == j);
	}
}

void test_mul() {
	typedef CanonicFactors<uint_fast64_t, 15> canonic_factors_type;
	canonic_factors_type::primes_array_type primes_array;
	canonic_factors_type a(primes_array), b(primes_array), c(primes_array), d(primes_array);
	canonic_factors_type::num_type i, j, k, l;
	for (i=1; i<=1024+1; ++i) {
		for (j=1; j<=1024+1; ++j) {
			a.assign(i);
			b.assign(j);
			//printf("a = "); a.dump(); printf("\n");
			//printf("b = "); b.dump(); printf("\n");
			c = a * b;
			d = a * j;
			k = c.value();
			l = d.value();
			assert(k == i * j);
			assert(l == i * j);
		}
	}
}

uint_fast64_t gcd(uint_fast64_t a, uint_fast64_t b) {
	if (a == 0) return b;
	return gcd(b % a, a);
}

uint_fast64_t phi(uint_fast64_t n) {
	uint_fast64_t result = 1;
	for (uint_fast64_t i=2; i<n; ++i) {
		if (gcd(i, n) == 1) ++result;
	}
	return result;
}

void test_eulers_phi() {
	typedef CanonicFactors<uint_fast64_t, 15> canonic_factors_type;
	canonic_factors_type::primes_array_type primes_array;
	canonic_factors_type a(primes_array), b(primes_array);
	unsigned int i, j;
	for (i=1; i<=1024*2+1; ++i) {
		a.assign(i);
		b = canonic_factors_type::eulers_phi(a);
		j = b.value();
		assert(j == phi(i));
	}
}

// checked by wolframalpha.com
uint_fast64_t phi_div_carmichael[] = {
	0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 4, 1, 1, 1, 2, 1, 2, 1,
	2, 2, 1, 2, 2, 1, 1, 2, 4, 1, 2, 1, 2, 2, 1, 1, 4, 1, 1, 2, 2, 1, 1, 2, 4, 2, 1, 1, 4, 1, 1, 6,
	2, 4, 2, 1, 2, 2, 2, 1, 4, 1, 1, 2, 2, 2, 2, 1, 8, 1, 1, 1, 4, 4, 1, 2, 4, 1, 2, 6, 2, 2, 1, 2,
	4, 1, 1, 2, 2, 1, 2, 1, 4, 4, 1, 1, 2, 1, 2, 2, 4, 1, 2, 2, 2, 6, 1, 2, 8, 1, 1, 2, 2, 1, 6, 1,
	2, 2, 4, 1, 4, 6, 1, 2, 4, 1, 2, 1, 4, 2, 1, 2, 4, 4, 1, 2, 2, 1, 2, 1, 4, 2, 2, 2, 4, 1, 1, 2,
	8, 2, 1, 1, 2, 4, 1, 1, 8, 1, 4, 6, 2, 1, 2, 2, 4, 2, 1, 1, 4, 1, 6, 2, 4, 4, 2, 2, 2, 6, 2, 1,
	4, 1, 1, 8, 2, 1, 2, 1, 4, 2, 1, 2, 4, 4, 1, 2, 8, 2, 4, 1, 2, 2, 1, 2, 4, 6, 1, 2, 4, 4, 2, 1,
	4, 2, 1, 1, 4, 1, 2, 4, 4, 1, 6, 2, 2, 2, 2, 1, 16, 1, 1, 1, 2, 2, 2, 6, 4, 2, 1, 1, 12, 2, 1, 8,
	2, 1, 2, 6, 8, 2, 1, 1, 8, 4, 6, 2, 2, 1, 2, 1, 8, 12, 1, 10, 4, 1, 1, 6, 8, 1, 2, 1, 2, 4, 2, 2,
	4, 1, 4, 2, 2, 1, 2, 2, 4, 2, 1, 2, 4, 6, 1, 2, 4, 4, 2, 1, 4, 2, 2, 1, 8, 1, 1, 12, 2, 1, 2, 2,
	8, 2, 2, 2, 2, 4, 1, 2, 4, 2, 4, 1, 2, 6, 1, 2, 8, 1, 1, 2, 8, 10, 6, 1, 4, 4, 1, 1, 4, 1, 2, 6,
	4, 1, 2, 2, 2, 4, 1, 1, 8, 1, 1, 2, 12, 4, 2, 1, 4, 2, 4, 2, 4, 1, 2, 2, 4, 4, 6, 1, 4, 2, 1, 1,
	4, 4, 1, 6, 2, 1, 8, 2, 4, 2, 1, 2, 4, 1, 1, 12, 8, 1, 2, 6, 2, 2, 2, 2, 8, 1, 4, 2, 2, 2, 2, 2,
	8, 2, 2, 1, 8, 1, 1, 2, 4, 4, 2, 6, 2, 4, 2, 1, 4, 1, 6, 8, 2, 2, 2, 1, 8, 6, 4, 1, 4, 4, 1, 2,
	4, 1, 2, 10, 2, 2, 1, 24, 8, 1, 1, 2, 4, 1, 4, 1, 8, 4, 1, 1, 12, 6, 2, 2, 4, 2, 2, 2, 4, 2, 1, 1,
	16, 12, 1, 4, 2, 4, 1, 1, 4, 2, 2, 1, 4, 4, 6, 4, 4, 2, 2, 1, 2, 2, 1, 1, 24, 4, 2, 2, 2, 1, 8, 6,
	2, 18, 1, 2, 4, 2, 6, 2, 16, 1, 2, 1, 2, 4, 1, 2, 8, 1, 4, 2, 12, 4, 2, 2, 4, 2, 1, 2, 4, 1, 1, 2,
	16, 4, 12, 1, 2, 6, 10, 2, 8, 6, 1, 8, 2, 1, 6, 6, 16, 4, 1, 1, 4, 4, 1, 6, 4, 1, 4, 1, 4, 2, 2, 2,
	4, 1, 1, 2, 8, 2, 2, 2, 4, 24, 1, 1, 4, 6, 2, 2, 8, 1, 2, 8, 2, 2, 2, 1, 8, 1, 6, 6, 2, 2, 2, 1,
	4, 4, 4, 2, 4, 1, 1, 8, 8, 1, 2, 1, 4, 2, 1, 2, 16, 1, 1, 4, 2, 4, 12, 1, 4, 2, 1, 2, 4, 6, 2, 2,
	8, 1, 2, 1, 4, 4, 2, 1, 4, 2, 4, 12, 2, 1, 2, 2, 8, 6, 2, 1, 8, 1, 1, 8, 4, 12, 6, 2, 2, 2, 2, 10,
	8, 1, 1, 2, 2, 1, 2, 6, 16, 2, 10, 1, 12, 4, 1, 2, 4, 4, 4, 1, 2, 12, 1, 2, 8, 8, 1, 2, 4, 1, 6, 18,
	4, 4, 1, 2, 4, 1, 2, 6, 4, 2, 4, 8, 2, 2, 1, 1, 16, 6, 1, 2, 2, 4, 2, 1, 24, 1, 4, 2, 4, 1, 1, 4,
	4, 2, 2, 1, 8, 12, 2, 1, 8, 4, 1, 2, 4, 2, 2, 1, 4, 2, 4, 2, 12, 1, 1, 4, 8, 1, 2, 6, 2, 8, 1, 2,
	4, 1, 4, 2, 2, 1, 6, 10, 4, 12, 1, 2, 16, 10, 2, 2, 4, 4, 2, 1, 2, 2, 2, 2, 8, 12, 1, 8, 2, 1, 12, 2,
	8, 2, 1, 2, 4, 4, 6, 2, 4, 1, 2, 1, 4, 2, 2, 2, 16, 6, 1, 36, 8, 1, 2, 1, 4, 20, 2, 1, 4, 1, 2, 2,
	8, 2, 2, 2, 4, 6, 1, 1, 16, 1, 1, 2, 2, 4, 2, 2, 8, 2, 4, 2, 4, 1, 6, 12, 4, 1, 4, 1, 4, 4, 1, 1,
	4, 4, 1, 2, 12, 2, 8, 6, 4, 6, 2, 2, 4, 1, 1, 2, 16, 1, 6, 1, 8, 4, 1, 1, 8, 6, 4, 2, 2, 2, 2, 2,
	4, 4, 1, 2, 4, 4, 10, 12, 4, 4, 2, 1, 2, 2, 24, 1, 8, 2, 1, 8, 2, 2, 2, 1, 8, 2, 1, 2, 8, 4, 1, 6,
	8, 1, 4, 6, 2, 2, 1, 8, 24, 1, 6, 2, 4, 1, 2, 2, 4, 12, 2, 1, 4, 12, 2, 2, 8, 1, 2, 2, 2, 4, 1, 2,
	16, 1, 12, 2, 2, 4, 4, 1, 4, 4, 4, 1, 2, 6, 1, 8, 8, 1, 2, 2, 4, 6, 1, 1, 8, 4, 4, 4, 12, 2, 4, 1,
	4, 2, 2, 2, 4, 1, 1, 18, 4, 12, 2, 2, 2, 4, 1, 2, 24, 1, 4, 2, 4, 1, 2, 8, 4, 2, 1, 1, 16, 1, 6, 20
};

void test_carmichael() {
	typedef CanonicFactors<uint_fast64_t, 15> canonic_factors_type;
	canonic_factors_type::primes_array_type primes_array;
	canonic_factors_type a(primes_array), b(primes_array), c(primes_array);
	unsigned int i, j, k;
	for (i=1; i<sizeof(phi_div_carmichael)/sizeof(phi_div_carmichael[0]); ++i) {
		a.assign(i);
		b = canonic_factors_type::carmichael(a);
		c = canonic_factors_type::eulers_phi(a);
		j = b.value();
		k = c.value();
		assert(k % j == 0);
		assert(k / j == phi_div_carmichael[i]);
	}
}

void tests_suite() {
	test_constructor_and_value();
	test_mul();
	test_eulers_phi();
	test_carmichael();
}

int main() {
	tests_suite();
	return 0;
}

