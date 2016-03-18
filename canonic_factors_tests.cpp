#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "canonic_factors.h"

void test_constructor_and_value() {
	num_type i, j;
	CanonicFactors a;
	for (i=1; i<1024; ++i) {
		a = i;
		//printf("%u = ", i); a.dump(); printf("\n");
		j = a.value();
		assert(i == j);
	}
	for (i=UINT32_MAX; i>=UINT32_MAX-1024; --i) {
		a = i;
		j = a.value();
		assert(i == j);
	}
}

void test_mul() {
	num_type i, j, k;
	CanonicFactors a, b, c;
	for (i=1; i<1024; ++i) {
		for (j=1; j<1024; ++j) {
			a = i;
			b = j;
			//printf("a = "); a.dump(); printf("\n");
			//printf("b = "); b.dump(); printf("\n");
			c = a * b;
			k = c.value();
			//printf("i = %u j = %u k = %u\n", i, j, k);
			assert(k == i * j);
		}
	}
}

void test_eulers_phi() {
	unsigned int i, j;
	CanonicFactors a, b;
	for (i=1; i<1024; ++i) {
		a = i;
		b = CanonicFactors::eulers_phi(a);
		j = b.value();
		printf("phi(%u) = %u\n", i, j);
	}
}

void test_carmichael() {
	unsigned int i, j;
	CanonicFactors a, b;
	for (i=1; i<1024; ++i) {
		a = i;
		b = CanonicFactors::carmichael(a);
		j = b.value();
		printf("carmichael(%u) = %u\n", i, j);
	}
}

void tests_suite() {
	test_constructor_and_value();
	test_mul();
	//test_eulers_phi();
	//test_carmichael();
}

int main() {
	tests_suite();
	return 0;
}

