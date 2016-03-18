#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "pow_mod.h"

void test_pow_mod() {
	typedef PowMod<uint_fast32_t, ((uint_fast32_t)1)<<31, uint_fast64_t> pow_mod_type;
	pow_mod_type::num_type mod, base, exp, pow;
	mod = 4294967291;
	base = 2;
	exp = (mod-1)/2;
	pow = pow_mod_type::pow_mod(mod, base, exp);
	printf("%u^%u mod %u = %u\n", base, exp, mod, pow);
}

void tests_suite() {
	test_pow_mod();
}

int main() {
	tests_suite();
	return 0;
}

