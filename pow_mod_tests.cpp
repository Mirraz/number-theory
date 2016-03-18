#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "pow_mod.h"

void test_pow_mod() {
	typedef PowMod<uint_fast32_t, ((uint_fast32_t)1)<<31, uint_fast64_t> pow_mod_type;
	pow_mod_type::num_type mod, base, exp, pow;
	for (mod = 2; mod<256; ++mod) {
		for (base = 0; base<256; ++base) {
			for (exp = 0; exp < 16; ++exp) {
				pow = pow_mod_type::pow_mod(mod, base, exp);
				printf(
					"mod = %" PRIuFAST32 " base = %" PRIuFAST32 " exp = %" PRIuFAST32 " pow = %" PRIuFAST32 "\n",
					mod, base, exp, pow
				);
			}
		}
	}
}

void tests_suite() {
	test_pow_mod();
}

int main() {
	tests_suite();
	return 0;
}

