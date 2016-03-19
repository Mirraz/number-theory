#ifndef MUL_GROUP_MOD_H
#define MUL_GROUP_MOD_H

#include <assert.h>
#include "canonic_factors.h"
#include "pow_mod.h"

template <typename NUM_TYPE, uint_fast8_t MAX_POW_COUNT, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class MulGroupMod {
public:
	typedef NUM_TYPE num_type;
private:
	typedef PowMod<num_type, NUM_TYPE_MAX_MASK, OPERATION_TYPE> pow_mod_type;
	typedef CanonicFactors<num_type, MAX_POW_COUNT> canonic_factors_type;
	typedef typename canonic_factors_type::PrimePow prime_pow_type;
	typedef typename canonic_factors_type::pow_count_type pow_count_type;
public:
	typedef typename canonic_factors_type::primes_array_type primes_array_type;
	
private:
	num_type modulo;
	canonic_factors_type group_exponent;
	
	MulGroupMod() = delete;
	MulGroupMod(const MulGroupMod &b) = delete;
	MulGroupMod& operator=(const MulGroupMod &b) = delete;
	
public:
	MulGroupMod(primes_array_type primes_array) : modulo(0), group_exponent(primes_array) {}
	
	MulGroupMod(primes_array_type primes_array, num_type b_modulo) : MulGroupMod(primes_array) {
		assign(b_modulo);
	}
	
	void assign(num_type b_modulo) {
		assert(b_modulo > 1);
		modulo = b_modulo;
		group_exponent.assign(modulo);
		group_exponent = canonic_factors_type::carmichael(group_exponent);
	}
	
	// gcd(modulo, element) == 1
	num_type element_order(num_type element) const {
		assert(modulo > 1);
		assert(pow_mod_type::pow_mod(modulo, element, group_exponent.value()) == 1);
		prime_pow_type exp_pows[MAX_POW_COUNT];
		pow_count_type exp_pow_count = group_exponent.copy(exp_pows, MAX_POW_COUNT);
		pow_count_type i;
		do {
			for (i=0; i<exp_pow_count; ++i) {
				if (exp_pows[i].exp == 0) continue;
				--exp_pows[i].exp;
				num_type exp_value = canonic_factors_type::value(exp_pows, exp_pow_count);
				if (pow_mod_type::pow_mod(modulo, element, exp_value) == 1) break;
				++exp_pows[i].exp;
			}
		} while (i<exp_pow_count);
		num_type exp_value = canonic_factors_type::value(exp_pows, exp_pow_count);
		assert(pow_mod_type::pow_mod(modulo, element, exp_value) == 1);
		return exp_value;
	}
};

#endif/*MUL_GROUP_MOD_H*/

