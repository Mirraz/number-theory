#ifndef PRIMITIVE_ROOTS_H
#define PRIMITIVE_ROOTS_H

#include <assert.h>
#include <stdint.h>
#include "canonic_factors.h"
#include "mul_mod.h"

template <typename NUM_TYPE, uint_fast8_t MAX_POW_COUNT, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class PrimitiveRoots {
public:
	typedef NUM_TYPE num_type;
private:
	typedef MulMod<num_type, NUM_TYPE_MAX_MASK, OPERATION_TYPE> mul_mod_type;
	typedef CanonicFactors<num_type, MAX_POW_COUNT> canonic_factors_type;
	typedef typename canonic_factors_type::pow_count_type pow_count_type;
	typedef typename canonic_factors_type::PrimePow prime_pow_type;
public:
	typedef typename canonic_factors_type::primes_array_type primes_array_type;
	
private:
	num_type exps[MAX_POW_COUNT];
	num_type modulo;
	pow_count_type exps_count;
	
	PrimitiveRoots() = delete;
	PrimitiveRoots(const PrimitiveRoots &b) = delete;
	PrimitiveRoots& operator=(const PrimitiveRoots &b) = delete;
	
public:
	// modulo must be prime
	PrimitiveRoots(primes_array_type b_primes_array, num_type b_modulo) : modulo(b_modulo) {
		assert(modulo >= 2);
		if (modulo <= 3) {
			exps_count = 0;
			return;
		}
		
		num_type modulo_1 = modulo - 1;
		canonic_factors_type phi(b_primes_array, modulo_1);
		prime_pow_type prime_factors[MAX_POW_COUNT];
		exps_count = phi.copy(prime_factors, MAX_POW_COUNT);
		
		for (pow_count_type i=0; i<exps_count; ++i) {
			assert(modulo_1 > prime_factors[i].prime);
			assert(modulo_1 % prime_factors[i].prime == 0);
			exps[exps_count-i-1] = modulo_1 / prime_factors[i].prime;
		}
	}
	
	// gcd(modulo, root) == 1
	bool is_primitive_root(num_type root) const {
		assert(mul_mod_type::pow_mod(modulo, root, modulo-1) == 1);
		for (pow_count_type i=0; i<exps_count; ++i) {
			num_type pow = mul_mod_type::pow_mod(modulo, root, exps[i]);
			if (pow == 1) return false;
		}
		return true;
	}
};

#endif/*PRIMITIVE_ROOTS_H*/

