#ifndef PRIMITIVE_ROOTS_H
#define PRIMITIVE_ROOTS_H

#include <assert.h>
#include <stdint.h>
#include "factorize.h"
#include "pow_mod.h"

template <typename NUM_TYPE, uint_fast8_t MAX_PRIME_FACTORS_COUNT, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class PrimitiveRoots {
public:
	typedef NUM_TYPE num_type;
private:
	typedef GetPrimeFactors<num_type> get_prime_factors_type;
	typedef PowMod<num_type, NUM_TYPE_MAX_MASK, OPERATION_TYPE> pow_mod_type;
public:
	typedef typename get_prime_factors_type::prime_factors_count_type prime_factors_count_type;
	typedef typename get_prime_factors_type::primes_array_type primes_array_type;
	
private:
	num_type modulo;
	num_type exps[MAX_PRIME_FACTORS_COUNT];
	prime_factors_count_type exps_count;
	
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
		
		num_type prime_factors[MAX_PRIME_FACTORS_COUNT];
		get_prime_factors_type get_prime_factors(b_primes_array);
		exps_count = get_prime_factors.get_prime_factors(modulo_1, prime_factors, MAX_PRIME_FACTORS_COUNT);
		
		for (prime_factors_count_type i=0; i<exps_count; ++i) {
			assert(modulo_1 > prime_factors[i]);
			assert(modulo_1 % prime_factors[i] == 0);
			exps[exps_count-i-1] = modulo_1 / prime_factors[i];
		}
	}
	
	// gcd(modulo, root) == 1
	bool is_primitive_root(num_type root) {
		assert(pow_mod_type::pow_mod(modulo, root, modulo-1) == 1);
		for (prime_factors_count_type i=0; i<exps_count; ++i) {
			num_type pow = pow_mod_type::pow_mod(modulo, root, exps[i]);
			if (pow == 1) return false;
		}
		return true;
	}
};

#endif/*PRIMITIVE_ROOTS_H*/

