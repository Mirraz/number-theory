#ifndef SQUARE_ROOT_MOD_H
#define SQUARE_ROOT_MOD_H

#include <assert.h>
#include "pow_mod.h"

template <typename NUM_TYPE, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class SquareRootMod {
public:
	typedef NUM_TYPE num_type;
private:
	typedef PowMod<num_type, NUM_TYPE_MAX_MASK, OPERATION_TYPE> pow_mod_type;
	
public:
	
	// legendre symbol (a / p)
	// p - odd prime
	static int legendre_symbol(num_type p, num_type a) {
		assert(p > 2);
		num_type pow = pow_mod_type::pow_mod(p, a, (p-1)>>1);
		if (pow == 1) {
			return 1;
		} else if (pow == p-1) {
			return -1;
		} else {
			return 0;
		}
	}
	
	// p - odd prime
	// primes_count enough to find least nonresidue: least nonresidue <= sqrt(p) * ln(p) + 1
	// returns 0 if least nonresidue is not found
	static num_type least_nonresidue(num_type primes[], size_t primes_count, num_type p) {
		assert(p > 2);
		for (size_t i=0; i < primes_count; ++i) {
			if (legendre_symbol(p, primes[i]) == -1) return primes[i];
		}
		return 0;
	}
};

#endif/*SQUARE_ROOT_MOD_H*/

