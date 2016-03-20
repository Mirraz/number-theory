#ifndef SQUARE_ROOT_MOD_H
#define SQUARE_ROOT_MOD_H

#include <assert.h>
#include "mul_mod.h"

template <typename NUM_TYPE, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class SquareRootMod {
public:
	typedef NUM_TYPE num_type;
private:
	typedef OPERATION_TYPE operation_type;
public:
	typedef MulMod<num_type, NUM_TYPE_MAX_MASK, operation_type> mul_mod_type;
private:
	static_assert(sizeof(num_type) <= 32, "NUM_TYPE is too big for pow_count_type");
	typedef uint_fast8_t pow_count_type;
	
public:
	// legendre symbol (a / p)
	// p - odd prime
	static int legendre_symbol(num_type p, num_type a) {
		assert(p > 2);
		num_type pow = mul_mod_type::pow_mod(p, a, (p-1)>>1);
		if (pow == 1) {
			return 1;
		} else if (pow == p-1) {
			return -1;
		} else {
			return 0;
		}
	}
	
	// p - odd prime
	// primes_count enough to find least nonresidue
	//     least nonresidue <= sqrt(p) * ln(p) + 1
	//     least nonresidue << O(log(p)^2) (on assumption of Generalised Riemann hypothesis)
	// returns 0 if least nonresidue is not found
	static num_type least_nonresidue(num_type primes[], size_t primes_count, num_type p) {
		assert(p > 2);
		for (size_t i=0; i < primes_count; ++i) {
			if (legendre_symbol(p, primes[i]) == -1) return primes[i];
		}
		return 0;
	}
	
	static num_type square_root_mod_algo_01(num_type p, num_type nr, num_type a) {
		assert(p > 2);
		assert(a > 0);
		if (mul_mod_type::mul_mod(p, a, (p-1)>>1) != 1) return 0;
		
		num_type k = p - 1;
		pow_count_type h = 0;
		while (!(k & 1)) {
			k >>= 1;
			++h;
		}
		k = (k - 1) >> 1;
		++h;
		assert(h >= 2);
		
		operation_type p_op = p;
		num_type z_count = 1 << (h-2);
		for (num_type z=0; z<z_count; ++z) {
			operation_type np = mul_mod_type::pow_mod(p, nr, z*(2*k+1));
			operation_type ap = mul_mod_type::pow_mod(p, a, k+1);
			operation_type x = (np * ap) % p_op;
			if ((x * x) % p_op == a) return x;
		}
		assert(0);
		return 0;
	}
	
	// p - odd prime
	// a - quadratic residue modulo p
	// nr - quadratic nonresidue modulo p
	// if p = 4k+3 then nr is not used and may be any number
	static num_type tonelli_shanks_algo(num_type p, num_type nr, num_type a) {
		assert(p > 2);
		assert(a > 0);
		assert(legendre_symbol(p, a) == 1);
		
		num_type q = p - 1;
		pow_count_type m = 0;
		while (!(q & 1)) {
			q >>= 1;
			++m;
		}
		
		if (m == 1) {
			assert((p+1) % 4 == 0);
			return mul_mod_type::pow_mod(p, a, (p+1)>>2);
		}
		assert(nr > 1);
		assert(legendre_symbol(p, nr) == -1);
		
		num_type r = mul_mod_type::pow_mod(p, a, (q+1)>>1);
		num_type t = mul_mod_type::pow_mod(p, a, q);
		num_type c = mul_mod_type::pow_mod(p, nr, q);
		while (t != 1) {
			pow_count_type i = 0;
			num_type tpow = t;
			do {
				tpow = mul_mod_type::square_mod(p, tpow);
				++i;
			} while (tpow != 1);
			assert(i < m);
			num_type b = mul_mod_type::pow_mod(p, c, ((num_type)1)<<(m-i-1));
			num_type bsq = mul_mod_type::square_mod(p, b);
			r = mul_mod_type::mul_mod(p, r, b);
			t = mul_mod_type::mul_mod(p, t, bsq);
			c = bsq;
			m = i;
		}
		return r;
	}
	
	// solve x^2 = a (mod p)
	// p - prime
	// nr - quadratic nonresidue for p (used only if a is residue and p != 4k+3)
	// if a is not quadratic residue returns 0
	// else returns any (of two) square root
	static num_type square_root_mod(num_type p, num_type nr, num_type a) {
		if (mul_mod_type::pow_mod(p, a, (p-1)>>1) != 1) return 0;
		if (p == 2) return 1;
		return tonelli_shanks_algo(p, nr, a);
	}
};

#endif/*SQUARE_ROOT_MOD_H*/

