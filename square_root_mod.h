#ifndef SQUARE_ROOT_MOD_H
#define SQUARE_ROOT_MOD_H

#include <assert.h>
#include <algorithm>
#include "mul_mod.h"

template <typename NUM_TYPE, uint_fast8_t NUM_TYPE_LEN, typename OPERATION_TYPE>
class SquareRootMod {
public:
	typedef NUM_TYPE num_type;
	static_assert(NUM_TYPE_LEN > 0, "NUM_TYPE_LEN can't be zero");
	static_assert(sizeof(num_type) <= 32, "NUM_TYPE is too big for NUM_TYPE_LEN");
	typedef MulMod<num_type, ((num_type)1)<<(NUM_TYPE_LEN-1), OPERATION_TYPE> mul_mod_type;
private:
	static_assert(sizeof(num_type) <= 32, "NUM_TYPE is too big for exp_type");
	typedef uint_fast8_t exp_type;
	
	// z == nr - quadratic nonresidue
	// zq = z^q mod p
	// zq_pows[i] = (zq)^(2^i) mod p
	num_type zq_pows[NUM_TYPE_LEN];
	// p-1 == q * 2^s
	// p_1d2 = (p-1)/2, q1d2 = (q+1)/2
	num_type p, p_1d2, q, q1d2;
	exp_type s;
	
	SquareRootMod() = delete;
	SquareRootMod(const SquareRootMod &b) = delete;
	SquareRootMod& operator=(const SquareRootMod &b) = delete;
	
public:
	// only for prime modulo = 4k+3
	SquareRootMod(num_type modulo) : p(modulo) {
		assert(p % 4 == 3);
		p_1d2 = q = (p-1) >> 1;
		q1d2 = (q+1) >> 1;
		s = 1;
	}
	
	// modulo - odd prime
	// nr - quadratic nonresidue modulo p (if needed)
	SquareRootMod(num_type modulo, num_type nr) : p(modulo) {
		assert(p > 2);
		assert((p & 1) == 1);
		p_1d2 = (p-1) >> 1;
		if ((p & 3) == 3) {
			q = p_1d2;
			s = 1;
		} else {
			assert(legendre_symbol(nr) == -1);
			q = p - 1;
			s = 0;
			do {
				q >>= 1;
				++s;
			} while (!(q & 1));
			std::fill(zq_pows, zq_pows+NUM_TYPE_LEN, 0);
			zq_pows[0] = mul_mod_type::pow_mod(p, nr, q);
		}
		q1d2 = (q+1) >> 1;
	}
	
	SquareRootMod(num_type modulo, num_type primes[], size_t primes_count) :
		SquareRootMod(modulo, ((modulo & 3) == 3 ? 0 : least_nonresidue(primes, primes_count, modulo))) {}
	
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
	
	int legendre_symbol(num_type a) const {
		num_type pow = mul_mod_type::pow_mod(p, a, p_1d2);
		if (pow == 1) {
			return 1;
		} else if (pow+1 == p) {
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
		if (legendre_symbol(p, a) != 1) return 0;
		
		num_type k = p - 1;
		exp_type h = 0;
		while (!(k & 1)) {
			k >>= 1;
			++h;
		}
		k = (k - 1) >> 1;
		++h;
		assert(h >= 2);
		
		num_type z_count = 1 << (h-2);
		for (num_type z=0; z<z_count; ++z) {
			num_type np = mul_mod_type::pow_mod(p, nr, z*(2*k+1));
			num_type ap = mul_mod_type::pow_mod(p, a, k+1);
			num_type x  = mul_mod_type::mul_mod(p, np, ap);
			if (mul_mod_type::square_mod(p, x) == a) return x;
		}
		assert(0);
		return 0;
	}
	
	// a - quadratic residue modulo p
	num_type tonelli_shanks_algo(num_type a) {
		assert(legendre_symbol(p, a) == 1);
		num_type r = mul_mod_type::pow_mod(p, a, q1d2);
		if (s == 1) return r;
		num_type t = mul_mod_type::pow_mod(p, a, q);
		num_type c = zq_pows[0];
		exp_type m = s;
		while (t != 1) {
			exp_type i = 0;
			num_type tpow = t;
			do {
				tpow = mul_mod_type::square_mod(p, tpow);
				++i;
			} while (tpow != 1);
			assert(i < m);
			
			num_type b = (
				zq_pows[s-i-1] == 0 ?
				zq_pows[s-i-1] = mul_mod_type::pow_mod(p, c, ((num_type)1)<<(m-i-1)) :
				zq_pows[s-i-1]
			);
			num_type bsq = (
				zq_pows[s-i] == 0 ?
				zq_pows[s-i] = mul_mod_type::square_mod(p, b) :
				zq_pows[s-i]
			);
			
			r = mul_mod_type::mul_mod(p, r, b);
			t = mul_mod_type::mul_mod(p, t, bsq);
			c = bsq;
			m = i;
		}
		return r;
	}
	
	// solve x^2 = a (mod p)
	// if a is not quadratic residue returns 0
	// else returns any (of two) square root
	num_type square_root_mod(num_type a) {
		if (legendre_symbol(a) != 1) return 0;
		return tonelli_shanks_algo(a);
	}
};

#endif/*SQUARE_ROOT_MOD_H*/

