#ifndef CANONIC_FACTORS_H
#define CANONIC_FACTORS_H

#include <assert.h>
#include <stdint.h>
#ifndef NDEBUG
#  include <stdio.h>
#endif
#include "factorize.h"

// primorial(4)  < 2^8  < primorial(5)
// primorial(6)  < 2^16 < primorial(7)
// primorial(9)  < 2^32 < primorial(10)
// primorial(15) < 2^64 < primorial(16)

template<typename NUM_TYPE, uint_fast8_t MAX_POW_COUNT>
class CanonicFactors {
public:
	typedef NUM_TYPE num_type;
	static_assert(sizeof(num_type) <= 32, "Too big num_type for exp_type");
	typedef uint_fast8_t exp_type;
	// prevent sum overflow
	static constexpr uint_fast8_t MAX_EXP = 127;
	// 256^286 < 2^2289 < primorial(256) < 2^2290 < 256^287
	static_assert(sizeof(num_type) <= 286, "Too big num_type for pow_count_type");
	typedef uint_fast8_t pow_count_type;
	static_assert(MAX_POW_COUNT > 0, "MAX_POW_COUNT can't be zero");
	
	struct PrimePow {
		num_type prime;
		exp_type exp;
		
		PrimePow() : prime(0), exp(0) {}
		PrimePow(num_type b_prime, exp_type b_exp) : prime(b_prime), exp(b_exp) {}
		PrimePow(const PrimePow &b) : prime(b.prime), exp(b.exp) {}
	};
	typedef PrimePow prime_pow_type;
	
	typedef Factorizer<num_type> factorizer_type;
	typedef typename factorizer_type::primes_array_type primes_array_type;
	
private:
	pow_count_type pow_count;
	PrimePow pows[MAX_POW_COUNT];
	
public:
	CanonicFactors() : pow_count(0) {}
	
	CanonicFactors(const CanonicFactors &b) {
		pow_count = b.pow_count;
		for (pow_count_type i=0; i<b.pow_count; ++i) pows[i] = b.pows[i];
	}
	
	CanonicFactors(const PrimePow &b) {
		pow_count = 1;
		pows[0] = b;
	}
	
	CanonicFactors(primes_array_type b_primes_array, num_type n) {
		pow_count = 0;
		typename factorizer_type::factorize_cb_type cb =
				[this] (typename factorizer_type::num_type prime, typename factorizer_type::exp_type exp) -> bool {
			assert(pow_count < MAX_POW_COUNT);
			pows[pow_count++] = PrimePow(prime, exp);
			return false;
		};
		factorizer_type factorizer(b_primes_array, cb);
		factorizer.factorize(n);
	}
	
	CanonicFactors(num_type n) : CanonicFactors(primes_array_type(), n) {}
	
	num_type value() const {
		num_type n = 1;
		for (pow_count_type i=0; i<pow_count; ++i) {
			// TODO fast pow
			for (exp_type j=0; j<pows[i].exp; ++j) {
				n *= pows[i].prime;
			}
		}
		return n;
	}
	
	pow_count_type copy(PrimePow result_pows[]) const {
		for (pow_count_type i=0; i<pow_count; ++i) result_pows[i] = pows[i];
		return pow_count;
	}
	
#ifndef NDEBUG
	void fdump(FILE *stream) const {
		fprintf(stream, "%u[", (unsigned int)pow_count);
		if (pow_count > 0) {
			for (pow_count_type i=0; i<pow_count; ++i) {
				fprintf(
					stream,
					"%llu^%u%s",
					(unsigned long long)pows[i].prime,
					(unsigned int)pows[i].exp,
					(i != pow_count-1 ? " " : "")
				);
			}
		}
		fprintf(stream, "]");
	}
	
	void dump() const {
		fdump(stdout);
	}
#endif
	
	// TODO maybe more optimal mul_assign_static
	static CanonicFactors mul_static(const CanonicFactors &a, const CanonicFactors &b) {
		CanonicFactors result;
		result.pow_count = 0;
		pow_count_type i = 0, j = 0;
		while (i < a.pow_count || j < b.pow_count) {
			if        (i < a.pow_count && (j == b.pow_count || a.pows[i].prime < b.pows[j].prime)) {
				assert(i < a.pow_count);
				assert(result.pow_count < MAX_POW_COUNT);
				result.pows[result.pow_count++] = a.pows[i++];
			} else if (j < b.pow_count && (i == a.pow_count || a.pows[i].prime > b.pows[j].prime)) {
				assert(j < b.pow_count);
				assert(result.pow_count < MAX_POW_COUNT);
				result.pows[result.pow_count++] = b.pows[j++];
			} else {
				assert(i < a.pow_count);
				assert(j < b.pow_count);
				assert(a.pows[i].prime == b.pows[j].prime);
				assert(result.pow_count < MAX_POW_COUNT);
				assert(a.pows[i].exp + b.pows[j].exp <= MAX_EXP);
				result.pows[result.pow_count++] = PrimePow(a.pows[i].prime, a.pows[i].exp + b.pows[j].exp);
				++i; ++j;
			}
		}
		return result;
	}
	
	CanonicFactors operator*(const CanonicFactors &b) const {
		return mul_static(*this, b);
	}
	
	CanonicFactors operator*(num_type b) const {
		return (*this) * CanonicFactors(b);
	}
	
	CanonicFactors& operator *=(const CanonicFactors &b) {
		CanonicFactors result = (*this) * b;
		*this = result;
		return *this;
	}
	
	CanonicFactors& operator *=(num_type b) {
		return (*this) *= CanonicFactors(b);
	}
	
	// maybe result === a, but it will be not efficient
	static void mul_pow_static(CanonicFactors &result, const CanonicFactors &a, const PrimePow &b) {
		pow_count_type i;
		for (i=0; i<a.pow_count && a.pows[i].prime < b.prime; ++i) {
			result.pows[i] = a.pows[i];
		}
		if (i < a.pow_count && a.pows[i].prime == b.prime) {
			assert(a.pows[i].exp < MAX_EXP);
			result.pows[i] = PrimePow(a.pows[i].prime, a.pows[i].exp+1);
			++i;
			for (; i<a.pow_count; ++i) {
				result.pows[i] = a.pows[i];
			}
			result.pow_count = a.pow_count;
		} else {
			assert(a.pow_count < MAX_POW_COUNT);
			for (pow_count_type j=a.pow_count; j>i ; --j) {
				result.pows[j] = a.pows[j-1];
			}
			result.pows[i] = b;
			result.pow_count = a.pow_count + 1;
		}
	}
	
	static void mul_pow_assign_static(CanonicFactors &result, const PrimePow &b) {
		pow_count_type i;
		for (i=0; i<result.pow_count && result.pows[i].prime < b.prime; ++i);
		if (i < result.pow_count && result.pows[i].prime == b.prime) {
			assert(result.pows[i].exp < MAX_EXP);
			++result.pows[i].exp;
		} else {
			assert(result.pow_count < MAX_POW_COUNT);
			for (pow_count_type j=result.pow_count; j>i ; --j) {
				result.pows[j] = result.pows[j-1];
			}
			result.pows[i] = b;
			++result.pow_count;
		}
	}
	
	CanonicFactors mul_pow(const PrimePow &b) const {
		CanonicFactors result;
		mul_pow_static(result, *this, b);
		return result;
	}
	
	void mul_pow_assign(const PrimePow &b) {
		mul_pow_assign_static(*this, b);
	}
	
private:
	template <typename T>
	static T max(T a, T b) {
		if (a < b) return b;
		return a;
	}
	
public:
	static CanonicFactors lcm(const CanonicFactors &a, const CanonicFactors &b) {
		CanonicFactors result;
		result.pow_count = 0;
		pow_count_type i = 0, j = 0;
		while (i < a.pow_count || j < b.pow_count) {
			if        (i < a.pow_count && (j == b.pow_count || a.pows[i].prime < b.pows[j].prime)) {
				assert(i < a.pow_count);
				assert(result.pow_count < MAX_POW_COUNT);
				result.pows[result.pow_count++] = a.pows[i++];
			} else if (j < b.pow_count && (i == a.pow_count || a.pows[i].prime > b.pows[j].prime)) {
				assert(j < b.pow_count);
				assert(result.pow_count < MAX_POW_COUNT);
				result.pows[result.pow_count++] = b.pows[j++];
			} else {
				assert(i < a.pow_count);
				assert(j < b.pow_count);
				assert(a.pows[i].prime == b.pows[j].prime);
				assert(result.pow_count < MAX_POW_COUNT);
				assert(a.pows[i].exp + b.pows[j].exp <= MAX_EXP);
				result.pows[result.pow_count++] = PrimePow(a.pows[i].prime, max(a.pows[i].exp, b.pows[j].exp));
				++i; ++j;
			}
		}
		return result;
	}
	
	static CanonicFactors eulers_phi(const CanonicFactors &b) {
		CanonicFactors result = 1;
		for (pow_count_type i=0; i<b.pow_count; ++i) {
			result *= CanonicFactors(b.pows[i].prime - 1);
			if (b.pows[i].exp > 1) result.mul_pow_assign(PrimePow(b.pows[i].prime, b.pows[i].exp - 1));
		}
		return result;
	}
	
	static CanonicFactors eulers_phi_pow(const PrimePow &b) {
		CanonicFactors result = CanonicFactors(b.prime - 1);
		if (b.exp > 1) result.mul_pow_assign(PrimePow(b.prime, b.exp - 1));
		return result;
	}
	
	static CanonicFactors carmichael(const CanonicFactors &b) {
		if (b.pow_count == 0) return 1;
		CanonicFactors result;
		if (b.pows[0].prime == 2 && b.pows[0].exp > 2) {
			result = PrimePow(2, b.pows[0].exp - 2);
		} else {
			result = eulers_phi_pow(b.pows[0]);
		}
		for (pow_count_type i=1; i<b.pow_count; ++i) {
			result = lcm(result, eulers_phi_pow(b.pows[i]));
		}
		return result;
	}
};

#endif/*CANONIC_FACTORS_H*/

