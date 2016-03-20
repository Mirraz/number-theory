#ifndef CANONIC_FACTORS_H
#define CANONIC_FACTORS_H

#include <assert.h>
#include <stdint.h>
#include <stddef.h>		// NULL
#ifndef NDEBUG
#  include <stdio.h>
#endif
#include <functional>
#include <algorithm>
#include "factorize.h"

// primorial(4)  < 2^8  < primorial(5)
// primorial(6)  < 2^16 < primorial(7)
// primorial(9)  < 2^32 < primorial(10)
// primorial(15) < 2^64 < primorial(16)

template<typename NUM_TYPE, uint_fast8_t MAX_POW_COUNT>
class CanonicFactorsTemplate {
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


class CanonicFactorizer {
private:
	typedef Factorizer<num_type> factorizer_type;
public:
	typedef typename factorizer_type::primes_array_type primes_array_type;

private:
	factorizer_type factorizer;
	PrimePow *m_pows;
	pow_count_type pow_count;
	
	CanonicFactorizer(const CanonicFactorizer &b) = delete;
	CanonicFactorizer& operator=(const CanonicFactorizer &b) = delete;

public:
	CanonicFactorizer(primes_array_type primes_array) :
		factorizer(
			primes_array,
			std::bind(
				&CanonicFactorizer::factorize_cb,
				this,
				std::placeholders::_1,
				std::placeholders::_2
			)
		) {}
	
	CanonicFactorizer() : CanonicFactorizer(primes_array_type()) {}

private:
	bool factorize_cb(typename factorizer_type::num_type prime, typename factorizer_type::exp_type exp) {
		assert(pow_count < MAX_POW_COUNT);
		m_pows[pow_count++] = PrimePow(prime, exp);
		return false;
	}

public:
	pow_count_type factorize(num_type n, PrimePow pows[]) {
		pow_count = 0;
		m_pows = pows;
		factorizer.factorize(n);
		pows = NULL;
		return pow_count;
	}
};


class CanonicFactors {
private:
	PrimePow pows[MAX_POW_COUNT];
	CanonicFactorizer &factorizer;
	pow_count_type pow_count;
	
public:
	CanonicFactors(CanonicFactorizer &b_factorizer) : factorizer(b_factorizer), pow_count(0) {}
	
	CanonicFactors(const CanonicFactors &b) : factorizer(b.factorizer), pow_count(b.pow_count) {
		std::copy(b.pows, b.pows+b.pow_count, pows);
	}
	
	CanonicFactors(CanonicFactors &&b) : factorizer(b.factorizer), pow_count(b.pow_count) {
		std::copy(b.pows, b.pows+b.pow_count, pows);
		b.pow_count = 0;
	}
	
private:
	void assign(const CanonicFactors &b) {
		pow_count = b.pow_count;
		std::copy(b.pows, b.pows+b.pow_count, pows);
	}
	
public:
	CanonicFactors& operator=(CanonicFactors &&b) {
		if (this != &b) {
			assign(b);
			b.pow_count = 0;
		}
		return *this;
	}
	
	CanonicFactors& operator=(const CanonicFactors &b) {
		if (this != &b) assign(b);
		return *this;
	}
	
	CanonicFactors(CanonicFactorizer &b_factorizer, const PrimePow &b) : factorizer(b_factorizer), pow_count(1) {
		pows[0] = b;
	}
	
	CanonicFactors(CanonicFactorizer &b_factorizer, num_type n) : factorizer(b_factorizer) {
		assign(n);
	}
	
	void assign(num_type n) {
		if (n == 1) {		// just for optimisation
			pow_count = 0;
			return;
		}
		pow_count = factorizer.factorize(n, pows);
	}
	
	// some pows[i].exp may be == 0
	static num_type value(const PrimePow pows[], pow_count_type pow_count) {
		num_type n = 1;
		for (pow_count_type i=0; i<pow_count; ++i) {
			// TODO fast pow
			for (exp_type j=0; j<pows[i].exp; ++j) {
				n *= pows[i].prime;
			}
		}
		return n;
	}
	
	num_type value() const {
		return value(pows, pow_count);
	}
	
	pow_count_type copy(PrimePow result_pows[], pow_count_type result_pows_size) const {
		assert(result_pows_size >= pow_count);
		std::copy(pows, pows+pow_count, result_pows);
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
	
private:
	// TODO maybe more optimal mul_assign_static
	static CanonicFactors mul_static(const CanonicFactors &a, const CanonicFactors &b) {
		CanonicFactors result(a.factorizer);
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
	
public:
	CanonicFactors operator*(const CanonicFactors &b) const {
		return mul_static(*this, b);
	}
	
	CanonicFactors operator*(num_type b) const {
		return (*this) * CanonicFactors(factorizer, b);
	}
	
	CanonicFactors& operator *=(const CanonicFactors &b) {
		CanonicFactors result = (*this) * b;
		*this = result;
		return *this;
	}
	
	CanonicFactors& operator *=(num_type b) {
		return (*this) *= CanonicFactors(factorizer, b);
	}
	
private:
	static CanonicFactors mul_pow_static(const CanonicFactors &a, const PrimePow &b) {
		CanonicFactors result(a.factorizer);
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
		return result;
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
	
public:
	CanonicFactors mul_pow(const PrimePow &b) const {
		return mul_pow_static(*this, b);
	}
	
	void mul_pow_assign(const PrimePow &b) {
		mul_pow_assign_static(*this, b);
	}
	
public:
	static CanonicFactors lcm(const CanonicFactors &a, const CanonicFactors &b) {
		CanonicFactors result(a.factorizer);
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
				result.pows[result.pow_count++] = PrimePow(a.pows[i].prime, std::max(a.pows[i].exp, b.pows[j].exp));
				++i; ++j;
			}
		}
		return result;
	}
	
	static CanonicFactors eulers_phi(const CanonicFactors &b) {
		CanonicFactors result(b.factorizer);
		for (pow_count_type i=0; i<b.pow_count; ++i) {
			result *= b.pows[i].prime - 1;
			if (b.pows[i].exp > 1) result.mul_pow_assign(PrimePow(b.pows[i].prime, b.pows[i].exp - 1));
		}
		return result;
	}
	
private:
	static CanonicFactors eulers_phi_pow(CanonicFactorizer &b_factorizer, const PrimePow &b) {
		CanonicFactors result(b_factorizer, b.prime - 1);
		if (b.exp > 1) result.mul_pow_assign(PrimePow(b.prime, b.exp - 1));
		return result;
	}
	
public:
	static CanonicFactors carmichael(const CanonicFactors &b) {
		if (b.pow_count == 0) return CanonicFactors(b.factorizer);
		CanonicFactors result(
			b.pows[0].prime == 2 && b.pows[0].exp > 2 ?
			CanonicFactors(b.factorizer, PrimePow(2, b.pows[0].exp - 2)) :
			eulers_phi_pow(b.factorizer, b.pows[0])
		);
		for (pow_count_type i=1; i<b.pow_count; ++i) {
			result = lcm(result, eulers_phi_pow(b.factorizer, b.pows[i]));
		}
		return result;
	}
};


};

#endif/*CANONIC_FACTORS_H*/

