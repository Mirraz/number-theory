#ifndef CANONIC_FACTORS_H
#define CANONIC_FACTORS_H

#include <assert.h>
#include <stdint.h>
#include <math.h>

typedef uint_fast8_t pow_count_type;
typedef uint_fast8_t exp_type;
typedef uint_fast32_t num_type;
#ifndef NDEBUG
#  define POW_COUNT_PRI PRIuFAST8
#  define EXP_PRI PRIuFAST8
#  define NUM_PRI PRIuFAST32
#endif
#define NUM_MAX_MASK (((num_type)1) << 31)

// 2^31 <= UINT32_MAX < 2^32
#define MAX_EXP 31
// primorial(9) <= UINT32_MAX < primorial(10)
#define MAX_POW_COUNT 9

struct PrimePow {
	num_type prime;
	exp_type exp;
	
	PrimePow() : prime(0), exp(0) {}
	PrimePow(num_type b_prime, exp_type b_exp) : prime(b_prime), exp(b_exp) {}
	PrimePow(const PrimePow &b) : prime(b.prime), exp(b.exp) {}
};

class CanonicFactors {
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
	
	CanonicFactors(num_type n) {
		assert(n > 0);
		pow_count = 0;
		if (n == 1) return;
		num_type n_sqrt = round(sqrt(n));
		// TODO prepare and use table of small primes
		for (num_type d=2; d<=n_sqrt; d+=2) {
			if (n % d == 0) {
				exp_type exp = 0;
				do {
					n /= d;
					assert(exp < MAX_EXP);
					++exp;
				} while (n % d == 0);
				assert(pow_count < MAX_POW_COUNT);
				pows[pow_count++] = PrimePow(d, exp);
				n_sqrt = round(sqrt(n));
			}
			if (d == 2) --d;
		}
		if (n > 1) {
			if (pows[pow_count-1].prime == n) {
				++pows[pow_count-1].exp;
			} else {
				pows[pow_count++] = PrimePow(n, 1);
			}
		}
	}
	
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
		fprintf(stream, "%" POW_COUNT_PRI "[", pow_count);
		if (pow_count > 0) {
			for (pow_count_type i=0; i<pow_count; ++i) {
				fprintf(
					stream,
					"%" NUM_PRI "^%" EXP_PRI "%s",
					pows[i].prime,
					pows[i].exp,
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

