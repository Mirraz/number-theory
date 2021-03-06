#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <assert.h>
#include <stddef.h>		// size_t, NULL
#include <stdint.h>
#include <math.h>
#include <functional>

template <typename NUM_TYPE> class Factorizer;

template <typename NUM_TYPE>
struct PrimesArray {
public:
	typedef NUM_TYPE num_type;

	num_type *primes;
	size_t count;
	
	PrimesArray() : primes(NULL), count(0) {}
	PrimesArray(num_type *b_primes, size_t b_count) : primes(b_primes), count(b_count) {}
	// use default copy constructor and assignment operator
	
	static size_t fill_primes(num_type primes[], size_t primes_size, num_type max_num) {
		return Factorizer<num_type>::fill_primes(primes, primes_size, max_num);
	}
};

template <typename NUM_TYPE>
class Factorizer {
public:
	typedef NUM_TYPE num_type;
	typedef PrimesArray<num_type> primes_array_type;
	static_assert(sizeof(num_type) <= 32, "Too big num_type for exp_type");
	typedef uint_fast8_t exp_type;
	typedef std::function<bool(num_type prime, exp_type exp)> factorize_cb_type;
	
private:
	primes_array_type primes_array;
	factorize_cb_type cb;
	
public:
	Factorizer(primes_array_type b_primes_array, factorize_cb_type b_cb) :
		primes_array(b_primes_array), cb(b_cb) {}
	
	// use default copy constructor and assignment operator
	
	primes_array_type get_primes_array() const {
		return primes_array;
	}
	
private:
	static_assert(sizeof(num_type) <= 8, "Too big num_type for round_sqrt");
	static inline num_type round_sqrt(num_type n) {
		return round(sqrt((double)n));
	}
	
public:
	#define NUM_MAX UINT64_MAX
	#define SQRT_MAX ((num_type)UINT32_MAX)
	static void test_round_sqrt() {
		const num_type max_sqrt = SQRT_MAX + (num_type)1;
		const num_type max_delta = 1024;
		for (num_type i=max_sqrt; i>=max_sqrt-1024; --i) {
			//printf("%llu\n", i);
			for (num_type j=0; j<=max_delta; ++j) {
				num_type s1, s2;
				if (i == max_sqrt) {
					if (j == 0) s1 = NUM_MAX;
					else s1 = NUM_MAX - j + 1;
					s2 = NUM_MAX;
				} else {
					num_type s = i * i;
					s1 = s - j;
					s2 = s + j;
				}
				//printf("%llu %llu\n", s1, s2);
				num_type r1 = round_sqrt(s1), r2 = round_sqrt(s2);
				assert(r1 == i && r2 == i);
			}
			//printf("\n");
		}
	}
	#undef NUM_MAX
	#undef SQRT_MAX
	
	// if cb returns true, factorize interrupts
	void factorize(num_type n) const {
		assert(n > 0);
		if (n == 1) return;
		
		if (!(n & 1)) {
			exp_type exp = 0;
			do {
				n >>= 1;
				++exp;
			} while (!(n & 1));
			if (cb(2, exp)) return;
		}
		
		num_type n_sqrt = round_sqrt(n);
		num_type p = 3;
		size_t idx = 2;
		
		while (p <= n_sqrt) {
			if (!(n % p)) {
				exp_type exp = 0;
				do {
					n /= p;
					++exp;
				} while (!(n % p));
				if (cb(p, exp)) return;
				n_sqrt = round_sqrt(n);
			}
			if (idx < primes_array.count) {
				p = primes_array.primes[idx++];
			} else {
				p += 2;
			}
		}
		if (n != 1) cb(n, 1);
	}
	
	static size_t fill_primes(num_type primes[], size_t primes_size, num_type max_num) {
		if (primes_size == 0 || max_num < 2) return 0;
		primes[0] = 2;
		size_t primes_count = 1;
		num_type n = 3;
		num_type n_sqrt = 1;
		while (primes_count < primes_size && n <= max_num) {
			bool is_prime = true;
			for (size_t i=0; primes[i]<=n_sqrt; ++i) {
				if (!(n % primes[i])) {
					is_prime = false;
					break;
				}
			}
			assert(primes_count < primes_size);
			if (is_prime) primes[primes_count++] = n;
			n += 2;
			n_sqrt = round_sqrt(n);
		}
		return primes_count;
	}
};

template< typename NUM_TYPE>
class PrimeChecker {
public:
	typedef NUM_TYPE num_type;
private:
	typedef Factorizer<num_type> factorizer_type;
public:
	typedef typename factorizer_type::primes_array_type primes_array_type;
	
private:
	factorizer_type factorizer;
	num_type m_n;
	bool result;
	
	PrimeChecker() = delete;
	PrimeChecker(const PrimeChecker &b) = delete;
	PrimeChecker& operator=(const PrimeChecker &b) = delete;
	
public:
	PrimeChecker(primes_array_type b_primes_array) :
		factorizer(
			b_primes_array,
			std::bind(
				&PrimeChecker::factorize_cb,
				this,
				std::placeholders::_1,
				std::placeholders::_2
			)
		) {}
	
private:
	bool factorize_cb(typename factorizer_type::num_type prime, typename factorizer_type::exp_type exp) {
			(void)exp;
			result = (prime == m_n);
			return true;
	}
	
public:
	bool is_prime(num_type n) {
		assert(n > 1);
		m_n = n;
		factorizer.factorize(n);
		return result;
	}
};

template< typename NUM_TYPE>
class DivisorsCounter {
public:
	typedef NUM_TYPE num_type;
private:
	typedef Factorizer<num_type> factorizer_type;
public:
	typedef typename factorizer_type::primes_array_type primes_array_type;
	
private:
	factorizer_type factorizer;
	num_type count;
	
	DivisorsCounter() = delete;
	DivisorsCounter(const DivisorsCounter &b) = delete;
	DivisorsCounter& operator=(const DivisorsCounter &b) = delete;
	
public:
	DivisorsCounter(primes_array_type b_primes_array) :
		factorizer(
			b_primes_array,
			std::bind(
				&DivisorsCounter::factorize_cb,
				this,
				std::placeholders::_1,
				std::placeholders::_2
			)
		) {}
	
private:
	bool factorize_cb(typename factorizer_type::num_type prime, typename factorizer_type::exp_type exp) {
			(void)prime;
			count *= exp + 1;
			return false;
	}
	
public:
	num_type divisors_count(num_type n) {
		assert(n > 0);
		if (n == 1) return 1;
		count = 1;
		factorizer.factorize(n);
		return count;
	}
};

template <typename NUM_TYPE>
class SumOfTwoSquaresChecker {
public:
	typedef NUM_TYPE num_type;
private:
	typedef Factorizer<num_type> factorizer_type;
public:
	typedef typename factorizer_type::primes_array_type primes_array_type;
	
private:
	factorizer_type factorizer;
	bool result;
	
	SumOfTwoSquaresChecker() = delete;
	SumOfTwoSquaresChecker(const SumOfTwoSquaresChecker &b) = delete;
	SumOfTwoSquaresChecker& operator=(const SumOfTwoSquaresChecker &b) = delete;
	
public:
	SumOfTwoSquaresChecker(primes_array_type b_primes_array) :
		factorizer(
			b_primes_array,
			std::bind(
				&SumOfTwoSquaresChecker::factorize_cb,
				this,
				std::placeholders::_1,
				std::placeholders::_2
			)
		) {}
	
private:
	// Theorem:
	//     A number n is a sum of two squares if and only if all prime factors of n
	//     of the form 4m+3 have even exponent in the prime factorization of n.
	bool factorize_cb(typename factorizer_type::num_type prime, typename factorizer_type::exp_type exp) {
		if ((prime & 3) == 3 && (exp & 1)) {
			result = false;
			return true;
		} else {
			return false;
		}
	}
	
public:
	// if n is a square then function returns true (it means one summand is 0)
	bool is_sum_of_two_squares(num_type n) {
		if (n <= 1) return true;
		result = true;
		factorizer.factorize(n);
		return result;
	}
};

#endif/*FACTORIZE_H*/

