#include <assert.h>
#include <stddef.h>		// size_t
#include <stdint.h>
#include <math.h>
#include <functional>

template <typename NUM_TYPE>
class Factorizer {
public:
	typedef NUM_TYPE num_type;
	static_assert(sizeof(num_type) <= 32, "Too big num_type for exp_type");
	typedef uint_fast8_t exp_type;
	typedef std::function<bool(num_type prime, exp_type exp)> factorize_cb_type;
	
private:
	num_type *primes;
	size_t primes_count;
	factorize_cb_type cb;
	
public:
	Factorizer(num_type *b_primes, size_t b_primes_count, factorize_cb_type b_cb) :
		primes(b_primes), primes_count(b_primes_count), cb(b_cb) {}
	
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
		num_type p = (primes_count == 0 ? 2 : primes[0]);
		size_t idx = 1;
		num_type n_sqrt = round_sqrt(n);
		 while (p <= n_sqrt) {
			if (!(n % p)) {
				exp_type exp = 0;
				do {
					n /= p;
					++exp;
				} while (!(n % p));
				if (cb(p, exp)) break;
				n_sqrt = round_sqrt(n);
			}
			if (idx < primes_count) {
				p = primes[idx];
				++idx;
			} else {
				if (p == 2) p += 1;
				else p += 2;
			}
		}
		if (n != 1) cb(n, 1);
	}
};

