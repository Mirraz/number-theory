#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "factorize.h"

void test_round_sqrt() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	fzr_type::test_round_sqrt();
}

struct MyPow {
	uint_fast64_t prime;
	uint_fast8_t exp;
};

struct MyFactors {
	uint_fast8_t pow_count;
	MyPow pows[16];
};

MyFactors my_factorize(uint_fast64_t n) {
	assert(n > 0);
	MyFactors factors;
	factors.pow_count = 0;
	if (n == 1) return factors;
	uint_fast64_t n_sqrt = round(sqrt(n));
	for (uint_fast64_t d=2; d<=n_sqrt; d+=2) {
		if (n % d == 0) {
			uint_fast8_t exp = 0;
			do {
				n /= d;
				++exp;
			} while (n % d == 0);
			assert(factors.pow_count < 16);
			factors.pows[factors.pow_count].prime = d;
			factors.pows[factors.pow_count].exp = exp;
			++factors.pow_count;
			n_sqrt = round(sqrt(n));
		}
		if (d == 2) --d;
	}
	if (n > 1) {
		assert(factors.pow_count == 0 || factors.pows[factors.pow_count-1].prime != n);
		factors.pows[factors.pow_count].prime = n;
		factors.pows[factors.pow_count].exp = 1;
		++factors.pow_count;
	}
	return factors;
}

void test_factorize() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	fzr_type::primes_array_type primes;
	
	MyFactors factors;
	fzr_type::factorize_cb_type cb = [&factors] (fzr_type::num_type prime, fzr_type::exp_type exp) -> bool {
		factors.pows[factors.pow_count].prime = prime;
		factors.pows[factors.pow_count].exp = exp;
		++factors.pow_count;
		return false;
	};
	fzr_type factorizer(primes, cb);
	
	for (fzr_type::num_type i=1; i<=1024*64+1; ++i) {
		factors.pow_count = 0;
		factorizer.factorize(i);
		MyFactors my_factors = my_factorize(i);
		assert(factors.pow_count == my_factors.pow_count);
		for (uint_fast8_t j=0; j<factors.pow_count; ++j) {
			assert(factors.pows[j].prime == my_factors.pows[j].prime);
			assert(factors.pows[j].exp == my_factors.pows[j].exp);
		}
	}
}

bool my_is_sum_of_two_squares(uint_fast64_t n) {
	if (n <= 1) return true;
	uint_fast64_t n_sqrt = floor(sqrt((double)n));
	while (!(n & 1)) n >>= 1;
	for (uint_fast64_t p=3; p<=n_sqrt; ++p) {
		if (n % p) continue;
		if ((p & 3) == 3) {
			uint_fast8_t p_pow = 0;
			do {
				n /= p;
				++p_pow;
			} while (!(n % p));
			if (p_pow & 1) return false;
		} else {
			do {n /= p;} while (!(n % p));
		}
		n_sqrt = floor(sqrt((double)n));
	}
	if (n != 1 && (n & 3) == 3) return false;
	return true;
}

void test_sum_of_two_squares() {
	typedef SumOfTwoSquaresChecker<uint_fast64_t> checker_type;
	checker_type::primes_array_type primes;
	checker_type checker(primes);
	for (checker_type::num_type i=0; i<1024*4+1; ++i) {
		bool my_res = my_is_sum_of_two_squares(i);
		bool res = checker.is_sum_of_two_squares(i);
		assert(my_res == res);
	}
}

size_t my_fill_primes(uint_fast64_t primes[], size_t primes_size, uint_fast64_t max_num) {
	uint_fast64_t n = 2;
	uint_fast64_t count = 0;
	while (true) {
		uint_fast64_t n_sqrt = round(sqrt(n));
		uint_fast64_t d;
		for (d=2; d<=n_sqrt; d+=2) {
			if (n % d == 0) break;
			if (d == 2) --d;
		}
		if (d > n_sqrt) {
			// n is prime
			assert(count < primes_size);
			primes[count++] = n;
			if (count == primes_size) break;
		}
		if (n >= max_num) break;
		if (n == 2) ++n;
		else n += 2;
	}
	return count;
}

void test_fill_primes() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	fzr_type::num_type primes[1024];
	fzr_type::primes_array_type primes_array(primes, 0);
	uint_fast64_t myprimes[1024];
	size_t count, mycount;
	
	count = fzr_type::fill_primes(primes, 1024, UINT64_MAX);
	assert(count == 1024);
	mycount = my_fill_primes(myprimes, 1024, UINT64_MAX);
	assert(mycount == 1024);
	for (size_t i=0; i<1024; ++i) assert(primes_array.primes[i] == myprimes[i]);
	
	for (size_t i=0; i<1024; ++i) primes_array.primes[i] = myprimes[i] = 0;
	const uint_fast64_t p_max = 8147;
	count = fzr_type::fill_primes(primes, 1024, p_max);
	mycount = my_fill_primes(myprimes, 1024, p_max);
	assert(count == mycount);
	for (size_t i=0; i<count; ++i) assert(primes_array.primes[i] == myprimes[i]);
}

void test_factorize_with_primes_array() {
	typedef Factorizer<uint_fast64_t> fzr_type;
	
	fzr_type::num_type primes[65536];
	size_t primes_count = fzr_type::fill_primes(primes, 65536, UINT64_MAX);
	fzr_type::primes_array_type primes_array(primes, primes_count);
	
	MyFactors factors;
	fzr_type::factorize_cb_type cb = [&factors] (fzr_type::num_type prime, fzr_type::exp_type exp) -> bool {
		factors.pows[factors.pow_count].prime = prime;
		factors.pows[factors.pow_count].exp = exp;
		++factors.pow_count;
		return false;
	};
	fzr_type factorizer(primes_array, cb);
	
	for (fzr_type::num_type i=1; i<=1024*64+1; ++i) {
		factors.pow_count = 0;
		factorizer.factorize(i);
		MyFactors my_factors = my_factorize(i);
		assert(factors.pow_count == my_factors.pow_count);
		for (uint_fast8_t j=0; j<factors.pow_count; ++j) {
			assert(factors.pows[j].prime == my_factors.pows[j].prime);
			assert(factors.pows[j].exp == my_factors.pows[j].exp);
		}
	}
}

void test_prime_checke() {
	typedef PrimeChecker<uint_fast64_t> prime_checker_type;
	// pi(2^16) = 6542
	prime_checker_type::num_type primes[6542];
	size_t primes_count = prime_checker_type::primes_array_type::fill_primes(primes, 6542, 65536);
	assert(primes_count == 6542);
	prime_checker_type::primes_array_type primes_array(primes, primes_count);
	prime_checker_type prime_checker(primes_array);
	
	size_t idx = 0;
	for (prime_checker_type::num_type i=2; i<65536; ++i) {
		bool my_is_prime = (idx < primes_count && i == primes[idx]);
		if (my_is_prime) ++idx;
		bool is_prime = prime_checker.is_prime(i);
		assert(is_prime == my_is_prime);
	}
}

void tests_suite() {
	test_round_sqrt();
	test_factorize();
	test_sum_of_two_squares();
	test_fill_primes();
	test_factorize_with_primes_array();
	test_prime_checke();
}

int main() {
	tests_suite();
	return 0;
}

