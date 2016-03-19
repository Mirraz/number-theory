#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include "mul_group_mod.h"

uint_fast64_t gcd(uint_fast64_t a, uint_fast64_t b) {
	if (a == 0) return b;
	return gcd(b % a, a);
}

void test_order() {
	typedef MulGroupMod<uint_fast32_t, 9, ((uint_fast32_t)1)<<31, uint_fast64_t> mgm_type;
	mgm_type::primes_array_type primes_array;
	mgm_type a(primes_array);
	mgm_type::num_type i, j, k;
	for (i=2; i<1024; ++i) {
		a.assign(i);
		for (j=0; j<i; ++j) {
			if (gcd(i,j) != 1) continue;
			printf("mod = %u element = %u ", (unsigned int)i, (unsigned int)j);
			k = a.element_order(j);
			printf("order = %u\n", (unsigned int)k);
		}
	}
}

// all primes less 65536
#define myprimes_size 6542
uint_fast16_t myprimes[myprimes_size];
uint_fast8_t max_roots[myprimes_size];

void fill_myprimes() {
	uint_fast16_t n = 2;
	uint_fast16_t count = 0;
	while (true) {
		uint_fast16_t n_sqrt = round(sqrt(n));
		uint_fast16_t d;
		for (d=2; d<=n_sqrt; d+=2) {
			if (n % d == 0) break;
			if (d == 2) --d;
		}
		if (d > n_sqrt) {
			// n is prime
			myprimes[count++] = n;
		}
		if (n == 65535) break;
		++n;
	}
	assert(count == myprimes_size);
}

// primorial(6) <= 65536 < primorial(7)
#define MAX_PRIME_DIVS 6

uint_fast8_t get_prime_divs(uint_fast16_t n, uint_fast16_t prime_divs[]) {
	uint_fast8_t count = 0;
	uint_fast16_t n_sqrt = round(sqrt(n));
	for (uint_fast8_t d_idx = 0; myprimes[d_idx] <= n_sqrt; ++d_idx) {
		uint_fast16_t d = myprimes[d_idx];
		if (n % d == 0) {
			assert(count < MAX_PRIME_DIVS);
			prime_divs[count++] = d;
			do {
				n /= d;
			} while (n % d == 0);
		}
		n_sqrt = round(sqrt(n));
	}
	if (n > 1) {
		assert(count < MAX_PRIME_DIVS);
		prime_divs[count++] = n;
	}
	return count;
}

uint_fast16_t pow_mod(uint_fast16_t mod, uint_fast16_t base, uint_fast16_t exp) {
	assert(mod > 1);
	assert(base > 0 || exp > 0);
	uint_fast32_t base_op = base, mod_op = mod, result = 1;
	uint_fast16_t mask = (uint_fast16_t)1 << 15;
	while (mask != 0 && !(exp & mask)) mask >>= 1;
	while (mask > 0) {
		result = (result * result) % mod_op;
		if (exp & mask) {
			result = (result * base_op) % mod_op;
		}
		mask >>= 1;
	}
	return result;
}

// max_roots[idx] = myprimes[idx] - (max primitive root modulo myprimes[idx])
void fill_max_roots() {
	max_roots[0] = 2 - 1; // for n = 2
	max_roots[1] = 3 - 2; // for n = 3

	// starting from 5
	for (uint_fast16_t n_idx = 2; n_idx<myprimes_size; ++n_idx) {
		uint_fast16_t n = myprimes[n_idx];
		uint_fast16_t n_1 = n - 1;
		// calculate prime divs of n-1
		uint_fast16_t prime_divs[MAX_PRIME_DIVS];
		uint_fast8_t prime_divs_count = get_prime_divs(n_1, prime_divs);
		uint_fast16_t exps[MAX_PRIME_DIVS];
		uint_fast8_t exp_count = prime_divs_count;
		for (uint_fast8_t i=0; i<prime_divs_count; ++i) {
			assert(n_1 > prime_divs[i]);
			assert(n_1 % prime_divs[i] == 0);
			exps[prime_divs_count-i-1] = n_1 / prime_divs[i];
		}
		
		uint_fast16_t root;
		for (root = n - 2; root > 1; --root) {
			uint_fast8_t i;
			for (i=0; i<exp_count; ++i) {
				uint_fast16_t pow = pow_mod(n, root, exps[i]);
				if (pow == 1) break;
			}
			if (i == exp_count) break;
		}
		assert(root > 1);
		assert(n - root < 256);
		max_roots[n_idx] = n - root;
	}
}

void test_max_primitive_root() {
	fill_myprimes();
	fill_max_roots();

	typedef MulGroupMod<uint_fast32_t, 9, ((uint_fast32_t)1)<<31, uint_fast64_t> mgm_type;
	// pi(2^16) = 6542
	mgm_type::num_type primes[6542];
	size_t primes_count = mgm_type::primes_array_type::fill_primes(
		primes,
		sizeof(primes) / sizeof(primes[0]),
		(mgm_type::num_type)UINT16_MAX + 1
	);
	assert(primes_count == sizeof(primes) / sizeof(primes[0]));
	mgm_type::primes_array_type primes_array(primes, primes_count);
	
	for (size_t idx=0; idx<sizeof(primes) / sizeof(primes[0]); ++idx) {
		mgm_type::num_type modulo = primes[idx];
		mgm_type mul_group_mod(primes_array, modulo);
		mgm_type::num_type max_root = 0;
		for (mgm_type::num_type i=modulo-1; i>=1; --i) {
			bool is_primitive_root = mul_group_mod.is_primitive_root(i);
			if (is_primitive_root) {
				assert(mul_group_mod.element_order(i) == modulo-1);
				max_root = i;
				break;
			}
		}
		assert(max_root != 0);
		
		assert(max_root == myprimes[idx] - max_roots[idx]);
	}
}

void tests_suite() {
	//test_order();
	test_max_primitive_root();
}

int main() {
	tests_suite();
	return 0;
}

