#ifndef POW_MOD_H
#define POW_MOD_H

#include <assert.h>

template <typename NUM_TYPE, NUM_TYPE NUM_TYPE_MAX_MASK, typename OPERATION_TYPE>
class PowMod {
public:
	typedef NUM_TYPE num_type;
	typedef OPERATION_TYPE operation_type;

	// on 0^0 returns 1
	static num_type pow_mod(num_type mod, num_type base, num_type exp) {
		assert(mod > 1);
		//assert(base > 0 || exp > 0);
		operation_type base_op = base, mod_op = mod, result = 1;
		num_type mask = NUM_TYPE_MAX_MASK;
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

};

#endif/*POW_MOD_H*/

