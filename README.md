# number-theory
Template classes with Number theory algorithms implementations.

**compile** all tests by `mkdir build; make`<br />
run `make clean` for clean build dir<br />
run `make clean_tests` for clean tests executables

### factorize
integer factorization by trial division

`factorize.h` - template classes<br />
`factorize_tests.cpp` - tests and usage examples, **compile** by `make factorize_tests`

##### `factorize.h` classes:
`PrimesArray` - holds array of primes<br />
`Factorizer` - integer factorization by trial division<br />
`PrimeChecker` - check whether number is a prime, uses Factorizer which uses trial division<br />
`DivisorsCounter` - calculate count of divisors of given number<br />
`SumOfTwoSquaresChecker` - check whether number is sum of two squares (including summand 0) by theorem about sum of two squares

### canonic_factors
canonical representation of integer

`canonic_factors.h` - template classes<br />
`canonic_factors_tests.cpp` - tests and usage examples, **compile** by `make canonic_factors_tests`

##### `canonic_factors.h` classes:
`CanonicFactorsTemplate` - surrounding template with typedefs<br />
`CanonicFactorsTemplate::PrimePow` - power of prime: stores prime and it's exponent<br />
`CanonicFactorsTemplate::CanonicFactorizer` - uses Factorizer<br />
`CanonicFactorsTemplate::CanonicFactors` - main class, canonical representation of integer

##### `CanonicFactors` methods and operators:
`CanonicFactors`, `=`, `assign` (empty, `PrimePow` or other object) - constructors and assign operators<br />
`CanonicFactors`, `=`, `assign` (basic integer) - constructors and assign operators which factorize given number using `CanonicFactorizer`<br />
`value` - compute value as product of powers of primes<br />
`*`, `*=` - multiplication with basic integer or other object<br />
`mul_pow`, `mul_pow_assign` - multiplication with `PrimePow`

##### `CanonicFactors` algorithms (static methods):
`mul_static` - multiplication of two objects<br />
`lcm` - calculate least common multiple of two objects<br />
`eulers_phi` - calculate Euler's phi function<br />
`carmichael` - calculate Carmichael function

### mul_mod
`mul_mod.h` - utility template class `MulMod` for multiplication by modulo

##### `MulMod` methods (all static):
`mul_mod` -  modular multiplication<br />
`square_mod` - modular squaring<br />
`pow_mod` - fast modular exponentiation by squaring

### mul_group_mod
multiplicative group modulo n

`mul_group_mod.h` - template class `MulGroupMod` for finding element order or primitive root<br />
`mul_group_mod_tests.cpp` - tests and usage examples, **compile** by `make mul_group_mod_tests`

##### `MulGroupMod` methods:
`MulGroupMod` - construct object from modulo n<br />
`element_order` - calculate order of element of multiplicative group modulo n<br />
`is_primitive_root` - check whether element is a primitive root modulo n

### primitive_roots
primitive root modulo n checker, maybe slightly more efficient then `MulGroupMod`

`primitive_roots.h` template class `PrimitiveRoots` for primitive root checking<br />
`primitive_roots_tests.cpp` - tests and usage examples, **compile** by `make primitive_roots_tests`

##### `PrimitiveRoots` methods:
`PrimitiveRoots` - construct object from modulo n<br />
`is_primitive_root` - check whether element is a primitive root modulo n

### square_root_mod
quadratic congruences modulo n

`square_root_mod.h` - template class `SquareRootMod` for solving quadratic congruences<br />
`square_root_mod_tests.cpp` - tests and usage examples, **compile** by `square_root_mod_tests`

##### `SquareRootMod` methods:
`SquareRootMod` - construct object from modulo n for storing and using already calculated data<br />
`legendre_symbol` - calculate Legendre symbol by Euler's criterion (not the most efficient)<br />
`least_nonresidue` - find Least quadratic non-residue modulo n<br />
`tonelli_shanks_algo` - Tonelli-Shanks algorithm implementation, optimization by storing and using already calculated data<br />
`square_root_mod` - wrapper for `tonelli_shanks_algo`

