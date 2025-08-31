# CFX

CFX is a C library for performing a host of arithmetic, combinatorial computations. Its name stands for a number of things: C 

# cfx

**cfx** is a small C library for fast and exact combinatorics using prime-factor exponents, among other arithmetic things

It provides:
- Factorization type (`cfx_factorization_t`) that stores integers as prime–exponent maps
- Utilities for factorials, binomial coefficients, gcd/lcm, etc. without overflow
- A minimal big-integer type (`cfx_big_t`) to materialize results
- prime generation (sieve of Eratosthenes)

## Build

```
cmake -B build
cmake --build build
```

## Install

```
cmake --install build
``` 