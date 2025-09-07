# cfx

cfx is a C library for performing a host of arithmetic computations with arbitrary precision integers.

cfx used to mean one thing, something around 'Factorization into prime eXponents in C', but it's constantly changing so forget about what it means - its just cfx.

It has:
- Factorization type (`cfx_fac_t`) that stores integers as prime–exponent maps
- Utilities for factorials, binomial coefficients, gcd/lcm, etc. without overflow
- A big-integer type (`cfx_big_t`) to materialize results, with fast prime factorization
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