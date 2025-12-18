# cfx - Rust Bindings

Safe Rust bindings for the cfx arbitrary-precision integer library.

## Usage

Add to your `Cargo.toml`:
```toml
[dependencies]
cfx = { path = "path/to/cfx/rust/cfx" }
```

### Basic Example
```rust
use cfx::Big;

fn main() {
    // Create from integers
    let a = Big::from(123u64);
    let b = Big::from(456u64);

    // Arithmetic operations
    let sum = &a + &b;
    let product = &a * &b;
    let quotient = &product / &a;

    println!("{} + {} = {}", a, b, sum);
    println!("{} * {} = {}", a, b, product);

    // Parse large numbers from strings
    let big = Big::from_dec("123456789012345678901234567890").unwrap();
    println!("Big number: {}", big);

    // Modular exponentiation
    let base = Big::from(2u64);
    let exp = Big::from(256u64);
    let modulus = Big::from_dec("1000000007").unwrap();
    let result = base.mod_exp(&exp, &modulus);
    println!("2^256 mod 1000000007 = {}", result);

    // Primality testing
    let prime = Big::from(104729u64);
    println!("{} is prime: {}", prime, prime.is_prime());
}
```

## Building & Testing

```bash
cd rust/cfx

# Build
cargo build

# Run tests
cargo test

# Build in release mode
cargo build --release
```

## API Overview

### Creating Big integers
- `Big::new()` - zero
- `Big::from(u64)` - from integer
- `Big::from_dec("123")` - from decimal string
- `Big::from_hex("abc")` - from hex string

### Arithmetic
- `+`, `-`, `*`, `/`, `%` operators (works with `&Big` references)
- `+=`, `-=`, `*=`, `/=` compound assignment
- `.pow(exp: u64)` - exponentiation
- `.mod_exp(&exp, &modulus)` - modular exponentiation
- `.square()` - in-place squaring

### Predicates
- `.is_zero()`, `.is_one()`
- `.is_even()`, `.is_odd()`
- `.is_prime()` - Miller-Rabin primality test

### Comparison
- Standard comparison operators: `==`, `<`, `>`, `<=`, `>=`
- Implements `Ord` for sorting

### Output
- `Display` trait for decimal output
- `.to_hex()` for hexadecimal

## Platform Notes

The C library is compiled automatically via `build.rs`. On Windows (MSVC), 32-bit limbs are used; on Unix with GCC/Clang, 64-bit limbs with 128-bit accumulators are used for better performance.
