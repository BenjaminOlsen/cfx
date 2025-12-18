//! find_prime - Generate an N-bit prime number
//!
//! This is a Rust port of the C `find_prime` example, simplified to single-threaded.
//! It demonstrates:
//! - Using the cfx Big type for arbitrary precision integers
//! - Command-line argument parsing in Rust
//! - The prime generation algorithm (trial division + Miller-Rabin)
//!
//! Run with: cargo run --example find_prime -- <bits> [--safe] [--top2] [--verbose]
//!
//! Example: cargo run --example find_prime --release -- 512

use cfx::Big;
use rand::Rng;
use std::env;
use std::process;

/// First 1024 small primes for trial division.
/// In C this comes from a generated file; here we generate them at compile time.
const SMALL_PRIME_CNT: usize = 1024;

/// Generate the first N primes using a simple sieve.
fn generate_small_primes() -> Vec<u32> {
    let mut primes = Vec::with_capacity(SMALL_PRIME_CNT);
    let mut candidate = 2u32;

    while primes.len() < SMALL_PRIME_CNT {
        let is_prime = primes.iter().take_while(|&&p| p * p <= candidate).all(|&p| candidate % p != 0);
        if is_prime {
            primes.push(candidate);
        }
        candidate += 1;
    }
    primes
}

/// Flags for prime generation.
#[derive(Default, Clone, Copy)]
struct PrimeFlags {
    /// Generate a safe prime p where (p-1)/2 is also prime (Sophie Germain prime)
    safe: bool,
    /// Set the top two bits (FIPS-style requirement for RSA moduli)
    top2: bool,
}

/// Generate a random N-bit odd number.
///
/// In C:
/// ```c
/// static int rand_nbit_odd(cfx_big_t* out, size_t bits, cfx_rng_fn rng, void* ctx, int flags)
/// ```
///
/// Here we use Rust's rand crate instead of the C RNG.
fn rand_nbit_odd(bits: usize, flags: PrimeFlags, rng: &mut impl Rng) -> Option<Big> {
    if bits < 2 {
        return None;
    }

    let nbytes = (bits + 7) / 8;
    let mut buf = vec![0u8; nbytes];

    // Fill with random bytes
    rng.fill(&mut buf[..]);

    // Mask the top partial byte so value < 2^bits
    let excess = 8 * nbytes - bits;
    if excess > 0 {
        buf[0] &= 0xFF >> excess;
    }

    // Force top bit to ensure exactly 'bits' bits
    buf[0] |= 1 << (7 - excess);

    // Force second highest bit if requested (for FIPS compliance)
    if flags.top2 && bits >= 2 {
        let pos = 7usize.saturating_sub(excess).saturating_sub(1);
        buf[0] |= 1 << pos;
    }

    // Force odd
    buf[nbytes - 1] |= 1;

    // For safe primes, force congruent to 3 mod 4
    // (every safe prime > 3 is 3 mod 4, because if p is prime and p = 1 mod 4,
    // then 2p+1 = 3 mod 4... wait that's backwards. Let me think again:
    // If q is an odd prime, then 2q+1 is either 3 mod 4 (if q = 1 mod 2, which is always true)
    // So safe prime p = 2q+1 means p mod 4 = 2*1 + 1 = 3 mod 4 when q = 1 mod 2.
    // Actually: q odd => q = 1 or 3 mod 4 => 2q+1 = 3 or 7 = 3 mod 4.
    // So safe primes are always 3 mod 4.)
    if flags.safe {
        let r = buf[nbytes - 1] & 3;
        if r != 3 {
            buf[nbytes - 1] = buf[nbytes - 1].wrapping_add(3 - r);
        }
    }

    Big::from_bytes_be(&buf)
}

/// Initialize remainder table: rem[i] = n % primes[i]
///
/// In C:
/// ```c
/// static void remainders_init(uint32_t* rem, const cfx_big_t* n,
///                             const uint32_t* primes, const size_t prime_cnt)
/// ```
fn remainders_init(n: &Big, primes: &[u32]) -> Vec<u32> {
    primes.iter().map(|&p| n.mod_sm(p as u64) as u32).collect()
}

/// Update remainders after adding `step` to the candidate.
///
/// In C:
/// ```c
/// static inline void remainders_add_step(uint32_t* rem, uint32_t step,
///                                        const uint32_t* primes, const size_t prime_cnt)
/// ```
fn remainders_add_step(rem: &mut [u32], step: u32, primes: &[u32]) {
    for (r, &p) in rem.iter_mut().zip(primes.iter()) {
        let new_r = *r + step;
        *r = if new_r >= p { new_r - p } else { new_r };
    }
}

/// Check if candidate passes small-prime trial division.
/// Returns false if divisible by any small prime.
///
/// In C:
/// ```c
/// static int passes_small_trial(const uint32_t* rem, const uint32_t* primes, const size_t prime_cnt)
/// ```
fn passes_small_trial(rem: &[u32]) -> bool {
    rem.iter().all(|&r| r != 0)
}

/// Check if a prime is a "safe prime" (Sophie Germain prime).
/// A safe prime p has the property that (p-1)/2 is also prime.
///
/// In C:
/// ```c
/// static int is_safe_prime(const cfx_big_t* p)
/// ```
fn is_safe_prime(p: &Big) -> bool {
    // q = (p - 1) / 2
    let mut q = p.clone();
    q -= 1u64;
    q.shr_assign(1);
    q.is_prime()
}

/// Generate an N-bit prime.
///
/// This is the core algorithm from the C version:
/// 1. Generate a random N-bit odd number
/// 2. Compute remainders mod small primes for fast trial division
/// 3. Loop: skip candidates divisible by small primes, test others with Miller-Rabin
/// 4. If we overflow N bits, redraw a new random candidate
///
/// In C:
/// ```c
/// int cfx_big_gen_nbit_prime(cfx_big_t* p, size_t bits,
///                            cfx_rng_fn rng, void* rngctx,
///                            unsigned flags, volatile int *stop_flag,
///                            const int threadid, int verbose)
/// ```
fn gen_nbit_prime(
    bits: usize,
    flags: PrimeFlags,
    primes: &[u32],
    verbose: bool,
) -> Big {
    let mut rng = rand::thread_rng();

    loop {
        // Draw initial random N-bit odd candidate
        let mut n = rand_nbit_odd(bits, flags, &mut rng)
            .expect("Failed to generate random number");

        // Prepare small prime remainders
        let mut rem = remainders_init(&n, primes);

        // Step by 2 for regular primes, 4 for safe primes (to stay 3 mod 4)
        let step: u32 = if flags.safe { 4 } else { 2 };

        let mut reject_cnt = 0u32;

        loop {
            // Small trial division
            if passes_small_trial(&rem) {
                // Probable-prime test (Miller-Rabin via cfx_big_is_prime)
                if n.is_prime() {
                    if flags.safe {
                        if is_safe_prime(&n) {
                            return n;
                        } else if verbose {
                            reject_cnt += 1;
                            println!("REJECTED prime #{} for not being safe: 0x{}", reject_cnt, n.to_hex());
                        }
                    } else {
                        return n;
                    }
                }
            }

            // Next odd candidate: n += step
            n.add_sm(step as u64);
            remainders_add_step(&mut rem, step, primes);

            // Keep within N bits: if we overflowed, redraw
            if n.bit_len() > bits {
                if verbose {
                    println!("Overflow, redrawing...");
                }
                break; // break inner loop to redraw
            }
        }
    }
}

fn print_usage(prog: &str) {
    eprintln!("Usage: {} <bits> [--safe] [--top2] [--verbose|-v]", prog);
    eprintln!();
    eprintln!("Generate an N-bit prime number.");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --safe     Generate a safe (Sophie Germain) prime");
    eprintln!("  --top2     Ensure top 2 bits are set (FIPS-style)");
    eprintln!("  -v, --verbose  Verbose output");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  {} 512           Generate a 512-bit prime", prog);
    eprintln!("  {} 2048 --safe   Generate a 2048-bit safe prime", prog);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let prog = &args[0];

    if args.len() < 2 {
        print_usage(prog);
        process::exit(1);
    }

    let mut bits: Option<usize> = None;
    let mut flags = PrimeFlags::default();
    let mut verbose = false;

    for arg in &args[1..] {
        match arg.as_str() {
            "--safe" => flags.safe = true,
            "--top2" => flags.top2 = true,
            "-v" | "--verbose" => verbose = true,
            "-h" | "--help" => {
                print_usage(prog);
                process::exit(0);
            }
            s => {
                // Try to parse as bits
                match s.parse::<usize>() {
                    Ok(b) => bits = Some(b),
                    Err(_) => {
                        eprintln!("Unknown argument: {}", s);
                        print_usage(prog);
                        process::exit(1);
                    }
                }
            }
        }
    }

    let bits = match bits {
        Some(b) if b >= 2 => b,
        Some(_) => {
            eprintln!("Error: bits must be >= 2");
            process::exit(1);
        }
        None => {
            eprintln!("Error: must specify number of bits");
            print_usage(prog);
            process::exit(1);
        }
    };

    if verbose {
        println!("Generating {}-bit {}prime...",
            bits,
            if flags.safe { "safe " } else { "" }
        );
    }

    // Generate small primes for trial division
    let small_primes = generate_small_primes();

    // Find the prime
    let prime = gen_nbit_prime(bits, flags, &small_primes, verbose);

    // Print result
    println!("---------------");
    println!("0x{}", prime.to_hex());
    println!("bits = {}", prime.bit_len());
}
