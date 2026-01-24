# cfx

cfx is a C library for performing a host of arithmetic computations with arbitrary precision integers. It's a work in progress for exploring algorithms in number theory, random number generators, cryptography, and algebra.

cfx used to mean one thing, something around 'Factorization into prime eXponents in C', but then the c started standing for 'cryptographic' and then the fx meant something else... for now it's just cfx.

## Requirements:
[cmake](https://cmake.org/download/), a c compiler, (optionally) [TestU01](https://simul.iro.umontreal.ca/testu01/tu01.html), (optionally, for benchmarks) [Google Benchmark](https://github.com/google/benchmark), (optionally, for benchmark comparisons) [OpenSSL](https://openssl-library.org/source/).

## Configure
Simplest default example : `cmake -S . -B build` 

choose your compiler, enable testu01, build benchmarks, release build, with 4096 primes in the static list: `CC=/usr/local/bin/clang cmake -B build -S . -DCFX_ENABLE_TESTU01=ON -DTESTU01_ROOT=$HOME/libs/TestU01 -DCFX_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release -DCFX_PRIMES_COUNT=4096`

force 32 bit limbs: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFX_BUILD_UTILS=ON -DCFX_BUILD_TESTS=ON -DCFX_FORCE_LIMB_32=ON`

build for armv7m: `cmake -B build-armv7m -S . -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-arm-none-eabi-gcc.cmake -DCFX_ARCH=armv7m -DCMAKE_BUILD_TYPE=Release -DCFX_PRIMES_COUNT=128 -DCFX_BUILD_UTILS=OFF -DCFX_BUILD_TESTS=ON`

build for ARM Cortex-M4 (optimized): `cmake -B build-m4 -DCFX_TARGET=arm_cortex_m4 -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-arm-none-eabi-gcc.cmake -DCFX_MEMORY_MODE=static`

### Security Options

**Paranoid Ed25519 verification:**

```bash
cmake -S . -B build -DCFX_ED25519_PARANOID=ON
```

This enables extra checks beyond RFC 8032 requirements: Small subgroup rejection; cofactored verification

If you accept Ed25519 public keys from untrusted sources. Without these checks, an attacker can provide a malicious "public key" (one of 8 special torsion points) and forge signatures for arbitrary messages.

**Cost:** ~6 extra point doublings per verification (~1% overhead; standard verification already does ~510 doublings).
## Compile

`cmake --build build -j` or `cd build && make` or `make VERBOSE=1`

other helpful commands: 

`cmake -L build `

### Windows (MSVC)

with 

```powershell
cmake -S . -B build
cmake --build build --config Release
```

Or for Debug builds: `cmake --build build --config Debug`

## Tests

The tests are divided into two categories, unit test and statistical tests:

### Unit tests

`ctest --test-dir build`

or for verbose output: `ctest --test-dir build -V`. Individual tests can be run `./build/test/unit/<testname>`

### Statistical tests with TestU01

To run statistical tests on the RNGs and hash functions in cfx using TestU01, install it following the instructions at https://github.com/umontreal-simul/TestU01-2009, then configure with:

```bash
cmake -B build -S . -DCFX_ENABLE_TESTU01=ON -DTESTU01_ROOT=/path/to/TestU01
cmake --build build
```

Example commands:

```bash
# See available RNGs and options
./build/test/stats/test_testu01 --help

# Run SmallCrush on ChaCha20 (fast, ~10 seconds)
./build/test/stats/test_testu01 --rng=cfx_chacha20 --smallcrush

# Run SmallCrush on SHA-256 in counter mode
./build/test/stats/test_testu01 --rng=cfx_sha256_ctr --smallcrush

# Run Crush on BLAKE2b (~30 minutes)
./build/test/stats/test_testu01 --rng=cfx_blake2b_ctr --crush

# Run BigCrush with custom seed (several hours)
./build/test/stats/test_testu01 --rng=cfx_xoshiro256starstar --bigcrush --seed=0x12345
```

### Code Coverage

Generate test coverage reports locally (requires GCC or Clang and gcovr):

```bash
pip install gcovr

# Configure with coverage enabled
cmake -S . -B build-cov -DCFX_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug

# Build and run tests
cmake --build build-cov -j
ctest --test-dir build-cov

# Generate HTML report
gcovr --root . --exclude 'test/*' --exclude 'utils/*' --exclude 'build-cov/*' --html-details coverage.html

# Open coverage.html in browser
```

Coverage is also run automatically on CI via GitHub Actions.

## Benchmarks

Enable benchmarks with `-DCFX_BUILD_BENCHMARKS=ON`. Requires [Google Benchmark](https://github.com/google/benchmark).

```bash
# Configure with benchmarks enabled
cmake -S . -B build -DCFX_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build --config Release

# Run a specific benchmark
./build/benchmark/Release/bench_ntt.exe      # Windows
./build/benchmark/bench_ntt                   # Linux/Mac
```

### Filtering Benchmarks

Use `--benchmark_filter` with a regex pattern to run specific benchmarks:

```bash
# Filter by argument size
./build/benchmark/bench_ntt --benchmark_filter=".*/(4096|8192)"

# Filter by name pattern
./build/benchmark/bench_poly1305 --benchmark_filter="BM_Poly1305.*"
```

### Other Useful Flags

```bash
--benchmark_format=json          # JSON output
--benchmark_out=results.json     # Save results to file
--benchmark_repetitions=5        # Run multiple times for statistics
--benchmark_list_tests           # List available benchmarks without running
```

## Examples

In utils, there are some interesting ways of using cfx, most have a `-h` or `--help` usage print.

### Installing as CLI Tools

The utils can be installed as command-line utilities for everyday use:

**Mac/Linux:**
```bash
# Build and install to ~/bin
cmake --build build && cmake --install build --prefix ~

# Add ~/bin to PATH (one-time, add to ~/.zshrc or ~/.bashrc)
export PATH="$HOME/bin:$PATH"
```

or to compile and install in one step:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=~
cmake --build build -j --target install
```  

**Windows:**
```powershell
# Build and install to %USERPROFILE%\bin
cmake --build build --config Release
cmake --install build --config Release --prefix %USERPROFILE%

# Add to PATH via System Properties > Environment Variables > User variables > Path
# Add: %USERPROFILE%\bin
```

or to compile and install in one step:

```powershell
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=%USERPROFILE%
cmake --build build --config Release -j --target install
```


After installation, run tools directly: `cfx_dc`, `cfx_factor 12345`, `cfx_primes 100`, etc.

To update after making changes, just re-run the build and install commands.

## Rust Bindings

cfx includes Rust bindings in `rust/cfx/`. The bindings provide a safe wrapper around the C library.

### Build and Test

```bash
cd rust/cfx
cargo test
```

The build script (`build.rs`) compiles the C library automatically using the `cc` crate.

## Cross-Platform / Embedded Targets

cfx includes optimized backends for ARM microcontrollers with Docker-based testing.

### ARM Cortex-M4

Optimized for embedded crypto using UMULL/UMLAL multiply-accumulate and barrel shifter rotations:

```bash
# Build and test with Docker (recommended - works on macOS/Linux/Windows)
docker build -t cfx-cortex-m4 docker/cortex-m4/
docker run --rm -v $(pwd):/cfx cfx-cortex-m4

# Interactive debugging
docker run --rm -it -v $(pwd):/cfx cfx-cortex-m4 shell
```

**Performance vs portable C:**
| Component | Speedup |
|-----------|---------|
| Big integer multiply | 3-4x |
| ChaCha20 | 1.5-2x |
| Poly1305 | 2-3x |

See `doc/CORTEX_M4.md` for design details and `doc/arm/` for optimization documentation.

### ARM NEON (ARMv7/AArch64)

```bash
docker build -t cfx-arm-neon docker/arm-neon/
docker run --rm -v $(pwd):/cfx cfx-arm-neon
```

### Available Targets

| Target | Backend | Key Optimizations |
|--------|---------|-------------------|
| `x86_64_avx2` | x86-64 with AVX2 | 8-lane SIMD for ChaCha20, BMI2 for bignum |
| `x86_64_bmi2` | x86-64 with BMI2 | MULX, ADCX, ADOX |
| `arm_cortex_m4` | ARM Cortex-M4 | UMULL/UMLAL, barrel shifter |
| `arm_neon` | ARMv7/v8 NEON | SIMD lanes |
| `portable` | Any | Standard C |

See `doc/EMULATION.md` for full cross-platform testing documentation.

### ChaCha20 Target-Specific Optimization

The ChaCha20 context size and implementation vary by target to eliminate per-call overhead:

| Target | ctx size | block8 implementation |
|--------|----------|----------------------|
| `x86_64_avx2` | 512 bytes | AVX2 8-lane SIMD, pre-broadcast state |
| others | 64 bytes | Scalar (portable) or target-optimized |

On AVX2, the context stores state in Structure-of-Arrays (SoA) layout with all values pre-broadcast to 8 lanes at init time. This eliminates ~20-25% overhead that would otherwise be spent broadcasting on every `block8` call.

```c
cfx_chacha20_ctx_t ctx;  // 512 bytes on AVX2, 64 bytes otherwise
cfx_chacha20_ctx_init(&ctx, key, nonce);
cfx_chacha20_block8(&ctx, counter, out);  // no broadcast overhead
```

The backend selection uses inheritance: if `x86_64_avx2/block4.c` doesn't exist, it falls back to `x86_64_bmi2`, then `x86_64`, then `portable`.

## License
The cfx library is dual-licensed.

You may choose to use it under either of the following licenses:

LGPL v3 (or later) — the GNU Lesser General Public License.

You can use cfx in proprietary or open-source projects.

If you modify cfx itself and distribute it, you must publish your changes under the LGPL.

Dynamic linking is always fine; static linking is permitted if you provide a way for users to relink with a modified cfx.

GPL v2 (or later) — the GNU General Public License.

You can use cfx in GPL-licensed software, including GPLv2-only projects.

If you distribute a project that includes cfx under the GPL path, the entire project must be licensed under the GPL as well.

In short:
- If you’re writing open source, pick GPL v2+ or GPL v3+.
- If you’re writing closed source, pick LGPL v3+.

See the LICENSE* texts for the full legal text
