# cfx

cfx is a C library for performing arithmetic computations with arbitrary precision integers. It's a work in progress for exploring algorithms in number theory, random number generators, cryptography, and algebra.

cfx used to mean one thing, something around 'Factorization into prime eXponents in C', but then the c started standing for 'cryptographic' and then the fx meant something else... for now it's just cfx.

## Requirements:
[cmake](https://cmake.org/download/), a c compiler, (optionally) [TestU01](https://simul.iro.umontreal.ca/testu01/tu01.html), (optionally, for benchmarks) [Google Benchmark](https://github.com/google/benchmark), (optionally, for benchmark comparisons) [OpenSSL](https://openssl-library.org/source/).

## Configure
Simplest default example : `cmake -S . -B build` 

example: choose your compiler, enable testu01, build benchmarks, release build, with 4096 primes in the static list: `CC=/usr/local/bin/clang cmake -B build -S . -DCFX_ENABLE_TESTU01=ON -DCFX_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release -DCFX_PRIMES_COUNT=4096`

example: force 32 bit limbs: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFX_BUILD_UTILS=ON -DCFX_BUILD_TESTS=ON -DCFX_FORCE_LIMB_32=ON`

example: build for ARM Cortex-M4: `cmake -B build-m4 -DCFX_TARGET=arm_cortex_m4 -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-arm-none-eabi-gcc.cmake -DCFX_MEMORY_MODE=dynamic`


## Compile

`cmake --build build -j` or `cd build && make` or `make VERBOSE=1`

list all cache variables: `cmake -L build`

### Windows (MSVC)

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
cmake -B build -S . -DCFX_ENABLE_TESTU01=ON
cmake --build build
```

If TestU01 is installed in a non-standard location, specify the root: `-DTESTU01_ROOT=/path/to/TestU01`

Example commands:

```bash
# See available RNGs and options
./build/test/stats/test_testu01 --help

# Run SmallCrush on ChaCha20
./build/test/stats/test_testu01 --rng=cfx_chacha20 --smallcrush
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
cmake -S . -B build -DCFX_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Benchmark binaries are in `build/benchmark/`. Use `--benchmark_filter=<regex>` to run specific benchmarks, `--benchmark_list_tests` to list them.

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


## Cross-Platform / Embedded Targets

cfx includes backends for ARM microcontrollers, with Docker & QEMU for testing.

### ARM Cortex-M4


```bash
# Build and test with Docker (recommended - works on macOS/Linux/Windows)
docker build -t cfx-cortex-m4 docker/cortex-m4/
docker run --rm -v $(pwd):/cfx cfx-cortex-m4

# Interactive debugging
docker run --rm -it -v $(pwd):/cfx cfx-cortex-m4 shell

# Inside the container, run individual tests with QEMU:
qemu-system-arm -M lm3s6965evb -cpu cortex-m4 -nographic -semihosting -kernel build-cortex-m4/test/unit/test_sha256
```

### ARM NEON (ARMv7/AArch64)

```bash
docker build -t cfx-arm-neon docker/arm-neon/
docker run --rm -v $(pwd):/cfx cfx-arm-neon
```

### Available Targets

Set with `-DCFX_TARGET=<target>`. See `cmake/cfx_target.cmake` for the full list. Targets include `portable`, `x86_64_avx2`, `x86_64_bmi2`, `arm_cortex_m4`, `arm_neon`, and `aarch64_neon`.

## License
The cfx library is dual-licensed.

You may choose to use it under either of the following licenses:

**LGPL v3 (or later), the GNU Lesser General Public License.**

You can use cfx in proprietary or open-source projects.

If you modify cfx itself and distribute it, you must publish your changes under the LGPL.

Dynamic linking is always fine; static linking is permitted if you provide a way for users to relink with a modified cfx.

**GPL v2 (or later), the GNU General Public License.**

You can use cfx in GPL-licensed software, including GPLv2-only projects.

If you distribute a project that includes cfx under the GPL path, the entire project must be licensed under the GPL as well.

In short:
- If you’re writing open source, pick GPL v2+ or GPL v3+.
- If you’re writing closed source, pick LGPL v3+.

See the LICENSE* texts for the full legal text

Thank you
