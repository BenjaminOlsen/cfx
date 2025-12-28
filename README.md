# cfx

cfx is a C library for performing a host of arithmetic computations with arbitrary precision integers. It's a work in progress for exploring algorithms in number theory, random number generators, cryptography, and algebra.

cfx used to mean one thing, something around 'Factorization into prime eXponents in C', but then the c started standing for 'cryptographic' and then the fx meant something else... for now it's just cfx.

## Requirements:
[cmake](https://cmake.org/download/), a c compiler, (optionally) [TestU01](https://simul.iro.umontreal.ca/testu01/tu01.html), (optionally, for benchmarks) [Google Benchmark](https://github.com/google/benchmark), (optionally, for benchmark comparisons) [OpenSSL](https://openssl-library.org/source/).

## Configure
Simplest default example : `cmake -S . -B build` 

choose your compiler, enable testu01, build benchmarks, release build, with 4096 primes in the static list: `CC=/usr/local/bin/clang cmake -B build -S . -DCFX_ENABLE_TESTU01=ON -DTESTU01_ROOT=$HOME/libs/TestU01 -DCFX_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release -DCFX_PRIMES_COUNT=4096`

force 32 bit limbs: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFX_BUILD_EXAMPLES=ON -DCFX_BUILD_TESTS=ON -DCFX_FORCE_LIMB_32=ON`

build for armv7m: `cmake -B build-armv7m -S . -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-arm-none-eabi-gcc.cmake -DCFX_ARCH=armv7m -DCMAKE_BUILD_TYPE=Release -DCFX_PRIMES_COUNT=128 -DCFX_BUILD_EXAMPLES=OFF -DCFX_BUILD_TESTS=ON`

## Compile

`cmake --build build -j` or `cd build && make` or `make VERBOSE=1`

### Windows (MSVC)

```
cmake -S . -B build
cmake --build build --config Release
```

Or for Debug builds: `cmake --build build --config Debug`

Note: Benchmarks require GCC/Clang intrinsics and are not available on MSVC. Some POSIX-only examples are also skipped on Windows.

## Tests

The tests are divided into two categories, unit test and statistical tests:

### Unit tests

`ctest --test-dir build`

or for verbose output: `ctest --test-dir build -V`. Individual tests can be run `./build/test/unit/<testname>`

### Statistical tests with TestU01

to run statistical tests on the zoo of RNGs in cfx using TestU01, you have to install it following their instructions: https://github.com/blep/TestU01, and then pass its install location to cmake configure:

`cmake -B build -S . -DCFX_ENABLE_TESTU01=ON -DTESTU01_ROOT=/path/to/TestU01`

You can see the options for running SmallCrush, Crush, and BigCrush: 

`./build/test/stats/test_testu01 --help`

Note, some of the RNGs are toy examples (like using poly1305 as an RNG), but others pass BigCrush.

## Examples

In examples, there are some interesting ways of using cfx, most have a `-h` or `--help` usage print.

### Installing as CLI Tools

The examples can be installed as command-line utilities for everyday use:

**Mac/Linux:**
```bash
# Build and install to ~/bin
cmake --build build && cmake --install build --prefix ~

# Add ~/bin to PATH (one-time, add to ~/.zshrc or ~/.bashrc)
export PATH="$HOME/bin:$PATH"
```

**Windows:**
```powershell
# Build and install to %USERPROFILE%\bin
cmake --build build --config Release
cmake --install build --config Release --prefix %USERPROFILE%

# Add to PATH via System Properties > Environment Variables > User variables > Path
# Add: %USERPROFILE%\bin
```

After installation, run tools directly: `cfx_dc`, `cfx_factor 12345`, `cfx_primes 100`, etc.

To update after making changes, just re-run the build and install commands.

## Rust Bindings

cfx includes Rust bindings in `rust/cfx/`. The bindings provide a safe wrapper around the C library.

### Build and Test

```
cd rust/cfx
cargo test
```

The Rust bindings require the C library to be built first. The build script (`build.rs`) will compile the C library automatically using cmake.

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
