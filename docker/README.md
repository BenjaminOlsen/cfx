# cfx Cross-Platform Docker Testing

This directory contains Dockerfiles for testing cfx on emulated architectures.

## Quick Start

```bash
# From the cfx root directory:

# ARM Cortex-M4 (bare-metal, 32-bit, DSP)
docker build -t cfx-cortex-m4 docker/cortex-m4/
docker run --rm -v $(pwd):/cfx cfx-cortex-m4

# ARM NEON (Linux, 32-bit and 64-bit)
docker build -t cfx-arm-neon docker/arm-neon/
docker run --rm -v $(pwd):/cfx cfx-arm-neon
```

## Available Targets

| Target | Dockerfile | CPU | Key Features |
|--------|------------|-----|--------------|
| Cortex-M4 | `cortex-m4/` | ARM Cortex-M4 | UMULL/UMLAL, barrel shifter, bare-metal |
| ARM NEON | `arm-neon/` | Cortex-A15/A72 | NEON SIMD, Linux userspace |

## Usage

### Run Tests (Default)

```bash
docker run --rm -v $(pwd):/cfx cfx-cortex-m4
```

### Build Only

```bash
docker run --rm -v $(pwd):/cfx cfx-cortex-m4 build
```

### Interactive Shell

```bash
docker run --rm -it -v $(pwd):/cfx cfx-cortex-m4 shell
```

### ARM NEON - Specific Architecture

```bash
# ARMv7 only
docker run --rm -v $(pwd):/cfx cfx-arm-neon armv7

# AArch64 only
docker run --rm -v $(pwd):/cfx cfx-arm-neon aarch64
```

## CI Integration

Add to GitHub Actions:

```yaml
jobs:
  cortex-m4:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: docker build -t cfx-m4 docker/cortex-m4/
      - run: docker run --rm -v ${{ github.workspace }}:/cfx cfx-m4

  arm-neon:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: docker build -t cfx-neon docker/arm-neon/
      - run: docker run --rm -v ${{ github.workspace }}:/cfx cfx-neon
```

## What's Inside

Each Docker image contains:

- **Cross-compiler**: GCC for the target architecture
- **QEMU user-mode**: CPU emulation for running binaries
- **CMake + Ninja**: Build system
- **cfx toolchain file**: Automatic compiler selection

## How It Works

1. Docker container starts with cross-compiler and QEMU installed
2. Source is mounted from host at `/cfx`
3. CMake configures with target-specific toolchain
4. Build produces target binaries
5. QEMU runs test binaries under emulation
6. Results are printed to stdout

## Troubleshooting

### "Permission denied" on entrypoint.sh

```bash
chmod +x docker/*/entrypoint.sh
```

### Tests timeout or hang

Some bare-metal tests may not have semihosting support. Check that:
- QEMU was started with `-semihosting` flag
- Test binary has proper exit mechanism

### "Illegal instruction"

The target CPU in QEMU may not support required instructions. Check:
- Correct `-cpu` flag in QEMU invocation
- Binary was compiled for correct target

## Adding New Targets

1. Create `docker/new-target/Dockerfile`
2. Create `docker/new-target/entrypoint.sh`
3. Add toolchain file in `cmake/toolchain-*.cmake`
4. Update this README

See `doc/EMULATION.md` for detailed documentation.
