# cfx ARM NEON Build and Test Environment
#
# Supports both ARMv7 (32-bit) and ARMv8 (64-bit) with NEON SIMD
#
# Build:   docker build -t cfx-arm-neon docker/arm-neon/
# Run:     docker run --rm -v $(pwd):/cfx cfx-arm-neon
# Debug:   docker run --rm -it -v $(pwd):/cfx cfx-arm-neon /bin/bash

FROM ubuntu:22.04

LABEL maintainer="cfx project"
LABEL description="ARM NEON cross-compilation and QEMU testing environment"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install toolchains and QEMU
RUN apt-get update && apt-get install -y \
    # Cross-compiler for ARM Linux (32-bit, hard float)
    gcc-arm-linux-gnueabihf \
    g++-arm-linux-gnueabihf \
    libc6-dev-armhf-cross \
    # Cross-compiler for AArch64 Linux (64-bit)
    gcc-aarch64-linux-gnu \
    g++-aarch64-linux-gnu \
    libc6-dev-arm64-cross \
    # QEMU user-mode emulators
    qemu-user \
    qemu-user-static \
    # Build tools
    cmake \
    ninja-build \
    make \
    # Utilities
    git \
    file \
    && rm -rf /var/lib/apt/lists/*

# Register QEMU with binfmt for transparent execution
# (allows running ARM binaries directly)
RUN update-binfmts --enable qemu-arm || true
RUN update-binfmts --enable qemu-aarch64 || true

# Working directory
WORKDIR /cfx

# Copy entrypoint script
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Default: run build and tests
ENTRYPOINT ["/entrypoint.sh"]
CMD ["test"]
