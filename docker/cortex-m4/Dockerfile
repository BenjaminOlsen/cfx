# cfx ARM Cortex-M4 Build and Test Environment
#
# Build:   docker build -t cfx-cortex-m4 docker/cortex-m4/
# Run:     docker run --rm -v $(pwd):/cfx cfx-cortex-m4
# Debug:   docker run --rm -it -v $(pwd):/cfx cfx-cortex-m4 /bin/bash

FROM ubuntu:22.04

LABEL maintainer="cfx project"
LABEL description="ARM Cortex-M4 cross-compilation and QEMU testing environment"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install toolchain and QEMU
RUN apt-get update && apt-get install -y \
    # Cross-compiler for bare-metal ARM
    gcc-arm-none-eabi \
    libnewlib-arm-none-eabi \
    # QEMU system emulation (bare-metal + semihosting)
    qemu-system-arm \
    # Build tools
    cmake \
    ninja-build \
    make \
    # Utilities
    git \
    file \
    binutils-arm-none-eabi \
    && rm -rf /var/lib/apt/lists/*

# Working directory (will be mounted from host)
WORKDIR /cfx

# Entrypoint runs from the mounted volume so changes take effect
# without rebuilding the image
ENTRYPOINT ["/cfx/docker/cortex-m4/entrypoint.sh"]
CMD ["test"]
