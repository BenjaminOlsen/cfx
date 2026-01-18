#!/bin/bash
# cfx ARM Cortex-M4 Build and Test Script
#
# Usage:
#   ./entrypoint.sh build    # Build only
#   ./entrypoint.sh test     # Build and test (default)
#   ./entrypoint.sh shell    # Drop into bash shell

set -e

BUILD_DIR="/cfx/build-cortex-m4"
TOOLCHAIN="/cfx/cmake/toolchain-arm-none-eabi-gcc.cmake"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

do_build() {
    log_info "Building cfx for ARM Cortex-M4..."

    # Create build directory
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    # Configure
    # Note: Tests disabled for bare-metal (qemu-arm user mode can't run arm-none-eabi binaries)
    # Use docker/arm-neon for testing ARM-optimized code paths
    log_info "Configuring with CMake..."
    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_cortex_m4 \
        -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN}" \
        -DCFX_MEMORY_MODE=static \
        -DCFX_BUILD_TESTS=OFF \
        -DCFX_BUILD_UTILS=OFF \
        ..

    # Build
    log_info "Compiling..."
    cmake --build . -j$(nproc)

    log_info "Build complete!"

    # Show binary info
    log_info "Library info:"
    file "${BUILD_DIR}/libcfx.a"
    arm-none-eabi-size "${BUILD_DIR}/libcfx.a" || true
}

do_test() {
    log_warn "Tests not available for bare-metal Cortex-M4 target."
    log_info "The arm-none-eabi toolchain produces bare-metal binaries that"
    log_info "cannot run under qemu-arm user-mode emulation (which expects Linux binaries)."
    log_info ""
    log_info "To test ARM-optimized code paths, use the ARM NEON docker instead:"
    log_info "  docker build -t cfx-arm-neon docker/arm-neon/"
    log_info "  docker run --rm -v \$(pwd):/cfx cfx-arm-neon"
    log_info ""
    log_info "The Cortex-M4 docker validates that the library compiles correctly"
    log_info "for bare-metal embedded targets."

    do_build
    log_info "Library built successfully for Cortex-M4!"
}

do_shell() {
    log_info "Dropping into interactive shell..."
    log_info "Build with: cmake -DCFX_TARGET=arm_cortex_m4 ..."
    log_info "Test with:  qemu-arm -semihosting -cpu cortex-m4 ./test_binary"
    exec /bin/bash
}

# Main entry point
case "${1:-test}" in
    build)
        do_build
        ;;
    test)
        do_test
        ;;
    shell|bash)
        do_shell
        ;;
    *)
        log_error "Unknown command: $1"
        echo "Usage: $0 {build|test|shell}"
        exit 1
        ;;
esac
