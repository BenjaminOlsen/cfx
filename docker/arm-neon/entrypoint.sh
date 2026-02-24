#!/bin/bash
# cfx ARM NEON Build and Test Script
#
# Supports both ARMv7 (32-bit) and ARMv8/AArch64 (64-bit)
#
# Usage:
#   ./entrypoint.sh build    # Build for both architectures
#   ./entrypoint.sh test     # Build and test (default)
#   ./entrypoint.sh armv7    # Build and test ARMv7 only
#   ./entrypoint.sh aarch64  # Build and test AArch64 only
#   ./entrypoint.sh shell    # Interactive shell

set -e

# Source common test runner
source /cfx/docker/common/test-runner.sh

# Build for ARMv7 with NEON (32-bit)
# Note: 32-bit ARM cannot use 128-bit types, so curve25519/ed25519 are excluded
build_armv7() {
    log_section "Building for ARMv7 (32-bit, NEON)"

    BUILD_DIR="/cfx/build-armv7-neon"
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_neon \
        -DCMAKE_TOOLCHAIN_FILE=/cfx/cmake/toolchain-arm-linux-gnueabihf.cmake \
        -DCFX_LIMB_BITS=32 \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        ..

    cmake --build . -j$(nproc)

    log_info "ARMv7 build complete!"
    file "${BUILD_DIR}/libcfx.a"
}

# Build for AArch64 with NEON (64-bit)
# NEON is always available on AArch64, no -mfpu flag needed
build_aarch64() {
    log_section "Building for AArch64 (64-bit, NEON)"

    BUILD_DIR="/cfx/build-aarch64-neon"
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=aarch64_neon \
        -DCMAKE_TOOLCHAIN_FILE=/cfx/cmake/toolchain-aarch64-linux-gnu.cmake \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        ..

    cmake --build . -j$(nproc)

    log_info "AArch64 build complete!"
    file "${BUILD_DIR}/libcfx.a"
}

test_armv7() {
    build_armv7
    run_all_tests "/cfx/build-armv7-neon" \
        "qemu-arm -cpu cortex-a15 -L /usr/arm-linux-gnueabihf" \
        "ARMv7 NEON"
}

test_aarch64() {
    build_aarch64
    run_all_tests "/cfx/build-aarch64-neon" \
        "qemu-aarch64 -cpu cortex-a72 -L /usr/aarch64-linux-gnu" \
        "AArch64 NEON"
}

do_build() {
    build_armv7
    build_aarch64
}

do_test() {
    FAILED=0

    test_armv7 || FAILED=1
    test_aarch64 || FAILED=1

    if [ $FAILED -eq 1 ]; then
        log_error "Some tests failed!"
        exit 1
    fi

    log_info "All tests passed!"
}

do_shell() {
    log_info "Interactive shell. Example commands:"
    echo "  # Build ARMv7:"
    echo "  arm-linux-gnueabihf-gcc -mcpu=cortex-a15 -mfpu=neon ..."
    echo ""
    echo "  # Build AArch64:"
    echo "  aarch64-linux-gnu-gcc ..."
    echo ""
    echo "  # Run ARMv7 binary:"
    echo "  qemu-arm -cpu cortex-a15 -L /usr/arm-linux-gnueabihf ./binary"
    echo ""
    echo "  # Run AArch64 binary:"
    echo "  qemu-aarch64 -cpu cortex-a72 -L /usr/aarch64-linux-gnu ./binary"
    exec /bin/bash
}

# Main
case "${1:-test}" in
    build)
        do_build
        ;;
    test)
        do_test
        ;;
    armv7)
        test_armv7
        ;;
    aarch64)
        test_aarch64
        ;;
    shell|bash)
        do_shell
        ;;
    *)
        log_error "Unknown command: $1"
        echo "Usage: $0 {build|test|armv7|aarch64|shell}"
        exit 1
        ;;
esac
