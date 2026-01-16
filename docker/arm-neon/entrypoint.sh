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

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }
log_section() { echo -e "\n${BLUE}=== $1 ===${NC}\n"; }

# Build for ARMv7 with NEON
build_armv7() {
    log_section "Building for ARMv7 (32-bit, NEON)"

    BUILD_DIR="/cfx/build-armv7-neon"
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_neon \
        -DCMAKE_C_COMPILER=arm-linux-gnueabihf-gcc \
        -DCMAKE_CXX_COMPILER=arm-linux-gnueabihf-g++ \
        -DCMAKE_SYSTEM_NAME=Linux \
        -DCMAKE_SYSTEM_PROCESSOR=arm \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        ..

    cmake --build . -j$(nproc)

    log_info "ARMv7 build complete!"
    file "${BUILD_DIR}/libcfx.a"
}

# Build for AArch64 with NEON
build_aarch64() {
    log_section "Building for AArch64 (64-bit, NEON)"

    BUILD_DIR="/cfx/build-aarch64-neon"
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_neon \
        -DCMAKE_C_COMPILER=aarch64-linux-gnu-gcc \
        -DCMAKE_CXX_COMPILER=aarch64-linux-gnu-g++ \
        -DCMAKE_SYSTEM_NAME=Linux \
        -DCMAKE_SYSTEM_PROCESSOR=aarch64 \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        ..

    cmake --build . -j$(nproc)

    log_info "AArch64 build complete!"
    file "${BUILD_DIR}/libcfx.a"
}

# Run tests for a given build directory and QEMU command
run_tests() {
    local BUILD_DIR=$1
    local QEMU_CMD=$2
    local ARCH_NAME=$3

    log_section "Testing ${ARCH_NAME}"

    cd "${BUILD_DIR}"

    local TEST_COUNT=0
    local PASS_COUNT=0
    local FAIL_COUNT=0

    for test_bin in test/unit/test_*; do
        if [ -x "$test_bin" ]; then
            test_name=$(basename "$test_bin")
            TEST_COUNT=$((TEST_COUNT + 1))

            echo -n "  $test_name: "

            if $QEMU_CMD "./$test_bin" >/dev/null 2>&1; then
                PASS_COUNT=$((PASS_COUNT + 1))
                echo -e "${GREEN}PASS${NC}"
            else
                FAIL_COUNT=$((FAIL_COUNT + 1))
                echo -e "${RED}FAIL${NC}"
            fi
        fi
    done

    echo ""
    echo "  Results: $PASS_COUNT/$TEST_COUNT passed"

    if [ $FAIL_COUNT -gt 0 ]; then
        return 1
    fi
    return 0
}

test_armv7() {
    build_armv7
    run_tests "/cfx/build-armv7-neon" "qemu-arm -L /usr/arm-linux-gnueabihf" "ARMv7 NEON"
}

test_aarch64() {
    build_aarch64
    run_tests "/cfx/build-aarch64-neon" "qemu-aarch64 -L /usr/aarch64-linux-gnu" "AArch64 NEON"
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
    echo "  qemu-arm -L /usr/arm-linux-gnueabihf ./binary"
    echo ""
    echo "  # Run AArch64 binary:"
    echo "  qemu-aarch64 -L /usr/aarch64-linux-gnu ./binary"
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
