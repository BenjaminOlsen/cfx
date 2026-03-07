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

    # Create build directory (clean stale test binaries from prior runs)
    mkdir -p "${BUILD_DIR}"
    rm -f "${BUILD_DIR}"/test/unit/test_* 2>/dev/null || true
    cd "${BUILD_DIR}"

    # Configure with tests enabled (semihosting support allows QEMU execution)
    log_info "Configuring with CMake..."
    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_cortex_m4 \
        -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN}" \
        -DCFX_MEMORY_MODE=dynamic \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        -DCFX_PRIMES_COUNT=1000 \
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
    do_build

    log_info "Running tests with QEMU system emulation..."

    local TEST_COUNT=0
    local PASS_COUNT=0
    local FAIL_COUNT=0

    for test_bin in ${BUILD_DIR}/test/unit/test_*; do
        if [ ! -x "$test_bin" ]; then
            continue
        fi

        test_name=$(basename "$test_bin")
        TEST_COUNT=$((TEST_COUNT + 1))

        printf "  %-40s" "$test_name"

        QEMU_OUT=$(timeout 30 qemu-system-arm \
            -M lm3s6965evb \
            -cpu cortex-m4 \
            -nographic \
            -semihosting \
            -kernel "$test_bin" 2>&1) && true
        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ]; then
            PASS_COUNT=$((PASS_COUNT + 1))
            echo -e "${GREEN}PASS${NC}"
        else
            FAIL_COUNT=$((FAIL_COUNT + 1))
            echo -e "${RED}FAIL${NC} (exit: $EXIT_CODE)"
            if [ -n "$QEMU_OUT" ]; then
                echo "$QEMU_OUT" | sed 's/^/    /'
            fi
        fi
    done

    echo ""
    echo "  Results: ${PASS_COUNT}/${TEST_COUNT} passed"

    if [ $FAIL_COUNT -gt 0 ]; then
        log_error "Some tests failed!"
        exit 1
    fi

    log_info "All Cortex-M4 tests passed!"
}

do_shell() {
    log_info "Dropping into interactive shell..."
    log_info "Build with: cmake -DCFX_TARGET=arm_cortex_m4 ..."
    log_info "Test with:  qemu-system-arm -M lm3s6965evb -cpu cortex-m4 -nographic -semihosting -kernel ./test_binary"
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
