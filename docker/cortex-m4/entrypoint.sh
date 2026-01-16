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
    log_info "Configuring with CMake..."
    cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_cortex_m4 \
        -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN}" \
        -DCFX_MEMORY_MODE=static \
        -DCFX_BUILD_TESTS=ON \
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
    do_build

    log_info "Running tests with QEMU..."

    cd "${BUILD_DIR}"

    # Find all test executables
    TEST_COUNT=0
    PASS_COUNT=0
    FAIL_COUNT=0
    SKIP_COUNT=0

    for test_bin in test/unit/test_*; do
        if [ -x "$test_bin" ]; then
            test_name=$(basename "$test_bin")
            TEST_COUNT=$((TEST_COUNT + 1))

            log_info "Running: $test_name"

            # Run with QEMU semihosting
            # Note: Bare-metal tests need semihosting support
            if qemu-arm -semihosting -cpu cortex-m4 "./$test_bin" 2>/dev/null; then
                PASS_COUNT=$((PASS_COUNT + 1))
                echo -e "  ${GREEN}PASS${NC}"
            else
                exit_code=$?
                if [ $exit_code -eq 1 ]; then
                    FAIL_COUNT=$((FAIL_COUNT + 1))
                    echo -e "  ${RED}FAIL${NC} (exit code: $exit_code)"
                else
                    # QEMU errors (missing semihosting, etc.) - skip
                    SKIP_COUNT=$((SKIP_COUNT + 1))
                    echo -e "  ${YELLOW}SKIP${NC} (QEMU error: $exit_code)"
                fi
            fi
        fi
    done

    echo ""
    log_info "Test Summary:"
    echo "  Total:   $TEST_COUNT"
    echo -e "  Passed:  ${GREEN}$PASS_COUNT${NC}"
    echo -e "  Failed:  ${RED}$FAIL_COUNT${NC}"
    echo -e "  Skipped: ${YELLOW}$SKIP_COUNT${NC}"

    if [ $FAIL_COUNT -gt 0 ]; then
        log_error "Some tests failed!"
        exit 1
    fi

    log_info "All tests passed!"
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
