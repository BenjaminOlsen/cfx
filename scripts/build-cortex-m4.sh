#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later
#
# Build and test cfx for ARM Cortex-M4 under QEMU system emulation.
#
# This script runs INSIDE the cfx-cortex-m4 docker container, which provides
# arm-none-eabi-gcc and qemu-system-arm. From the host:
#
#   docker build -t cfx-cortex-m4 docker/cortex-m4/
#   docker run --rm -v "$(pwd):/cfx" cfx-cortex-m4 build
#   docker run --rm -v "$(pwd):/cfx" cfx-cortex-m4 test
#   docker run --rm -it -v "$(pwd):/cfx" cfx-cortex-m4 shell
#
# Usage (in container):
#   ./scripts/build-cortex-m4.sh [build|test|test-single <name>|clean|shell]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build-cortex-m4"
TOOLCHAIN="${PROJECT_ROOT}/cmake/toolchain-arm-none-eabi-gcc.cmake"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

check_build_prerequisites() {
    if ! command -v arm-none-eabi-gcc &> /dev/null; then
        error "arm-none-eabi-gcc not found.
Run this script inside the cfx-cortex-m4 container:
  docker run --rm -v \"\$(pwd):/cfx\" cfx-cortex-m4 build"
    fi
    info "Found: $(arm-none-eabi-gcc --version | head -1)"
}

check_test_prerequisites() {
    if ! command -v qemu-system-arm &> /dev/null; then
        error "qemu-system-arm not found (tests run under QEMU system emulation).
Run this script inside the cfx-cortex-m4 container:
  docker run --rm -v \"\$(pwd):/cfx\" cfx-cortex-m4 test"
    fi
    info "Found: $(qemu-system-arm --version | head -1)"
}

run_qemu() {
    timeout 30 qemu-system-arm \
        -M lm3s6965evb \
        -cpu cortex-m4 \
        -nographic \
        -semihosting \
        -kernel "$1" 2>&1
}

do_build() {
    check_build_prerequisites

    info "Configuring Cortex-M4 build..."

    # A cache configured from a different source path (host vs container
    # mount) poisons the build dir; wipe it rather than fail cryptically.
    if [ -f "${BUILD_DIR}/CMakeCache.txt" ] && \
       ! grep -q "CMAKE_HOME_DIRECTORY:INTERNAL=${PROJECT_ROOT}$" "${BUILD_DIR}/CMakeCache.txt"; then
        warn "Build dir was configured from a different path, wiping ${BUILD_DIR}"
        rm -rf "${BUILD_DIR}"
    fi

    # Clean stale test binaries from prior runs
    rm -f "${BUILD_DIR}"/test/unit/test_* 2>/dev/null || true
    cmake -B "${BUILD_DIR}" -S "${PROJECT_ROOT}" \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCFX_TARGET=arm_cortex_m4 \
        -DCMAKE_TOOLCHAIN_FILE="${TOOLCHAIN}" \
        -DCFX_MEMORY_MODE=dynamic \
        -DCFX_BUILD_TESTS=ON \
        -DCFX_BUILD_UTILS=OFF \
        -DCFX_PRIMES_COUNT=1000

    info "Building..."
    cmake --build "${BUILD_DIR}" -j"$(nproc)"

    info "Build complete: ${BUILD_DIR}"
    file "${BUILD_DIR}/libcfx.a"
    arm-none-eabi-size "${BUILD_DIR}/libcfx.a" || true
}

do_test() {
    check_test_prerequisites
    do_build

    info "Running tests under QEMU system emulation..."

    local total=0 passed=0 failed=0

    for test_bin in "${BUILD_DIR}"/test/unit/test_*; do
        if [ ! -x "$test_bin" ]; then
            continue
        fi

        total=$((total + 1))
        printf "  %-40s" "$(basename "$test_bin")"

        QEMU_OUT=$(run_qemu "$test_bin") && true
        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ]; then
            passed=$((passed + 1))
            echo -e "${GREEN}PASS${NC}"
        else
            failed=$((failed + 1))
            echo -e "${RED}FAIL${NC} (exit: $EXIT_CODE)"
            if [ -n "$QEMU_OUT" ]; then
                echo "$QEMU_OUT" | sed 's/^/    /'
            fi
        fi
    done

    echo ""
    echo "  Results: ${passed}/${total} passed"

    if [ $failed -gt 0 ]; then
        error "Some tests failed!"
    fi
    info "All Cortex-M4 tests passed!"
}

do_test_single() {
    check_test_prerequisites

    local test_name="$1"
    if [ -z "${test_name}" ]; then
        error "Usage: $0 test-single <test_name>"
    fi

    TEST_BIN="${BUILD_DIR}/test/unit/${test_name}"
    if [ ! -f "${TEST_BIN}" ]; then
        error "Test not found: ${TEST_BIN}"
    fi

    info "Running ${test_name} under QEMU..."
    run_qemu "${TEST_BIN}"
}

do_clean() {
    info "Cleaning build directory..."
    rm -rf "${BUILD_DIR}"
    info "Clean complete"
}

do_shell() {
    info "Build with: ./scripts/build-cortex-m4.sh build"
    info "Test with:  ./scripts/build-cortex-m4.sh test"
    exec /bin/bash
}

usage() {
    echo "Usage: $0 [command]"
    echo ""
    echo "Commands:"
    echo "  build        Configure and build for Cortex-M4"
    echo "  test         Build and run all tests under QEMU (default)"
    echo "  test-single  Run a single test: $0 test-single test_big_mul"
    echo "  clean        Remove build directory"
    echo "  shell        Drop into an interactive shell (for docker)"
    echo "  help         Show this help"
}

case "${1:-test}" in
    build)
        do_build
        ;;
    test)
        do_test
        ;;
    test-single)
        do_test_single "$2"
        ;;
    clean)
        do_clean
        ;;
    shell|bash)
        do_shell
        ;;
    help|--help|-h)
        usage
        ;;
    *)
        error "Unknown command: $1"
        ;;
esac
