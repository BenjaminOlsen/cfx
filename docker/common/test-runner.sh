#!/bin/bash
# Common test runner for cfx cross-platform testing
#
# This script provides shared functionality for running tests
# across different emulated architectures.
#
# Usage:
#   source test-runner.sh
#   run_all_tests "/path/to/build" "qemu-arm" "ARMv7"

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Logging functions
log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }
log_section() { echo -e "\n${BLUE}=== $1 ===${NC}\n"; }

# Run all test binaries in a build directory
# Args:
#   $1 - Build directory path
#   $2 - QEMU command (e.g., "qemu-arm -L /usr/arm-linux-gnueabihf")
#   $3 - Architecture name for display
#   $4 - (optional) Additional QEMU flags
run_all_tests() {
    local BUILD_DIR="$1"
    local QEMU_CMD="$2"
    local ARCH_NAME="$3"
    local QEMU_EXTRA="${4:-}"

    log_section "Testing ${ARCH_NAME}"

    if [ ! -d "${BUILD_DIR}" ]; then
        log_error "Build directory not found: ${BUILD_DIR}"
        return 1
    fi

    cd "${BUILD_DIR}"

    local TEST_COUNT=0
    local PASS_COUNT=0
    local FAIL_COUNT=0
    local SKIP_COUNT=0

    # Find all test executables
    for test_bin in test/unit/test_*; do
        if [ ! -x "$test_bin" ]; then
            continue
        fi

        local test_name=$(basename "$test_bin")
        TEST_COUNT=$((TEST_COUNT + 1))

        printf "  %-40s" "$test_name"

        # Run with timeout to prevent hangs
        local OUTPUT
        local EXIT_CODE=0

        OUTPUT=$(timeout 60 $QEMU_CMD $QEMU_EXTRA "./$test_bin" 2>&1) || EXIT_CODE=$?

        case $EXIT_CODE in
            0)
                PASS_COUNT=$((PASS_COUNT + 1))
                echo -e "${GREEN}PASS${NC}"
                ;;
            124)
                SKIP_COUNT=$((SKIP_COUNT + 1))
                echo -e "${YELLOW}TIMEOUT${NC}"
                ;;
            *)
                FAIL_COUNT=$((FAIL_COUNT + 1))
                echo -e "${RED}FAIL${NC} (exit: $EXIT_CODE)"
                # Show first few lines of output for debugging
                echo "$OUTPUT" | head -5 | sed 's/^/    /'
                ;;
        esac
    done

    # Summary
    echo ""
    echo "  Summary for ${ARCH_NAME}:"
    echo "    Total:   $TEST_COUNT"
    echo -e "    Passed:  ${GREEN}$PASS_COUNT${NC}"
    echo -e "    Failed:  ${RED}$FAIL_COUNT${NC}"
    echo -e "    Skipped: ${YELLOW}$SKIP_COUNT${NC}"
    echo ""

    # Return failure if any tests failed
    [ $FAIL_COUNT -eq 0 ]
}

# Run a single test with verbose output
# Args:
#   $1 - Test binary path
#   $2 - QEMU command
run_single_test() {
    local TEST_BIN="$1"
    local QEMU_CMD="$2"

    log_info "Running: $TEST_BIN"
    echo "Command: $QEMU_CMD $TEST_BIN"
    echo "---"

    $QEMU_CMD "$TEST_BIN"
    local EXIT_CODE=$?

    echo "---"
    if [ $EXIT_CODE -eq 0 ]; then
        log_info "Test passed"
    else
        log_error "Test failed with exit code $EXIT_CODE"
    fi

    return $EXIT_CODE
}

# Check QEMU availability
# Args:
#   $1 - QEMU binary name (e.g., "qemu-arm")
check_qemu() {
    local QEMU="$1"

    if ! command -v "$QEMU" &>/dev/null; then
        log_error "QEMU not found: $QEMU"
        log_info "Install with: apt-get install qemu-user"
        return 1
    fi

    log_info "Found: $($QEMU --version | head -1)"
    return 0
}

# Check cross-compiler availability
# Args:
#   $1 - Compiler binary name (e.g., "arm-linux-gnueabihf-gcc")
check_compiler() {
    local CC="$1"

    if ! command -v "$CC" &>/dev/null; then
        log_error "Cross-compiler not found: $CC"
        return 1
    fi

    log_info "Found: $($CC --version | head -1)"
    return 0
}

# Generate JUnit XML test report
# Args:
#   $1 - Output file path
#   $2 - Test suite name
#   $3 - Total tests
#   $4 - Failures
#   $5 - Skipped
generate_junit_report() {
    local OUTPUT="$1"
    local SUITE="$2"
    local TOTAL="$3"
    local FAILURES="$4"
    local SKIPPED="$5"

    cat > "$OUTPUT" << EOF
<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="${SUITE}" tests="${TOTAL}" failures="${FAILURES}" skipped="${SKIPPED}">
</testsuite>
EOF

    log_info "JUnit report written to: $OUTPUT"
}
