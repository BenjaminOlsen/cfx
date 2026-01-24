# convenience Makefile - wraps CMake (C) and Cargo (Rust)

BUILD_DIR = build
BUILD_DIR_COV = build-cov
RUST_DIR = rust/cfx

.PHONY: all release debug test install utils coverage benchmark clean help
.PHONY: rust rust-release rust-test rust-clean all-test all-clean

all: release

# ============================================================================
# C Library (CMake)
# ============================================================================

release:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~
	cmake --build $(BUILD_DIR) -j --config Release

debug:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) -j --config Debug

test: release
	ctest --test-dir $(BUILD_DIR) -C Release --output-on-failure -j

install: release
	cmake --install $(BUILD_DIR)

utils:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCFX_BUILD_TESTS=OFF -DCFX_BUILD_UTILS=ON
	cmake --build $(BUILD_DIR) -j --config Release --target install

coverage:
	cmake -S . -B $(BUILD_DIR_COV) -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCFX_COVERAGE=ON -DCFX_BUILD_UTILS=OFF
	cmake --build $(BUILD_DIR_COV) -j
	ctest --test-dir $(BUILD_DIR_COV) --output-on-failure -j
	rm -rf coverage
	mkdir -p coverage
	gcovr --root . --filter src/ --html --html-details -o coverage/index.html
	gcovr --root . --filter src/ --print-summary

benchmark:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCFX_BUILD_BENCHMARKS=ON
	cmake --build $(BUILD_DIR) -j --config Release
	@echo ""
	@echo "Benchmarks built. Run with:"
	@echo "  ./$(BUILD_DIR)/benchmark/bench_chacha20"
	@echo "  ./$(BUILD_DIR)/benchmark/bench_poly1305"
	@echo "  ./$(BUILD_DIR)/benchmark/bench_ntt"
	@echo "  etc."

clean:
	rm -rf $(BUILD_DIR) $(BUILD_DIR_COV)

# ============================================================================
# Rust Bindings (Cargo)
# ============================================================================

rust:
	cd $(RUST_DIR) && cargo build

rust-release:
	cd $(RUST_DIR) && cargo build --release

rust-test:
	cd $(RUST_DIR) && cargo test

rust-clean:
	cd $(RUST_DIR) && cargo clean

# ============================================================================
# Combined
# ============================================================================

all-test: test rust-test

all-clean: clean rust-clean

# ============================================================================
# Help
# ============================================================================

help:
	@echo "C Library (CMake):"
	@echo "  release      Build release (default)"
	@echo "  debug        Build debug"
	@echo "  test         Run C tests"
	@echo "  install      Install to CMAKE_INSTALL_PREFIX"
	@echo "  utils        Build and install utilities only"
	@echo "  coverage     Build with coverage and generate report"
	@echo "  benchmark    Build benchmarks (requires Google Benchmark)"
	@echo "  clean        Remove build directories"
	@echo ""
	@echo "Rust Bindings (Cargo):"
	@echo "  rust         Build debug"
	@echo "  rust-release Build release"
	@echo "  rust-test    Run Rust tests"
	@echo "  rust-clean   Clean Rust build"
	@echo ""
	@echo "Combined:"
	@echo "  all-test     Run all tests (C + Rust)"
	@echo "  all-clean    Clean everything"
