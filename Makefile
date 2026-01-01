# convenience Makefile - wraps common CMake commands

BUILD_DIR = build
BUILD_DIR_COV = build-cov

.PHONY: all release debug test install utils coverage clean help

all: release

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
	cmake -S . -B $(BUILD_DIR_COV)  -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCFX_COVERAGE=ON -DCFX_BUILD_UTILS=OFF
	cmake --build $(BUILD_DIR_COV) -j
	ctest --test-dir $(BUILD_DIR_COV) --output-on-failure -j
	rm -rf coverage
	mkdir -p coverage
	gcovr --root . --filter src/ --html --html-details -o coverage/index.html
	gcovr --root . --filter src/ --print-summary

clean:
	rm -rf $(BUILD_DIR) $(BUILD_DIR_COV)

help:
	@echo "Available targets:"
	@echo "  release   - Build release configuration (default)"
	@echo "  debug     - Build debug configuration"
	@echo "  test      - Build and run tests"
	@echo "  install   - Install to CMAKE_INSTALL_PREFIX"
	@echo "  coverage  - Build with coverage and generate report"
	@echo "  clean     - Remove build directories"
