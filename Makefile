CXX ?= c++
CPPFLAGS += -Isrc
CXXFLAGS ?= -O2
# The numerical guards require IEEE finite-value semantics; do not add fast-math.
CXXFLAGS += -std=c++17 -pthread
HEADERS := $(wildcard src/*.hpp src/*.h)

.PHONY: all test smoke clean
all: build/fractal-render build/fractal-self-test

build:
	mkdir -p build

build/fractal-render: tools/render.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) -o $@

build/fractal-self-test: tests/cpu_self_tests.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) -o $@

test: all
	./build/fractal-self-test
	./build/fractal-render --list

# Optional adapter tests use Python's standard library, not the renderer.
smoke: build/fractal-render
	python3 tools/test_render.py build/fractal-render

clean:
	rm -rf build
