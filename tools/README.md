# Portable rendering tool

`render.cpp` exposes the 32 registered v1.8.6 renderers for scripts and future
visualization work. It calls the preserved engine's catalogue, state, scene,
pixel renderer, worker pool, presets, BMP writer and view serialization directly.
It adds no mathematical model or replacement colouring. The native Windows
application remains available through the root `run.cmd`.

Build from the repository root using `make`, or directly with an existing C++17
compiler:

```sh
mkdir -p build
g++ -std=c++17 -O2 -pthread tools/render.cpp -o build/fractal-render
./build/fractal-render --list
./build/fractal-render --renderer extra_30 --width 1200 --height 900 --output kernel.bmp
./build/fractal-render --renderer extra_31 --preset 4 --output feedback.bmp --save-view feedback.view
./build/fractal-render --view feedback.view --output feedback-large.bmp --width 1600 --height 1200
```

On Windows, use the same source and compiler flags with MinGW-w64, choose an
`.exe` output, and invoke the resulting executable from CMD. This portable
wrapper does not include Windows or D3D11 headers.

Run `--help` for every option. Dimensions are limited to 8192 per axis and
16,777,216 pixels total; iteration overrides accept 32–1024. CPU work can be
expensive at large sizes. `--samples 2` performs four samples per pixel.
The output is a 32-bit BMP using the engine's existing writer. Current supported
targets are little-endian x86/x64; the inherited writer stores native words.

The default state matches the native catalogue: RH entries use the first
**Currents** preset; additional entries use their engine defaults. `--preset`
selects 1–4 in that Currents bank for RH entries and the native preset array for
additional entries. It is applied before iteration and palette overrides.
Native `.view` files retain camera coordinates and model parameters; only
views belonging to one of the current 32 registered entries are accepted.
Save dimensions, sample count and compiler settings separately when exact
image reproduction matters. Floating-point image hashes may vary by compiler
or hardware. This tool renders CPU double-precision visualizations; its output
does not establish GPU parity or prove RH/QQRH.

For optional adapter regression checks with Python 3 installed:

```sh
python3 tools/test_render.py build/fractal-render
```

These checks render all 32 entries, compare worker counts, round-trip native
view files, exercise presets and reject invalid CLI inputs. Python is not
required to build or use the rendering executable.
