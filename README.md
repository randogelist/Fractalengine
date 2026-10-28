# Fractalengine

Reusable FRAKTALFAKTUR visualization engine, based on
**FRAKTALFAKTUR_RH_LAB_GPU_v1_8_6**: native C++17 / Direct3D11, a CPU
double-precision renderer, 32 registered catalogue entries, saved views and
image export.

The original engine source and Windows launchers are preserved. A portable
command-line adapter makes the same CPU renderer available to scripts and
future visualization work.

## Generate an image

On Linux, use an existing C++17 compiler and GNU Make:

```sh
make -j2
./build/fractal-render --list
./build/fractal-render --renderer extra_30 --width 1200 --height 900 --output kernel.bmp
./build/fractal-render --renderer extra_31 --preset 4 --output feedback.bmp --save-view feedback.view
./build/fractal-render --view feedback.view --width 1600 --height 1200 --output feedback-large.bmp
```

No Python or graphics driver is required for these CPU renders. The output is
the engine's native BMP format. The [rendering tool guide](tools/README.md)
documents all options and direct compiler commands, including MinGW-w64.

## Open the native Windows application

From Windows CMD, with an existing x64 MinGW-w64 compiler:

```bat
call run.cmd
call run.cmd --renderer extra_30
```

The upstream launcher builds the application and runs CPU/GPU checks. See
[the original usage guide](README.txt) for zoom, pan, palettes, inspector,
image export and saved-view controls.

## Reuse the engine

| Need | Entry point |
| --- | --- |
| List available renderers | `build/fractal-render --list` or `catalogue::builtins()` |
| Reproduce a visualization | `--renderer`, `--preset`, `--save-view`, `--view` |
| Initialize parameters | `amp::defaults(family, generation)` in `src/amp_core.hpp` |
| Render a CPU pixel | `amp::scene(state)` and `amp::pixel(...)` |
| Save an image or parameters | `amp::writeBmp(...)`, `amp::writeView(...)` |
| Native GPU/UI | `src/amp_native.cpp`, `src/rh_lab_gpu.hpp` |

For future work, refer to this repository and a commit, then choose a renderer
ID and saved view. Keep new mathematics or sibling-engine merges explicit.

## Why this base

v1.8.6 preserves the accepted RH_LAB_GPU architecture and all 18 original
candidates, adds 14 integrated drivers, and includes the Hilbert point
inspector. FRONTIER_TRACE v1.9.0 is a separate 23-entry sibling. The older
Native Catalog has a broader collection of descriptors but less developed
engine infrastructure. See [source comparison and provenance](docs/ENGINE_SELECTION.md).

## Validation

```sh
make test
# Optional CLI adapter tests (Python 3 standard library):
make smoke
```

Fresh Linux/GCC 13.3 validation passed all 11 portable CPU suites. All 14
extension rendering hashes reproduce the archived v1.8.6 receipt, including
serial/four-worker agreement. All 97 entries in the original source manifest
verify. The [fresh test log](validation/portable_cpu_2026-10-28.log) records the
executed checks; original receipts remain alongside the preserved source.
The portable adapter passed 137 smoke invocations covering all 32 renderers,
parallel equality, saved-view round trips, presets and invalid inputs.

Native Windows compilation, GUI interaction, HLSL compilation and hardware
D3D11 execution remain unverified for this import. Numerical visualizations
do not certify RH or full QQRH positivity.
