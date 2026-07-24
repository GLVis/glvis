# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is GLVis

GLVis is an OpenGL visualization tool and library for finite element meshes and functions, built on MFEM. It runs as a local GUI app, a socket server (default port 19916) receiving data from remote simulations, or in headless/batch mode. It also compiles to WebAssembly via Emscripten for browser use.

## Build Systems

GLVis must use the same build system as MFEM.

### GNU Make

MFEM must be built at `../mfem` (or set `MFEM_DIR`).

```bash
make            # Build glvis binary
make opt        # Optimized build (GLVIS_DEBUG=NO)
make debug      # Debug build (GLVIS_DEBUG=YES)
make -j4        # Parallel build
make style      # Format source with astyle (required before PRs)
make status     # Show current config
make install [PREFIX=<dir>]
make clean
make distclean
```

Key variables: `MFEM_DIR`, `GLVIS_DEBUG`, `GLVIS_USE_LIBPNG` (default YES), `GLVIS_USE_EGL` (default NO), `PREFIX`.

### CMake

In-source builds are prohibited. Requires CMake >= 3.10.

```bash
mkdir build && cd build
cmake -G Ninja \
    -D CMAKE_BUILD_TYPE=Release \
    -D mfem_DIR=/path/to/mfem/build \
    -D GLVIS_USE_LIBPNG=ON \
    -D ENABLE_TESTS=TRUE \
    ..
cmake --build . --parallel 4
cmake --install .
```

Key variables: `mfem_DIR`, `CMAKE_BUILD_TYPE`, `GLVIS_USE_EGL` (ON for headless Linux), `ENABLE_TESTS`, `GLVIS_BASELINE_SYS`.

## Testing

Tests are screenshot-based visual regression tests (SSIM comparison), only available via CMake. The `tests/data` directory is a git submodule from https://github.com/GLVis/data.git.

```bash
# Install Python deps first
pip install -r tests/requirements.txt

# Run all tests
xvfb-run -a ctest --verbose   # Linux
ctest --verbose                # macOS

# Rebaseline (create new local baselines)
cmake --build . --target rebaseline

# Run a single test manually
python3 tests/glvis_driver.py \
    -e ./glvis \
    -s tests/data/streams/<test_name>.saved \
    -b tests/data/baselines/<GLVIS_BASELINE_SYS> \
    -a "-lw 1 -mslw 1 -nohidpi"
```

Pass threshold is SSIM >= 0.999. Test names: `ex1`, `ex2`, `ex3`, `ex5`, `ex6-dofs-numbering`, `ex9`, `ex22-scalar`, `ex22-vector`, `ex27`, `capacitor`, `distance`, `klein-bottle`, `laghos`, `mesh-explorer`, `mfem-logo`, `minimal-surface`, `mobius-strip`, `navier`, `quad`, `quadrature-1D`, `quadrature-lor`, `remhos`, `shaper`, `shifted`, `snake`.

## Code Style

```bash
make style   # Runs astyle 3.1 with ../mfem/config/mfem.astylerc
```

This must pass before any PR is merged. Excludes `lib/gl2ps.c` and `lib/gl2ps.h` (external).

## Architecture

### Key Source Directories

- **`glvis.cpp`** — main executable entry point; parses args, starts server or file reader
- **`lib/`** — core library built as `libglvis.a`
  - **`gl/`** — OpenGL abstraction: `renderer_core` (modern GL core profile), `renderer_ff` (legacy fixed-function), shader management, GLSL sources in `gl/shaders/`
  - **`sdl/`** — SDL2 windowing/event platform layer with platform-specific files (`sdl_mac.mm`, `sdl_x11.cpp`, `sdl_windows.cpp`)
  - **`egl/`** — EGL headless rendering (Linux only)
  - **`openglvis.cpp/hpp`** — `VisualizationScene` base class and camera
  - **`vsdata.cpp/hpp`** — `VisualizationSceneScalarData`: base for all data visualization
  - **`vssolution.cpp/hpp`** / **`vssolution3d.cpp/hpp`** — 2D/3D scalar field visualization
  - **`vsvector.cpp/hpp`** / **`vsvector3d.cpp/hpp`** — 2D/3D vector field visualization
  - **`aux_vis.cpp/hpp`** — main visualization thread initialization API
  - **`stream_reader.cpp/hpp`** — socket server mode (receives `.saved` stream data)
  - **`file_reader.cpp/hpp`** — file-based mesh/solution input
  - **`script_controller.cpp/hpp`** — batch/scripting command processing
  - **`palettes.cpp/hpp`** — color palette management
  - **`font.cpp/hpp`** — FreeType font rendering
  - **`gltf.cpp/hpp`** — glTF 3D export
  - **`gl2ps.c/h`** — PDF/EPS/SVG export (external library, do not style-check)
  - **`aux_js.cpp`** — Emscripten/WebAssembly bindings (web builds only)

### Data Flow

GLVis receives MFEM `Mesh` + `GridFunction` data either from a file or over a TCP socket (port 19916). The data is parsed in `stream_reader` or `file_reader`, stored in `data_state`, then dispatched to the appropriate `VisualizationScene` subclass (`vssolution`, `vssolution3d`, `vsvector`, `vsvector3d`) based on mesh dimension and field type. The scene renders via the GL renderer (`renderer_core` or `renderer_ff`), mediated by the SDL or EGL platform layer.

## Dependencies

**Required:** MFEM, OpenGL, GLEW >=2.1.0, GLM >=0.9.9.8, SDL2, FreeType 2, Fontconfig, `xxd` (build tool, not needed on Windows).

**Optional:** libpng (default ON), libtiff (default OFF), EGL (headless Linux), GnuTLS (secure sockets via MFEM).

Install on macOS: `brew install fontconfig freetype sdl2 glew glm libpng`

Install on Ubuntu: `apt-get install libfontconfig1-dev libfreetype-dev libsdl2-dev libglew-dev libglm-dev libpng-dev`

## Contribution Guidelines

- C++17 required (as of v4.5)
- Feature branches use `-dev` suffix; PRs target `master`
- PR checklist: builds cleanly, passes `make style`, updates `CHANGELOG`, updates `INSTALL` if dependencies changed
- New external library code must be `#ifdef`-guarded and off by default
- New classes/methods need Doxygen-style comments
