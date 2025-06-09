                           GLVis visualization tool

                     _/_/_/  _/      _/      _/  _/
                  _/        _/      _/      _/        _/_/_/
                 _/  _/_/  _/      _/      _/  _/  _/_/
                _/    _/  _/        _/  _/    _/      _/_/
                 _/_/_/  _/_/_/_/    _/      _/  _/_/_/

                             https://glvis.org

<a href="https://github.com/GLVis/glvis/releases/latest"><img alt="Release" src="https://img.shields.io/badge/release-v4.4-success.svg"></a>
<a href="https://github.com/GLVis/glvis/actions/workflows/builds.yml"><img alt="Build" src="https://github.com/GLVis/glvis/actions/workflows/builds.yml/badge.svg"></a>
<a href="https://github.com/glvis/glvis/blob/master/LICENSE"><img alt="License" src="https://img.shields.io/badge/License-BSD-success.svg"></a>
<a href="https://glvis.github.io/doxygen/html/index.html"><img alt="Doxygen" src="https://img.shields.io/badge/code-documented-success.svg"></a>
<a href="https://glvis.github.io/releases/glvis-macOS.dmg"><img alt="License" src="https://img.shields.io/badge/Download-Mac-success.svg"></a>
<a href="https://glvis.github.io/releases/glvis-Windows.zip"><img alt="License" src="https://img.shields.io/badge/Download-Windows-success.svg"></a>

[GLVis](https://glvis.org) is an OpenGL tool for visualization of finite element
meshes and functions. It is a multiplatform application that can be built on
Linux/Unix systems, including macOS, and under Windows. It can also be used in a
Jupyter notebook, or in a web browser, see https://glvis.org/live.

- For building instructions, see [INSTALL](INSTALL).

- The GLVis [key commands](#key-commands) and [mouse functions](#mouse-functions)
  are documented below.

- GLVis is distributed under the terms of the BSD-3 license. All new contributions
  must be made under this license. See [LICENSE](LICENSE) and [NOTICE](NOTICE) for details.

We welcome contributions and feedback from the community. Please see the file
[CONTRIBUTING.md](CONTRIBUTING.md) for additional details about our development
process.

When started without any options, glvis starts a server which waits for a socket
connections (on port `19916` by default) and visualizes any received data. This
way the results of simulations on a remote (parallel) machine can be visualized
on the local user desktop.

GLVis can also be used to visualize a mesh with or without a finite element
function (solution), as in

```
glvis -m cube.mesh3d
```

For parallel computations, GLVis supports input from several parallel socket
connections as well as the visualization of parallel meshes and grid functions
saved in separate files from the command line as in

```
glvis -np 4 -m mesh -g solution
```

When given parallel input, GLVis will stitch the results to show the global mesh
and solution. GLVis can also run a batch sequence of commands (GLVis scripts),
or display previously saved socket streams.

For a complete list of command line options, type

```
glvis -h
```

Depending on the data type, variety of manipulations can be performed by using
the mouse and by typing (case sensitive) keystrokes in the GLVis window. Below
is a partial list of the available functionality. Some of these keys can also be
provided as input, using the `-k` command-line option and the `keys` script
command.

For high-order meshes and/or solution data, GLVis performs element subdivision
to try to represent the data more accurately. However, for highly-varying data
or large meshes, the auto-selected subdivision factor (see the
[Auto-refinement](#auto-refinement) section below) may not be sufficient -- use
the keys <kbd>o</kbd> / <kbd>O</kbd>, described below, to manually adjust the
subdivision factor.

SPDX-License-Identifier: BSD-3-Clause <br>
LLNL Release Number: LLNL-CODE-443271 <br>
DOI: 10.11578/dc.20171025.1249


Mouse functions
===============

## Basic

- <kbd>Left</kbd> – Rotate the viewpoint
- <kbd>Right</kbd> – Zoom in (up) / Zoom out (down)
- <kbd>Middle</kbd> – Translate the viewpoint
- <kbd>Left</kbd> + <kbd>Shift</kbd> – Start spinning the viewpoint (according to the dragging vector)

## Advanced

- <kbd>Left</kbd> + <kbd>Alt</kbd> – Tilt
- <kbd>Left</kbd> + <kbd>Ctrl</kbd> – Spherical rotation
- <kbd>Left</kbd> + <kbd>Ctrl</kbd> + <kbd>Shift</kbd> – `z`-spinning
- <kbd>Right</kbd> + <kbd>Ctrl</kbd> – Object scaling (see also <kbd>[</kbd> and <kbd>]</kbd>)
- <kbd>Right</kbd> + <kbd>Shift</kbd>  – Change light source position (see <kbd>\\</kbd>)
- <kbd>Middle</kbd> + <kbd>Ctrl</kbd> – Object translation (moves the camera left/right/up/down)
- <kbd>Middle</kbd> + <kbd>Ctrl</kbd> + <kbd>Alt</kbd> – Moves the camera forward/backward (vertical mouse motion) and tilts the camera left/right (horizontal mouse motion)
- <kbd>Middle</kbd> + <kbd>Ctrl</kbd> + <kbd>Shift</kbd> – Object translation (turns the camera left/right/up/down)


Key commands
============

## Basic

- <kbd>h</kbd> – Print a short help message in the terminal
- <kbd>r</kbd> – Reset the plot to 3D view
- <kbd>R</kbd> – Cycle through the six 2D projections (camera looking above/below in `x`/`y`/`z` directions)
- <kbd>j</kbd> – Turn on/off perspective
- <kbd>s</kbd> – Turn on/off unit cube scaling
- <kbd>A</kbd> – Turn on/off the use of anti-aliasing/multi-sampling
- <kbd>L</kbd> – Turn on/off logarithmic scale
- <kbd>c</kbd> – Toggle the colorbar and caption display state
- <kbd>C</kbd> – Change the main plot caption
- <kbd>p</kbd> / <kbd>P</kbd> – Cycle forward/backwards through color palettes (lots of options, use <kbd>F6</kbd> for a menu)
- <kbd>t</kbd> – Cycle materials and lights (5 states)
- <kbd>i</kbd> – Toggle cutting plane (different options in 2D and 3D, see below)
- <kbd>o</kbd> / <kbd>O</kbd> – Control element subdivisions (different options in 2D and 3D, see below)
- <kbd>l</kbd> – Turn on/off the light
- <kbd>g</kbd> – Toggle background color (white/black)
- <kbd>a</kbd> – Toggle the bounding box *axes*. The options are:
  - none
  - bounding box with coordinates of the corners
  - bounding box without coordinates
  - red, green, blue colored main `x`, `y`, `z` axes + dashed axes
- <kbd>m</kbd> – Toggle the *mesh* state. The options are:
  - no mesh or level lines
  - draw the element edges (i.e. the mesh)
  - draw the level lines (use <kbd>F5 </kbd> to modify the level lines)
- <kbd>e</kbd> – Toggle the *elements* state (see below for vector functions). The options are:
  - show surface elements (corresponding to the function)
  - show no surface elements
  - (2D only) show element attributes
  - (2D only) show element `det(J)`
  - (2D only) show element `1/det(J)`
  - (2D only) show element `\kappa`
  - (2D only) show element `\kappa + 1/\kappa`
- <kbd>S</kbd> – Take an image snapshot or record a movie (in spinning mode). By
  default, the screenshots are taken in `png` format, using libpng. When GLVis
  is compiled with `libtiff` support (see [INSTALL](INSTALL)) then the
  screenshots are taken internally and saved in TIFF format (`.tif`
  extension). If both of these options are disabled during the build process,
  GLVis will use `SDL` to take screenshots in `bmp` format, which it will then
  convert to `png` if ImageMagick's `convert` tool is available.
- <kbd>G</kbd> – 3D scene export to [glTF format](https://www.khronos.org/gltf)
- <kbd>Ctrl</kbd> + <kbd>p</kbd> – Print to a PDF file using `gl2ps`. Other
  vector formats (SVG, EPS) are also possible, but keep in mind that the
  printing takes a while and the generated files are big.
- <kbd>Q</kbd> – Cycle between representations of the visualized *quadrature data*. The options are:
  - piece-wise constant refined (LOR)
  - L2 element dof collocation (interpolation)
  - L2 element projection (L2 projection)
- <kbd>q</kbd> – Exit

## Advanced

- <kbd>f</kbd> – Change the shading type (the way the elements and mesh are drawn). The options are:
  - one triangle / quad per element with a constant normal
  - one triangle / quad per element with normals averaged at the vertices
  - multiple triangles / quads per element, also allowing for the visualization of discontinuous
    functions and curvilinear elements (use <kbd>o</kbd> / <kbd>O</kbd> to control subdivisions)
- <kbd>\\</kbd> – Set light source position (see <kbd>Right</kbd> + <kbd>Shift</kbd>)
- <kbd>*</kbd> / <kbd>/</kbd> – Zoom in/out
- <kbd>+</kbd> / <kbd>-</kbd> – Stretch/compress in `z`-direction
- <kbd>[</kbd> / <kbd>]</kbd> – Shrink/enlarge the bounding box (relative to the colorbar)
- <kbd>(</kbd> / <kbd>)</kbd> – Shrink/enlarge the visualization window
- <kbd>.</kbd> – Start/stop `z`-spinning (speed/direction can be controlled with <kbd>0</kbd> / <kbd>Enter</kbd>)
- <kbd>←</kbd>, <kbd>→</kbd>, <kbd>↑</kbd>, <kbd>↓</kbd> – Manual rotation
- <kbd>1</kbd>, <kbd>2</kbd>, <kbd>3</kbd>, <kbd>4</kbd>, <kbd>5</kbd>, <kbd>6</kbd>, <kbd>7</kbd>, <kbd>8</kbd>, <kbd>9</kbd> – Manual rotation along coordinate axes
- <kbd>Alt</kbd> + <kbd>a</kbd> – Set axes number format
- <kbd>Alt</kbd> + <kbd>c</kbd> – Set colorbar number format
- <kbd>Ctrl</kbd> + <kbd>←</kbd>, <kbd>→</kbd>, <kbd>↑</kbd>, <kbd>↓</kbd> – Translate the viewpoint
- <kbd>Ctrl</kbd> + <kbd>o</kbd> – Toggle an element ordering curve
- <kbd>n</kbd> / <kbd>N</kbd> – Cycle through numberings. The options are:
  - none
  - show elements numbering
  - show edges numbering
  - show vertices numbering
  - show DOFs numbering
- <kbd>`</kbd> – Toggle a ruler, with initial origin at the center of the bounding box. The origin can be later moved with <kbd>~</kbd>. The options are:
  - none
  - coordinate axes lines
  - coordinate axes planes
- <kbd>~</kbd> - Enter new ruler origin
- <kbd>k</kbd> / <kbd>K</kbd> - Adjust the transparency level. The balance of
  transparency can be further adjusted with <kbd>,</kbd> and <kbd><</kbd>.
- <kbd>!</kbd> - Toggle the use of (1D) texture (smooth interpolation of colors). The options are:
  - use discrete texture, the number of colors used depends on the current palette
  - use smooth texture (interpolated from current palette)
- <kbd>F5</kbd> – Change the range and number of the level lines
- <kbd>F6</kbd> – Palette menu (negative repeat number flips the palette)
- <kbd>F7</kbd> – Change the minimum and maximum values
- <kbd>Shift</kbd> + <kbd>F7</kbd> – Set the bounding box from the terminal
- <kbd>Ctrl</kbd> + <kbd>a</kbd> – Cycle through the auto-scaling options:
  - `off`: do not change the bounding box and the value range
  - `on` (default): recompute both the bounding box and the value range
  - `value`: recompute only the value range
  - `mesh`: recompute only the bounding box

## 2D scalar data

- <kbd>i</kbd> – Toggle cutting (clipping) plane in 2D
- <kbd>y</kbd> / <kbd>Y</kbd> – Rotate cutting plane (`\theta`) in 2D
- <kbd>z</kbd> / <kbd>Z</kbd> – Translate cutting plane in 2D
- <kbd>o</kbd> / <kbd>O</kbd> – Control element subdivisions in 2D
  - there are two subdivision factors in this case: element (`s1`) and boundary (`s2`).
  - <kbd>O</kbd> cycles through the following "subdivision functions", (prints a message in the terminal when changed):
    - Increase element subdivision factor: `s1 += s2`
    - Decrease element subdivision factor: `s1 -= s2`
    - Increase boundary subdivision factor: `s2++`
    - Decrease boundary subdivision factor: `s2--`
  - <kbd>o</kbd> – performs the currently selected function
- <kbd>b</kbd> – Toggle the boundary in 2D scalar mode. The options are:
  - no boundary
  - black boundary
  - boundary colored with the boundary attribute
  - Use <kbd>Shift</kbd> + <kbd>F9</kbd> / <kbd>F10</kbd> to cycle through the boundary attributes.
- <kbd>F3</kbd> / <kbd>F4</kbd> – Shrink/Zoom each element towards its center, in order to better visualize the different element shapes
- <kbd>F8</kbd> – List of material subdomains to show
- <kbd>F9</kbd> / <kbd>F10</kbd> – Walk through material subdomains
- <kbd>F11</kbd> / <kbd>F12</kbd> – Shrink/Zoom material subdomains (to the centers of the attributes)

## 3D scalar data

- <kbd>i</kbd> – Toggle cutting (clipping) plane in 3D. The options are:
  - none
  - cut through the elements
  - show only elements behind the cutting plane
- <kbd>I</kbd> – Toggle the cutting plane algorithm used when the option *cut through the elements* is selected. The two algorithms are:
  - slower, more accurate algorithm for curved meshes (default)
  - faster algorithm suitable for meshes with planar faces
- <kbd>x</kbd> / <kbd>X</kbd> – Rotate cutting plane (`\phi`) in 3D
- <kbd>y</kbd> / <kbd>Y</kbd> – Rotate cutting plane (`\theta`) in 3D
- <kbd>z</kbd> / <kbd>Z</kbd> – Translate cutting plane in 3D
- <kbd>E</kbd> – Display/Hide the elements in the cutting plane
- <kbd>M</kbd> – Display/Hide the mesh in the cutting plane
- <kbd>o</kbd> / <kbd>O</kbd> – Refine/de-refine elements in 3D
- <kbd>u</kbd> / <kbd>U</kbd> – Move level surface value up/down
- <kbd>v</kbd> / <kbd>V</kbd> – Add/Delete a level surface value
- <kbd>w</kbd> / <kbd>W</kbd> – Move boundary elements up/down in direction of
  the normal (i.e. "plot" the boundary values in normal direction)
- <kbd>F3</kbd> / <kbd>F4</kbd> – Shrink/Zoom boundary elements (to the centers of the attributes)
- <kbd>F8</kbd> – List of boundary subdomains to show
- <kbd>F9</kbd> / <kbd>F10</kbd> – Walk through boundary subdomains
- <kbd>F11</kbd> / <kbd>F12</kbd> – Shrink/Zoom material subdomains (to the centers of the attributes)

## 2D vector data

- <kbd>v</kbd> – Toggle the *vector* state (uses vector subdivision factor, accept <kbd>u</kbd> / <kbd>U</kbd>). The options are:
  - do not show vectors
  - show vectors as displacement
  - show vector field; vectors are uniformly scaled; the color varies with the
    magnitude (or the current *vector-to-scalar function*, see keys <kbd>u</kbd> / <kbd>U</kbd>)
  - show vector field as above, but the vectors are scaled proportionally to their magnitude
  - show vector field; vectors are uniformly scaled; the color is gray; arrows are above the surface
  - show vector field as above, but the vectors are scaled proportionally to their magnitude
- <kbd>V</kbd> – Change the scaling of the vectors relative to the default
- <kbd>d</kbd> – Toggle the *displaced mesh* state: (see also keys <kbd>n</kbd> / <kbd>b</kbd>). The options are:
  - do not show displaced mesh
  - show displaced mesh
  - assuming displacement field show deformation using Cartesian lines
  - assuming displacement field show deformation using polar lines
- <kbd>n</kbd> – Increase the displacement amount in 10% steps, wraps around from 100% to 0%
- <kbd>b</kbd> – Decrease the displacement amount in 10% steps, wraps around from 0% to 100%
- <kbd>B</kbd> – Toggle the boundary in 2D vector mode
- <kbd>e</kbd> – Toggle the *elements* state (vector data version). The options are:
  - show surface elements corresponding to the current *vector-to-scalar function*
  - do not show surface elements
  - assuming a displacement field show `det(J)/det(J_d)`
  - assuming a displacement field show `det(J_d)/det(J)`
- <kbd>u</kbd> / <kbd>U</kbd> – Change the *vector-to-scalar function* and the vector subdivision factor
- <kbd>U</kbd> – Toggle the functionality of <kbd>u</kbd> (prints a message in the terminal when changed). The options are:
  - Increase the vector subdivision factor
  - Decrease the vector subdivision factor
  - Cycle through *vector-to-scalar functions* choices:
     - magnitude: `\sqrt{v_x^2+v_y^2}`
     - direction from `-\pi` to `\pi`: `atan2(v_y,v_x)`
     - `x`-component: `v_x`
     - `y`-component: `v_y`
     - divergence: `div(v)`
     - curl: `curl(v)` [skipped for H(div) elements]
     - anisotropy in `grad(v)`  [skipped for H(div) elements]

## 3D vector data

- <kbd>v</kbd> – Toggle the *vector* state. The options are:
   - do not show vectors
   - show vectors as displacement
   - show vector field; vectors are uniformly scaled; the color varies with the
     magnitude (or the current *vector-to-scalar function*, see key <kbd>F</kbd>)
   - show vector field as above, but the vectors are scaled proportionally to
     their magnitude
   - show the subset of the vector field with scalar function around a given
     value (adjusted with keys <kbd>u</kbd> / <kbd>U</kbd> and <kbd>w</kbd> / <kbd>W</kbd>)
   - show the vector field restricted to the boundary of the domain
- <kbd>V</kbd> – Cycle the *vector* state in the opposite direction of <kbd>v</kbd>
- <kbd>u</kbd> / <kbd>U</kbd> – Move the level field vectors (in the appropriate *vector* state)
- <kbd>w</kbd> / <kbd>W</kbd> – Add/Delete level field vector (in the appropriate *vector* state)
- <kbd>d</kbd> – Toggle the *displaced mesh* state (see also keys <kbd>n</kbd> / <kbd>b</kbd>). The options are:
  - do not show displaced mesh
  - show displaced mesh
- <kbd>n</kbd> – Increase the displacement amount in 10% steps, wraps around from 100% to 0%
- <kbd>b</kbd> – Decrease the displacement amount in 10% steps, wraps around from 0% to 100%
- <kbd>F</kbd> – Change the *vector-to-scalar function*. The options are:
  - magnitude: `\sqrt{v_x^2+v_y^2+v_z^2}`
  - `x`-component: `v_x`
  - `y`-component: `v_y`
  - `z`-component: `v_z`

## Auto-refinement

The GLVis auto-refinement algorithm selects a subdivision factor trying to
achieve an accurate representation of high-order meshes and solution data while
keeping the initial time to visualize the data reasonable. The algorithm can be
summarized as follows:
- GLVis draws surface elements; the number of drawn elements, `ne`, is either:
  - the number of elements in the mesh for 2D meshes (including surface meshes,
    i.e. 2D meshes embedded in 3D space), or
  - the number of boundary mesh elements described by the mesh in 3D.
- A tentative upper limit on the number of vertices to be drawn is defined based
  on the maximum order of the mesh and the solution, `max_order`:
  ```
  max_vert = ne * (max_order + 1) * (max_order + 1)
  ```
- To allow more accurate representation for small meshes, this number is
  potentially increased:
  ```
  max_vert = max(max_vert, 100 000)
  ```
- To keep the time to initially visualize the data reasonable, this number is
  potentially reduced:
  ```
  max_vert = min(max_vert, 2 000 000)
  ```
- Finally, the subdivision factor `ref` is chosen to be the largest number such
  that:
  - the number of vertices needed to draw the `ne` surface elements with `ref`
    subdivisions does not exceed `max_vert`:
    ```
    ne * (ref + 1) * (ref + 1) <= max_vert
    ```
  - for large meshes where the above limit cannot be satisfied, set `ref = 1`
  - for small meshes, avoid excessive refinements:
    ```
    ref <= 16
    ```

Note that, for highly-varying data or large meshes, this auto-selected
subdivision factor may not be sufficient for accurate representation. In such
cases the subdivision can be manually adjusted using the keys <kbd>o</kbd> /
<kbd>O</kbd>, described above.
