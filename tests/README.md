                           GLVis visualization tool

                     _/_/_/  _/      _/      _/  _/
                  _/        _/      _/      _/        _/_/_/
                 _/  _/_/  _/      _/      _/  _/  _/_/
                _/    _/  _/        _/  _/    _/      _/_/
                 _/_/_/  _/_/_/_/    _/      _/  _/_/_/

                             https://glvis.org


Automated testing for GLVis
===========================
This directory contains scripts and data for regression testing of GLVis:

- `glvis_driver.py` is a Python script for running test cases and performing
  comparisons on generated screenshots. Currently, the only supported mode is
  *stream-based testing*, where given a saved input stream, commands are added
  to the end of the stream to generate a screenshot and exit. This screenshot
  can then be compared against a baseline image.

- Test streams and baselines used in automated testing are stored in another
  repository, [GLVis/data](https://github.com/GLVis/data), which is linked
  into this repository as a submodule under the directory data.

  To fetch the test data submodule, call `git submodule update --init` in an
  existing clone of this repository.

Dependencies
------------
Running tests with the test driver requires an installation of Python, with the
`scikit-image` package. Install it using `pip` or your distribution's package
manager.

CMake is required to run the automated tests as they are used in Github Actions.
The tests also expect `.png` images, so GLVis should be built with PNG support.

CMake-based testing
-------------------
The test cases for automated testing are defined using `CTest`, and can be found
in the `CMakeLists.txt` file in this directory. This is the method that is used
by the Github Actions tests.

To run all the tests:
- Ensure that you have pulled in the data submodule.
- Build GLVis using CMake, passing `-DENABLE_TESTS=ON` as a configuration
  option.
- Run `ctest` or `make test` to run all the test cases.

Note that a set of baselines need to be generated manually for your specific
system, due to differences in rendering between operating systems and graphics
hardware. After the first run of `ctest` on master, run `make rebaseline` in
the build directory. This will copy the generated images from the test run to a
system-specific baseline directory.

Standalone testing
------------------
You can also run an individual test case outside of `CTest`, with the following
command:

```sh
python3 glvis_driver.py -e [path to glvis] \
                        -s [path to stream file] \
                        -b [path to stream baseline *directory*] \
                        -a [additional commands to pass to glvis]
```

for example

```sh
python3 glvis_driver.py -e ../glvis -s data/streams/shaper.saved -b data/baselines/local -a "-lw 1 -mslw 1 -nohidpi"
```

A screenshot named "test.[stream file name].png" will be generated. If no path
to a stream baseline directory is given, the test will exit; otherwise it will
attempt to find and open the corresponding screenshot in the baseline directory
and compare the two images.

Details: Screenshot-based testing
---------------------------------
As implemented in the `glvis_driver.py` test,
[structural similarity (SSIM)](https://en.wikipedia.org/wiki/Structural_similarity)
is used to determine similarity between the generated test screenshot and the
baseline screenshot.

The current cutoff for an image and its baseline is set to 0.999; this allows
for minor differences in rendering output.

Details: Baselines
------------------
Two sets of public baselines are maintained: `ubuntu-20.04` and `macos-10.15`,
each representing their corresponding Github Actions virtual environment. These
are selected based on the value of the `GLVIS_BASELINE_SYS` CMake variable (set
to `local` by default).

You should leave this variable unset for development purposes, even if your
development platform shares the same operating system; baselines are usually
generated in these virtual environments with software rendering, and may not
match the output of your specific system/graphics hardware.

Reminder: after running `ctest` on master, developers can run `make rebaseline`
to populate the local baseline with reference figures.
