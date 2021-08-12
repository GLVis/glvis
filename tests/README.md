                           GLVis visualization tool

                     _/_/_/  _/      _/      _/  _/
                  _/        _/      _/      _/        _/_/_/
                 _/  _/_/  _/      _/      _/  _/  _/_/
                _/    _/  _/        _/  _/    _/      _/_/
                 _/_/_/  _/_/_/_/    _/      _/  _/_/_/

                              http://glvis.org


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
python3 glvis_driver.py -e ../glvis -s data/streams/shaper.saved -b data/baselines -a "-lw 1 -mslw 1 -nohidpi"
```

A screenshot named "test.[stream file name].png" will be generated. If no path
to a stream baseline directory is given, the test will exit; otherwise it will
attempt to find and open the corresponding screenshot in the baseline directory
and compare the two images.

Details: Screenshot-based testing
---------------------------------
As implemented in the `glvis_driver.py` test,
[peak signal-to-noise ratio (PSNR)](https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio)
is used to determine similarity between the generated test screenshot and the
baseline screenshot.

The current cutoff for an image and its baseline is that PSNR must exceed 12dB;
this is much lower than would normally be expected, but is necessary to handle
differences in line rasterization between different OpenGL implementations.
