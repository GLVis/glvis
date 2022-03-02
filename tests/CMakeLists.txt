# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

if(NOT EMSCRIPTEN)

# Used for automated testing
set(GLVIS_BASELINE_SYS local CACHE STRING
    "The location to use for testing baselines.")

add_custom_target(
  glvis-test-driver ALL
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/glvis_driver.py
    ${CMAKE_CURRENT_BINARY_DIR}/glvis_driver.py
  DEPENDS glvis_driver.py
  COMMENT "Copying glvis_driver.py")

set(stream_tests
    "capacitor"
    "distance"
    "ex1"
    "ex2"
    "ex27"
    "ex3"
    "ex5"
    "ex9"
    "klein-bottle"
    "laghos"
    "mesh-explorer"
    "mfem-logo"
    "minimal-surface"
    "mobius-strip"
    "navier"
    "quad"
    "remhos"
    "shaper"
    "shifted"
    "snake")

add_custom_target(rebaseline)

foreach(test_name IN LISTS stream_tests)

  add_test(NAME stream_${test_name}
           COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/glvis_driver.py
           -e ${CMAKE_CURRENT_BINARY_DIR}/../glvis
           -s ${CMAKE_CURRENT_SOURCE_DIR}/data/streams/${test_name}.saved
           -b ${CMAKE_CURRENT_SOURCE_DIR}/data/baselines/${GLVIS_BASELINE_SYS}
           -a "-lw 1 -mslw 1 -nohidpi")

  add_custom_target(_rebaseline_stream_${test_name}
      COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/baselines/local
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/test.${test_name}.saved.png
        ${CMAKE_CURRENT_SOURCE_DIR}/data/baselines/local
      DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/test.${test_name}.saved.png
      VERBATIM)

  add_dependencies(rebaseline _rebaseline_stream_${test_name})

endforeach()

endif(NOT EMSCRIPTEN)
