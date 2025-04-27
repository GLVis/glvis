# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

"""
Code snippets to test glvis in the command-line. None of the code
contained here is currently being used.
"""

# Globals
screenshot_keys = "Sq"
screenshot_file = "GLVis_s01.png"

# Below are key commands that are passed to the -keys command-line argument for
# glvis in order to perform testing on raw mesh/grid function data (i.e. non-
# streams).
test_cases = {
    "magnify": "*****",
    "axes1": "a",
    "axes2": "aa",
    "mesh1": "m",
    "mesh2": "mm",
    "cut_plane": "i",
    "cut_plane_rotate": "iyyyy",
    "cut_plane_rotate_back": "iyyyyYYYY",
    "cut_plane_transl": "izzzz",
    "cut_plane_transl_back": "izzzzZZZZ",
    "orient2d_1": "R",
    "orient2d_2": "RR",
    "orient2d_3": "RRR",
    "orient2d_4": "RRRR",
    "orient2d_5": "RRRRR",
    "orient2d_6": "RRRRRR",
    "orient3d": "Rr",
    "perspective": "j",
}

# Function to test a given glvis command with a variety of key-based commands.
def test_case(exec_path, exec_args, baseline, t_group, t_name, cmd):
    print(f"Testing {t_group}:{t_name}...")
    full_screenshot_cmd = cmd + screenshot_keys
    cmd = f"{exec_path} {exec_args} -k \"{full_screenshot_cmd}\""
    print(f"Exec: {cmd}")
    ret = os.system(cmd + " > /dev/null 2>&1")
    if ret != 0:
        print(f"[FAIL] GLVis exited with error code {ret}.")
        return False
    if not os.path.exists(t_group):
        os.mkdir(t_group)
    output_name = f"{t_group}/{t_name}.png"

    ret = os.system(f"mv {screenshot_file} {output_name}")
    if ret != 0:
        print(f"[FAIL] Could not move output image: exit code {ret}.")
        return False

    if baseline:
        baseline_name = f"{baseline}/test.{t_name}.png"
        return compare_images(baseline_name, output_name)
    else:
        print("[IGNORE] No baseline exists to compare against.")
        return True

def test_cmd(exec_path, exec_args, tgroup, baseline):
    try:
        os.remove(screenshot_file)
    except OSError:
        pass
    all_tests_passed = True
    for testname, cmds in test_cases.items():
        result = test_case(exec_path, exec_args, baseline, tgroup, testname, cmds)
        all_tests_passed = all_tests_passed and result

    if all_tests_passed:
        print("All tests passed.")
    else:
        sys.exit(1)
