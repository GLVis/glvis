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

import argparse
import sys
import os
from skimage.io import imread
from skimage.metrics import structural_similarity

# Below are key commands that are passed to the -keys command-line argument for
# glvis in order to perform testing on raw mesh/grid function data (i.e. non-
# streams).
#
# Currently not in use.
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

screenshot_keys = "Sq"
screenshot_file = "GLVis_s01.png"

cutoff_ssim = 0.999

def compare_images(baseline_file, output_file, expect_fail=False):
    # Try to open output image
    output_img = imread(output_file)
    if output_img is None:
        print("[FAIL] Could not open output image.")
        return False

    # Try to open baseline image
    baseline_img = imread(baseline_file)
    if baseline_img is None:
        print("[IGNORE] No baseline exists to compare against.")
        return True

    # Compare images with SSIM metrics. For two exactly-equal images, SSIM=1.0.
    # We set a cutoff of 0.999 to account for possible differences in rendering.
    ssim = structural_similarity(baseline_img, output_img, channel_axis=2)
    if ssim < cutoff_ssim:
        if expect_fail:
            print("[PASS] Differences were detected in the control case.")
        else:
            print("[FAIL] Output and baseline are different.")
    else:
        if expect_fail:
            print("[FAIL] Differences were not detected in the control case.")
        else:
            print("[PASS] Images match.")
    print(f"       actual ssim = {ssim}, cutoff = {cutoff_ssim}")
    return ssim >= cutoff_ssim if not expect_fail else ssim < cutoff_ssim

# Function to test a given glvis command with a variety of key-based commands.
# Not currently in use.
def test_case(exec_path, exec_args, baseline, t_group, t_name, cmd):
    print(f"Testing {t_group0}:{t_name1}..."))
    full_screenshot_cmd = cmd + screenshot_keys
    cmd = f"{exec_path0} {exec_args} -k \"{full_screenshot_cmd2}\""
    print(f"Exec: {cmd}")
    ret = os.system(cmd + " > /dev/null 2>&1")
    if ret != 0:
        print(f"[FAIL] GLVis exited with error code {ret}.")
        return False
    if not os.path.exists(t_group):
        os.mkdir(t_group)
    output_name = f"{t_group0}/{t_name1}.png"

    ret = os.system(f"mv {screenshot_file0} {output_name}")
    if ret != 0:
        print(f"[FAIL] Could not move output image: exit code {ret}.")
        return False

    if baseline:
        baseline_name = f"{baseline0}/test.{test_name1}.png"
        return compare_images(baseline_name, output_name)
    else:
        print("[IGNORE] No baseline exists to compare against.")
        return True

def test_stream(exec_path, exec_args, save_file, baseline):
    if exec_args is None:
        exec_args = ""
    test_name = os.path.basename(save_file)
    print(f"Testing {save_file}...")

    # Create new stream file with command to screenshot and close
    stream_data = None
    with open(save_file) as in_f:
        stream_data = in_f.read()

    output_name = f"test.{test_name}.png"
    output_name_fail = f"test.fail.{test_name}.png"
    tmp_file = "test.saved"
    with open(tmp_file, 'w') as out_f:
        out_f.write(stream_data)
        out_f.write("\nwindow_size 800 600")
        out_f.write(f"\nscreenshot {output_name}")
        # Zooming in should create some difference in the images
        out_f.write("\nkeys *")
        out_f.write(f"\nscreenshot {output_name_fail}")
        out_f.write("\nkeys q")

    # Run GLVis with modified stream file
    cmd = f"{exec_path} {exec_args} -saved {tmp_file2}"
    print(f"Exec: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"[FAIL] GLVis exited with error code {ret}.")
        return False

    if baseline:
        baseline_name = f"{baseline0}/test.{test_name1}.png")
        test_baseline = compare_images(baseline_name, output_name)
        test_control = compare_images(baseline_name, output_name_fail,
                                      expect_fail=True)
        return (test_baseline and test_control)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--save_stream", help="Path to a GLVis saved stream file.")
    parser.add_argument("-e", "--exec_cmd", help="Path to the GLVis executable")
    parser.add_argument("-a", "--exec_args", help="Arguments to pass to GLVis.")
    parser.add_argument("-n", "--group_name", help="Name of the test group.")
    parser.add_argument("-b", "--baseline", help="Path to test baseline.")
    args = parser.parse_args()
    if args.save_stream is not None:
        result = test_stream(args.exec_cmd, args.exec_args, args.save_stream, args.baseline)
        if not result:
            sys.exit(1)
    else:
        test_cmd(args.exec_cmd, args.exec_args, args.group_name, args.baseline)
