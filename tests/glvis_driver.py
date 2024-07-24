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
import numpy as np
from skimage.io import imread
from skimage.metrics import structural_similarity
from skimage.color import rgb2gray, gray2rgb
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def compare_images(baseline_file, output_file, expect_fail=False, CUTOFF_SSIM=0.999):
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
    if ssim < CUTOFF_SSIM:
        if expect_fail:
            print("[PASS] Differences were detected in the control case.")
        else:
            print("[FAIL] Output and baseline are different.")
    else:
        if expect_fail:
            print("[FAIL] Differences were not detected in the control case.")
        else:
            print("[PASS] Images match.")
    print(f"       actual ssim = {ssim}, cutoff = {CUTOFF_SSIM}")
    return ssim >= CUTOFF_SSIM if not expect_fail else ssim < CUTOFF_SSIM

def color_distance(I1, I2):
    """
    L2-norm in rgb space. There are better ways but this is probably good enough.
    """
    NORM_CONSTANT = (3*(255**2))**0.5 # max distance
    l2norm = lambda x: np.linalg.norm(x, ord=2, axis=2)
    delta = l2norm(I2.astype(int)-I1.astype(int)) / NORM_CONSTANT # output is NxM [0,1]
    # now we scale to [0,255] and cast as uint8 so it is a "proper" image
    Idiff_abs = (delta * 255).astype(np.uint8)
    # get relative version
    Idiff_rel = (Idiff_abs / Idiff_abs.max() * 255).astype(np.uint8)
    return {'abs': Idiff_abs,
            'rel': Idiff_rel,}

def generate_image_diff(
    image1_filename: str,
    image2_filename: str,
    imagediff_filename: str,
) -> None:
    # Images are read as NxMx3 [uint8] arrays from [0,255]
    I1 = imread(image1_filename)
    I2 = imread(image2_filename)
    # Take absolute "diff"
    Idiff = color_distance(I1, I2) # output is NxM [0,1]
    # Illustrate results
    fig = make_subplots(rows=1, cols=3, subplot_titles=('I1', 'I2', 'ΔI (normalized)'))
    fig.add_trace(go.Image(z=I1), 1, 1)
    fig.add_trace(go.Image(z=I2), 1, 2)
    fig.add_trace(go.Heatmap(z=Idiff['rel'], colorscale='inferno'), 1, 3)
    fig.update_yaxes(autorange='reversed', scaleanchor='x', constrain='domain', row=1, col=3)
    fig.update_xaxes(constrain='domain', row=1, col=3)
    fig.write_html(imagediff_filename)

def test_stream(exec_path, exec_args, save_file, baseline):
    if exec_args is None:
        exec_args = ""
    test_name = os.path.basename(save_file).replace(".saved", "")
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
    cmd = f"{exec_path} {exec_args} -saved {tmp_file}"
    print(f"Exec: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"[FAIL] GLVis exited with error code {ret}.")
        return False

    if baseline:
        baseline_name = f"{baseline}/test.{test_name}.saved.png"
        test_baseline = compare_images(baseline_name, output_name)
        test_control = compare_images(baseline_name, output_name_fail,
                                      expect_fail=True)
        return (test_baseline and test_control)
    else:
        print("[IGNORE] No baseline exists to compare against.")
        return True

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
        raise Exception("test_cmd() is unused. Import from `test_cmd.py`")
