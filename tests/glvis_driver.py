# Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
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
from base64 import b64encode
from skimage.io import imread, imsave
from skimage.metrics import structural_similarity
from skimage.color import rgb2gray, gray2rgb
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def compare_images(
    baseline_file: str,
    output_file: str,
    expect_fail: bool = False,
    CUTOFF_SSIM: float = 0.999
) -> bool:

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

def color_distance(I1: np.array, I2: np.array) -> dict[str, np.array]:
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

def generate_image_diffs(
    image1_filename: str,
    image2_filename: str,
    absdiff_filename: str,
    reldiff_filename: str,
) -> None:
    # Images are read as NxMx3 [uint8] arrays from [0,255]
    I1 = imread(image1_filename)
    I2 = imread(image2_filename)
    # Get the image diffs (abs and rel)
    Idiffs = color_distance(I1, I2) # output is NxM [0,1]
    # Save 3-channel image to file
    imsave(absdiff_filename, gray2rgb(Idiffs['abs']))
    imsave(reldiff_filename, gray2rgb(Idiffs['rel']))

# For the source= argument in plotly
def _get_image_src(filename):
    with open(filename, "rb") as f:
        image_bytes = b64encode(f.read()).decode()
        return f"data:image/png;base64,{image_bytes}"

def image_comparison_plot(
    image_filenames: list[str],
    image_names: list[str],  # for subtitles
    output_filename: str,
):
    """
    Illustrate results as an interactive plotly figure (html)
    """
    assert len(image_filenames) == len(image_names)
    n = len(image_filenames)
    fig = make_subplots(rows=1, cols=n,
                        shared_xaxes=True,
                        shared_yaxes=True,
                        subplot_titles=image_filenames)
    for idx, filename in enumerate(image_filenames):
        fig.add_trace(go.Image(source=_get_image_src(filename)), 1, idx)
    fig.update_xaxes(matches='x', showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(matches='y', showticklabels=False, showgrid=False, zeroline=False)
    fig.write_html(output_filename)

def test_stream(
    exec_path: str,
    exec_args: str,
    save_file: str,
    baseline: str
) -> bool:

    if exec_args is None:
        exec_args = ""
    print(f"Testing {save_file}...")
    test_name = os.path.basename(save_file).replace(".saved", "") # e.g. "ex3"
    output_dir = f"outputs/{test_name}"
    os.makedirs(output_dir, exist_ok=True)

    # Create new stream file with command to screenshot and close
    stream_data = None
    with open(save_file) as in_f:
        stream_data = in_f.read()

    output_name = f"{output_dir}/test.nominal.{test_name}.png"
    output_name_fail = f"{output_dir}/test.zoom.{test_name}.png"
    absdiff_name = f"{output_dir}/test.nominal.absdiff.{test_name}.png"
    reldiff_name = f"{output_dir}/test.nominal.reldiff.{test_name}.png"
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
        generate_image_diffs(baseline_name, output_name, absdiff_name, reldiff_name)
        # Generate an interactive html plot, only if the test fails
        if not test_baseline:
            image_comparison_plot([baseline_name, output_name, reldiff_name],
                                  ["Baseline", "Test Output", "Normalized Diff"],
                                  reldiff_name.replace(".png", ".html"))
        test_control = compare_images(baseline_name, output_name_fail, expect_fail=True)
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

    # Make a directory for storing test outputs
    os.makedirs("outputs", exist_ok=True)
    # Run tests
    if args.save_stream is not None:
        result = test_stream(args.exec_cmd, args.exec_args, args.save_stream, args.baseline)
        if not result:
            sys.exit(1)
    else:
        raise Exception("--save_stream must be specified. test_cmd() is unused. Import from `test_cmd.py`")
