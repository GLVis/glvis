import argparse
import sys
import os
import cv2

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

def test_case(exec_path, exec_args, baseline, t_group, t_name, cmd):
    print("Testing {0}:{1}...".format(t_group, t_name))
    full_screenshot_cmd = cmd + screenshot_keys
    cmd = "{0} {1} -k \"{2}\"".format(exec_path, exec_args, full_screenshot_cmd)
    print("Test: {}".format(cmd))
    os.system(cmd + " > /dev/null 2>&1")
    if not os.path.exists(t_group):
        os.mkdir(t_group)
    output_name = "{0}/{1}.png".format(t_group, t_name)

    ret = os.system("mv {0} {1}".format(screenshot_file, output_name))
    if ret != 0:
        print("[FAIL] GLVis exited with error code {}.".format(ret))
        return False

    # Try to open output image
    output_img = cv2.imread(output_name)
    if output_img is None:
        print("[FAIL] Could not open output image.")
        return False

    # Try to open baseline image
    baseline_img = None
    if baseline:
        baseline_name = "{0}/{1}.png".format(baseline, t_name)
        baseline_img = cv2.imread(baseline_name)
    if baseline_img is None:
        print("[IGNORE] No baseline exists to compare against.")
        return True
    diff_img = cv2.subtract(output_img, baseline_img)
    diff = cv2.norm(diff_img)
    if diff > 0.001:
        print("[FAIL] Output and baseline differ by more than 0.1%")
    else:
        print("[PASS] Images match.")
    return diff <= 0.001

def main(exec_path, exec_args, tgroup, baseline):
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
    parser.add_argument("-e", "--exec_cmd", help="Path to the GLVis executable")
    parser.add_argument("-a", "--exec_args", help="Arguments to pass to GLVis.")
    parser.add_argument("-n", "--group_name", help="Name of the test group.")
    parser.add_argument("-b", "--baseline", help="Path to test baseline.")
    args = parser.parse_args()
    main(args.exec_cmd, args.exec_args, args.group_name, args.baseline)
