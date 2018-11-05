R"(
uniform bool useClipPlane;
varying float fClipVal;

void fragmentClipPlane() {
    if (useClipPlane && fClipVal < 0.0) {
        discard;
    }
}
)"
