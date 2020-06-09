R"(
uniform sampler2D alphaTex;
uniform sampler2D colorTex;
 
varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord;

uniform bool useClipPlane;
varying float fClipVal;

void fragmentClipPlane() {
    if (useClipPlane && fClipVal < 0.0) {
        discard;
    }
}

void main() 
{
    fragmentClipPlane();
    vec4 color = fColor * texture2D(colorTex, vec2(fTexCoord));
    color = blinnPhong(fPosition, fNormal, color);
#ifdef GL_ES
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).a;
#else
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).r;
#endif
    gl_FragColor = color;
}
)"
