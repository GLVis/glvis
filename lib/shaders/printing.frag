R"(
uniform bool containsText; 
uniform bool useColorTex;
 
uniform sampler2D fontTex; 
uniform sampler2D colorTex;

varying vec4 fColor;
varying float fClipCoord;

uniform bool useClipPlane;
 
void main()
{
    gl_FragColor = fColor;
}
)"
