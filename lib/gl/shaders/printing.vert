R"(
attribute vec3 vertex;
attribute vec2 textVertex;
attribute vec4 color;
attribute vec3 normal;
attribute vec2 texCoord0;

uniform bool containsText;
uniform bool useColorTex;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 textProjMatrix;
uniform mat3 normalMatrix; 
 
uniform vec4 clipPlane;

varying vec4 fColor;
varying float fClipCoord;

uniform sampler2D colorTex;
uniform sampler2D alphaTex;

void main() 
{ 
    vec4 pos = modelViewMatrix * vec4(vertex, 1.0);
    vec3 eye_normal = normalize(normalMatrix * normal);
    fColor = color * texture2DLod(colorTex, vec2(texCoord0), 0.0);
    fColor.a = texture2DLod(alphaTex, vec2(texCoord0), 0.0).r;
    fColor = blinnPhong(pos.xyz, eye_normal, fColor);
    //colors normally get clamped after fragment shader stage
    fColor = clamp(fColor, 0.0, 1.0);
    fClipCoord = dot(vec4(pos.xyz, 1.0), clipPlane);
    gl_Position = projectionMatrix * pos;
})"

