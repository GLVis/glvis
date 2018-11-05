R"(
attribute vec3 vertex;
attribute vec2 textVertex;
attribute vec4 color;
attribute vec3 normal;
attribute vec2 texCoord0;
attribute vec2 texCoord1;

uniform bool containsText;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 textProjMatrix;
uniform mat3 normalMatrix; 
 
uniform vec4 clipPlane;

varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord; 

void setupClipPlane(in float dist);

void main() 
{ 
    vec4 pos = modelViewMatrix * vec4(vertex, 1.0);
    fPosition = pos.xyz; 
    fNormal = normalize(normalMatrix * normal); 
    fColor = color; 
    fTexCoord = texCoord0.xy;
    setupClipPlane(dot(vec4(pos.xyz, 1.0), clipPlane));
    pos = projectionMatrix * pos;
    gl_Position = pos;
    if (containsText) {
        vec4 textOffset = textProjMatrix * vec4(textVertex, 0.0, 0.0);
        fTexCoord = texCoord1.xy;
        gl_Position += vec4((textOffset.xy * pos.w), -0.005, 0.0);
    }
})"

