R"(
uniform bool containsText; 
uniform bool useColorTex;
 
uniform sampler2D fontTex; 
uniform sampler2D colorTex;
 
varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord;

void fragmentClipPlane();
vec4 blinnPhong(in vec3 pos, in vec3 norm, in vec4 color);

void main() 
{
    fragmentClipPlane();
    vec4 color;
    if (containsText) {

#ifdef GL_ES
        color = fColor * vec4(1.0, 1.0, 1.0, texture2D(fontTex, fTexCoord).a);
#else
        color = fColor * vec4(1.0, 1.0, 1.0, texture2D(fontTex, fTexCoord).r);
#endif
    } else {
        if (useColorTex) {
            color.xyz = texture2D(colorTex, vec2(fTexCoord.x, 0.0)).xyz;
            color.w = fTexCoord.y;
        } else {
            color = fColor; 
        }
        color = blinnPhong(fPosition, fNormal, color);
    }
    gl_FragColor = color;
}
)"
