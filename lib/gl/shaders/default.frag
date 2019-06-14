R"(
uniform sampler2D alphaTex;
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
    vec4 color = fColor * texture2D(colorTex, vec2(fTexCoord));
#ifdef GL_ES
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).a;
#else
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).r;
#endif
    color = blinnPhong(fPosition, fNormal, color);
    gl_FragColor = color;
}
)"
