#include "glstate.hpp"
#include <string>
#include <iostream>

using std::cerr;
using std::endl;

#ifdef __EMSCRIPTEN__
const std::string _glsl_add = "precision mediump float;\n";
#else
const std::string _glsl_add = "#version 120\n";
#endif

const std::string vertex_shader_file = _glsl_add +
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

uniform float matAlpha;
uniform float matAlphaCenter;

varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord; 
varying vec2 fFontTexCoord;
varying float fClipVal;
varying float fAlpha;
 
void main() 
{ 
    fNormal = normalize(normalMatrix * normal); 
    vec4 pos = modelViewMatrix * vec4(vertex, 1.0);
    fPosition = pos.xyz; 
    fColor = color; 
    fTexCoord = texCoord0.xy;
    fFontTexCoord = texCoord1.xy;
    fClipVal = dot(vec4(pos.xyz, 1.0), clipPlane);
    vec4 textOffset = textProjMatrix * vec4(textVertex, 0.0, 0.0);
    pos = projectionMatrix * pos;
    gl_Position = pos;
    fAlpha = matAlpha;
    if (matAlpha < 1.0) {
        if (matAlphaCenter > 1.0) {
            fAlpha *= exp(-(matAlphaCenter) * abs(texCoord0.x - 1.0));
        } else if (matAlphaCenter < 0.0) {
            fAlpha *= exp((matAlphaCenter - 1.0) * abs(texCoord0.x));
        } else {
            fAlpha *= exp(-abs(texCoord0.x - matAlphaCenter));
        }
    }
    if (containsText) {
        gl_Position += vec4((textOffset.xy * pos.w), -0.005, 0.0);
    }
})";

const std::string fragment_shader_file = _glsl_add +
R"(
uniform bool containsText; 
uniform bool useColorTex;
uniform bool useClipPlane;
 
uniform sampler2D fontTex; 
uniform sampler2D colorTex;
 
varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord; 
varying vec2 fFontTexCoord; 
varying float fClipVal;
varying float fAlpha;

struct PointLight { 
    vec3 position; 
    vec4 diffuse; 
    vec4 specular; 
}; 
 
uniform int numLights; 
uniform PointLight lights[3]; 
uniform vec4 g_ambient;

struct Material {
    vec4 specular; 
    float shininess; 
}; 
 
uniform Material material;
 
void main() 
{
    if (useClipPlane && fClipVal < 0.0) {
        discard;
    }
    if (containsText) { 
        vec4 colorOut = vec4(0.0, 0.0, 0.0, texture2D(fontTex, fFontTexCoord).a);
        gl_FragColor = colorOut;
    } else { 
        vec4 color = fColor; 
        if (useColorTex) { 
            color.xyz = texture2D(colorTex, fTexCoord).xyz;
            color.w = fAlpha;
        }
        if (numLights == 0) {
            gl_FragColor = color;
        } else {
            float normSgn = float(int(gl_FrontFacing) * 2 - 1);
            vec4 ambient_light = g_ambient * color; 
            vec4 diffuse_light = vec4(0.0, 0.0, 0.0, 0.0); 
            vec4 specular_light = vec4(0.0, 0.0, 0.0, 0.0); 
            for (int i = 0; i < 3; i++) {
                if (i == numLights)
                    break;
                vec3 light_dir = normalize(lights[i].position - fPosition); 
                diffuse_light += lights[i].diffuse * color * max(dot(fNormal * normSgn, light_dir), 0.0); 
     
                vec3 eye_to_vert = normalize(-fPosition);
                vec3 half_v = normalize(eye_to_vert + light_dir);
                float specular_factor = max(dot(half_v, normSgn * fNormal), 0.0); 
                specular_light += lights[i].specular * material.specular * pow(specular_factor, material.shininess); 
            } 
            gl_FragColor.xyz = ambient_light.xyz + diffuse_light.xyz + specular_light.xyz;
            gl_FragColor.w = color.w;
        }
    } 
})";

bool GlState::compileShaders() {
    GLuint vtx_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    GLint success = 0;
    int vertex_shader_len = vertex_shader_file.length();
    const char * vert_str = vertex_shader_file.c_str();
    glShaderSource(vtx_shader, 1, &vert_str, &vertex_shader_len);
    glCompileShader(vtx_shader);
    glGetShaderiv(vtx_shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "FATAL: Vertex shader compilation failed." << endl;
        int err_len;
        glGetShaderiv(vtx_shader, GL_INFO_LOG_LENGTH, &err_len);
        char * error_text = new char[err_len];
        glGetShaderInfoLog(vtx_shader, err_len, &err_len, error_text);
        cerr << error_text << endl;
        delete [] error_text;
        return false;
    }
    int frag_shader_len = fragment_shader_file.length();
    const char * frag_str = fragment_shader_file.c_str();
    glShaderSource(frag_shader, 1, &frag_str, &frag_shader_len);
    glCompileShader(frag_shader);
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "FATAL: Fragment shader compilation failed." << endl;
        int err_len;
        glGetShaderiv(frag_shader, GL_INFO_LOG_LENGTH, &err_len);
        char * error_text = new char[err_len];
        glGetShaderInfoLog(frag_shader, err_len, &err_len, error_text);
        cerr << error_text << endl;
        delete [] error_text;
        return false;
    }
    program = glCreateProgram();
    //for OSX, attrib 0 must be bound to render an object
    glBindAttribLocation(program, 0, "vertex");
    glAttachShader(program, vtx_shader);
    glAttachShader(program, frag_shader);

    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "FATAL: Shader linking failed." << endl;
        glDeleteProgram(program);
        program = 0;
    } else {
        glUseProgram(program);
    }

    glDetachShader(program, vtx_shader);
    glDetachShader(program, frag_shader);
    glDeleteShader(vtx_shader);
    glDeleteShader(frag_shader);
    return (success != GL_FALSE);
}

void GlState::initShaderState() {
    _attr_locs[ATTR_VERTEX] = glGetAttribLocation(program, "vertex");
    _attr_locs[ATTR_TEXT_VERTEX] = glGetAttribLocation(program, "textVertex");
    _attr_locs[ATTR_COLOR] = glGetAttribLocation(program, "color");
    _attr_locs[ATTR_NORMAL] = glGetAttribLocation(program, "normal");
    _attr_locs[ATTR_TEXCOORD0] = glGetAttribLocation(program, "texCoord0");
    _attr_locs[ATTR_TEXCOORD1] = glGetAttribLocation(program, "texCoord1");

    locUseClipPlane = glGetUniformLocation(program, "useClipPlane");
    locClipPlane = glGetUniformLocation(program, "clipPlane");
    
    locModelView = glGetUniformLocation(program, "modelViewMatrix");
    locProject = glGetUniformLocation(program, "projectionMatrix");
    locProjectText = glGetUniformLocation(program, "textProjMatrix");
    locNormal = glGetUniformLocation(program, "normalMatrix");

    locSpec = glGetUniformLocation(program, "material.specular");
    locShin = glGetUniformLocation(program, "material.shininess");

    locMatAlpha = glGetUniformLocation(program, "matAlpha");
    locMatAlphaCenter = glGetUniformLocation(program, "matAlphaCenter");
    glUniform1f(locMatAlpha, 1.f);

    //Texture unit 0: color palettes
    //Texture unit 1: font atlas
    GLuint locColorTex = glGetUniformLocation(program, "colorTex");
    GLuint locFontTex = glGetUniformLocation(program, "fontTex");
    glUniform1i(locColorTex, 0);
    glUniform1i(locFontTex, 1);
    GLuint locContainsText = glGetUniformLocation(program, "containsText");
    glUniform1i(locContainsText, GL_FALSE);
    GLuint locUseColorTex = glGetUniformLocation(program, "useColorTex");
    glUniform1i(locUseColorTex, GL_FALSE);
    _shaderMode = RENDER_COLOR;
    modelView.identity();
    projection.identity();
    loadMatrixUniforms();
}
