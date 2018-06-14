#include "glstate.hpp"
#include <string>
#include <iostream>

using std::cerr;
using std::endl;

const std::string vertex_shader_file = 
R"(#version 120 
 
//attribute vec4 gl_Vertex; 
//attribute vec4 gl_Color; 
//attribute vec4 gl_Normal; 
//attribute vec4 gl_MultiTexCoord0; 
 
uniform mat4 modelViewMatrix; 
uniform mat4 projectionMatrix; 
uniform mat3 normalMatrix; 
 
varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord; 
varying vec2 fFontTexCoord; 
 
void main() 
{ 
    fNormal = normalize(normalMatrix * gl_Normal); 
    vec4 pos = modelViewMatrix * gl_Vertex; 
    fPosition = pos.xyz; 
    fColor = gl_Color; 
    fTexCoord = gl_MultiTexCoord0.xy; 
    fFontTexCoord = gl_MultiTexCoord1.xy; 
    gl_Position = projectionMatrix * pos; 
})";

const std::string fragment_shader_file =
R"(#version 120 
 
uniform bool containsText; 
uniform bool useColorTex; 
 
uniform sampler2D fontTex; 
uniform sampler2D colorTex; 
 
varying vec3 fNormal; 
varying vec3 fPosition; 
varying vec4 fColor; 
varying vec2 fTexCoord; 
varying vec2 fFontTexCoord; 
 
struct PointLight { 
    vec3 position; 
    vec4 diffuse; 
    vec4 specular; 
}; 
 
uniform int numLights; 
uniform PointLight lights[3]; 
uniform vec4 g_ambient; 
 
struct Material { 
    vec4 ambient; 
    vec4 diffuse; 
    vec4 specular; 
    float shininess; 
}; 
 
uniform Material material; 
 
void main() 
{ 
    if (containsText) { 
        gl_FragColor = texture2D(fontTex, fFontTexCoord); 
    } else { 
        vec4 color = fColor; 
        if (useColorTex) { 
            color = texture2D(colorTex, fTexCoord); 
        }
        if (numLights == 0) {
            gl_FragColor = color;
        } else {
            float normSgn = float(int(gl_FrontFacing) * 2 - 1);
            vec4 ambient_light = g_ambient * material.ambient; 
            vec4 diffuse_light = vec4(0.0, 0.0, 0.0, 0.0); 
            vec4 specular_light = vec4(0.0, 0.0, 0.0, 0.0); 
            for (int i = 0; i < numLights; i++) { 
                vec3 light_dir = normalize(lights[i].position - fPosition); 
                diffuse_light += lights[i].diffuse * material.diffuse * max(dot(fNormal * normSgn, light_dir), 0.0); 
     
                //vec3 eye_to_vert = normalize(-fPosition); 
                vec3 half_v = normalize(vec3(0,0,1) + light_dir);
                float specular_factor = max(dot(half_v, fNormal * normSgn), 0.0); 
                specular_light += lights[i].specular * material.specular * pow(specular_factor, material.shininess); 
            } 
            gl_FragColor = color * (vec4(ambient_light.xyz, 0.0) + diffuse_light + vec4(specular_light.xyz, 0.0));
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
        return false;
    }
    int frag_shader_len = fragment_shader_file.length();
    const char * frag_str = fragment_shader_file.c_str();
    glShaderSource(frag_shader, 1, &frag_str, &frag_shader_len);
    glCompileShader(frag_shader);
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "FATAL: Fragment shader compilation failed." << endl;
        return false;
    }
    program = glCreateProgram();
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
