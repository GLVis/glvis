#include "glstate.hpp"
#include <string>
#include <regex>
#include <iostream>

using std::cerr;
using std::endl;

#ifdef __EMSCRIPTEN__
const std::string _glsl_add = "precision mediump float;\n";
#else
const std::string _glsl_add = "#version GLSL_VER\n";
#endif

enum ShaderFile {
    VS_CLIP_PLANE = 0,
    FS_CLIP_PLANE,
    VS_LIGHTING,
    FS_LIGHTING,
    VS_DEFAULT,
    FS_DEFAULT,
    VS_PRINTING,
    FS_PRINTING,
    NUM_SHADERS
};

const std::string shader_files[] = {
#include "shaders/clip_plane.vert"
,
#include "shaders/clip_plane.frag"
,
#include "shaders/lighting.glsl"
,
#include "shaders/lighting.glsl"
,
#include "shaders/default.vert"
,
#include "shaders/default.frag"
,
#include "shaders/printing.vert"
,
#include "shaders/printing.frag"
};

GLuint compileShaderFile(GLenum shaderType, const std::string& shaderText, int glslVersion) {
    GLuint shader = glCreateShader(shaderType);
    GLint success = 0;
    std::string fmt_shader = shaderText;
    if (glslVersion >= 130) {
        if (shaderType == GL_VERTEX_SHADER) {
            fmt_shader = std::regex_replace(fmt_shader, std::regex("attribute"), "in");
            fmt_shader = std::regex_replace(fmt_shader, std::regex("varying"), "out");
        } else { // FRAGMENT_SHADER
            fmt_shader = std::regex_replace(fmt_shader, std::regex("varying"), "in");
            if (glslVersion >= 140) {
                fmt_shader = "layout(location = 0) out vec4 fragColor;\n" + fmt_shader;
                fmt_shader = std::regex_replace(fmt_shader, std::regex("gl_FragColor"), "fragColor");
            }
        }
        fmt_shader = std::regex_replace(fmt_shader, std::regex("texture2D"), "texture");
    }
    fmt_shader = _glsl_add + fmt_shader;
    fmt_shader = std::regex_replace(fmt_shader, std::regex("GLSL_VER"), std::to_string(glslVersion));
    int shader_len = fmt_shader.length();
    const char * shader_cstr = fmt_shader.c_str();
    glShaderSource(shader, 1, &shader_cstr, &shader_len);
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "Shader compilation failed." << endl;
        int err_len; 
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &err_len);
        char * error_text = new char[err_len];
        glGetShaderInfoLog(shader, err_len, &err_len, error_text);
        cerr << error_text << endl;
        delete [] error_text;
        return 0;
    }
    return shader;
}

template<int Count>
bool linkShaders(GLuint prgm, const GLuint (&shaders)[Count]) {
    // explicitly specify attrib positions so we don't have to reset VAO
    // bindings when switching programs

    // for OSX, attrib 0 must be bound to render an object
    glBindAttribLocation(prgm, GlState::ATTR_VERTEX, "vertex");
    glBindAttribLocation(prgm, GlState::ATTR_TEXT_VERTEX, "textVertex");
    glBindAttribLocation(prgm, GlState::ATTR_NORMAL, "normal");
    glBindAttribLocation(prgm, GlState::ATTR_COLOR, "color");
    glBindAttribLocation(prgm, GlState::ATTR_TEXCOORD0, "texCoord0");
    glBindAttribLocation(prgm, GlState::ATTR_TEXCOORD1, "texCoord1");
    for (int i = 0; i < Count; i++) {
        glAttachShader(prgm, shaders[i]);
    }
    glLinkProgram(prgm);
    GLint success = 0;
    glGetProgramiv(prgm, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
        cerr << "FATAL: Shader linking failed." << endl;
    }
    return (success == GL_TRUE);
}

bool GlState::compileShaders() {
    GLuint refshaders[NUM_SHADERS];
    int glsl_ver = 100;
#ifndef __EMSCRIPTEN__
    if (GLEW_VERSION_3_0) {
        //GLSL 1.30-1.50, 3.30+
        int ver_major, ver_minor;
        glGetIntegerv(GL_MAJOR_VERSION, &ver_major);
        glGetIntegerv(GL_MINOR_VERSION, &ver_minor);
        glsl_ver = ver_major * 100 + ver_minor * 10;
        if (glsl_ver < 330) {
            // OpenGL 3.2 -> GLSL 1.50
            // OpenGL 3.1 -> GLSL 1.40
            glsl_ver -= 170;
        }
    } else if (GLEW_VERSION_2_0) {
        glsl_ver = 110;
    }
#endif
    for (int i = 0; i < NUM_SHADERS; i++) {
        GLenum shader_type = (i % 2 == 0) ? GL_VERTEX_SHADER
                                          : GL_FRAGMENT_SHADER;
        refshaders[i] = compileShaderFile(shader_type, shader_files[i], glsl_ver);
        if (refshaders[i] == 0)
            return false;
    }
    default_program = glCreateProgram();
    GLuint default_pipeline[] = {
        refshaders[VS_CLIP_PLANE],
        refshaders[VS_DEFAULT],
        refshaders[FS_LIGHTING],
        refshaders[FS_CLIP_PLANE],
        refshaders[FS_DEFAULT]
    };
    if (!linkShaders(default_program, default_pipeline)) {
        glDeleteProgram(default_program);
        default_program = 0;
        return false;
    }
    //TODO: enable a legacy path for opengl2.1 without ext_tranform_feedback?
    if (GLEW_EXT_transform_feedback || GLEW_VERSION_3_0) {
        feedback_program = glCreateProgram();
        const char * xfrm_varyings[] = {
            "gl_Position",
            "fColor",
            "fClipCoord",
        };
        glTransformFeedbackVaryings(feedback_program, 3, xfrm_varyings,
                                    GL_INTERLEAVED_ATTRIBS);
        GLuint print_pipeline[] = {
            refshaders[VS_LIGHTING],
            refshaders[VS_PRINTING],
            refshaders[FS_PRINTING]
        };
        if (!linkShaders(feedback_program, print_pipeline)) {
            glDeleteProgram(feedback_program);
            feedback_program = 0;
        }
    }
    glUseProgram(default_program);
    initShaderState(default_program);
    if (GLEW_VERSION_3_0) {
        if (global_vao == 0) {
            glGenVertexArrays(1, &global_vao);
        }
        glBindVertexArray(global_vao);
    }
    return true;
}

/**
 * Loads uniform locations and reloads previously-set uniform values for a
 * program.
 */
void GlState::initShaderState(GLuint program) {
    locUseClipPlane = glGetUniformLocation(program, "useClipPlane");
    locClipPlane = glGetUniformLocation(program, "clipPlane");
    
    locModelView = glGetUniformLocation(program, "modelViewMatrix");
    locProject = glGetUniformLocation(program, "projectionMatrix");
    locProjectText = glGetUniformLocation(program, "textProjMatrix");
    locNormal = glGetUniformLocation(program, "normalMatrix");

    locNumLights = glGetUniformLocation(program, "num_lights");
    locGlobalAmb = glGetUniformLocation(program, "g_ambient");

    locSpec = glGetUniformLocation(program, "material.specular");
    locShin = glGetUniformLocation(program, "material.shininess");

    for (int i = 0; i < MAX_LIGHTS; i++) {
        std::string location = "lights[" + std::to_string(i) + "].";
        locPosition[i] = glGetUniformLocation(program, (location + "position").c_str());
        locDiffuse[i] = glGetUniformLocation(program, (location + "diffuse").c_str());
        locSpecular[i] = glGetUniformLocation(program, (location + "specular").c_str());
    }

    //Texture unit 0: color palettes
    //Texture unit 1: font atlas
    GLuint locColorTex = glGetUniformLocation(program, "colorTex");
    GLuint locFontTex = glGetUniformLocation(program, "fontTex");
    glUniform1i(locColorTex, 0);
    glUniform1i(locFontTex, 1);
    // Set render type uniforms
    locContainsText = glGetUniformLocation(program, "containsText");
    locUseColorTex = glGetUniformLocation(program, "useColorTex");
    if (_shaderMode == RENDER_COLOR) {
        glUniform1i(locContainsText, GL_FALSE);
        glUniform1i(locUseColorTex, GL_FALSE);
    } else if (_shaderMode == RENDER_COLOR_TEX) {
        glUniform1i(locContainsText, GL_FALSE);
        glUniform1i(locUseColorTex, GL_TRUE);
    } else { //_shaderMode == RENDER_TEXT
        glUniform1i(locContainsText, GL_TRUE);
        glUniform1i(locUseColorTex, GL_FALSE);
    }
    // Set lighting uniforms
    glUniform1i(locNumLights, gl_lighting ? _num_lights : 0);
    glUniform4fv(locGlobalAmb, 1, _ambient);
    glUniform4fv(locSpec, 1, _mat.specular);
    glUniform1f(locShin, _mat.shininess);
    for (int i = 0; i < MAX_LIGHTS; i++) {
        glUniform3fv(locPosition[i], 1, _pt_lights[i].position);
        glUniform4fv(locDiffuse[i], 1, _pt_lights[i].diffuse);
        glUniform4fv(locSpecular[i], 1, _pt_lights[i].specular);
    }
    // Set clip plane uniforms
    glUniform1i(locUseClipPlane, gl_clip_plane);
    glUniform4fv(locClipPlane, 1, glm::value_ptr(_clip_plane));
    // Set transform matrix uniforms
    loadMatrixUniforms(true);
}
