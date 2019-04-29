#include "glstate.hpp"
#include "visual.hpp"
#include <string>
#include <regex>
#include <iostream>

using std::cerr;
using std::endl;

#ifdef __EMSCRIPTEN__
// TODO: for webgl glsl the version ends in "es": #version GLSL_VER es
const std::string GLSL_HEADER = "precision mediump float;\n";
#else
const std::string GLSL_HEADER = "#version GLSL_VER\n";
#endif


// weird but loads them inline
const std::string CLIP_PLANE_VS =
#include "shaders/clip_plane.vert"
;

const std::string CLIP_PLANE_FS =
#include "shaders/clip_plane.frag"
;

const std::string BLINN_PHONG_FS =
#include "shaders/lighting.glsl"
;

const std::string DEFAULT_VS =
#include "shaders/default.vert"
;
const std::string DEFAULT_FS =
#include "shaders/default.frag"
;
const std::string PRINTING_VS =
#include "shaders/printing.vert"
;
const std::string PRINTING_FS =
#include "shaders/printing.frag"
;


std::string formatShaderString(const std::string & shader_string, GLenum shader_type, int glsl_version) {
  // webgl does not allow name resolution across shaders, we have to substitute them in here
  // add lighting
  std::string formatted = std::regex_replace(shader_string,
      std::regex(R"(vec4 blinnPhong\(in vec3 pos, in vec3 norm, in vec4 color\);)"),
      BLINN_PHONG_FS);

  // add clip plane
  formatted = std::regex_replace(formatted, std::regex(R"(void fragmentClipPlane\(\);)"),
      CLIP_PLANE_FS);
  formatted = std::regex_replace(formatted, std::regex(R"(void setupClipPlane\(in float dist\);)"),
      CLIP_PLANE_VS);

  // replace some identifiers depending on the verions of glsl we're using
  if (glsl_version >= 130) {
    if (shader_type == GL_VERTEX_SHADER) {
      formatted = std::regex_replace(formatted, std::regex("attribute"), "in");
      formatted = std::regex_replace(formatted, std::regex("varying"), "out");
    }

    else if (shader_type == GL_FRAGMENT_SHADER) {
      formatted = std::regex_replace(formatted, std::regex("varying"), "in");

      // requires GL_ARB_explicit_attrib_location extension or GLSL 3.30
      // although gl_FragColor was depricated in GLSL 1.3
      if (glsl_version > 130 && glsl_version < 330) {
          formatted = "out vec4 fragColor;\n" + formatted;
          formatted = std::regex_replace(formatted, std::regex("gl_FragColor"), "fragColor");
      }
      else if (glsl_version >= 330) {
          formatted = "layout(location = 0) out vec4 fragColor;\n" + formatted;
          formatted = std::regex_replace(formatted, std::regex("gl_FragColor"), "fragColor");
      }
    }

    else {
      std::cerr << "buildShaderString: unknown shader type" << std::endl;
      return {};
    }

    formatted = std::regex_replace(formatted, std::regex("texture2D"), "texture");
  }

  // add the header
  formatted = std::regex_replace(GLSL_HEADER, std::regex("GLSL_VER"), std::to_string(glsl_version)) + formatted;

  return formatted;
}


GLuint compileShaderString(const std::string & shader_str, GLenum shader_type, int glsl_version) {
    int shader_len = shader_str.length();
    const char * shader_cstr = shader_str.c_str();

    GLuint shader = glCreateShader(shader_type);
    glShaderSource(shader, 1, &shader_cstr, &shader_len);
    glCompileShader(shader);

    GLint stat;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &stat);
    // glGetObjectParameteriv
    if (stat == GL_FALSE) {
      std::cerr << "failed to compile shader" << std::endl;
      int err_len;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &err_len);
      char * error_text = new char[err_len];
      glGetShaderInfoLog(shader, err_len, &err_len, error_text);
      std::cerr << error_text << std::endl;
      delete [] error_text;
      return 0;
    }
    return shader;
}


GLuint compileRawShaderString(const std::string & shader_string, GLenum shader_type, int glsl_version) {
  std::string formatted = formatShaderString(shader_string, shader_type, glsl_version);
  //std::cout << "compiling '''" << formatted << "...";
  GLuint shader_ref = compileShaderString(formatted, shader_type, glsl_version);
  //if (shader_ref != 0) { std::cerr << "successfully compiled shader" << std::endl; }
  return shader_ref;
}


bool linkShaders(GLuint prgm, const std::vector<GLuint> & shaders) {
    // explicitly specify attrib positions so we don't have to reset VAO
    // bindings when switching programs

    // for OSX, attrib 0 must be bound to render an object
    glBindAttribLocation(prgm, GlState::ATTR_VERTEX, "vertex");
    glBindAttribLocation(prgm, GlState::ATTR_TEXT_VERTEX, "textVertex");
    glBindAttribLocation(prgm, GlState::ATTR_NORMAL, "normal");
    glBindAttribLocation(prgm, GlState::ATTR_COLOR, "color");
    glBindAttribLocation(prgm, GlState::ATTR_TEXCOORD0, "texCoord0");
    glBindAttribLocation(prgm, GlState::ATTR_TEXCOORD1, "texCoord1");

    for (GLuint i : shaders) {
        glAttachShader(prgm, i);
    }
    glLinkProgram(prgm);

    GLint stat;
    glGetProgramiv(prgm, GL_LINK_STATUS, &stat);
    if (stat == GL_FALSE) {
        cerr << "fatal: Shader linking failed" << endl;
    }
    return (stat == GL_TRUE);
}


bool GlState::compileShaders() {
    int ver_major, ver_minor;
    SDL_GL_GetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, &ver_major);
    SDL_GL_GetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, &ver_minor);
    int opengl_ver = ver_major * 100 + ver_minor * 10;
    int glsl_ver = -1;

#ifndef __EMSCRIPTEN__
    // The GLSL verion is the same as the OpenGL version
    // when the OpenGL version is >= 3.30, otherwise
    // it is:
    //
    // GL Version | GLSL Version
    // -------------------------
    //    2.0     |   1.10
    //    2.1     |   1.20
    //    3.0     |   1.30
    //    3.1     |   1.40
    //    3.2     |   1.50

    if (opengl_ver < 330) {
      if (ver_major == 2) {
        glsl_ver = opengl_ver - 90;
      }
      else if (ver_major == 3) {
        glsl_ver = opengl_ver - 170;
      }
      else {
        std::cerr << "fatal: unsupported OpenGL version " << opengl_ver << std::endl;
        return false;
      }
    }
    else {
      glsl_ver = opengl_ver;
    }
#else
    // GL Version | WebGL Version | GLSL Version
    //    2.0     |      1.0      |   1.00 ES
    //    3.0     |      2.0      |   3.00 ES
    //    3.1     |               |   3.10 ES
    if (opengl_ver < 300) {
      glsl_ver = 100;
    }
    else {
      glsl_ver = 300;
    }
#endif

    std::cerr << "Using GLSL " << glsl_ver << std::endl;

    GLuint default_vs = compileRawShaderString(DEFAULT_VS, GL_VERTEX_SHADER, glsl_ver);
    GLuint default_fs = compileRawShaderString(DEFAULT_FS, GL_FRAGMENT_SHADER, glsl_ver);

    default_program = glCreateProgram();
    if (!linkShaders(default_program, {default_vs, default_fs})) {
      std::cerr << "failed to link shaders for the default program" << std::endl;
      glDeleteProgram(default_program);
      default_program = 0;
      return false;
    }


#ifndef __ENSCRIPTEN__
    //TODO: enable a legacy path for opengl2.1 without ext_tranform_feedback?
    // This program is for rendering to pdf/ps
    if (GLEW_EXT_transform_feedback || GLEW_VERSION_3_0) {
        feedback_program = glCreateProgram();
        const char * xfrm_varyings[] = {
            "gl_Position",
            "fColor",
            "fClipCoord",
        };
        glTransformFeedbackVaryings(feedback_program, 3, xfrm_varyings, GL_INTERLEAVED_ATTRIBS);

        GLuint printing_vs = compileRawShaderString(PRINTING_VS, GL_VERTEX_SHADER, glsl_ver);
        GLuint printing_fs = compileRawShaderString(PRINTING_FS, GL_FRAGMENT_SHADER, glsl_ver);

        if (!linkShaders(feedback_program, {printing_vs, printing_fs})) {
          std::cerr << "failed to link shaders for the printing program" << std::endl;
          glDeleteProgram(feedback_program);
          feedback_program = 0;
          return false;
        }
    }
#endif

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
    //Texture unit 1: transparency palette
    //Texture unit 2: font atlas
    GLuint locColorTex = glGetUniformLocation(program, "colorTex");
    GLuint locAlphaTex = glGetUniformLocation(program, "alphaTex");
    glUniform1i(locColorTex, 0);
    glUniform1i(locAlphaTex, 1);
    // Set render type uniforms
    locContainsText = glGetUniformLocation(program, "containsText");
    if (_shaderMode == RENDER_COLOR) {
        glUniform1i(locContainsText, GL_FALSE);
    } else if (_shaderMode == RENDER_COLOR_TEX) {
        glUniform1i(locContainsText, GL_FALSE);
    } else { //_shaderMode == RENDER_TEXT
        glUniform1i(locContainsText, GL_TRUE);
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

void GlState::setModeRenderText()
{
   if (_shaderMode != RENDER_TEXT)
   {
      glDepthMask(GL_FALSE);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glUniform1i(locContainsText, GL_TRUE);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, 0);
      GetFont()->bindFontTex();
      _shaderMode = RENDER_TEXT;
   }
}
