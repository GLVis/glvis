#include "renderer.hpp"
#include "aux_vis.hpp"

#include <regex>

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

namespace gl3
{

std::string formatShaderString(const std::string &shader_string, GLenum shader_type, int glsl_version)
{
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
    if (glsl_version >= 130)
    {
        if (shader_type == GL_VERTEX_SHADER)
        {
            formatted = std::regex_replace(formatted, std::regex("attribute"), "in");
            formatted = std::regex_replace(formatted, std::regex("varying"), "out");
        }

        else if (shader_type == GL_FRAGMENT_SHADER)
        {
            formatted = std::regex_replace(formatted, std::regex("varying"), "in");

            // requires GL_ARB_explicit_attrib_location extension or GLSL 3.30
            // although gl_FragColor was depricated in GLSL 1.3
            if (glsl_version > 130 && glsl_version < 330)
            {
                formatted = "out vec4 fragColor;\n" + formatted;
                formatted = std::regex_replace(formatted, std::regex("gl_FragColor"), "fragColor");
            }
            else if (glsl_version >= 330)
            {
                formatted = "layout(location = 0) out vec4 fragColor;\n" + formatted;
                formatted = std::regex_replace(formatted, std::regex("gl_FragColor"), "fragColor");
            }
        }

        else
        {
            std::cerr << "buildShaderString: unknown shader type" << std::endl;
            return {};
        }

        formatted = std::regex_replace(formatted, std::regex("texture2D"), "texture");
    }

    // add the header
    formatted = std::regex_replace(GLSL_HEADER, std::regex("GLSL_VER"), std::to_string(glsl_version)) + formatted;

    return formatted;
}

GLuint compileShaderString(const std::string &shader_str, GLenum shader_type, int glsl_version)
{
    int shader_len = shader_str.length();
    const char *shader_cstr = shader_str.c_str();

    GLuint shader = glCreateShader(shader_type);
    glShaderSource(shader, 1, &shader_cstr, &shader_len);
    glCompileShader(shader);

    GLint stat;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &stat);
    // glGetObjectParameteriv
    if (stat == GL_FALSE)
    {
        std::cerr << "failed to compile shader" << std::endl;
        int err_len;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &err_len);
        char *error_text = new char[err_len];
        glGetShaderInfoLog(shader, err_len, &err_len, error_text);
        std::cerr << error_text << std::endl;
        delete[] error_text;
        return 0;
    }
    return shader;
}

GLuint compileRawShaderString(const std::string &shader_string, GLenum shader_type, int glsl_version)
{
    std::string formatted = formatShaderString(shader_string, shader_type, glsl_version);
    //std::cout << "compiling '''" << formatted << "...";
    GLuint shader_ref = compileShaderString(formatted, shader_type, glsl_version);
    //if (shader_ref != 0) { std::cerr << "successfully compiled shader" << std::endl; }
    return shader_ref;
}

bool linkShaders(GLuint prgm, const std::vector<GLuint> &shaders)
{
    // explicitly specify attrib positions so we don't have to reset VAO
    // bindings when switching programs

    // for OSX, attrib 0 must be bound to render an object
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_VERTEX, "vertex");
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_TEXT_VERTEX, "textVertex");
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_NORMAL, "normal");
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_COLOR, "color");
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_TEXCOORD0, "texCoord0");
    glBindAttribLocation(prgm, CoreGLDevice::ATTR_TEXCOORD1, "texCoord1");

    for (GLuint i : shaders)
    {
        glAttachShader(prgm, i);
    }
    glLinkProgram(prgm);

    GLint stat;
    glGetProgramiv(prgm, GL_LINK_STATUS, &stat);
    if (stat == GL_FALSE)
    {
        cerr << "fatal: Shader linking failed" << endl;
    }
    return (stat == GL_TRUE);
}

bool CoreGLDevice::compileShaders()
{
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

    if (opengl_ver < 330)
    {
        if (ver_major == 2)
        {
            glsl_ver = opengl_ver - 90;
        }
        else if (ver_major == 3)
        {
            glsl_ver = opengl_ver - 170;
        }
        else
        {
            std::cerr << "fatal: unsupported OpenGL version " << opengl_ver << std::endl;
            return false;
        }
    }
    else
    {
        glsl_ver = opengl_ver;
    }
#else
    // GL Version | WebGL Version | GLSL Version
    //    2.0     |      1.0      |   1.00 ES
    //    3.0     |      2.0      |   3.00 ES
    //    3.1     |               |   3.10 ES
    if (opengl_ver < 300)
    {
        glsl_ver = 100;
    }
    else
    {
        glsl_ver = 300;
    }
#endif
    std::cerr << "Using GLSL " << glsl_ver << std::endl;

    GLuint default_vs = compileRawShaderString(DEFAULT_VS, GL_VERTEX_SHADER, glsl_ver);
    GLuint default_fs = compileRawShaderString(DEFAULT_FS, GL_FRAGMENT_SHADER, glsl_ver);

    _default_prgm = glCreateProgram();
    if (!linkShaders(_default_prgm, {default_vs, default_fs}))
    {
        std::cerr << "Failed to link shaders for the default shader program" << std::endl;
        glDeleteProgram(_default_prgm);
        _default_prgm = 0;
        return false;
    }
}

void CoreGLDevice::initializeShaderState()
{
#ifdef GLVIS_DEBUG
    // verify that uniform map is consisted with current shaders
    int num_attribs;
    glGetProgramiv(_default_prgm, GL_ACTIVE_UNIFORMS, &num_attribs);
    if (num_attribs != _uniforms.size())
    {
        std::cerr << "Warning: Unexpected number of uniforms in shader.\n"
                  << "Expected " << _uniforms.size() << " uniforms, got "
                  << num_attribs << std::endl;
    }

#endif
    for (auto &uf : _uniforms)
    {
        uf.second = glGetUniformLocation(_default_prgm, uf.first.c_str());
    }
    glUniform1i(_uniforms["colorTex"], 0);
    glUniform1i(_uniforms["alphaTex"], 1);
    glUniform1i(_uniforms["fontTex"], 1);
}

void CoreGLDevice::init()
{
    if (!this->compileShaders())
    {
        std::cerr << "Unable to initialize CoreGLDevice." << std::endl;
        return;
    }
    glUseProgram(_default_prgm);
    this->initializeShaderState();
    if (GLEW_VERSION_3_0)
    {
        if (_global_vao == 0)
        {
            glGenVertexArrays(1, &_global_vao);
        }
        glBindVertexArray(_global_vao);
    }
}

void CoreGLDevice::setTransformMatrices(glm::mat4 model_view, glm::mat4 projection)
{
    _model_view_mtx = model_view;
    glm::mat4 proj_text = glm::ortho<float>(0, _vp_width, 0, _vp_height, -5.0, 5.0);
    glm::mat3 inv_normal = glm::inverseTranspose(glm::mat3(model_view));
    glUniformMatrix4fv(_uniforms["modelViewMatrix"], 1, GL_FALSE, glm::value_ptr(model_view));
    glUniformMatrix4fv(_uniforms["projectionMatrix"], 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(_uniforms["textProjMatrix"], 1, GL_FALSE, glm::value_ptr(proj_text));
    glUniformMatrix3fv(_uniforms["normalMatrix"], 1, GL_FALSE, glm::value_ptr(inv_normal));
}

void CoreGLDevice::setNumLights(int i)
{
    if (i > LIGHTS_MAX)
    {
        return;
    }
    glUniform1i(_uniforms["num_lights"], i);
}

void CoreGLDevice::setMaterial(Material mat)
{
    glUniform4fv(_uniforms["material.specular"], 1, mat.specular);
    glUniform1f(_uniforms["material.shininess"], mat.shininess);
}

void CoreGLDevice::setPointLight(int i, Light lt)
{
    if (i > LIGHTS_MAX)
    {
        return;
    }
    std::string lt_index = "lights[" + std::to_string(i) + "]";
    glUniform3fv(_uniforms[lt_index + ".position"], 1, lt.position);
    glUniform4fv(_uniforms[lt_index + ".diffuse"], 1, lt.diffuse);
    glUniform4fv(_uniforms[lt_index + ".specular"], 1, lt.specular);
}

void CoreGLDevice::setAmbientLight(const std::array<float, 4> &amb)
{
    glUniform4fv(_uniforms["g_ambient"], 1, amb.data());
}

void CoreGLDevice::setClipPlaneUse(bool enable)
{
    glUniform1i(_uniforms["useClipPlane"], enable);
}

void CoreGLDevice::setClipPlaneEqn(const std::array<double, 4> &eqn)
{
    glm::vec4 clip_plane(eqn[0], eqn[1], eqn[2], eqn[3]);
    clip_plane = glm::inverseTranspose(_model_view_mtx) * clip_plane;
    glUniform4fv(_uniforms["clipPlane"], 1, glm::value_ptr(clip_plane));
}

void CoreGLDevice::bufferToDevice(array_layout layout, IVertexBuffer &buf)
{
    if (buf.count() == 0)
    {
        return;
    }
    if (buf.get_handle() == 0)
    {
        GLuint handle;
        glGenBuffers(1, &handle);
        buf.set_handle(handle);
    }
    glBindBuffer(GL_ARRAY_BUFFER, buf.get_handle());
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, buf.count() * buf.get_stride(),
                 buf.data_begin(), GL_STATIC_DRAW);
}

void CoreGLDevice::bufferToDevice(TextBuffer &t_buf)
{
    std::vector<float> buf_data;
    float tex_w = GetFont()->getAtlasWidth();
    float tex_h = GetFont()->getAtlasHeight();
    for (auto &e : t_buf)
    {
        float x = 0.f, y = 0.f;
        for (char c : e.text)
        {
            GlVisFont::glyph g = GetFont()->GetTexChar(c);
            float cur_x = x + g.bear_x;
            float cur_y = -y - g.bear_y;
            x += g.adv_x;
            y += g.adv_y;
            if (!g.w || !g.h)
            {
                continue;
            }
            float tris[] = {
                e.rx, e.ry, e.rz, cur_x, -cur_y, g.tex_x, 0, 0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y, g.tex_x + g.w / tex_w, 0, 0,
                e.rx, e.ry, e.rz, cur_x, -cur_y - g.h, g.tex_x, g.h / tex_h, 0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y, g.tex_x + g.w / tex_w, 0, 0,
                e.rx, e.ry, e.rz, cur_x, -cur_y - g.h, g.tex_x, g.h / tex_h, 0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y - g.h, g.tex_x + g.w / tex_w, g.h / tex_h, 0};
            buf_data.insert(buf_data.end(), tris, tris + 8 * 6);
        }
    }
    if (t_buf.get_handle() == 0)
    {
        GLuint handle;
        glGenBuffers(1, &handle);
        t_buf.set_handle(handle);
    }
    glBindBuffer(GL_ARRAY_BUFFER, t_buf.get_handle());
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * buf_data.size(), buf_data.data(), GL_STATIC_DRAW);
}

void CoreGLDevice::drawDeviceBuffer(array_layout layout, const IVertexBuffer &buf)
{
  if (buf.get_handle() == 0) { return; }
  if (buf.count() == 0) { return; }
  switch (layout) {
    case Vertex::layout:
    {
      Vertex::setupAttribLayout();
    }
    break;
    case VertexNorm::layout:
    {
      VertexNorm::
    }
  }
  if (layout == Vertex::layout || layout == VertexNorm::layout) {
    glVertexAttrib4fv(ATTR_COLOR, _static_color.data());
  }
  if (!(layout == VertexNorm::layout
      || layout == VertexNormColor::layout
      || layout == VertexNormTex::layout)) {
    glVertexAttrib3f(ATTR_NORMAL, 0.f, 0.f, 1.f);
  }
  glBindBuffer(GL_ARRAY_BUFFER, buf.get_handle());
  glDrawArrays(buf.get_shape(), 0, buf.count());
}
void CoreGLDevice::drawDeviceBuffer(const TextBuffer& t_buf)
{
}

}