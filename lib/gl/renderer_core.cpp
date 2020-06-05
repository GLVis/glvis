#include "attr_traits.hpp"
#include "renderer_core.hpp"
#include "../aux_vis.hpp"

#include <regex>
#include <type_traits>

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

const std::vector<std::string> CoreGLDevice::_unif_list =
{
   "useClipPlane",
   "clipPlane",
   "containsText",
   "modelViewMatrix",
   "projectionMatrix",
   "textProjMatrix",
   "normalMatrix",
   "num_lights",
   "g_ambient",
   "material.specular",
   "material.shininess",
   "lights[0].position",
   "lights[0].diffuse",
   "lights[0].specular",
   "lights[1].position",
   "lights[1].diffuse",
   "lights[1].specular",
   "lights[2].position",
   "lights[2].diffuse",
   "lights[2].specular",
   "colorTex",
   "alphaTex"
};

template<typename TVtx>
void setupVtxAttrLayout()
{
   static_assert(AttrCoord<TVtx>::exists,
                 "Invalid vertex type, requires at least TVtx::coord to be present.");
   AttrCoord<TVtx>::setup();
   AttrNormal<TVtx>::setup();
   AttrColor<TVtx>::setup();
   AttrTexcoord<TVtx>::setup();
}

template<typename TVtx>
void clearVtxAttrLayout()
{
   static_assert(AttrCoord<TVtx>::exists,
                 "Invalid vertex type, requires at least TVtx::coord to be present.");
   AttrCoord<TVtx>::clear();
   AttrNormal<TVtx>::clear();
   AttrColor<TVtx>::clear();
   AttrTexcoord<TVtx>::clear();
}

std::string formatShaderString(const std::string &shader_string,
                               GLenum shader_type, int glsl_version)
{
   // webgl does not allow name resolution across shaders, we have to substitute them in here
   // add lighting
   std::string formatted = std::regex_replace(shader_string,
                                              std::regex(R"(vec4 blinnPhong\(in vec3 pos, in vec3 norm, in vec4 color\);)"),
                                              BLINN_PHONG_FS);

   // add clip plane
   formatted = std::regex_replace(formatted,
                                  std::regex(R"(void fragmentClipPlane\(\);)"),
                                  CLIP_PLANE_FS);
   formatted = std::regex_replace(formatted,
                                  std::regex(R"(void setupClipPlane\(in float dist\);)"),
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
            formatted = std::regex_replace(formatted, std::regex("gl_FragColor"),
                                           "fragColor");
         }
         else if (glsl_version >= 330)
         {
            formatted = "layout(location = 0) out vec4 fragColor;\n" + formatted;
            formatted = std::regex_replace(formatted, std::regex("gl_FragColor"),
                                           "fragColor");
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
   formatted = std::regex_replace(GLSL_HEADER, std::regex("GLSL_VER"),
                                  std::to_string(glsl_version)) + formatted;

   return formatted;
}

GLuint compileShaderString(const std::string &shader_str, GLenum shader_type)
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

GLuint compileRawShaderString(const std::string &shader_string,
                              GLenum shader_type, int glsl_version)
{
   std::string formatted = formatShaderString(shader_string, shader_type,
                                              glsl_version);
   //std::cout << "compiling '''" << formatted << "...";
   GLuint shader_ref = compileShaderString(formatted, shader_type);
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
   std::string verStr = (char*)glGetString(GL_VERSION);
   int ver_major, ver_minor;
   int vs_idx = verStr.find_first_of(".");
   ver_major = std::stoi(verStr.substr(vs_idx - 1, vs_idx));
   ver_minor = std::stoi(verStr.substr(vs_idx + 1, 1));
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

   GLuint default_vs = compileRawShaderString(DEFAULT_VS, GL_VERTEX_SHADER,
                                              glsl_ver);
   GLuint default_fs = compileRawShaderString(DEFAULT_FS, GL_FRAGMENT_SHADER,
                                              glsl_ver);

   _default_prgm = glCreateProgram();
   if (!linkShaders(_default_prgm, {default_vs, default_fs}))
   {
      std::cerr << "Failed to link shaders for the default shader program" <<
                std::endl;
      glDeleteProgram(_default_prgm);
      _default_prgm = 0;
      return false;
   }

#ifndef __ENSCRIPTEN__
   //TODO: enable a legacy path for opengl2.1 without ext_tranform_feedback?
   // This program is for rendering to pdf/ps
   if (GLEW_EXT_transform_feedback || GLEW_VERSION_3_0)
   {
      _feedback_prgm = glCreateProgram();
      const char * xfrm_varyings[] =
      {
         "gl_Position",
         "fColor",
         "fClipCoord",
      };
      glTransformFeedbackVaryings(_feedback_prgm, 3, xfrm_varyings,
                                  GL_INTERLEAVED_ATTRIBS);

      GLuint printing_vs = compileRawShaderString(PRINTING_VS, GL_VERTEX_SHADER,
                                                  glsl_ver);
      GLuint printing_fs = compileRawShaderString(PRINTING_FS, GL_FRAGMENT_SHADER,
                                                  glsl_ver);

      if (!linkShaders(_feedback_prgm, {printing_vs, printing_fs}))
      {
         std::cerr << "failed to link shaders for the printing program" << std::endl;
         glDeleteProgram(_feedback_prgm);
         _feedback_prgm = 0;
         return false;
      }
   }
#endif

   return true;
}

void CoreGLDevice::initializeShaderState(CoreGLDevice::RenderMode mode)
{
   GLuint curr_prgm = 0;
   if (mode == RenderMode::Default)
   {
      glUseProgram(_default_prgm);
      curr_prgm = _default_prgm;
   }
   else if (mode == RenderMode::Feedback)
   {
      glUseProgram(_feedback_prgm);
      curr_prgm = _feedback_prgm;
   }
   _uniforms.clear();
   // initialize uniform map
   int num_unifs;
   glGetProgramiv(curr_prgm, GL_ACTIVE_UNIFORMS, &num_unifs);
   std::unordered_map<std::string, bool> uniform_query{};
   for (const auto& uf : _unif_list)
   {
      uniform_query.emplace(uf, false);
   }
   for (int i = 0; i < num_unifs; i++)
   {
      const size_t max_length = 100;
      char unif_name[max_length];
      GLsizei name_length;
      GLint gl_size;
      GLenum gl_type;
      glGetActiveUniform(curr_prgm, i, max_length, &name_length, &gl_size, &gl_type,
                         unif_name);
      auto it = uniform_query.find(std::string(unif_name, name_length));
      if (it == uniform_query.end())
      {
#ifdef GLVIS_DEBUG
         std::cerr << "Warning: unexpected uniform \""
                   << std::string(unif_name, name_length)
                   << "\" found in shader." << std::endl;
#endif
      }
      else
      {
         it->second = true;
         GLuint uf_idx = glGetUniformLocation(curr_prgm, it->first.c_str());
         _uniforms.emplace(it->first, uf_idx);
      }
   }
   for (const auto& uf : uniform_query)
   {
      if (!uf.second)
      {
#ifdef GLVIS_DEBUG
         std::cerr << "Uniform \"" << uf.first
                   << "\" missing in shader, ignoring." << std::endl;
#endif
         // set uniform index to -1 so glUniform ignores data
         _uniforms.emplace(uf.first, -1);
      }
   }
   glUniform1i(_uniforms["colorTex"], 0);
   glUniform1i(_uniforms["alphaTex"], 1);
   _use_clip_plane = false;
}

void CoreGLDevice::init()
{
   GLDevice::init();
   if (!this->compileShaders())
   {
      std::cerr << "Unable to initialize CoreGLDevice." << std::endl;
      return;
   }
   this->initializeShaderState(RenderMode::Default);
   if (GLEW_VERSION_3_0 || GLEW_ARB_vertex_array_object)
   {
      GLuint hnd_vao;
      glGenVertexArrays(1, &hnd_vao);
      _global_vao = VtxArrayHandle_(hnd_vao);
      glBindVertexArray(_global_vao);
   }
   GLuint hnd_fb_buf;
   glGenBuffers(1, &hnd_fb_buf);
   _feedback_vbo = BufObjHandle_(hnd_fb_buf);
}

void CoreGLDevice::setTransformMatrices(glm::mat4 model_view,
                                        glm::mat4 projection)
{
   GLDevice::setTransformMatrices(model_view, projection);
   glm::mat4 proj_text = glm::ortho<float>(0, _vp_width, 0, _vp_height, -5.0, 5.0);
   glm::mat3 inv_normal = glm::inverseTranspose(glm::mat3(model_view));
   glUniformMatrix4fv(_uniforms["modelViewMatrix"], 1, GL_FALSE,
                      glm::value_ptr(model_view));
   glUniformMatrix4fv(_uniforms["projectionMatrix"], 1, GL_FALSE,
                      glm::value_ptr(projection));
   glUniformMatrix4fv(_uniforms["textProjMatrix"], 1, GL_FALSE,
                      glm::value_ptr(proj_text));
   glUniformMatrix3fv(_uniforms["normalMatrix"], 1, GL_FALSE,
                      glm::value_ptr(inv_normal));
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
   glUniform4fv(_uniforms["material.specular"], 1, mat.specular.data());
   glUniform1f(_uniforms["material.shininess"], mat.shininess);
}

void CoreGLDevice::setPointLight(int i, Light lt)
{
   if (i > LIGHTS_MAX)
   {
      return;
   }
   std::string lt_index = "lights[" + std::to_string(i) + "]";
   glUniform3fv(_uniforms[lt_index + ".position"], 1, lt.position.data());
   glUniform4fv(_uniforms[lt_index + ".diffuse"], 1, lt.diffuse.data());
   glUniform4fv(_uniforms[lt_index + ".specular"], 1, lt.specular.data());
}

void CoreGLDevice::setAmbientLight(const std::array<float, 4> &amb)
{
   glUniform4fv(_uniforms["g_ambient"], 1, amb.data());
}

void CoreGLDevice::setClipPlaneUse(bool enable)
{
   _use_clip_plane = enable;
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
   if (buf.getHandle() == 0)
   {
      if (buf.count() == 0) { return; }
      GLuint handle;
      glGenBuffers(1, &handle);
      buf.setHandle(_vbos.size());
      _vbos.emplace_back(VBOData_{handle, 0, buf.getShape(), buf.count(), layout});
   }
   else
   {
      _vbos[buf.getHandle()].count = buf.count();
   }
   glBindBuffer(GL_ARRAY_BUFFER, _vbos[buf.getHandle()].vert_buf);
   glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
   glBufferData(GL_ARRAY_BUFFER, buf.count() * buf.getStride(),
                buf.getData(), GL_STATIC_DRAW);
}

void CoreGLDevice::bufferToDevice(array_layout layout, IIndexedBuffer& buf)
{
   if (buf.getHandle() == 0)
   {
      if (buf.count() == 0) { return; }
      GLuint handle[2];
      glGenBuffers(2, &handle[0]);
      buf.setHandle(_vbos.size());
      _vbos.emplace_back(VBOData_{handle[0], handle[1], buf.getShape(), buf.getIndices().size(), layout});
   }
   else
   {
      _vbos[buf.getHandle()].count = buf.getIndices().size();
   }
   // Buffer vertex array
   glBindBuffer(GL_ARRAY_BUFFER, _vbos[buf.getHandle()].vert_buf);
   glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
   glBufferData(GL_ARRAY_BUFFER, buf.count() * buf.getStride(),
                buf.getData(), GL_STATIC_DRAW);
   // Buffer element array
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbos[buf.getHandle()].elem_buf);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, buf.getIndices().size() * sizeof(int),
                buf.getIndices().data(), GL_STATIC_DRAW);
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
         float tris[] =
         {
            e.rx, e.ry, e.rz, cur_x, -cur_y, g.tex_x, 0, 0,
            e.rx, e.ry, e.rz, cur_x + g.w, -cur_y, g.tex_x + g.w / tex_w, 0, 0,
            e.rx, e.ry, e.rz, cur_x, -cur_y - g.h, g.tex_x, g.h / tex_h, 0,
            e.rx, e.ry, e.rz, cur_x + g.w, -cur_y, g.tex_x + g.w / tex_w, 0, 0,
            e.rx, e.ry, e.rz, cur_x, -cur_y - g.h, g.tex_x, g.h / tex_h, 0,
            e.rx, e.ry, e.rz, cur_x + g.w, -cur_y - g.h, g.tex_x + g.w / tex_w, g.h / tex_h, 0
         };
         buf_data.insert(buf_data.end(), tris, tris + 8 * 6);
      }
   }
   if (buf_data.size() == 0) { return; }
   if (t_buf.getHandle() == 0)
   {
      GLuint handle;
      glGenBuffers(1, &handle);
      t_buf.setHandle(handle);
   }
   glBindBuffer(GL_ARRAY_BUFFER, t_buf.getHandle());
   glBufferData(GL_ARRAY_BUFFER, sizeof(float) * buf_data.size(), buf_data.data(),
                GL_STATIC_DRAW);
}

template<typename T>
void CoreGLDevice::drawDeviceBufferImpl(GLenum shape, int count, bool indexed)
{
   if (!AttrNormal<T>::exists)
   {
      glVertexAttrib3f(CoreGLDevice::ATTR_NORMAL, 0.f, 0.f, 1.f);
   }
   if (!AttrColor<T>::exists && AttrTexcoord<T>::exists)
   {
      glVertexAttrib4f(CoreGLDevice::ATTR_COLOR, 1.f, 1.f, 1.f, 1.f);
   }
   setupVtxAttrLayout<T>();
   if (indexed)
   {
      glDrawElements(shape, count, GL_UNSIGNED_INT, (void*)0);
   }
   else
   {
      glDrawArrays(shape, 0, count);
   }
   clearVtxAttrLayout<T>();
}

void CoreGLDevice::drawDeviceBuffer(int hnd)
{
   if (hnd == 0) { return; }
   if (_vbos[hnd].count == 0) { return; }
   glBindBuffer(GL_ARRAY_BUFFER, _vbos[hnd].vert_buf);
   bool indexed = false;
   if (_vbos[hnd].elem_buf != 0)
   {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbos[hnd].elem_buf);
      indexed = true;
   }
   if (_vbos[hnd].layout == Vertex::layout
       || _vbos[hnd].layout == VertexNorm::layout)
   {
      glVertexAttrib4fv(ATTR_COLOR, _static_color.data());
   }
   GLenum shape = _vbos[hnd].shape;
   int count = _vbos[hnd].count;
   switch (_vbos[hnd].layout)
   {
      case Vertex::layout:
         drawDeviceBufferImpl<Vertex>(shape, count, indexed);
         break;
      case VertexColor::layout:
         drawDeviceBufferImpl<VertexColor>(shape, count, indexed);
         break;
      case VertexTex::layout:
         drawDeviceBufferImpl<VertexTex>(shape, count, indexed);
         break;
      case VertexNorm::layout:
         drawDeviceBufferImpl<VertexNorm>(shape, count, indexed);
         break;
      case VertexNormColor::layout:
         drawDeviceBufferImpl<VertexNormColor>(shape, count, indexed);
         break;
      case VertexNormTex::layout:
         drawDeviceBufferImpl<VertexNormTex>(shape, count, indexed);
         break;
      default:
         cerr << "WARNING: Unhandled vertex layout " << _vbos[hnd].layout << endl;
   }
}
void CoreGLDevice::drawDeviceBuffer(const TextBuffer& t_buf)
{
   if (t_buf.getHandle() == 0) { return; }
   if (t_buf.count() == 0) { return; }
   glUniform1i(_uniforms["containsText"], GL_TRUE);
   glEnableVertexAttribArray(ATTR_VERTEX);
   glEnableVertexAttribArray(ATTR_TEXT_VERTEX);
   glEnableVertexAttribArray(ATTR_TEXCOORD0);
   glBindBuffer(GL_ARRAY_BUFFER, t_buf.getHandle());

   glVertexAttrib4fv(ATTR_COLOR, _static_color.data());
   glVertexAttribPointer(ATTR_VERTEX, 3, GL_FLOAT, GL_FALSE, t_buf.getStride(), 0);
   glVertexAttribPointer(ATTR_TEXT_VERTEX, 2, GL_FLOAT, GL_FALSE,
                         t_buf.getStride(), (void*)(sizeof(float) * 3));
   glVertexAttribPointer(ATTR_TEXCOORD0, 2, GL_FLOAT, GL_FALSE, t_buf.getStride(),
                         (void*)(sizeof(float) * 5));
   glDrawArrays(GL_TRIANGLES, 0, t_buf.count());

   glDisableVertexAttribArray(ATTR_TEXT_VERTEX);
   glDisableVertexAttribArray(ATTR_TEXCOORD0);
   glUniform1i(_uniforms["containsText"], GL_FALSE);
}

#ifndef __EMSCRIPTEN__
inline FeedbackVertex XFBPostTransform(CoreGLDevice::ShaderXfbVertex v,
                                       float half_w, float half_h)
{
   glm::vec3 coord = glm::make_vec3(v.pos);
   glm::vec4 color = glm::make_vec4(v.color);
   // clip coords -> ndc
   coord /= v.pos[3];
   // ndc -> device coords
   coord.x = half_w * coord.x + half_w;
   coord.y = half_h * coord.y + half_h;
   return { coord, color };
}

void CoreGLDevice::processTriangleXfbBuffer(CaptureBuffer& cbuf,
                                            const vector<ShaderXfbVertex>& verts)
{
   float half_w = _vp_width * 0.5f;
   float half_h = _vp_height * 0.5f;
   if (!_use_clip_plane)
   {
      // all triangles into capture buf
      for (size_t i = 0; i < verts.size(); i++)
      {
         cbuf.triangles.emplace_back(XFBPostTransform(verts[i], half_w, half_h));
      }
      return;
   }
   // clipping needed
   for (size_t t_i = 0; t_i < verts.size() / 3; t_i++)
   {
      if (verts[t_i*3].clipCoord >= 0.f
          && verts[t_i*3+1].clipCoord >= 0.f
          && verts[t_i*3+2].clipCoord >= 0.f)
      {
         // triangle fully in the unclipped region
         for (int vert_i = 0; vert_i < 3; vert_i++)
         {
            cbuf.triangles.emplace_back(
               XFBPostTransform(verts[t_i*3 + vert_i], half_w, half_h));
         }
      }
      else if (verts[3*t_i].clipCoord < 0.f
               && verts[3*t_i+1].clipCoord < 0.f
               && verts[3*t_i+2].clipCoord < 0.f)
      {
         //triangle fully in clipped region
         continue;
      }
      else
      {
         //clip through middle of triangle
         for (int vert_i = 0; vert_i < 3; vert_i++)
         {
            int i_a = 3*t_i+vert_i;
            int i_b = 3*t_i+((vert_i+1) % 3);
            int i_c = 3*t_i+((vert_i+2) % 3);
            // find two points on the same side of clip plane
            if ((verts[i_a].clipCoord < 0.f) == (verts[i_b].clipCoord < 0.f))
            {
               //pts a, b are on same side of clip plane, c on other side
               //perspective-correct interpolation factors, needed for colors
               float c_w_a = verts[i_a].clipCoord / verts[i_a].pos[3];
               float c_w_b = verts[i_b].clipCoord / verts[i_b].pos[3];
               float c_w_c = verts[i_c].clipCoord / verts[i_c].pos[3];
               //compute clip pts
               glm::vec4 pos[2], color[2];
               // a --- n_0 --- c
               pos[0] = (glm::make_vec4(verts[i_a].pos) * verts[i_c].clipCoord
                         - glm::make_vec4(verts[i_c].pos) * verts[i_a].clipCoord);
               color[0] = (glm::make_vec4(verts[i_a].color) * c_w_c
                           - glm::make_vec4(verts[i_c].color) * c_w_a)
                          / (c_w_c - c_w_a);
               // b --- n_1 --- c
               pos[1] = (glm::make_vec4(verts[i_b].pos) * verts[i_c].clipCoord
                         - glm::make_vec4(verts[i_c].pos) * verts[i_b].clipCoord);
               color[1] = (glm::make_vec4(verts[i_b].color) * c_w_c
                           - glm::make_vec4(verts[i_c].color) * c_w_b)
                          / (c_w_c - c_w_b);
               for (int i = 0; i < 2; i++)
               {
                  // perform transform to device coords
                  pos[i] /= pos[i].w;
                  pos[i].x *= half_w; pos[i].x += half_w;
                  pos[i].y *= half_h; pos[i].y += half_h;
               }

               if (verts[i_c].clipCoord < 0.f)
               {
                  //pts a, b are in clip plane
                  //create quadrilateral a -- n_0 -- n_1 -- b
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_a], half_w, half_h));
                  cbuf.triangles.emplace_back(pos[0], color[0]);
                  cbuf.triangles.emplace_back(pos[1], color[1]);
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_a], half_w, half_h));
                  cbuf.triangles.emplace_back(pos[1], color[1]);
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_b], half_w, half_h));
               }
               else
               {
                  //pt c is in clip plane
                  //add triangle c -- n_0 -- n_1
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_c], half_w, half_h));
                  cbuf.triangles.emplace_back(pos[0], color[0]);
                  cbuf.triangles.emplace_back(pos[1], color[1]);
               }
               break;
            }
         }
      }
   }
}

void CoreGLDevice::processLineXfbBuffer(CaptureBuffer& cbuf,
                                        const vector<ShaderXfbVertex>& verts)
{
   float half_w = _vp_width * 0.5f;
   float half_h = _vp_height * 0.5f;
   for (size_t i = 0; i < verts.size(); i += 2)
   {
      if (!_use_clip_plane ||
          (verts[i].clipCoord >= 0.f && verts[i+1].clipCoord >= 0.f))
      {
         cbuf.lines.emplace_back(XFBPostTransform(verts[i], half_w, half_h));
         cbuf.lines.emplace_back(XFBPostTransform(verts[i+1], half_w, half_h));
      }
      else if (verts[i].clipCoord < 0.f && verts[i+1].clipCoord < 0.f)
      {
         //outside of clip plane;
         continue;
      }
      else
      {
         int i_a, i_b;
         if (verts[i].clipCoord < 0.f)
         {
            //vertex i lies in clipped region
            i_a = i+1;
            i_b = i;
         }
         else     //verts[i+1].clipCoord < 0.f
         {
            //vertex i+1 lies in clipped region
            i_a = i;
            i_b = i+1;
         }
         //compute new vertex (CbVa - CaVb), where Vb lies in the clipped region
         ShaderXfbVertex clip_vert;
         //perspective-correct interpolation factors for color
         float c_w_a = verts[i_a].clipCoord / verts[i_a].pos[3];
         float c_w_b = verts[i_b].clipCoord / verts[i_b].pos[3];
         for (int j = 0; j < 4; j++)
         {
            clip_vert.pos[j] = verts[i_a].pos[j] * verts[i_b].clipCoord
                               - verts[i_b].pos[j] * verts[i_a].clipCoord;
            clip_vert.color[j] = (verts[i_a].color[j] * c_w_b
                                  - verts[i_b].color[j] * c_w_a)
                                 / (c_w_b - c_w_a);
         }
         cbuf.lines.emplace_back(XFBPostTransform(clip_vert, half_w, half_h));
         cbuf.lines.emplace_back(XFBPostTransform(verts[i_a], half_w, half_h));
      }
   }
}


void CoreGLDevice::captureXfbBuffer(
   CaptureBuffer& cbuf, int hnd)
{
   // allocate feedback buffer
   int buf_size = _vbos[hnd].count * sizeof(ShaderXfbVertex);
   glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER,
                buf_size, nullptr, GL_STATIC_READ);
   glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, _feedback_vbo);
   // Draw objects in feedback-only mode
   glBeginTransformFeedback(_vbos[hnd].shape);
   drawDeviceBuffer(hnd);
   glEndTransformFeedback();
   // Read back feedback buffer
   vector<ShaderXfbVertex> xfbBuf(buf_size);
   glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, buf_size, xfbBuf.data());
   if (_vbos[hnd].shape == GL_TRIANGLES)
   {
      processTriangleXfbBuffer(cbuf, xfbBuf);
   }
   else if (_vbos[hnd].shape == GL_LINES)
   {
      processLineXfbBuffer(cbuf, xfbBuf);
   }
   else
   {
      std::cerr << "Warning: GL_POINTS handling not implemented in transform "
                << "feedback processing" << std::endl;
   }
}
#else
void CoreGLDevice::captureXfbBuffer(
   CaptureBuffer& cbuf, int hnd)
{
   std::cerr << "CoreGLDevice::captureXfbBuffer: "
             << "Not implemented for WebGL." << std::endl;
}
#endif // __EMSCRIPTEN__
}
