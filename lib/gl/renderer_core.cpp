// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include <unordered_map>

#include "attr_traits.hpp"
#include "renderer_core.hpp"
#include "../aux_vis.hpp"

#ifdef GLVIS_DEBUG
#include <unordered_set>
#endif

// weird but loads them inline

const std::string BLINN_PHONG_FS =
#include "shaders/lighting.glsl"
   ;
const std::string DEFAULT_VS =
#include "shaders/default.vert"
   ;
const std::string DEFAULT_FS =
   BLINN_PHONG_FS +
#include "shaders/default.frag"
   ;
const std::string PRINTING_VS =
   BLINN_PHONG_FS +
#include "shaders/printing.vert"
   ;
const std::string PRINTING_FS =
#include "shaders/printing.frag"
   ;

namespace gl3
{

const std::vector<std::string> CoreGLDevice::unif_list =
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

bool CoreGLDevice::compileShaders()
{
   std::unordered_map<int, std::string> attribMap =
   {
      { CoreGLDevice::ATTR_VERTEX, "vertex"},
      { CoreGLDevice::ATTR_TEXT_VERTEX, "textVertex"},
      { CoreGLDevice::ATTR_NORMAL, "normal"},
      { CoreGLDevice::ATTR_COLOR, "color"},
      { CoreGLDevice::ATTR_TEXCOORD0, "texCoord0"}
   };

   if (!default_prgm.create(DEFAULT_VS, DEFAULT_FS, attribMap, 1))
   {
      std::cerr << "Failed to create the default shader program." <<
                std::endl;
      return false;
   }

#ifndef __EMSCRIPTEN__
   if (GLEW_EXT_transform_feedback || GLEW_VERSION_3_0)
   {
      const char * xfrm_varyings[] =
      {
         "gl_Position",
         "fColor",
         "fClipCoord",
      };
      glTransformFeedbackVaryings(feedback_prgm.getProgramId(), 3, xfrm_varyings,
                                  GL_INTERLEAVED_ATTRIBS);

      if (!feedback_prgm.create(PRINTING_VS, PRINTING_FS, attribMap, 1))
      {
         std::cerr << "Failed to create the printing capture program." << std::endl;
         return false;
      }
   }
#endif

   return true;
}

void CoreGLDevice::initializeShaderState(const ShaderProgram& prog)
{
   prog.bind();
   uniforms = prog.getUniformMap();
   for (const auto& uf : unif_list)
   {
      if (uniforms.find(uf) == uniforms.end())
      {
#ifdef GLVIS_DEBUG
         std::cerr << "Uniform \"" << uf
                   << "\" missing in shader, ignoring." << std::endl;
#endif
         // set uniform index to -1 so glUniform ignores data
         uniforms.emplace(uf, -1);
      }
   }
#ifdef GLVIS_DEBUG
   std::unordered_set<string> expectedUnifs(unif_list.begin(), unif_list.end());
   for (const auto& pairunif : uniforms)
   {
      if (expectedUnifs.find(pairunif.first) == expectedUnifs.end())
      {
         std::cerr << "Warning: unexpected uniform \""
                   << pairunif.first
                   << "\" found in shader." << std::endl;
      }
   }
#endif
   glUniform1i(uniforms["colorTex"], 0);
   glUniform1i(uniforms["alphaTex"], 1);
   use_clip_plane = false;
}

void CoreGLDevice::init()
{
   GLDevice::init();
   if (!this->compileShaders())
   {
      std::cerr << "Unable to initialize CoreGLDevice." << std::endl;
      return;
   }
   this->initializeShaderState(default_prgm);
   if (GLEW_VERSION_3_0 || GLEW_ARB_vertex_array_object)
   {
      GLuint hnd_vao;
      glGenVertexArrays(1, &hnd_vao);
      global_vao = resource::VtxArrayHandle(hnd_vao);
      glBindVertexArray(global_vao);
   }
   GLuint hnd_fb_buf;
   glGenBuffers(1, &hnd_fb_buf);
   feedback_vbo = resource::BufObjHandle(hnd_fb_buf);
}

void CoreGLDevice::setTransformMatrices(glm::mat4 model_view,
                                        glm::mat4 projection)
{
   GLDevice::setTransformMatrices(model_view, projection);
   glm::mat4 proj_text = glm::ortho<float>(0, vp_width, 0, vp_height, -5.0, 5.0);
   glm::mat3 inv_normal = glm::inverseTranspose(glm::mat3(model_view));
   glUniformMatrix4fv(uniforms["modelViewMatrix"], 1, GL_FALSE,
                      glm::value_ptr(model_view));
   glUniformMatrix4fv(uniforms["projectionMatrix"], 1, GL_FALSE,
                      glm::value_ptr(projection));
   glUniformMatrix4fv(uniforms["textProjMatrix"], 1, GL_FALSE,
                      glm::value_ptr(proj_text));
   glUniformMatrix3fv(uniforms["normalMatrix"], 1, GL_FALSE,
                      glm::value_ptr(inv_normal));
}

void CoreGLDevice::setNumLights(int i)
{
   if (i > LIGHTS_MAX)
   {
      return;
   }
   glUniform1i(uniforms["num_lights"], i);
}

void CoreGLDevice::setMaterial(Material mat)
{
   glUniform4fv(uniforms["material.specular"], 1, mat.specular.data());
   glUniform1f(uniforms["material.shininess"], mat.shininess);
}

void CoreGLDevice::setPointLight(int i, Light lt)
{
   if (i > LIGHTS_MAX)
   {
      return;
   }
   std::string lt_index = "lights[" + std::to_string(i) + "]";
   glUniform4fv(uniforms[lt_index + ".position"], 1, lt.position.data());
   glUniform4fv(uniforms[lt_index + ".diffuse"], 1, lt.diffuse.data());
   glUniform4fv(uniforms[lt_index + ".specular"], 1, lt.specular.data());
}

void CoreGLDevice::setAmbientLight(const std::array<float, 4> &amb)
{
   glUniform4fv(uniforms["g_ambient"], 1, amb.data());
}

void CoreGLDevice::setClipPlaneUse(bool enable)
{
   use_clip_plane = enable;
   glUniform1i(uniforms["useClipPlane"], enable);
}

void CoreGLDevice::setClipPlaneEqn(const std::array<double, 4> &eqn)
{
   glm::vec4 clip_plane(eqn[0], eqn[1], eqn[2], eqn[3]);
   clip_plane = glm::inverseTranspose(model_view_mtx) * clip_plane;
   glUniform4fv(uniforms["clipPlane"], 1, glm::value_ptr(clip_plane));
}

void CoreGLDevice::bufferToDevice(array_layout layout, IVertexBuffer &buf)
{
   if (buf.getHandle() == 0)
   {
      if (buf.count() == 0) { return; }
      GLuint handle;
      glGenBuffers(1, &handle);
      buf.setHandle(vbos.size());
      vbos.emplace_back(VBOData{handle, 0, buf.getShape(), buf.count(), layout});
   }
   else
   {
      vbos[buf.getHandle()].count = buf.count();
   }
   glBindBuffer(GL_ARRAY_BUFFER, vbos[buf.getHandle()].vert_buf);
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
      buf.setHandle(vbos.size());
      vbos.emplace_back(VBOData{handle[0], handle[1], buf.getShape(), buf.getIndices().size(), layout});
   }
   else
   {
      vbos[buf.getHandle()].count = buf.getIndices().size();
   }
   // Buffer vertex array
   glBindBuffer(GL_ARRAY_BUFFER, vbos[buf.getHandle()].vert_buf);
   glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
   glBufferData(GL_ARRAY_BUFFER, buf.count() * buf.getStride(),
                buf.getData(), GL_STATIC_DRAW);
   // Buffer element array
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos[buf.getHandle()].elem_buf);
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
      float pen_x = e.ox, pen_y = e.oy;
      char prev_c = '\0';
      for (char c : e.text)
      {
         const GlVisFont::glyph &g = GetFont()->GetTexChar(c);
         pen_x += GetFont()->GetKerning(prev_c, c);
         // note: subtract 1 to account for the padding in the texture glyphs
         float cur_x = pen_x + g.bear_x - 1;
         float cur_y = -pen_y - g.bear_y - 1;
         pen_x += g.adv_x;
         pen_y += g.adv_y;
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
         prev_c = c;
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
   if (vbos[hnd].count == 0) { return; }
   glBindBuffer(GL_ARRAY_BUFFER, vbos[hnd].vert_buf);
   bool indexed = false;
   if (vbos[hnd].elem_buf != 0)
   {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos[hnd].elem_buf);
      indexed = true;
   }
   if (vbos[hnd].layout == Vertex::layout
       || vbos[hnd].layout == VertexNorm::layout)
   {
      glVertexAttrib4fv(ATTR_COLOR, static_color.data());
   }
   GLenum shape = vbos[hnd].shape;
   int count = vbos[hnd].count;
   switch (vbos[hnd].layout)
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
         cerr << "WARNING: Unhandled vertex layout " << vbos[hnd].layout << endl;
   }
}
void CoreGLDevice::drawDeviceBuffer(const TextBuffer& t_buf)
{
   if (t_buf.getHandle() == 0) { return; }
   if (t_buf.count() == 0) { return; }
   glUniform1i(uniforms["containsText"], GL_TRUE);
   glEnableVertexAttribArray(ATTR_VERTEX);
   glEnableVertexAttribArray(ATTR_TEXT_VERTEX);
   glEnableVertexAttribArray(ATTR_TEXCOORD0);
   glBindBuffer(GL_ARRAY_BUFFER, t_buf.getHandle());

   glVertexAttrib4fv(ATTR_COLOR, static_color.data());
   glVertexAttribPointer(ATTR_VERTEX, 3, GL_FLOAT, GL_FALSE, t_buf.getStride(), 0);
   glVertexAttribPointer(ATTR_TEXT_VERTEX, 2, GL_FLOAT, GL_FALSE,
                         t_buf.getStride(), (void*)(sizeof(float) * 3));
   glVertexAttribPointer(ATTR_TEXCOORD0, 2, GL_FLOAT, GL_FALSE, t_buf.getStride(),
                         (void*)(sizeof(float) * 5));
   glDrawArrays(GL_TRIANGLES, 0, t_buf.count());

   glDisableVertexAttribArray(ATTR_TEXT_VERTEX);
   glDisableVertexAttribArray(ATTR_TEXCOORD0);
   glUniform1i(uniforms["containsText"], GL_FALSE);
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
   float half_w = vp_width * 0.5f;
   float half_h = vp_height * 0.5f;
   if (!use_clip_plane)
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
         // triangle fully in clipped region
         continue;
      }
      else
      {
         // clip through middle of triangle
         for (int vert_i = 0; vert_i < 3; vert_i++)
         {
            int i_a = 3*t_i+vert_i;
            int i_b = 3*t_i+((vert_i+1) % 3);
            int i_c = 3*t_i+((vert_i+2) % 3);
            // find two points on the same side of clip plane
            if ((verts[i_a].clipCoord < 0.f) == (verts[i_b].clipCoord < 0.f))
            {
               // pts a, b are on same side of clip plane, c on other side
               // perspective-correct interpolation factors, needed for colors
               float c_w_a = verts[i_a].clipCoord / verts[i_a].pos[3];
               float c_w_b = verts[i_b].clipCoord / verts[i_b].pos[3];
               float c_w_c = verts[i_c].clipCoord / verts[i_c].pos[3];
               // compute clip pts
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
                  // pts a, b are in clip plane
                  // create quadrilateral a -- n_0 -- n_1 -- b
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_a], half_w, half_h));
                  cbuf.triangles.emplace_back(pos[0], color[0]);
                  cbuf.triangles.emplace_back(pos[1], color[1]);
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_a], half_w, half_h));
                  cbuf.triangles.emplace_back(pos[1], color[1]);
                  cbuf.triangles.emplace_back(XFBPostTransform(verts[i_b], half_w, half_h));
               }
               else
               {
                  // pt c is in clip plane
                  // add triangle c -- n_0 -- n_1
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
   float half_w = vp_width * 0.5f;
   float half_h = vp_height * 0.5f;
   for (size_t i = 0; i < verts.size(); i += 2)
   {
      if (!use_clip_plane ||
          (verts[i].clipCoord >= 0.f && verts[i+1].clipCoord >= 0.f))
      {
         cbuf.lines.emplace_back(XFBPostTransform(verts[i], half_w, half_h));
         cbuf.lines.emplace_back(XFBPostTransform(verts[i+1], half_w, half_h));
      }
      else if (verts[i].clipCoord < 0.f && verts[i+1].clipCoord < 0.f)
      {
         // outside of clip plane;
         continue;
      }
      else
      {
         int i_a, i_b;
         if (verts[i].clipCoord < 0.f)
         {
            // vertex i lies in clipped region
            i_a = i+1;
            i_b = i;
         }
         else     // verts[i+1].clipCoord < 0.f
         {
            // vertex i+1 lies in clipped region
            i_a = i;
            i_b = i+1;
         }
         // compute new vertex (CbVa - CaVb), where Vb lies in the clipped region
         ShaderXfbVertex clip_vert;
         // perspective-correct interpolation factors for color
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
   PaletteState& pal, CaptureBuffer& cbuf, int hnd)
{
   if (hnd == 0) { return; }
   if (vbos[hnd].count == 0) { return; }
   // allocate feedback buffer
   int buf_size = vbos[hnd].count * sizeof(ShaderXfbVertex);
   glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER,
                buf_size, nullptr, GL_STATIC_READ);
   glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, feedback_vbo);
   // Draw objects in feedback-only mode
   glBeginTransformFeedback(vbos[hnd].shape);
   drawDeviceBuffer(hnd);
   glEndTransformFeedback();
   // Read back feedback buffer
   vector<ShaderXfbVertex> xfbBuf(buf_size);
   glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, buf_size, xfbBuf.data());
   if (vbos[hnd].shape == GL_TRIANGLES)
   {
      processTriangleXfbBuffer(cbuf, xfbBuf);
   }
   else if (vbos[hnd].shape == GL_LINES)
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

void CoreGLDevice::captureXfbBuffer(PaletteState & /*unused*/,
                                    CaptureBuffer& /*unused*/, int /*unused*/)
{
   std::cerr << "CoreGLDevice::captureXfbBuffer: "
             << "Not implemented for WebGL." << std::endl;
}

#endif // __EMSCRIPTEN__

}
