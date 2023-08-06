// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_TYPES_HPP
#define GLVIS_TYPES_HPP

#include <vector>
#include <array>
#include <iostream>
#include <memory>
#include <iterator>
#include <algorithm>

#ifdef __INTEL_COMPILER
// for some reason, icpc is choking on GLM's constexpr stuff - force it to use
// only C++98 features as a workaround.
#define GLM_FORCE_CXX98
#endif
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "platform_gl.hpp"

using namespace std;

class VisualizationScene;

namespace gl3
{

namespace resource
{
// RAII scope guard for OpenGL handles.
template<void(*GLFinalizer)(GLuint)>
class Handle
{
   GLuint hnd;
public:
   Handle() : hnd{0} {}
   Handle(GLuint h) : hnd{h} {}
   ~Handle() { if (hnd) { GLFinalizer(hnd); } }
   Handle(Handle&& other)
      : hnd{other.hnd} { other.hnd = 0; }
   Handle& operator = (Handle&& other) noexcept
   {
      if (this != &other)
      {
         hnd = other.hnd;
         other.hnd = 0;
      }
      return *this;
   }
   operator GLuint() const { return hnd; }
};

inline void boCleanup(GLuint vbo_hnd)
{
   glDeleteBuffers(1, &vbo_hnd);
}

inline void dspListCleanup(GLuint dlist)
{
   glDeleteLists(dlist, 1);
}

inline void prgmCleanup(GLuint prgm)
{
   glDeleteProgram(prgm);
}

inline void shdrCleanup(GLuint shdr)
{
   glDeleteShader(shdr);
}

inline void vaoCleanup(GLuint vao)
{
   glDeleteVertexArrays(1, &vao);
}

inline void texCleanup(GLuint tex)
{
   glDeleteTextures(1, &tex);
}

inline void fboCleanup(GLuint fbo)
{
   glDeleteFramebuffers(1, &fbo);
}

inline void rboCleanup(GLuint rbo)
{
   glDeleteRenderbuffers(1, &rbo);
}

using BufObjHandle = Handle<boCleanup>;
using DispListHandle = Handle<dspListCleanup>;
using VtxArrayHandle = Handle<vaoCleanup>;
using ShaderPrgmHandle = Handle<prgmCleanup>;
using ShaderHandle = Handle<shdrCleanup>;
using TextureHandle = Handle<texCleanup>;
using FBOHandle = Handle<fboCleanup>;
using RenderBufHandle = Handle<rboCleanup>;

} // end namespace resource

struct GlMatrix
{
   glm::mat4 mtx;

   /// Applies a rotation transform to the matrix.
   void rotate(float angle, double x, double y, double z)
   {
      mtx = glm::rotate(mtx, glm::radians(angle), glm::vec3(x,y,z));
   }

   void mult(glm::mat4 rhs)
   {
      mtx = mtx * rhs;
   }

   /// Applies a translation transform to the matrix.
   void translate(double x, double y, double z)
   {
      mtx = glm::translate(mtx, glm::vec3(x, y, z));
   }

   /// Applies a scale transform to the matrix.
   void scale(double x, double y, double z)
   {
      mtx = glm::scale(mtx, glm::vec3(x, y, z));
   }

   /// Sets the matrix to an orthographic projection.
   void ortho(double left,
              double right,
              double bottom,
              double top,
              double z_near,
              double z_far)
   {
      mtx = glm::ortho(left, right, bottom, top, z_near, z_far);
   }

   /// Sets the matrix to a perspective projection.
   void perspective(double fov, double aspect, double z_near, double z_far)
   {
      mtx = glm::perspective(glm::radians(fov), aspect, z_near, z_far);
   }

   /// Sets the matrix to the identity matrix.
   void identity()
   {
      mtx = glm::mat4(1.0);
   }
};

enum array_layout
{
   LAYOUT_VTX = 0,
   LAYOUT_VTX_NORMAL,
   LAYOUT_VTX_COLOR,
   LAYOUT_VTX_TEXTURE0,
   LAYOUT_VTX_NORMAL_COLOR,
   LAYOUT_VTX_NORMAL_TEXTURE0,
   NUM_LAYOUTS
};

inline std::array<uint8_t, 4> ColorU8(float rgba[])
{
   return
   {
      (rgba[0] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[0] * 256.),
      (rgba[1] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[1] * 256.),
      (rgba[2] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[2] * 256.),
      (rgba[3] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[3] * 256.),
   };
}

inline std::array<uint8_t, 4> ColorU8(float r, float g, float b, float a)
{
   return
   {
      (r >= 1.0) ? (uint8_t) 255 : (uint8_t)(r * 256.),
      (g >= 1.0) ? (uint8_t) 255 : (uint8_t)(g * 256.),
      (b >= 1.0) ? (uint8_t) 255 : (uint8_t)(b * 256.),
      (a >= 1.0) ? (uint8_t) 255 : (uint8_t)(a * 256.),
   };
}

struct alignas(16) Vertex
{
   std::array<float, 3> coord;

   static Vertex create(double * d)
   {
      return Vertex {(float) d[0], (float) d[1], (float) d[2]};
   }
   static const int layout = LAYOUT_VTX;
};

struct alignas(16) VertexColor
{
   std::array<float, 3> coord;
   std::array<uint8_t, 4> color;

   static const int layout = LAYOUT_VTX_COLOR;
};

struct alignas(16) VertexTex
{
   std::array<float, 3> coord;
   std::array<float, 2> texCoord;

   static const int layout = LAYOUT_VTX_TEXTURE0;
};

struct alignas(16) VertexNorm
{
   std::array<float, 3> coord;
   std::array<float, 3> norm;

   static const int layout = LAYOUT_VTX_NORMAL;
};

struct alignas(16) VertexNormColor
{
   std::array<float, 3> coord;
   std::array<float, 3> norm;
   std::array<uint8_t, 4> color;

   static const int layout = LAYOUT_VTX_NORMAL_COLOR;
};

struct alignas(16) VertexNormTex
{
   std::array<float, 3> coord;
   std::array<float, 3> norm;
   std::array<float, 2> texCoord;

   static const int layout = LAYOUT_VTX_NORMAL_TEXTURE0;
};



class GlDrawable;

/// Crude fixed-function OpenGL emulation helper.
class GlBuilder
{
   GlDrawable * parent_buf;
   GLenum render_as;
   int count;

   bool is_line;

   bool use_norm;
   bool use_color;
   bool use_tex;

   struct FFState
   {
      std::array<float, 3> coords;
      std::array<float, 3> norm;
      std::array<uint8_t, 4> color;
      std::array<float, 2> texcoord;
   };

   FFState saved[3];
   FFState curr;

   void saveVertex(const FFState& v);

public:
   GlBuilder(GlDrawable * buf)
      : parent_buf(buf)
      , count(0)
      , is_line(false)
      , use_norm(false)
      , use_color(false)
      , use_tex(false) { }

   void glBegin(GLenum e)
   {
      if (e == GL_LINES || e == GL_LINE_STRIP || e == GL_LINE_LOOP)
      {
         is_line = true;
      }
      else
      {
         is_line = false;
      }
      render_as = e;
      count = 0;
   }

   void glEnd()
   {
      // create degenerate primitives if necessary
      if (render_as == GL_LINES && count % 2 != 0)
      {
         saveVertex(curr);
      }
      if (render_as == GL_TRIANGLES)
      {
         for (int i = 0; i < count % 3; i++)
         {
            saveVertex(curr);
         }
      }
      if (render_as == GL_LINE_LOOP && count > 2)
      {
         // link first and last vertex
         saveVertex(saved[0]);
         saveVertex(saved[1]);
      }
      count = 0;
   }

   void glVertex3d(double x, double y, double z)
   {
      curr.coords = { (float) x, (float) y, (float) z };
      if (render_as == GL_LINES || render_as == GL_TRIANGLES)
      {
         // Lines and triangles are stored as-is
         saveVertex(curr);
      }
      else if (is_line)
      {
         // LINE_LOOP and LINE_STRIP: each vertex creates a new line with the
         // previous vertex
         if (count == 0)
         {
            saved[0] = curr;
         }
         else
         {
            saveVertex(saved[1]);
            saveVertex(curr);
         }
         saved[1] = curr;
      }
      else if (render_as == GL_QUADS)
      {
         if (count % 4 == 3)
         {
            // split on 0-2 diagonal
            saveVertex(saved[0]);
            saveVertex(saved[1]);
            saveVertex(saved[2]);
            saveVertex(saved[0]);
            saveVertex(saved[2]);
            saveVertex(curr);
            count -= 4;
         }
         else
         {
            saved[count] = curr;
         }
      }
      else if (render_as != GL_NONE)
      {
         // TriangleStrip, TriangleFan, Polygon
         if (count >= 2)
         {
            saveVertex(saved[0]);
            saveVertex(saved[1]);
            saveVertex(curr);
            if (render_as == GL_TRIANGLE_STRIP)
            {
               // pop earliest vertex
               saved[0] = saved[1];
            }
            saved[1] = curr;
         }
         else
         {
            saved[count] = curr;
         }
      }
      count++;
   }

   void glVertex3dv(const double * d)
   {
      glVertex3d(d[0], d[1], d[2]);
   }

   void glNormal3d(double nx, double ny, double nz)
   {
      use_norm = true;
      curr.norm = { (float) nx, (float) ny, (float) nz };
   }

   void glNormal3dv(const double * d) { glNormal3d(d[0], d[1], d[2]); }

   void glColor4f(float r, float g, float b, float a)
   {
      if (count == 0)
      {
         use_color = true;
         use_tex = false;
      }
      curr.color =
      {
         (r >= 1.0) ? (uint8_t) 255 : (uint8_t)(r * 256.),
         (g >= 1.0) ? (uint8_t) 255 : (uint8_t)(g * 256.),
         (b >= 1.0) ? (uint8_t) 255 : (uint8_t)(b * 256.),
         (a >= 1.0) ? (uint8_t) 255 : (uint8_t)(a * 256.),
      };
   }

   void glColor3f(float r, float g, float b) { glColor4f(r, g, b, 1.f); }

   void glColor4fv(float * cv) { glColor4f(cv[0], cv[1], cv[2], cv[3]); }

   void glTexCoord2f(float coord_u, float coord_v)
   {
      if (count == 0)
      {
         use_tex = true;
         use_color = false;
      }
      curr.texcoord = { coord_u, coord_v };
   }
};

class IVertexBuffer
{
private:
   int handle;
public:
   IVertexBuffer() : handle(0) { }
   virtual ~IVertexBuffer() { }

   int getHandle() const { return handle; }
   void setHandle(int dev_hnd) { handle = dev_hnd; }
   /// Clears the data stored in the vertex buffer.
   virtual void clear() = 0;
   /// Gets the number of vertices contained in the buffer.
   virtual size_t count() const = 0;
   /// Gets the primitive type contained by the vertex buffer.
   virtual GLenum getShape() const = 0;
   /// Gets the stride between vertices.
   virtual size_t getStride() const = 0;

   virtual const void* getData() const = 0;
};

template<typename T>
class VertexBuffer : public IVertexBuffer
{
private:
   GLenum primitive;
   std::vector<T> vertex_data;

public:
   typedef typename std::vector<T>::const_iterator ConstIterator;

   VertexBuffer(GLenum shape) : primitive(shape) { }
   ~VertexBuffer() { }

   virtual void clear() { vertex_data.clear(); }
   virtual size_t count() const { return vertex_data.size(); }

   virtual GLenum getShape() const { return primitive; }
   virtual size_t getStride() const { return sizeof(T); }

   ConstIterator begin() const { return vertex_data.begin(); };
   ConstIterator end() const { return vertex_data.end(); };

   virtual const void* getData() const { return vertex_data.data(); }

   /// Add a vertex to the buffer.
   void addVertex(const T& vertex)
   {
      vertex_data.emplace_back(vertex);
   }

   /// Add vertices to a buffer.
   void addVertices(const std::vector<T>& verts)
   {
      vertex_data.insert(vertex_data.end(), verts.begin(), verts.end());
   }
};

class IIndexedBuffer : public IVertexBuffer
{
public:
   virtual const std::vector<int>& getIndices() const = 0;
};

template<typename T>
class IndexedVertexBuffer : public IIndexedBuffer
{
private:
   GLenum primitive;
   std::vector<T> vertex_data;
   std::vector<int> vertex_indices;

public:
   typedef typename std::vector<T>::const_iterator ConstIterator;
   IndexedVertexBuffer(GLenum shape) : primitive(shape) { }
   ~IndexedVertexBuffer() { }
   virtual void clear()
   {
      vertex_data.clear();
      vertex_indices.clear();
   }

   virtual size_t count() const { return vertex_data.size(); }

   virtual GLenum getShape() const { return primitive; }
   virtual size_t getStride() const { return sizeof(T); }

   ConstIterator begin() const { return vertex_data.begin(); };
   ConstIterator end() const { return vertex_data.end(); };

   virtual const void* getData() const { return vertex_data.data(); }

   const std::vector<int>& getIndices() const { return vertex_indices; }

   void addVertices(const std::vector<T>& verts, const std::vector<int>& ids)
   {
      int index_offset = vertex_data.size();
      std::copy(verts.begin(), verts.end(), std::back_inserter(vertex_data));
      int old_end = vertex_indices.size();
      std::copy(ids.begin(), ids.end(), std::back_inserter(vertex_indices));
      for (size_t i = old_end; i < vertex_indices.size(); i++)
      {
         vertex_indices[i] += index_offset;
      }
   }
};

class TextBuffer : public IVertexBuffer
{
public:
   struct Entry
   {
      float rx, ry, rz;
      int ox, oy;
      std::string text;
      Entry() = default;
      Entry(float x, float y, float z, int ox, int oy, const std::string& txt)
         : rx(x), ry(y), rz(z), ox(ox), oy(oy), text(txt) { }
   };
   typedef std::vector<Entry>::iterator Iterator;
   typedef std::vector<Entry>::const_iterator ConstIterator;
private:
   std::vector<Entry> vertex_data;
   size_t num_chars;

public:
   TextBuffer() : num_chars(0) { }
   ~TextBuffer() { }

   /// Adds a text element at the specified local space (pre-transform)
   /// coordinates.
   void addText(float x, float y, float z, int ox, int oy,
                const std::string& text)
   {
      vertex_data.emplace_back(x, y, z, ox, oy, text);
      num_chars += text.length();
   }

   /// Gets an iterator referencing the first text entry.
   Iterator begin() { return vertex_data.begin(); }
   ConstIterator begin() const { return vertex_data.begin(); }

   /// Gets an iterator pointing to the last text entry.
   Iterator end() { return vertex_data.end(); }
   ConstIterator end() const { return vertex_data.end(); }

   virtual void clear() { vertex_data.clear(); num_chars = 0; }
   virtual size_t count() const { return num_chars * 6; };
   virtual GLenum getShape() const { return GL_TRIANGLES; };
   virtual size_t getStride() const { return sizeof(float) * 8; };
   virtual const void* getData() const { return nullptr; }
};

class GlDrawable
{
private:
   const static size_t NUM_SHAPES = 2;
   std::unique_ptr<IVertexBuffer> buffers[NUM_LAYOUTS][NUM_SHAPES];
   std::unique_ptr<IIndexedBuffer> indexed_buffers[NUM_LAYOUTS][NUM_SHAPES];
   TextBuffer text_buffer;

   friend class GlBuilder;
   friend class MeshRenderer;
   friend class ::VisualizationScene; // needed for glTF export

   template<typename Vert>
   VertexBuffer<Vert> * getBuffer(GLenum shape)
   {
      int idx = -1;
      if (shape == GL_LINES)
      {
         idx = 0;
      }
      else if (shape == GL_TRIANGLES)
      {
         idx = 1;
      }
      else
      {
         return nullptr;
      }
      if (!buffers[Vert::layout][idx])
      {
         buffers[Vert::layout][idx].reset(new VertexBuffer<Vert>(shape));
      }
      VertexBuffer<Vert> * buf = static_cast<VertexBuffer<Vert>*>
                                 (buffers[Vert::layout][idx].get());
      return buf;
   }

   template<typename Vert>
   IndexedVertexBuffer<Vert> * getIndexedBuffer(GLenum shape)
   {
      int idx = -1;
      if (shape == GL_LINES)
      {
         idx = 0;
      }
      else if (shape == GL_TRIANGLES)
      {
         idx = 1;
      }
      else
      {
         return nullptr;
      }
      if (!indexed_buffers[Vert::layout][idx])
      {
         indexed_buffers[Vert::layout][idx]
         .reset(new IndexedVertexBuffer<Vert>(shape));
      }
      return static_cast<IndexedVertexBuffer<Vert>*>
             (indexed_buffers[Vert::layout][idx].get());
   }
public:

   /// Adds a string at the given position in object coordinates.
   void addText(float x, float y, float z, const std::string& text)
   {
      text_buffer.addText(x, y, z, 0, 0, text);
   }

   /** @brief Adds a string at the given position in object coordinates with an
       offset in pixels. */
   void addText(float x, float y, float z, int ox, int oy,
                const std::string& text)
   {
      text_buffer.addText(x, y, z, ox, oy, text);
   }

   template<typename Vert>
   void addLine(const Vert& v1, const Vert& v2)
   {
      getBuffer<Vert>(GL_LINES)->addVertex(v1);
      getBuffer<Vert>(GL_LINES)->addVertex(v2);
   }

   template<typename Vert>
   void addLines(const std::vector<Vert>& verts)
   {
      getBuffer<Vert>(GL_LINES)->addVertices(verts);
   }

   template<typename Vert>
   void addTriangle(const Vert& v1, const Vert& v2, const Vert& v3)
   {
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v2);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
   }

   template<typename Vert>
   void addQuad(const Vert& v1, const Vert& v2, const Vert& v3, const Vert& v4)
   {
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v2);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
      getBuffer<Vert>(GL_TRIANGLES)->addVertex(v4);
   }

   template<typename Vert>
   void addTriangleIndexed(const std::vector<Vert>& verts,
                           const std::vector<int>& inds)
   {
      getIndexedBuffer<Vert>(GL_TRIANGLES)->addVertices(verts, inds);
   }

   template<typename Vert>
   void addQuadIndexed(const std::vector<Vert>& verts,
                       const std::vector<int>& inds)
   {
      std::vector<int> new_inds;
      for (size_t i = 0; i < inds.size() / 4; i++)
      {
         new_inds.emplace_back(inds[4*i]);
         new_inds.emplace_back(inds[4*i + 1]);
         new_inds.emplace_back(inds[4*i + 2]);
         new_inds.emplace_back(inds[4*i]);
         new_inds.emplace_back(inds[4*i + 2]);
         new_inds.emplace_back(inds[4*i + 3]);
      }
      getIndexedBuffer<Vert>(GL_TRIANGLES)->addVertices(verts, new_inds);
   }

   void addCone(float x, float y, float z,
                float vx, float vy, float vz,
                float cone_scale = 0.075);

   GlBuilder createBuilder()
   {
      return GlBuilder(this);
   }

   /// Clears the drawable object.
   void clear()
   {
      for (int i = 0; i < NUM_LAYOUTS; i++)
      {
         for (size_t j = 0; j < NUM_SHAPES; j++)
         {
            if (buffers[i][j])
            {
               buffers[i][j]->clear();
            }
            if (indexed_buffers[i][j])
            {
               indexed_buffers[i][j]->clear();
            }
         }
      }
      text_buffer.clear();
   }
};

}

#endif // GLVIS_TYPES_HPP
