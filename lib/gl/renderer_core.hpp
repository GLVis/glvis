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

#ifndef GLVIS_RENDERER_CORE_HPP
#define GLVIS_RENDERER_CORE_HPP

#include <unordered_map>

#include "renderer.hpp"
#include "shader.hpp"

namespace gl3
{

// Renderer for OpenGL versions with access to the programmable pipeline
class CoreGLDevice : public GLDevice
{
public:
   enum ShaderAttrib
   {
      ATTR_VERTEX = 0,
      ATTR_TEXT_VERTEX,
      ATTR_NORMAL,
      ATTR_COLOR,
      ATTR_TEXCOORD0,
      NUM_ATTRS
   };

   struct ShaderXfbVertex
   {
      float pos[4];
      float color[4];
      float clipCoord;
   };

private:
   ShaderProgram default_prgm;
   ShaderProgram feedback_prgm;
   resource::VtxArrayHandle global_vao;

   resource::BufObjHandle feedback_vbo;

   const static std::vector<std::string> unif_list;

   std::unordered_map<std::string, GLuint> uniforms;

   bool use_clip_plane;

   struct VBOData
   {
      resource::BufObjHandle vert_buf;
      resource::BufObjHandle elem_buf;
      GLenum shape;
      size_t count;
      array_layout layout;
   };

   std::vector<VBOData> vbos;

   bool compileShaders();
   void initializeShaderState(const ShaderProgram& prog);

   template<typename T>
   void drawDeviceBufferImpl(GLenum shape, int count, bool indexed);

   void processTriangleXfbBuffer(CaptureBuffer& cbuf,
                                 const vector<ShaderXfbVertex>& verts);
   void processLineXfbBuffer(CaptureBuffer& cbuf,
                             const vector<ShaderXfbVertex>& verts);

public:
   CoreGLDevice()
      : global_vao(0)
   {
      vbos.emplace_back(VBOData{}); // dummy for index 0
   }

   DeviceType getType() override { return GLDevice::CORE_DEVICE; }

   void init() override;
   void setTransformMatrices(glm::mat4 model_view, glm::mat4 projection) override;
   void setNumLights(int i) override;
   void setMaterial(Material mat) override;
   void setPointLight(int i, Light lt) override;
   void setAmbientLight(const std::array<float, 4>& amb) override;
   void setClipPlaneUse(bool enable) override;
   void setClipPlaneEqn(const std::array<double, 4>& eqn) override;

   void bufferToDevice(array_layout layout, IVertexBuffer& buf) override;
   void bufferToDevice(array_layout layout, IIndexedBuffer& buf) override;
   void bufferToDevice(TextBuffer& t_buf) override;
   void drawDeviceBuffer(int hnd) override;
   void drawDeviceBuffer(const TextBuffer& t_buf) override;

   void initXfbMode() override
   {
      glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, feedback_vbo);
      initializeShaderState(feedback_prgm);
      glEnable(GL_RASTERIZER_DISCARD);
   }
   void exitXfbMode() override
   {
      glDisable(GL_RASTERIZER_DISCARD);
      initializeShaderState(default_prgm);
      glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);
   }
   void bindExternalProgram(const ShaderProgram& prog)
   {
      glDisable(GL_RASTERIZER_DISCARD);
      initializeShaderState(prog);
      glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);
   }
   void captureXfbBuffer(PaletteState& pal, CaptureBuffer& cbuf, int hnd) override;
};

}

#endif // GLVIS_RENDERER_CORE_HPP
