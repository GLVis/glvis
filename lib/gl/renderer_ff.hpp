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

#ifndef GLVIS_RENDERER_FF_HPP
#define GLVIS_RENDERER_FF_HPP

#include "renderer.hpp"

namespace gl3
{

// Render for legacy OpenGL systems with access to only the fixed-function
// pipeline
class FFGLDevice : public GLDevice
{
   struct DispListData_
   {
      DispListHandle list;
      GLenum shape;
      size_t count;
      array_layout layout;
   };
   std::vector<DispListData_> disp_lists;

   template<typename TVtx>
   void bufferFFDeviceImpl(const VertexBuffer<TVtx>& buf);

   template<typename TVtx>
   void bufferFFDeviceImpl(const IndexedVertexBuffer<TVtx>& buf);
public:
   FFGLDevice()
   {
      disp_lists.emplace_back(DispListData_{}); // dummy for index 0
   }

   DeviceType getType() override { return GLDevice::FF_DEVICE; }

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
   void captureXfbBuffer(PaletteState& pal, CaptureBuffer& cbuf, int hnd) override;
};

}

#endif // GLVIS_RENDERER_FF_HPP
