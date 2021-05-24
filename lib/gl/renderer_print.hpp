// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_RENDERER_PRINT_HPP
#define GLVIS_RENDERER_PRINT_HPP

#include "renderer.hpp"

namespace gl3
{

struct FeedbackVertex
{
   glm::vec3 position;
   glm::vec4 color;

   FeedbackVertex() = default;
   FeedbackVertex(glm::vec3 pos, glm::vec4 c)
      : position(pos), color(c) { }
   FeedbackVertex(glm::vec4 pos, glm::vec4 c)
      : position(pos), color(c) { }
};

struct FeedbackText
{
   glm::vec3 offset;
   glm::vec4 color;
   std::string text;

   FeedbackText() = default;
   FeedbackText(glm::vec3 t_off, glm::vec4 t_color, std::string txt)
      : offset(t_off), color(t_color), text(std::move(txt)) { }
};

struct CaptureBuffer
{
   vector<FeedbackVertex> lines;
   vector<FeedbackVertex> triangles;
   vector<FeedbackText> text;
};

class CapturePass : public IMainRenderPass
{
public:
   CapturePass() { }
   virtual bool Filter(const RenderParams& param) { return true; }

   virtual void PreRender() { device->initXfbMode(); }
   virtual void Render(const RenderQueue& queued);
   virtual void PostRender() { device->initRenderMode(); }

   CaptureBuffer GetLastCaptureBuffer()
   {
      CaptureBuffer b_mov = std::move(cbuf);
      cbuf = {};
      return b_mov;
   }

private:
   CaptureBuffer cbuf;
};

}

#endif // GLVIS_RENDERER_PRINT_HPP
