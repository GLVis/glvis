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

#ifndef GLVIS_RENDERER_MSAA_HPP
#define GLVIS_RENDERER_MSAA_HPP

#include "renderer.hpp"

namespace gl3
{

class MultisamplePass : public IRenderPass
{
public:
   MultisamplePass() { }
   virtual void SetGLDevice(GLDevice* dev);
   virtual void PreRender();
   virtual void PostRender();

   void SetAntialiasing(bool aa_status);
   bool GetAntialiasing() { return msaa_enable; }
   void SetNumSamples(int samples)
   {
      msaa_samples = samples;
   }
   int GetNumSamples() { return msaa_samples; }

   void SetLineWidth(float w);
   float GetLineWidth() { return line_w; }
   void SetLineWidthMS(float w);
   float GetLineWidthMS() { return line_w_aa; }
private:
   bool feat_use_fbo_antialias;
   bool msaa_enable{false};
   int msaa_samples;
   float line_w, line_w_aa;

   RenderBufHandle renderBufs[2];
   FBOHandle msaaFb;
};

}

#endif // GLVIS_RENDERER_MSAA_HPP
