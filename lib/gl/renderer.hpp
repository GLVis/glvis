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

#ifndef GLVIS_RENDERER_HPP
#define GLVIS_RENDERER_HPP

#include <memory>
#include <vector>
#include <set>
#include <unordered_map>

#include "platform_gl.hpp"
#include "types.hpp"
#include "../material.hpp"
#include "../palettes.hpp"

namespace gl3
{

using namespace resource;

const int LIGHTS_MAX = 3;
#ifdef GLVIS_MS_LINEWIDTH
const float LINE_WIDTH_AA = GLVIS_MS_LINEWIDTH;
#else
const float LINE_WIDTH_AA = 1.4;
#endif

struct RenderParams
{
   // Transformation matrices
   GlMatrix model_view;
   GlMatrix projection;

   // Lighting settings
   Material mesh_material;
   int num_pt_lights;
   std::array<Light, LIGHTS_MAX> lights;
   std::array<float, 4> light_amb_scene;
   std::array<float, 4> static_color;

   // Clip plane parameters
   bool use_clip_plane;
   std::array<double, 4> clip_plane_eqn;

   // If true, batch contains translucent drawables
   bool contains_translucent;
};

typedef vector<pair<RenderParams, GlDrawable*>> RenderQueue;

struct CaptureBuffer;

struct SceneInfo
{
   vector<GlDrawable*> needs_buffering;
   RenderQueue queue;
};

// OpenGL device interface representing rendering capabilities
class GLDevice
{
protected:
   int vp_width;
   int vp_height;
   glm::mat4 model_view_mtx;
   glm::mat4 proj_mtx;

   std::array<float, 4> static_color;
   std::array<float, 4> clear_color { 1.f, 1.f, 1.f, 1.f };

   bool blend_enabled = false;

protected:
   TextureHandle passthrough_texture;

public:

   enum DeviceType
   {
      NO_DEVICE,
      FF_DEVICE,
      CORE_DEVICE
   };

   virtual ~GLDevice() = default;
   const static int SAMPLER_COLOR = 0;
   const static int SAMPLER_ALPHA = 1;

   void detachTexture(int tex_unit)
   {
      glActiveTexture(GL_TEXTURE0 + tex_unit);
      glBindTexture(GL_TEXTURE_2D, passthrough_texture);
   }
   void attachTexture(int tex_unit, int tex_id)
   {
      glActiveTexture(GL_TEXTURE0 + tex_unit);
      glBindTexture(GL_TEXTURE_2D, tex_id);
   };

   void attachFramebuffer(const FBOHandle& fbo);

   // If true, use unsized internal formats and GL_ALPHA for single-channel
   // data. Otherwise, use the newer sized internal formats and GL_RED.
   static bool useLegacyTextureFmts();

   void enableBlend() { blend_enabled = true; glEnable(GL_BLEND); }
   void disableBlend() { blend_enabled = false; glDisable(GL_BLEND); }
   bool isBlendEnabled() { return blend_enabled; }
   void enableDepthWrite() { glDepthMask(GL_TRUE); }
   void disableDepthWrite() { glDepthMask(GL_FALSE); }
   void setLineWidth(float w) { glLineWidth(w); }

   virtual void init();
   virtual DeviceType getType() = 0;
   // Sets the window viewport.
   void setViewport(GLsizei w, GLsizei h);
   // Gets the current window viewport.
   void getViewport(GLint (&vp)[4]);
   // Set the color to use, if a color attribute is not provided.
   void setStaticColor(const std::array<float, 4>& rgba) { static_color = rgba; }
   // Set the background clear color to use.
   void setClearColor(const std::array<float, 4>& rgba) { clear_color = rgba; }
   // Gets the current background color.
   std::array<float, 4> getClearColor() const { return clear_color; }

   // === Render pipeline functions ===

   // Set the current transform matrices.
   virtual void setTransformMatrices(glm::mat4 model_view, glm::mat4 projection);
   // Set the number of lights to use. Setting number of lights to 0 disables
   // lighting.
   virtual void setNumLights(int i) = 0;
   // Set the parameters to use for the mesh material.
   virtual void setMaterial(Material mat) = 0;
   // Set the array of parameters to use for each light.
   virtual void setPointLight(int i, Light lt) = 0;
   // Set the color value of the global ambient light.
   virtual void setAmbientLight(const std::array<float, 4>& amb) = 0;
   // Set whether to enable or disable the clip plane.
   virtual void setClipPlaneUse(bool enable) = 0;
   // Set the equation to use for the clip plane.
   virtual void setClipPlaneEqn(const std::array<double, 4>& eqn) = 0;

   // === Buffer management functions ===

   // Load a client-side vertex buffer into a device buffer.
   virtual void bufferToDevice(ArrayLayout layout, IVertexBuffer& buf) = 0;
   virtual void bufferToDevice(ArrayLayout layout, IIndexedBuffer& buf) = 0;
   virtual void bufferToDevice(TextBuffer& t_buf) = 0;
   // Draw the data loaded in a device buffer.
   virtual void drawDeviceBuffer(int hnd) = 0;
   virtual void drawDeviceBuffer(const TextBuffer& t_buf) = 0;

   // === Transform feedback functions ===

   // Initializes state needed for transform feedback.
   virtual void initXfbMode() {}
   // Prepares state when exiting transform feedback.
   virtual void exitXfbMode() {}
   // Capture the next drawn vertex buffer to a feedback buffer instead of
   // drawing to screen.
   virtual void captureXfbBuffer(PaletteState& pal, CaptureBuffer& capture,
                                 int hnd) = 0;
   // Capture the next text buffer instead of drawing to screen.
   void captureXfbBuffer(CaptureBuffer& capture, const TextBuffer& t_buf);

};

class IRenderPass
{
public:
   IRenderPass() { }
   virtual void SetGLDevice(GLDevice* dev) { device = dev; }
   virtual void PreRender() = 0;
   virtual void PostRender() = 0;

   virtual const FBOHandle& GetSourceFramebuffer() const
   { return default_target; }

   virtual void SetTargetFramebuffer(const FBOHandle& fbo) { target = &fbo; }
protected:
   GLDevice* device;
   const FBOHandle* target = nullptr;
   const FBOHandle default_target{0};
};

// Generic interface for an action on a queue of renderable objects.
class IMainRenderPass : public IRenderPass
{
public:
   IMainRenderPass() { }
   // If the class will handle objects with the given set of render
   // parameters, returns true.
   virtual bool Filter(const RenderParams& param) = 0;
   // Renders objects in the queue to the setup output surfaces.
   virtual void Render(const RenderQueue& queue) = 0;
   // Sets the palette state for this render pass.
   void setPalette(PaletteState& pal) { palette = &pal; }
   // Sets the texture handle of the font atlas.
   void setFontTexture(GLuint tex_h) { font_tex = tex_h; }
protected:
   bool use_default_output = true;
   PaletteState* palette = nullptr;
   GLuint font_tex = 0;
};

class DefaultPass : public IMainRenderPass
{
public:
   DefaultPass() { }
   virtual bool Filter(const RenderParams& param) { return true; }

   virtual void PreRender() {}
   virtual void Render(const RenderQueue& queued);
   virtual void PostRender() {}
};

class MeshRenderer
{
   unique_ptr<GLDevice> device;
   float line_w;

public:
   MeshRenderer()
      : line_w(1.f) { }

   template<typename TDevice>
   void setDevice()
   {
      device.reset(new TDevice());
      device->setLineWidth(line_w);
      device->init();
   }

   template<typename TDevice>
   void setDevice(TDevice&& device)
   {
      device.reset(new TDevice(device));
   }

   void SetLineWidth(float w)
   {
      line_w = w;
      if (device) { device->setLineWidth(w); }
   }
   float GetLineWidth() { return line_w; }

   void setClearColor(float r, float g, float b, float a)
   { device->setClearColor({r,g,b,a}); }
   void setViewport(GLsizei w, GLsizei h) { device->setViewport(w, h); }

   void render(const std::vector<IMainRenderPass*>& main_passes,
               const std::vector<IRenderPass*>& extra_passes,
               const RenderQueue& queued);

   void buffer(GlDrawable* buf);
};

}

#endif // GLVIS_RENDERER_HPP
