#ifndef __RENDERER_CORE_HPP__
#define __RENDERER_CORE_HPP__
#include "renderer.hpp"

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
   ShaderPrgmHandle default_prgm;
   ShaderPrgmHandle feedback_prgm;
   VtxArrayHandle global_vao;

   BufObjHandle feedback_vbo;

   enum class RenderMode
   {
      Default,
      Feedback
   };

   const static std::vector<std::string> unif_list;

   std::unordered_map<std::string, GLuint> uniforms;

   bool use_clip_plane;

   struct VBOData
   {
      BufObjHandle vert_buf;
      BufObjHandle elem_buf;
      GLenum shape;
      size_t count;
      array_layout layout;
   };

   std::vector<VBOData> vbos;

   bool compileShaders();
   void initializeShaderState(RenderMode mode);

   template<typename T>
   void drawDeviceBufferImpl(GLenum shape, int count, bool indexed);

   void processTriangleXfbBuffer(CaptureBuffer& cbuf,
                                 const vector<ShaderXfbVertex>& verts);
   void processLineXfbBuffer(CaptureBuffer& cbuf,
                             const vector<ShaderXfbVertex>& verts);
public:
   CoreGLDevice()
      : default_prgm(0), feedback_prgm(0), global_vao(0)
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
      initializeShaderState(RenderMode::Feedback);
      glEnable(GL_RASTERIZER_DISCARD);
   }
   void exitXfbMode() override
   {
      glDisable(GL_RASTERIZER_DISCARD);
      initializeShaderState(RenderMode::Default);
      glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);
   }
   void captureXfbBuffer(CaptureBuffer& cbuf, int hnd) override;
};

}
#endif // __RENDERER_CORE_HPP__
