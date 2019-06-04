#ifndef __RENDERER_FF_HPP__
#define __RENDERER_FF_HPP__

#include "renderer.hpp"

namespace gl3
{
// Render for legacy OpenGL systems with access to only the fixed-function pipeline
class FFGLDevice : public GLDevice
{
public:

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
    void bufferToDevice(TextBuffer& t_buf) override;
    void drawDeviceBuffer(array_layout layout, const IVertexBuffer& buf) override;
    void drawDeviceBuffer(const TextBuffer& t_buf) override;
    void captureXfbBuffer(CaptureBuffer& cbuf,
                          array_layout layout,
                          const IVertexBuffer& buf) override;
};

}
#endif // __RENDERER_FF_HPP__
