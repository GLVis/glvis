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


private:
    GLuint _default_prgm;
    GLuint _global_vao;

    std::unordered_map<std::string, GLuint> _uniforms = {
        { "useClipPlane", 0 },
        { "clipPlane", 0 },
        { "containsText", 0 },
        { "modelViewMatrix", 0 },
        { "projectionMatrix", 0 },
        { "textProjMatrix", 0 },
        { "normalMatrix", 0 },
        { "num_lights", 0 },
        { "g_ambient", 0 },
        { "material.specular", 0 },
        { "material.shininess", 0 },
        { "lights[0].position", 0 },
        { "lights[0].diffuse", 0 },
        { "lights[0].specular", 0 },
        { "lights[1].position", 0 },
        { "lights[1].diffuse", 0 },
        { "lights[1].specular", 0 },
        { "lights[2].position", 0 },
        { "lights[2].diffuse", 0 },
        { "lights[2].specular", 0 },
        { "colorTex", 0 },
        { "alphaTex", 0 }
    };

    bool compileShaders();
    void initializeShaderState();
public:
    virtual DeviceType getType() { return GLDevice::CORE_DEVICE; }

    virtual void init();
    virtual void setTransformMatrices(glm::mat4 model_view, glm::mat4 projection);
    virtual void setNumLights(int i);
    virtual void setMaterial(Material mat);
    virtual void setPointLight(int i, Light lt);
    virtual void setAmbientLight(const std::array<float, 4>& amb);
    virtual void setClipPlaneUse(bool enable);
    virtual void setClipPlaneEqn(const std::array<double, 4>& eqn);

    virtual void bufferToDevice(array_layout layout, IVertexBuffer& buf);
    virtual void bufferToDevice(TextBuffer& t_buf);
    virtual void drawDeviceBuffer(array_layout layout, const IVertexBuffer& buf);
    virtual void drawDeviceBuffer(const TextBuffer& t_buf);
};

}
#endif // __RENDERER_CORE_HPP__