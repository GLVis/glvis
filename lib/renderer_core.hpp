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
    GLuint _default_prgm;
    GLuint _feedback_prgm;
    GLuint _global_vao;

    GLuint _feedback_vbo;

    enum class RenderMode
    {
        Default,
        Feedback
    };

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
    
    bool _use_clip_plane;

    bool compileShaders();
    void initializeShaderState(RenderMode mode);
    void processTriangleXfbBuffer(CaptureBuffer& cbuf, const vector<ShaderXfbVertex>& verts);
    void processLineXfbBuffer(CaptureBuffer& cbuf, const vector<ShaderXfbVertex>& verts);
public:
	CoreGLDevice()
		: _default_prgm(0), _feedback_prgm(0), _global_vao(0) { }

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
    void bufferToDevice(TextBuffer& t_buf) override;
    void drawDeviceBuffer(array_layout layout, const IVertexBuffer& buf) override;
    void drawDeviceBuffer(const TextBuffer& t_buf) override;

    void initXfbMode() override
    {
	    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _feedback_vbo);
        initializeShaderState(RenderMode::Feedback);
        glEnable(GL_RASTERIZER_DISCARD);
    }
    void exitXfbMode() override
    {
        glDisable(GL_RASTERIZER_DISCARD);
        initializeShaderState(RenderMode::Default);
	    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);
    }
    void captureXfbBuffer(CaptureBuffer& cbuf,
                          array_layout layout,
                          const IVertexBuffer& buf) override;
};

}
#endif // __RENDERER_CORE_HPP__
