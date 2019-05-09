#ifndef __RENDERER_HPP__
#define __RENDERER_HPP__

#include <memory>
#include <set>
#include <unordered_map>

#include "platform_gl.hpp"
#include "aux_gl3.hpp"
#include "material.hpp"

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

struct GlMatrix
{
   glm::mat4 mtx;

   /**
    * Applies a rotation transform to the matrix.
    */
   void rotate(float angle, double x, double y, double z)
   {
      mtx = glm::rotate(mtx, glm::radians(angle), glm::vec3(x,y,z));
   }

   void mult(glm::mat4 rhs)
   {
      mtx = mtx * rhs;
   }

   /**
    * Applies a translation transform to the matrix.
    */
   void translate(double x, double y, double z)
   {
      mtx = glm::translate(mtx, glm::vec3(x, y, z));
   }

   /**
    * Applies a scale transform to the matrix.
    */
   void scale(double x, double y, double z)
   {
      mtx = glm::scale(mtx, glm::vec3(x, y, z));
   }

   /**
    * Sets the matrix to an orthographic projection.
    */
   void ortho(double left,
              double right,
              double bottom,
              double top,
              double z_near,
              double z_far)
   {
      mtx = glm::ortho(left, right, bottom, top, z_near, z_far);
   }

   /**
    * Sets the matrix to a perspective projection.
    */
   void perspective(double fov, double aspect, double z_near, double z_far)
   {
      mtx = glm::perspective(glm::radians(fov), aspect, z_near, z_far);
   }

   /**
    * Sets the matrix to the identity matrix.
    */
   void identity()
   {
      mtx = glm::mat4(1.0);
   }
};
namespace gl3
{

const int LIGHTS_MAX = 3;
#ifdef GLVIS_MS_LINEWIDTH
const float LINE_WIDTH_AA = GLVIS_MS_LINEWIDTH
#else
const float LINE_WIDTH_AA = 1.4;
#endif

struct RenderParams
{
    //Transformation matrices
    GlMatrix model_view;
    GlMatrix projection;

    //Lighting settings
    Material mesh_material;
    int num_pt_lights;
    std::array<Light, LIGHTS_MAX> lights;
    std::array<float, 4> light_amb_scene;
    std::array<float, 4> static_color;

    //Clip plane params
    bool use_clip_plane;
    std::array<double, 4> clip_plane_eqn;

    //If true, batch contains translucent drawables
    bool contains_translucent;
};

typedef vector<pair<RenderParams, GlDrawable*>> RenderQueue;

struct SceneInfo
{
    vector<GlDrawable*> needs_buffering;
    RenderQueue queue;
};

// OpenGL device interface representing rendering capabilities
class GLDevice
{
protected:
    int _vp_width;
    int _vp_height;
    glm::mat4 _model_view_mtx;
    glm::mat4 _proj_mtx;

    std::array<float, 4> _static_color;
public:

    void detachTexture(int tex_unit) {
        glActiveTexture(GL_TEXTURE0 + tex_unit);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    void attachTexture(int tex_unit, int tex_id) {
        glActiveTexture(GL_TEXTURE0 + tex_unit);
        glBindTexture(GL_TEXTURE_2D, tex_id);
    };

    void enableBlend() { glEnable(GL_BLEND); }
    void disableBlend() { glDisable(GL_BLEND); }
    void enableDepthWrite() { glDepthMask(GL_TRUE); }
    void disableDepthWrite() { glDepthMask(GL_FALSE); }
    void enableMultisample() { glEnable(GL_MULTISAMPLE); }
    void disableMultisample() { glDisable(GL_MULTISAMPLE); }
    void setLineWidth(float w) { glLineWidth(w); }

    virtual void init();
    // Sets the window viewport.
    void setViewport(GLsizei w, GLsizei h);
    // Gets the current window viewport.
    void getViewport(GLint (&vp)[4]);
    // Set the color to use, if a color attribute is not provided.
    void setStaticColor(const std::array<float, 4>& rgba) { _static_color = rgba; }

    // === Render pipeline functions ===
    // Set the current transform matrices.
    virtual void setTransformMatrices(glm::mat4 model_view, glm::mat4 projection);
    // Set the number of lights to use. Setting number of lights to 0 disables lighting.
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
    virtual void bufferToDevice(array_layout layout, IVertexBuffer& buf) = 0;
    virtual void bufferToDevice(TextBuffer& t_buf) = 0;
    // Draw the data loaded in a device buffer.
    virtual void drawDeviceBuffer(array_layout layout, const IVertexBuffer& buf) = 0;
    virtual void drawDeviceBuffer(const TextBuffer& t_buf) = 0;

    virtual void preRender() = 0;
};

class MeshRenderer
{
    unique_ptr<GLDevice> _device;
    bool _msaa_enable;
    int _color_tex, _alpha_tex, _font_tex;
    float _line_w, _line_w_aa;
public:
    template<typename TDevice>
    void setDevice() {
        _device.reset(new TDevice());
        _device->init();
    }

    template<typename TDevice>
    void setDevice(TDevice&& device) {
        _device.reset(new TDevice(device));
    }

    void setAntialiasing(bool aa_status);
    bool getAntialiasing() { return _msaa_enable; }

    void setLineWidth(float w);
    float getLineWidth() { return _line_w; }
    void setLineWidthMS(float w);
    float getLineWidthMS() { return _line_w_aa; }

    void setClearColor(float r, float g, float b, float a) { glClearColor(r, g, b, a); }
	void setViewport(GLsizei w, GLsizei h) { glViewport(0, 0, w, h); }

    void init();
    void render(const RenderQueue& queued);

    void buffer(GlDrawable* buf);
};


// Render for legacy OpenGL systems with access to only the fixed-function pipeline
class FFGLDevice : public GLDevice
{
public:
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

    virtual void preRender();
};

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
        ATTR_TEXCOORD1,
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
        { "alphaTex", 0 },
        { "fontTex", 0 }
    };

    bool compileShaders();
    void initializeShaderState();
public:
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

    virtual void preRender();
};

}

#endif // __RENDERER_HPP__
