#include "renderer.hpp"
#include "aux_vis.hpp"

namespace gl3
{

void MeshRenderer::render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for (auto& q_elem : _queued) {
        const MeshRenderingParams& params = q_elem.first;
        _device->setTransformMatrices(params.model_view.mtx, params.projection.mtx);
        _device->setMaterial(params.mesh_material);
        _device->setNumLights(params.num_pt_lights);
        for (int i = 0; i < LIGHTS_MAX; i++) {
            _device->setPointLight(i, params.lights[i]);
        }
        _device->setAmbientLight(params.light_amb_scene);
        _device->setStaticColor(params.static_color);
        _device->setClipPlaneUse(params.use_clip_plane);
        _device->setClipPlaneEqn(params.clip_plane_eqn);
        //aggregate buffers with common parameters
        std::vector<pair<array_layout, IVertexBuffer*>> texture_bufs, no_texture_bufs;
        std::vector<TextBuffer*> text_bufs;
        for (GlDrawable* batch_elem : q_elem.second) {
            for (int i = 0; i < NUM_LAYOUTS; i++) {
                for (int j = 0; j < GlDrawable::NUM_SHAPES; j++) {
                    if (!batch_elem->buffers[i][j])
                        continue;
                    if (i == LAYOUT_VTX_TEXTURE0
                        || i == LAYOUT_VTX_NORMAL_TEXTURE0) {
                        texture_bufs.emplace_back(i, batch_elem->buffers[i][j].get());
                    } else {
                        no_texture_bufs.emplace_back(i, batch_elem->buffers[i][j].get());
                    }
                }
            }
            text_bufs.emplace_back(&batch_elem->text_buffer);
        }
        if (params.contains_translucent) {
            _device->enableBlend();
            _device->disableDepthWrite();
        }
        _device->attachTexture(0, _color_tex);
        _device->attachTexture(1, _alpha_tex);
        for (auto& buf : texture_bufs) {
            _device->drawDeviceBuffer(buf.first, *buf.second);
        }
        _device->detachTexture(0);
        _device->detachTexture(1);
        for (auto& buf : no_texture_bufs) {
            _device->drawDeviceBuffer(buf.first, *buf.second);
        }
        if (!params.contains_translucent) {
            _device->enableBlend();
            _device->disableDepthWrite();
        }
        _device->attachTexture(1, _font_tex);
        _device->setNumLights(0);
        for (TextBuffer* buf : text_bufs) {
            _device->drawDeviceBuffer(*buf);
        }
        _device->enableDepthWrite();
        if (!_msaa_enable) { _device->disableBlend(); }
    }
}

void MeshRenderer::buffer(GlDrawable* buf)
{
    for (int i = 0; i < NUM_LAYOUTS; i++) {
        for (int j = 0; j < GlDrawable::NUM_SHAPES; j++) {
            if (!buf->buffers[i][j])
                continue;
            _device->bufferToDevice((array_layout) i, *(buf->buffers[i][j].get()));
        }
    }
    _device->bufferToDevice(buf->text_buffer);
}

void MeshRenderer::queueDraw(MeshRenderingParams params,
                             const std::vector<GlDrawable*>& queue)
{
    _queued[params] = queue;
}

void GLDevice::init()
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glLineWidth(Get_LineWidth());
    glPolygonOffset(1,1);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Generate the default texture: a 1x1 texel with color (1, 1, 1, 1)
    glBindTexture(GL_TEXTURE_2D, 0);
    int blk = 0xFFFFFFFF;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &blk);
}

void GLDevice::setViewport(GLsizei w, GLsizei h)
{
    _vp_width = w;
    _vp_height = h;
    glViewport(0, 0, w, h);
}

void GLDevice::getViewport(GLint (&vp)[4])
{
    vp[0] = vp[1] = 0;
    vp[2] = _vp_width;
    vp[3] = _vp_height;
}

void GLDevice::setTransformMatrices(glm::mat4 model_view, glm::mat4 projection)
{
    _model_view_mtx = model_view;
    _proj_mtx = projection;
}

}