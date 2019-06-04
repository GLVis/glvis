#include "renderer.hpp"

namespace gl3
{

GLenum MeshRenderer::getDeviceAlphaChannel()
{
    if (!_device) return GL_NONE;
#ifdef __EMSCRIPTEN__
    return GL_ALPHA;
#else
    if (_device->getType() == GLDevice::FF_DEVICE) {
        return GL_ALPHA;
    } else {
        return GL_RED;
    }
#endif
}

void MeshRenderer::setAntialiasing(bool aa_status)
{
    if (_msaa_enable != aa_status) {
        _msaa_enable = aa_status;
        if (_msaa_enable) {
            _device->enableMultisample();
            _device->enableBlend();
            _device->setLineWidth(_line_w_aa);
        } else {
            _device->disableMultisample();
            _device->disableBlend();
            _device->setLineWidth(_line_w);
        }
    }
}

void MeshRenderer::setLineWidth(float w)
{
    _line_w = w;
    if (_device && !_msaa_enable) {
        _device->setLineWidth(_line_w);
    }
}

void MeshRenderer::setLineWidthMS(float w)
{
    _line_w_aa = w;
    if (_device && _msaa_enable) {
        _device->setLineWidth(_line_w_aa);
    }
}

void MeshRenderer::init()
{
}

void MeshRenderer::render(const RenderQueue& queue)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for (auto& q_elem : queue) {
        const RenderParams& params = q_elem.first;
        _device->setTransformMatrices(params.model_view.mtx, params.projection.mtx);
        _device->setMaterial(params.mesh_material);
        _device->setNumLights(params.num_pt_lights);
        for (int i = 0; i < params.num_pt_lights; i++) {
            _device->setPointLight(i, params.lights[i]);
        }
        _device->setAmbientLight(params.light_amb_scene);
        _device->setStaticColor(params.static_color);
        _device->setClipPlaneUse(params.use_clip_plane);
        _device->setClipPlaneEqn(params.clip_plane_eqn);
        //aggregate buffers with common parameters
        std::vector<pair<array_layout, IVertexBuffer*>> texture_bufs, no_texture_bufs;
        std::vector<TextBuffer*> text_bufs;
        GlDrawable* curr_drawable = q_elem.second;
        for (int i = 0; i < NUM_LAYOUTS; i++) {
            for (int j = 0; j < GlDrawable::NUM_SHAPES; j++) {
                if (!curr_drawable->buffers[i][j])
                    continue;
                if (i == LAYOUT_VTX_TEXTURE0
                    || i == LAYOUT_VTX_NORMAL_TEXTURE0) {
                    texture_bufs.emplace_back((array_layout) i, curr_drawable->buffers[i][j].get());
                } else {
                    no_texture_bufs.emplace_back((array_layout) i, curr_drawable->buffers[i][j].get());
                }
            }
        }
        text_bufs.emplace_back(&curr_drawable->text_buffer);
        if (params.contains_translucent) {
            _device->enableBlend();
        } else {
            _device->enableDepthWrite();
        }
        _device->attachTexture(GLDevice::SAMPLER_COLOR, _color_tex);
        _device->attachTexture(GLDevice::SAMPLER_ALPHA, _alpha_tex);
        for (auto& buf : texture_bufs) {
            _device->drawDeviceBuffer(buf.first, *buf.second);
        }
        _device->detachTexture(GLDevice::SAMPLER_COLOR);
        _device->detachTexture(GLDevice::SAMPLER_ALPHA);
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

CaptureBuffer MeshRenderer::capture(const RenderQueue& queue)
{
    CaptureBuffer cbuf;
    _device->initXfbMode();
    for (auto& q_elem : queue) {
        const RenderParams& params = q_elem.first;
        _device->setTransformMatrices(params.model_view.mtx, params.projection.mtx);
        _device->setMaterial(params.mesh_material);
        _device->setNumLights(params.num_pt_lights);
        for (int i = 0; i < params.num_pt_lights; i++) {
            _device->setPointLight(i, params.lights[i]);
        }
        _device->setAmbientLight(params.light_amb_scene);
        _device->setStaticColor(params.static_color);
        _device->setClipPlaneUse(params.use_clip_plane);
        _device->setClipPlaneEqn(params.clip_plane_eqn);
        //aggregate buffers with common parameters
        std::vector<pair<array_layout, IVertexBuffer*>> texture_bufs, no_texture_bufs;
        std::vector<TextBuffer*> text_bufs;
        GlDrawable* curr_drawable = q_elem.second;
        for (int i = 0; i < NUM_LAYOUTS; i++) {
            for (int j = 0; j < GlDrawable::NUM_SHAPES; j++) {
                if (!curr_drawable->buffers[i][j])
                    continue;
                if (i == LAYOUT_VTX_TEXTURE0
                    || i == LAYOUT_VTX_NORMAL_TEXTURE0) {
                    texture_bufs.emplace_back((array_layout) i, curr_drawable->buffers[i][j].get());
                } else {
                    no_texture_bufs.emplace_back((array_layout) i, curr_drawable->buffers[i][j].get());
                }
            }
        }
        text_bufs.emplace_back(&curr_drawable->text_buffer);

        _device->attachTexture(GLDevice::SAMPLER_COLOR, _color_tex);
        _device->attachTexture(GLDevice::SAMPLER_ALPHA, _alpha_tex);
        for (auto& buf : texture_bufs) {
			_device->captureXfbBuffer(cbuf, buf.first, *buf.second);
        }
        _device->detachTexture(GLDevice::SAMPLER_COLOR);
        _device->detachTexture(GLDevice::SAMPLER_ALPHA);
        for (auto& buf : no_texture_bufs) {
            _device->captureXfbBuffer(cbuf, buf.first, *buf.second);
        }
        if (!params.contains_translucent) {
            _device->enableBlend();
            _device->disableDepthWrite();
        }
        _device->attachTexture(1, _font_tex);
        _device->setNumLights(0);
        for (TextBuffer* buf : text_bufs) {
            _device->captureXfbBuffer(cbuf, *buf);
        }
    }
    _device->exitXfbMode();
    return cbuf;
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

void GLDevice::init()
{
    // enable depth testing
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    // enable polygon offset to expose mesh lines
    glPolygonOffset(1,1);
    glEnable(GL_POLYGON_OFFSET_FILL);
    // use "over" blending equation
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //generate a white default texture
    //modulation with default texture will just pass through input color
    glBindTexture(GL_TEXTURE_2D, 0);
    int black_color = 0xFFFFFFFF;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &black_color);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
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

void GLDevice::captureXfbBuffer(CaptureBuffer& capture, const TextBuffer& t_buf)
{
	for (const auto& entry : t_buf)
	{
		glm::vec3 raster = glm::project(
				glm::vec3(entry.rx, entry.ry, entry.rz),
				_model_view_mtx,
				_proj_mtx,
				glm::vec4(0, 0, _vp_width, _vp_height));
		capture.text.emplace_back(raster, glm::make_vec4(_static_color.data()), entry.text); 
	}
}

}
