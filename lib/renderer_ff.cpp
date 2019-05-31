#include "renderer_ff.hpp"
#include "aux_vis.hpp"
#include "palettes.hpp"
#include "gl2ps.h"

namespace gl3
{

void FFDrawVertex(Vertex v)
{
    glVertex3fv(v.coord.data());
}

void FFDrawVertexColor(VertexColor vc)
{
    glColor4ubv(vc.color.data());
    glVertex3fv(vc.coord.data());
}

void FFDrawVertexTex(VertexTex vt)
{
    glMultiTexCoord2fv(GL_TEXTURE0, vt.texCoord.data());
    glMultiTexCoord2fv(GL_TEXTURE1, vt.texCoord.data());
    glVertex3fv(vt.coord.data());
}

void FFDrawVertexNorm(VertexNorm vn)
{
    glNormal3fv(vn.norm.data());
    glVertex3fv(vn.coord.data());
}

void FFDrawVertexNormColor(VertexNormColor vnc)
{
    glColor4ubv(vnc.color.data());
    glNormal3fv(vnc.norm.data());
    glVertex3fv(vnc.coord.data());
}

void FFDrawVertexNormTex(VertexNormTex vnc)
{
    glMultiTexCoord2fv(GL_TEXTURE0, vnc.texCoord.data());
    glMultiTexCoord2fv(GL_TEXTURE1, vnc.texCoord.data());
    glNormal3fv(vnc.norm.data());
    glVertex3fv(vnc.coord.data());
}

void FFGLDevice::init()
{
    GLDevice::init();
    // Fixed-function pipeline parameters
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    // Texturing is set up such that the output color is computed as follows:
    // - color_out.rgb = color_in.rgb * tex0.rgb
    // - color_out.a = tex1.a
    // Texture unit 0 should contain the color palette, while texture unit 1
    // contains either the transparency alpha channel, or the font texture
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);
    glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_REPLACE);
    glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
    glTexEnvi(GL_TEXTURE_ENV, GL_SRC0_RGB, GL_PREVIOUS);
    glTexEnvi(GL_TEXTURE_ENV, GL_SRC0_ALPHA, GL_TEXTURE);
}

void FFGLDevice::setTransformMatrices(glm::mat4 model_view, glm::mat4 projection)
{
    GLDevice::setTransformMatrices(model_view, projection);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(glm::value_ptr(model_view));
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(glm::value_ptr(projection));
}

void FFGLDevice::setNumLights(int i) {
    if (i == 0) {
        glDisable(GL_LIGHTING);
        return;
    }
    glEnable(GL_LIGHTING);
    for (int light_id = 0; light_id < i; light_id++) {
        glEnable(GL_LIGHT0 + light_id);
    }
    for (int light_id = i; light_id < LIGHTS_MAX; light_id++) {
        glDisable(GL_LIGHT0 + light_id);
    }
}

void FFGLDevice::setMaterial(Material mat)
{
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat.diffuse.data());
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat.ambient.data());
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat.specular.data());
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat.shininess);

    GLfloat memis[] = { 0.0, 0.0, 0.0, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, memis);
}

void FFGLDevice::setPointLight(int i, Light lt)
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glLightfv(GL_LIGHT0 + i, GL_POSITION, lt.position.data());

    GLfloat lambi[] = { 0.0, 0.0, 0.0, 1.0 };
    glLightfv(GL_LIGHT0 + i, GL_AMBIENT, lambi);

    glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, lt.diffuse.data());
    glLightfv(GL_LIGHT0 + i, GL_SPECULAR, lt.specular.data());
    glPopMatrix();
}

void FFGLDevice::setAmbientLight(const std::array<float, 4>& amb)
{
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, &amb[0]);
}

void FFGLDevice::setClipPlaneUse(bool enable)
{
    if (enable) { glEnable(GL_CLIP_PLANE0); }
    else { glDisable(GL_CLIP_PLANE0); }
}

void FFGLDevice::setClipPlaneEqn(const std::array<double, 4>& eqn)
{
    glClipPlane(GL_CLIP_PLANE0, eqn.data());
}

void FFGLDevice::bufferToDevice(array_layout layout, IVertexBuffer& buf)
{
    if (buf.count() == 0) { return; }
    if (buf.get_handle() == 0) {
        int new_hnd = glGenLists(1);
        buf.set_handle(new_hnd);
    }
    glNewList(buf.get_handle(), GL_COMPILE);
    glBegin(buf.get_shape());

    switch(layout) {
        case LAYOUT_VTX:
        {
            for (const Vertex& v : static_cast<VertexBuffer<Vertex>&>(buf)) {
                FFDrawVertex(v);
            }
        }
        break;
        case LAYOUT_VTX_COLOR:
        {
            for (const VertexColor& v : static_cast<VertexBuffer<VertexColor>&>(buf)) {
                FFDrawVertexColor(v);
            }
        }
        break;
        case LAYOUT_VTX_TEXTURE0:
        {
            for (const VertexTex& v : static_cast<VertexBuffer<VertexTex>&>(buf)) {
                FFDrawVertexTex(v);
            }
        }
        break;
        case LAYOUT_VTX_NORMAL:
        {
            for (const VertexNorm& v : static_cast<VertexBuffer<VertexNorm>&>(buf)) {
                FFDrawVertexNorm(v);
            }
        }
        break;
        case LAYOUT_VTX_NORMAL_COLOR:
        {
            for (const VertexNormColor& v : static_cast<VertexBuffer<VertexNormColor>&>(buf)) {
                FFDrawVertexNormColor(v);
            }
        }
        break;
        case LAYOUT_VTX_NORMAL_TEXTURE0:
        {
            for (const VertexNormTex& v : static_cast<VertexBuffer<VertexNormTex>&>(buf)) {
                FFDrawVertexNormTex(v);
            }
        }
        break;
    }
    glEnd();
    glEndList();
}

void FFGLDevice::bufferToDevice(TextBuffer& buf)
{
    // we can't really do anything here
    // can only compute offset matrix at draw
}

void FFGLDevice::drawDeviceBuffer(array_layout layout, const IVertexBuffer& buf)
{
    if (buf.get_handle() == 0)
        return;
    if (buf.count() == 0)
        return;
    if (layout == Vertex::layout || layout == VertexNorm::layout) {
        glColor4fv(_static_color.data());
    } else {
        glColor4f(1.f, 1.f, 1.f, 1.f);
    }
    if (!(layout == VertexNorm::layout
        || layout == VertexNormColor::layout
        || layout == VertexNormTex::layout)) {
        glNormal3f(0.f, 0.f, 1.f);
    }
    glCallList(buf.get_handle());
    // reset texturing parameters
    //glMultiTexCoord2f(GL_TEXTURE0, 0.f, 0.f);
    //glMultiTexCoord2f(GL_TEXTURE1, 0.f, 0.f);
}

void FFGLDevice::drawDeviceBuffer(const TextBuffer& buf)
{
    glColor4fv(_static_color.data());
    glNormal3f(0.f, 0.f, 1.f);
    glMultiTexCoord2f(GL_TEXTURE0, 0.f, 0.f);
    float tex_w = GetFont()->getAtlasWidth();
    float tex_h = GetFont()->getAtlasHeight();
    // Model-view transform:
    // - scale bounding boxes to relative clip-space/NDC coords
    // - add z-offset of -0.005 to reduce text hiding by mesh
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glScalef(2.f / _vp_width, 2.f / _vp_height, 0.f);
    glTranslatef(0.f, 0.f, -0.005);
    glMatrixMode(GL_PROJECTION);
    for (const TextBuffer::Entry& e : buf) {
        glm::vec4 pos(e.rx, e.ry, e.rz, 1.0);
        // transform text starting position into NDC
        pos = _model_view_mtx * pos;
        pos = _proj_mtx * pos;
        pos = pos / pos.w;
        // Projection transform:
        // - add starting offset in NDC 
        glPushMatrix();
        glLoadIdentity();
        glTranslatef(pos.x, pos.y, pos.z);
        float x = 0.f, y = 0.f;
        for (char c : e.text) {
            GlVisFont::glyph g = GetFont()->GetTexChar(c);
            float cur_x = x + g.bear_x;
            float cur_y = -y - g.bear_y;
            x += g.adv_x;
            y += g.adv_y;
            if (!g.w || !g.h)
                continue;
            glBegin(GL_TRIANGLE_STRIP);
                glMultiTexCoord2f(GL_TEXTURE1, g.tex_x, 0);
                glVertex2f(cur_x, -cur_y);
                glMultiTexCoord2f(GL_TEXTURE1, g.tex_x + g.w / tex_w, 0);
                glVertex2f(cur_x + g.w, -cur_y);
                glMultiTexCoord2f(GL_TEXTURE1, g.tex_x, g.h / tex_h);
                glVertex2f(cur_x, -cur_y - g.h);
                glMultiTexCoord2f(GL_TEXTURE1, g.tex_x + g.w / tex_w, g.h / tex_h);
                glVertex2f(cur_x + g.w, -cur_y - g.h);
            glEnd();
        }
        glPopMatrix();
    }
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

unique_ptr<GLDevice::XfbVertexCapture>
FFGLDevice::captureXfbBuffer(array_layout layout, const IVertexBuffer& buf)
{
    GLenum fbType;
    int fbStride;
    unique_ptr<GLDevice::XfbVertexCapture> capture;
    if (layout == VertexTex::layout || layout == VertexNormTex::layout)
    {
        //capture texture values too
        // [ X Y Z ] [ R G B A ] [ U V - - ]
        fbType = GL_3D_COLOR_TEXTURE;
        fbStride = 11;
    } else {
        // only capture pos and color
        // [ X Y Z ] [ R G B A ]
        fbType = GL_3D_COLOR;
        fbStride = 7;
    }
    // compute feedback buffer size
    int sizebuf = 0;
    if (buf.get_shape() == GL_LINES) {
        // for each line: LINE_TOKEN [Vtx] [Vtx]
        sizebuf = (buf.count() / 2) + buf.count() * fbStride;
    } else if (buf.get_shape() == GL_TRIANGLES) {
        // for each tri: POLY_TOKEN 3 [Vtx] [Vtx] [Vtx]
        // NOTE: when clip plane is enabled, we might get two triangles
        // or a quad for an input triangle. However, the other clipped
        // triangles get discarded, so this *should* be enough space.
        sizebuf = (buf.count() / 3) * (2 + fbStride * 4);
    } else {
        std::cerr << "Warning: unhandled primitive type in FFPrinter::preDraw()" << std::endl;

        capture.reset(new VertexCapture(vector<float>(), fbStride));
        return capture;
    }
    // allocate feedback buffer
    vector<float> xfb_buf;
    xfb_buf.reserve(sizebuf);
    glFeedbackBuffer(sizebuf, fbType, xfb_buf.data());
    // draw with feedback capture
    glRenderMode(GL_FEEDBACK);
    drawDeviceBuffer(layout, buf);
    int result = glRenderMode(GL_RENDER);
#ifdef GLVIS_DEBUG
    if (result < 0) {
        std::cerr << "Warning: feedback data exceeded available buffer size" << std::endl;
    }
#endif
    capture.reset(new VertexCapture(std::move(xfb_buf), fbStride));
    return capture;
}

bool FFGLDevice::VertexCapture::next()
{
    if (_tok_idx > _xfb_buf.size())
        return false;
    switch ((GLuint)_xfb_buf[_tok_idx])
    {
        case GL_POINT_TOKEN:
            //unused
            _tok_idx += 1 + _stride;
            break;
        case GL_LINE_TOKEN:
            _tok_idx += 1 + _stride * 2;
            break;
        case GL_POLYGON_TOKEN:
        {
            int num_vtx = _xfb_buf[_tok_idx + 1];
            _poly_idx++;
            if (_poly_idx > num_vtx - 3) {
                _tok_idx += 2 + _stride * num_vtx;
                _poly_idx = 0;
            }
        }
        break;
        default:
            //ignore unused token
            _tok_idx++;
        break;
    }
}

vector<GLDevice::FeedbackVertex> FFGLDevice::VertexCapture::getShape()
{
    vector<FeedbackVertex> shape;
    switch((GLuint)_xfb_buf[_tok_idx])
    {
        case GL_LINE_TOKEN:
        {
            glm::vec3 coord0 = glm::make_vec3(&_xfb_buf[_tok_idx + 1]),
                      coord1 = glm::make_vec3(&_xfb_buf[_tok_idx + 1 + _stride]);
            glm::vec4 color0 = glm::make_vec4(&_xfb_buf[_tok_idx + 4]),
                      color1 = glm::make_vec4(&_xfb_buf[_tok_idx + 4 + _stride]);
            if (_stride == 11) {
                // get texture
                GetColorFromVal(_xfb_buf[_tok_idx + 7], glm::value_ptr(color0));
                GetColorFromVal(_xfb_buf[_tok_idx + 7 + _stride], glm::value_ptr(color1));
            }
            shape = {
                { coord0, color0 },
                { coord1, color1 }
            };
        }
        break;
        case GL_POLYGON_TOKEN:
        {
            int vtx_ids[3] = { 0, _poly_idx + 1, _poly_idx + 2 };
            for (int i = 0; i < 3; i++)
            {
                int tri_idx = _tok_idx + 1 + vtx_ids[i] * _stride;
                glm::vec3 coord = glm::make_vec3(&_xfb_buf[tri_idx]);
                glm::vec4 color = glm::make_vec4(&_xfb_buf[tri_idx + 3]);
                if (_stride == 11) {
                    // get texture
                    GetColorFromVal(_xfb_buf[tri_idx + 7], glm::value_ptr(color));
                }
                shape.emplace_back(coord, color);
            }
        }
        break;
    }
    return shape;
}

}
