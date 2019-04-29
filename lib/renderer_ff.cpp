#include "renderer.hpp"
#include "aux_vis.hpp"

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
    glTexCoord2fv(vt.texCoord.data());
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
    glTexCoord2fv(vnc.texCoord.data());
    glNormal3fv(vnc.norm.data());
    glVertex3fv(vnc.coord.data());
}

void FFGLDevice::init()
{
    GLDevice::init();
    // Fixed-function pipeline parameters
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    // Texturing is set up such that the output color is computed as follows:
    // - color_out.rgb = color_in.rgb * tex0.rgb
    // - color_out.a = tex1.r
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
    glTexEnvi(GL_TEXTURE_ENV, GL_SRC0_RGB, GL_TEXTURE0);
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
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat.diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat.ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat.specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat.shininess);
}

void FFGLDevice::setPointLight(int i, Light lt)
{
    glLightfv(GL_LIGHT0 + i, GL_POSITION, lt.position);
    glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, lt.diffuse);
    glLightfv(GL_LIGHT0 + i, GL_SPECULAR, lt.specular);
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
}

void FFGLDevice::drawDeviceBuffer(const TextBuffer& buf)
{
    glColor4fv(_static_color.data());
    glNormal3f(0.f, 0.f, 1.f);
    float tex_w = GetFont()->getAtlasWidth();
    float tex_h = GetFont()->getAtlasHeight();
    glm::mat4 projTextMtx = glm::ortho<float>(0.f, _vp_width,
                                              0.f, _vp_height,
                                              -5.f, 5.f);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadMatrixf(glm::value_ptr(projTextMtx));
    glMatrixMode(GL_MODELVIEW);
    for (const TextBuffer::Entry& e : buf) {
        glm::vec4 pos(e.rx, e.ry, e.rz, 1.0);
        pos = _model_view_mtx * pos;
        pos = _proj_mtx * pos;
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
                glVertex2f(cur_x, -cur_y);
                glTexCoord2f(g.tex_x, 0);
                glVertex2f(cur_x + g.w, -cur_y);
                glTexCoord2f(g.tex_x + g.w / tex_w, 0);
                glVertex2f(cur_x + g.w, -cur_y - g.h);
                glTexCoord2f(g.tex_x + g.w / tex_w, g.h / tex_h);
                glVertex2f(cur_x, -cur_y - g.h);
                glTexCoord2f(g.tex_x, g.h / tex_h);
            glEnd();
        }
        glPopMatrix();
    }
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}

}