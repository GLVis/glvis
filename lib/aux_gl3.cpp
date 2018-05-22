#include "aux_gl3.hpp"

using namespace gl3;

void VertexBuffer::BufferData(GLenum renderAs) {
    if (!handles_created) { init(); }
    glBindBufferARB(GL_ARRAY_BUFFER, vbo_handles[0]);
    glBufferDataARB(GL_ARRAY_BUFFER, sizeof(GlVertex) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
    glVertexPointer(3, GL_FLOAT, sizeof(GlVertex), 0);
    glNormalPointer(GL_FLOAT, sizeof(GlVertex), (void*) offsetof(GlVertex, norm));
    this->render_type = renderAs;
    if (this->is_textured) {
        glBindBufferARB(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferDataARB(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        glTexCoordPointer(1, GL_FLOAT, 0, 0);
    } else {
        glBindBufferARB(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferDataARB(GL_ARRAY_BUFFER, sizeof(float) * color_data.size(), color_data.data(), GL_STATIC_DRAW);
        glColorPointer(4, GL_FLOAT, 0, 0);
    }
}

void VertexBuffer::DrawObject(bool drawNow) {
    if (!handles_created) { init(); }
    glBindBufferARB(GL_ARRAY_BUFFER, vbo_handles[0]);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    if (this->is_textured) {
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    } else {
        glEnableClientState(GL_COLOR_ARRAY);
    }
    glDrawArrays(render_type, 0, size);
    if (this->is_textured) {
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    } else {
        glDisableClientState(GL_COLOR_ARRAY);
    }
    if (drawNow) { glFlush(); }
}
