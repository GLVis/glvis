#include "aux_gl3.hpp"
#include <iostream>
using namespace gl3;

void VertexBuffer::BufferData() {
    if (!handles_created) { init(); }
    if (size == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GlVertex) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
    glVertexPointer(3, GL_FLOAT, sizeof(GlVertex), 0);
    glNormalPointer(GL_FLOAT, sizeof(GlVertex), (void*) offsetof(GlVertex, norm));
    if (this->is_textured) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        glTexCoordPointer(1, GL_FLOAT, 0, 0);
    } else {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * color_data.size(), color_data.data(), GL_STATIC_DRAW);
        glColorPointer(4, GL_FLOAT, 0, 0);
    }
}

void VertexBuffer::DrawObject(GLenum renderAs, bool drawNow) {
    if (!handles_created) { init(); }
    if (size == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    if (this->is_textured) {
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    } else {
        glEnableClientState(GL_COLOR_ARRAY);
    }
    glDrawArrays(renderAs, 0, size);
    if (this->is_textured) {
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    } else {
        glDisableClientState(GL_COLOR_ARRAY);
    }
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
//    if (drawNow) { glFlush(); }
}

void LineLoopBuffer::BufferData() {
    if (!handles_created) { init(); }
    if (size == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
    glVertexPointer(3, GL_FLOAT, 0, 0);
}

void LineLoopBuffer::DrawObject(GLenum renderAs, bool drawNow) {
    if (!handles_created) { init(); }
    if (size == 0) return;
    if (renderAs != GL_LINE_LOOP || renderAs != GL_LINES) {
        std::cout << "WARNING: LineLoopBuffer::DrawObject:"
                  << "renderAs parameter not line type, ignored" << std::endl;
    }
    size_t n = 0;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glEnableClientState(GL_VERTEX_ARRAY);
    for (size_t count : loop_strides) {
        glDrawArrays(renderAs, n, count);
        n += count;
    }
    std::cout << "num_vertices: " << n << std::endl;
    glDisableClientState(GL_VERTEX_ARRAY);
//    if (drawNow) { glFlush(); }
}
