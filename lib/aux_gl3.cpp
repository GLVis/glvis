#include "aux_gl3.hpp"
#include <iostream>
#include <utility>
using namespace gl3;

void LineBuilder::glVertex3d(double x, double y, double z) {
#ifdef GLVIS_OGL3
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(z);
#else
    glVertex3d(x, y, z);
#endif
}

void LineBuilder::glEnd() {
#ifdef GLVIS_OGL3
    if (render_as == GL_LINES) {
        std::move(pts.begin(), pts.end(), std::back_inserter(parent_buf->pt_data));
    } else {
        bool first = true;
        for (auto v : pts) {
            parent_buf->pt_data.push_back(v);
            if (!first) {
                parent_buf->pt_data.push_back(v);
            }
            first = false;
        }
        if (render_as == GL_LINE_STRIP) {
            parent_buf->pt_data.pop_back();
        } else {
            parent_buf->pt_data.push_back(*pts.begin());
        }
    }
    //if we've std::moved the data, pts is junked
    pts.clear();
#else
    glEnd();
#endif
}

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

void LineBuffer::BufferData() {
    if (!handles_created) { init(); }
    if (size == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
    glVertexPointer(3, GL_FLOAT, 0, 0);
}

void LineBuffer::DrawObject(GLenum renderAs, bool drawNow) {
    renderHint = renderAs;
    DrawObject();
}

void LineBuffer::DrawObject() {
    if (!handles_created) { init(); }
    if (size == 0) return;
    if (renderHint != GL_LINES) {
        std::cout << "WARNING: LineLoopBuffer::DrawObject:"
                  << "renderAs parameter not line type" << std::endl;
        return;
    }
    size_t n = 0;
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDrawArrays(GL_LINES, 0, size);
    std::cout << "num_vertices: " << n << std::endl;
    glDisableClientState(GL_VERTEX_ARRAY);
//    if (drawNow) { glFlush(); }
}
