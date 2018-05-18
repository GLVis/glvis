#include "aux_gl3.hpp"

using namespace gl3;

void VertexBuffer::BufferData(GLenum renderAs, std::vector<Vertex>& vertex_data) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vertex_data.size(), vertex_data.data(), GL_STATIC_DRAW);
    this->render_type = renderAs;
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(Vertex), 0);
    glNormalPointer(GL_FLOAT, sizeof(Vertex), sizeof(float) * 3);
    //TODO: this ""should"" assume that color is in RGBA format, is there a more robust method?
    glColorPointer(4, GL_FLOAT, sizeof(Vertex), sizeof(float) * 6);
    /*
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    */
}

void VertexBuffer::DrawObject() {
    glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
    glDrawArrays(render_type, 0, size);
}
