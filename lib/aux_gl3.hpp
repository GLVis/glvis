// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. [Insert Release # here] All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef GLVIS_AUX_GL3
#define GLVIS_AUX_GL3
#include <GL/glew.h>

namespace gl3 {

struct Vertex {
    float pos[3];
    float norm[3];
    float rgba[4];
    
    Vertex(const double pos[]) {
        for (int i = 0; i < 3; i++) {
            this->pos[i] = pos[i];
        }
   }

    Vertex(const double pos[], const double norm[]) {
        for (int i = 0; i < 3; i++) {
            this->pos[i] = pos[i];
            this->norm[i] = norm[i];
        }
    }
};

/* *
 * Class to manage vertex buffers
 */
class VertexBuffer {
private:
    GLenum render_type;
    GLuint vbo_handle;
    size_t size;
public:
    /**
     * Constructs a new Vertex buffer object.
     * Prior to creating any VBOs a Vertex Array Object must be binded.
     */
    VertexBuffer()
        : size(0) {
        glGenBuffers(1, &vbo_handle);
    }

    ~VertexBuffer() {
        glDeleteBuffers(1, &vbo_handle);
    }
    /**
     * Buffers vertex data onto the VBO.
     * The VAO associated with this VBO must be binded before calling this.
     */
    void BufferData(GLenum renderAs, std::vector<Vertex>& vertex_data);
    /**
     * Draws the VBO.
     */
    void DrawObject();

    GLenum getRenderType() { return render_type; }
    GLuint getGlHandle() { return vbo_handle; }
    size_t getNumVertices() { return size; } 
};

}
#endif

