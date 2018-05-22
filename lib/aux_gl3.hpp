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
#include <vector>
namespace gl3 {

struct GlVertex
{
    float pos[3];
    float norm[3];

    GlVertex() {}

    GlVertex(const double pos[]) {
        for (int i = 0; i < 3; i++) {
            this->pos[i] = pos[i];
        }
   }

    GlVertex(const double pos[], const double norm[]) {
        for (int i = 0; i < 3; i++) {
            this->pos[i] = pos[i];
            this->norm[i] = norm[i];
        }
    }
};

/* *
 * Class to manage vertex buffers
 */
class VertexBuffer
{
private:
    GLenum render_type;

    bool handles_created;
    GLuint vbo_handles[2];
    bool is_textured;
    size_t size;

    std::vector<GlVertex> pt_data;
    std::vector<float> color_data;
    std::vector<float> texcoord_data;
    
    void init() {
        glGenBuffers(2, vbo_handles);
        handles_created = true;
    }
public:
    /**
     * Constructs a new Vertex buffer object.
     */
    VertexBuffer()
        : size(0), handles_created(false) {
    }

    ~VertexBuffer() {
        glDeleteBuffers(2, vbo_handles);
    }

    void clear() {
        pt_data.clear();
        color_data.clear();
        texcoord_data.clear();
        size = 0;
    }

    void addVertex(GlVertex gv, float texCoord) {
        is_textured = true;
        pt_data.push_back(gv);
        texcoord_data.push_back(texCoord);
        size++;
    }
    
    void addVertex(GlVertex gv, float (&rgba)[4]) {
        is_textured = false;
        pt_data.push_back(gv); 
        color_data.insert(color_data.end(), rgba, rgba+4);
        size++;
    }

    void BufferData(GLenum renderAs);

    /**
     * Draws the VBO.
     */
    void DrawObject(bool drawNow = true);

    GLenum getRenderType() { return render_type; }
    size_t getNumVertices() { return size; } 
};

}
#endif

