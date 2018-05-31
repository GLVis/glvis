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
#include <map>
#include <iostream>

using namespace std;

namespace gl3 {

struct GlVertex
{
    float pos[3];
    float norm[3];

    GlVertex() = default;

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

struct PolyBuilder
{
    float norm[3];
    float color[3];
    float uv[2];
};

class LineBuffer;

class LineBuilder
{
    LineBuffer * parent_buf;
    GLenum render_as;
    std::vector<float> pts;

public:
    LineBuilder(LineBuffer * buf)
        : parent_buf(buf) { }
    
    void glBegin(GLenum e) {
#ifdef GLVIS_OGL3
        render_as = e;
#else
        glBegin(e);
#endif
    }

    void glVertex3d(double x, double y, double z);

    void glEnd();
};

/* *
 * Class to manage vertex buffers
 */
class VertexBuffer
{
protected:

    bool handles_created;
    GLuint vbo_handles[2];
    bool is_textured;
    size_t size;

    std::vector<GlVertex> pt_data;
    std::vector<float> color_data;
    std::vector<float> texcoord_data;
    
    void init() {
        glewInit(); //just in case
        glGenBuffers(2, vbo_handles);
        if (vbo_handles[0] != 0 && vbo_handles[1] != 0) {
            cout << "Handles created" << endl;
            handles_created = true;
        } else {
            const GLubyte* string = gluErrorString(glGetError());
            cout << string << endl;
        }
    }
public:
    /**
     * Constructs a new Vertex buffer object.
     */
    VertexBuffer()
        : handles_created(false), size(0) {
    }

    ~VertexBuffer() {
        cout << "Handles destroyed" << endl;
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

    virtual void BufferData();

    /**
     * Draws the VBO.
     */
    virtual void DrawObject(GLenum renderAs, bool drawNow = true);

    bool isEmpty() { return size == 0; }
};

class LineBuffer : public VertexBuffer {
private:
    std::vector<float> pt_data;
    size_t curr_count;
    GLenum renderHint;
    
    friend void LineBuilder::glEnd();

public:
    LineBuffer()
        : curr_count(0),
          renderHint(GL_LINE_LOOP) {
    }

    ~LineBuffer() {
    }

    void clear() {
        pt_data.clear();
        size = 0;
    }

    void addVertex(float x, float y, float z) {
        pt_data.push_back(x);
        pt_data.push_back(y);
        pt_data.push_back(z);
        curr_count++;
        size++;
    }

    void renderAs(GLenum type) {
        renderHint = type;
    }

    LineBuilder createBuilder() {
        return LineBuilder(this);
    }

    virtual void BufferData();

    /**
     * Draws the VBO.
     */
    virtual void DrawObject(GLenum renderAs, bool drawNow = true);
    
    void DrawObject();
};

}
#endif
