// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
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
#include <vector>
#include <map>
#include <iostream>
#include <memory>

#include "platform_gl.hpp"

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
    int count;
    bool has_color;
    bool has_stipple;
    float color[4];
public:
    LineBuilder(LineBuffer * buf, bool save_color)
        : parent_buf(buf)
        , has_color(save_color)
        , has_stipple(false) { }

    void setUseColor(bool use = true) {
        has_color = use;
    }

    void glBegin(GLenum e) {
#ifdef GLVIS_OGL3
        render_as = e;
        count = 0;
#else
        ::glBegin(e);
#endif
    }
    
    void glEnable(GLenum e) {
#ifdef GLVIS_OGL3
        if (pts.size() != 0) {
            cerr << "WARNING: Disabling stipple in the middle of a glEnable/glDisable block" << endl;
        }
        if (e == GL_LINE_STIPPLE) {
            has_stipple = true;
        }
#else
        ::glEnable(e);
#endif
    }

    void glDisable(GLenum e) {
#ifdef GLVIS_OGL3
        if (pts.size() != 0) {
            cerr << "WARNING: Disabling stipple in the middle of a glEnable/glDisable block" << endl;
        }
        if (e == GL_LINE_STIPPLE) {
            has_stipple = false;
        }
#else
        ::glDisable(e);
#endif
    }
    void glVertex3d(double x, double y, double z);
    void glColor3f(float r, float g, float b);
    void glColor4fv(float * color);

    void glEnd();
};

/* *
 * Class to manage vertex buffers
 */
class VertexBuffer
{
protected:

    struct _handle {
        GLuint vbo_handles[3];
        _handle() {
            glGenBuffers(3, vbo_handles);
        }
        ~_handle() {
            glDeleteBuffers(3, vbo_handles);
        }

        GLuint get(int i) {
            return vbo_handles[i];
        }
    };
    std::shared_ptr<_handle> vbo;
    
    // contains only point data
    std::vector<float> pt_data;
    int pt_cnt;
    // contains interleaved point and color data
    // default layout: V(3f)|N(3f)|C(4f)
    std::vector<float> color_data;
    int color_cnt;
    // contains interleaved point and texture coord data
    // default layout: V(3f)|N(3f)|T(1f)
    std::vector<float> texcoord_data;
    int texcoord_cnt;
    
public:
    /**
     * Constructs a new Vertex buffer object.
     */
    VertexBuffer()
        : vbo(std::make_shared<_handle>()) {
    }

    ~VertexBuffer() { }

    virtual void clear() {
        pt_data.clear();
        color_data.clear();
        texcoord_data.clear();
        pt_cnt = 0;
        color_cnt = 0;
        texcoord_cnt = 0;
    }

    void addVertex(GlVertex gv, float texCoord) {
        std::move(gv.pos, gv.pos + 3, std::back_inserter(texcoord_data));
        std::move(gv.norm, gv.norm + 3, std::back_inserter(texcoord_data));
        texcoord_data.push_back(texCoord);
        texcoord_data.push_back(0);
    }
    
    void addVertex(GlVertex gv, float (&rgba)[4]) {
        std::move(gv.pos, gv.pos + 3, std::back_inserter(color_data));
        std::move(gv.norm, gv.norm + 3, std::back_inserter(color_data));
        color_data.insert(color_data.end(), rgba, rgba+4);
    }

    virtual void BufferData();

    /**
     * Draws the VBO.
     */
    virtual void DrawObject(GLenum renderAs);
};

class LineBuffer : public VertexBuffer {
    //texcoord_data is used to store stipple lines
    friend void LineBuilder::glEnd();
public:

    void addVertex(float x, float y, float z) {
        pt_data.push_back(x);
        pt_data.push_back(y);
        pt_data.push_back(z);
    }

    LineBuilder createBuilder(bool save_color = false) {
        return LineBuilder(this, save_color);
    }

    virtual void BufferData();

    /**
     * Draws the VBO.
     */
    virtual void DrawObject(GLenum renderAs);
    
    virtual void DrawObject();
};

class TextBuffer : public LineBuffer {
    struct _text_entry {
        float x, y, z;
        std::string str;
        _text_entry(float x, float y, float z, std::string str)
            : x(x), y(y), z(z), str(std::move(str)) { }
    };
    std::vector<_text_entry> entries;
public:
    TextBuffer() { }
    ~TextBuffer() { }

    virtual void clear() {
        entries.clear();
        LineBuffer::clear();
    }

    void SetText(float x, float y, float z, std::string text);
    virtual void BufferData();
    virtual void DrawObject(GLenum renderAs);
    virtual void DrawObject();
};

}
#endif

