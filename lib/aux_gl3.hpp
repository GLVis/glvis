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
#include <unordered_map>
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

class GlDrawable;

class LineBuilder
{
    GlDrawable * parent_buf;
    GLenum render_as;
    std::vector<float> pts;
    int count;
    bool has_color;
    bool has_stipple;
    float color[4];
public:
    LineBuilder(GlDrawable * buf, bool save_color)
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
public:
    enum array_layout {
        LAYOUT_NONE = 0,
        LAYOUT_VTX,
        LAYOUT_VTX_COLOR,
        LAYOUT_VTX_TEXTURE0,
        LAYOUT_VTX_NORMAL_COLOR,
        LAYOUT_VTX_NORMAL_TEXTURE0
    };
private:
    array_layout _layout;
    std::unique_ptr<GLuint> _handle;
    std::vector<float> _pt_data;
    int _buffered_size;

    friend void LineBuilder::glEnd();
public:

    VertexBuffer(array_layout layout)
        : _layout(layout)
        , _handle(new GLuint)
        , _buffered_size(0) {
        glGenBuffers(1, _handle.get());
    }

    ~VertexBuffer() {
        glDeleteBuffers(1, _handle.get());
    }

    void clear() {
        _pt_data.clear();
        _layout = LAYOUT_NONE;
    }

    array_layout getArrayLayout() { return _layout; }

    void addVertex(float (&vtx)[3], float (&rgba)[4]) {
        if (_pt_data.empty()) {
            _layout = LAYOUT_VTX_COLOR;
        } else if (_layout != LAYOUT_VTX_COLOR) {
            cerr << "Unexpected vertex of layout VTX_COLOR" << endl;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(rgba, rgba+4, std::back_inserter(_pt_data));
    }

    void addVertex(float (&vtx)[3], float colorTexCoord) {
        if (_pt_data.empty()) {
            _layout = LAYOUT_VTX_TEXTURE0;
        } else if (_layout != LAYOUT_VTX_TEXTURE0) {
            cerr << "Unexpected vertex of layout VTX_TEXTURE0" << endl;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        _pt_data.emplace_back(colorTexCoord);
        _pt_data.emplace_back(0);
    }

    void addVertex(float (&vtx)[3], float (&norm)[3], float (&rgba)[4]) {
        if (_pt_data.empty()) {
            _layout = LAYOUT_VTX_NORMAL_COLOR;
        } else if (_layout != LAYOUT_VTX_NORMAL_COLOR) {
            cerr << "Unexpected vertex of layout LAYOUT_VTX_NORMAL_COLOR" << endl;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(norm, norm+3, std::back_inserter(_pt_data));
        std::copy(rgba, rgba+4, std::back_inserter(_pt_data));
    }

    void addVertex(float (&vtx)[3], float (&norm)[3], float colorTexCoord) {
        if (_pt_data.empty()) {
            _layout = LAYOUT_VTX_NORMAL_TEXTURE0;
        } else if (_layout != LAYOUT_VTX_NORMAL_TEXTURE0) {
            cerr << "Unexpected vertex of layout VTX_NORMAL_TEXTURE0" << endl;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(norm, norm+3, std::back_inserter(_pt_data));
        _pt_data.emplace_back(colorTexCoord);
        _pt_data.emplace_back(0);
    }

    /**
     * Buffers the vertex data onto the GPU.
     */
    virtual void bufferData();

    /**
     * Draws the vertex data.
     */
    virtual void drawObject(GLenum renderAs);
};

class GlDrawable {
private:
    std::unordered_map<GLenum, VertexBuffer> buffers[6];

    VertexBuffer& getBuffer(VertexBuffer::array_layout layout, GLenum shape) {
        if (buffers[layout].find(shape) == buffers[layout].end()) {
            buffers[layout].emplace(shape, VertexBuffer(layout));
        }
        return buffers[layout][shape];
    }

    friend void LineBuilder::glEnd();
public:

    void addText(float x, float y, float z, std::string text);
    
    /**
     * Creates a new LineBuilder associated with the current drawable object.
     */
    LineBuilder createBuilder(bool save_color = false) {
        return LineBuilder(this, save_color);
    }

    /**
     * Buffers the drawable object onto the GPU.
     */
    void buffer();

    /**
     * Draws the object.
     */
    void draw();
};

}
#endif

