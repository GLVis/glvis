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

struct GlColor 
{
    float rgba[4];
    float texcoord;
    bool use_texture;

    GlColor() = default;
    
    GlColor(float val)
        : texcoord(val)
        , use_texture(true) { }

    GlColor(float r, float g, float b, float a)
        : rgba{r,g,b,a}
        , use_texture(false) { }
};

class GlDrawable;

class GlBuilder
{
    GlDrawable * parent_buf;
    GLenum render_as;
    std::vector<float> pts;
    int count;

    bool is_line;
    bool use_color;
    bool use_color_tex;

    float norm[3];
    float color[4];
    float texcoord;
public:
    GlBuilder(GlDrawable * buf)
        : parent_buf(buf)
        , count(0) 
        , is_line(false) 
        , use_color(false)
        , use_color_tex(false) { }

    void glBegin(GLenum e) {
#ifdef GLVIS_OGL3
        if (e == GL_LINES || e == GL_LINE_STRIP || e == GL_LINE_LOOP) {
            is_line = true;
        } else {
            is_line = false;
        }
        render_as = e;
        count = 0;
#else
        ::glBegin(e);
#endif
    }
    
    void glEnd();

    void glVertex3d(double x, double y, double z) {
#ifdef GLVIS_OGL3
        if (count >= 2 && (render_as == GL_LINE_STRIP
                        || render_as == GL_LINE_LOOP)) {
            int offset = 4;
            if (use_color) offset = 8;
            else if (use_color_tex) offset = 5;
            std::copy_n(pts.end() - offset, offset, std::back_inserter(pts));
        }
        pts.emplace_back(x);
        pts.emplace_back(y);
        pts.emplace_back(z);
        if (!is_line) {
            std::copy(norm, norm+3, std::back_inserter(pts));
        }
        if (use_color) {
            std::copy(color, color+4, std::back_inserter(pts));
        } else if (use_color_tex) {
            pts.emplace_back(texcoord);
        }
        if (is_line || use_color_tex) {
            pts.emplace_back(0);
        }
        count++;
#else
        ::glVertex3d(x, y, z);
#endif
    }

    void glVertex3dv(const double * d) {
        glVertex3d(d[0], d[1], d[2]);
    }

    void glNormal3d(double nx, double ny, double nz) {
#ifdef GLVIS_OGL3
        norm[0] = (float) nx;
        norm[1] = (float) ny;
        norm[2] = (float) nz;
#else
        ::glNormal3d(nx, ny, nz);
#endif
    }

    void glNormal3dv(const double * d) {
#ifdef GLVIS_OGL3
        norm[0] = (float) d[0];
        norm[1] = (float) d[1];
        norm[2] = (float) d[2];
#else
        ::glNormal3dv(d);
#endif
    }

    void glColor3f(float r, float g, float b) {
#ifdef GLVIS_OGL3
        this->color[0] = r;
        this->color[1] = g;
        this->color[2] = b;
        this->color[3] = 1.0;
#else
        ::glColor3f(r,g,b);
#endif
    }

    void glColor4fv(float (&rgba)[4]) {
#ifdef GLVIS_OGL3
        if (pts.empty()) {
            use_color = true;
            use_color_tex = false;
        }
        std::copy(rgba, rgba + 4, color);
#else
        ::glColor4fv(rgba);
#endif
    }
    
    void glTexCoord1f(float coord) {
#ifdef GLVIS_OGL3
        if (pts.empty()) {
            use_color_tex = true;
            use_color = false;
        }
        texcoord = coord;
#else
        ::glTexCoord1f(coord);
#endif
    }
};

/* *
 * Class to manage vertex buffers
 */
class VertexBuffer
{
public:
    enum array_layout {
        LAYOUT_VTX,
        LAYOUT_VTX_NORMAL,
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

    friend void GlBuilder::glEnd();
public:

    VertexBuffer(array_layout layout)
        : _layout(layout)
        , _handle(new GLuint)
        , _buffered_size(0) {
        glGenBuffers(1, _handle.get());
    }

    ~VertexBuffer() {
        if (_handle)
            glDeleteBuffers(1, _handle.get());
    }

    VertexBuffer(VertexBuffer&&) = default;
    VertexBuffer& operator = (VertexBuffer&&) = default;

    void clear() {
        _pt_data.clear();
    }

    array_layout getArrayLayout() { return _layout; }

    void addVertex(float (&vtx)[3]) {
        if (_layout != LAYOUT_VTX) {
            cerr << "Unexpected vertex of layout VTX" << endl;
            return;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
    }

    void addVertex(float (&vtx)[3], float (&rgba)[4]) {
        if (_layout != LAYOUT_VTX_COLOR) {
            cerr << "Unexpected vertex of layout VTX_COLOR" << endl;
            return;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(rgba, rgba+4, std::back_inserter(_pt_data));
    }

    void addVertex(float (&vtx)[3], float colorTexCoord) {
        if (_layout != LAYOUT_VTX_TEXTURE0) {
            cerr << "Unexpected vertex of layout VTX_TEXTURE0" << endl;
            return;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        _pt_data.emplace_back(colorTexCoord);
        _pt_data.emplace_back(0);
    }

    void addVertex(float (&vtx)[3], float (&norm)[3], float (&rgba)[4]) {
        if (_layout != LAYOUT_VTX_NORMAL_COLOR) {
            cerr << "Unexpected vertex of layout LAYOUT_VTX_NORMAL_COLOR" << endl;
            return;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(norm, norm+3, std::back_inserter(_pt_data));
        std::copy(rgba, rgba+4, std::back_inserter(_pt_data));
    }

    void addVertex(float (&vtx)[3], float (&norm)[3], float colorTexCoord) {
        if (_layout != LAYOUT_VTX_NORMAL_TEXTURE0) {
            cerr << "Unexpected vertex of layout VTX_NORMAL_TEXTURE0" << endl;
            return;
        }
        std::copy(vtx, vtx+3, std::back_inserter(_pt_data));
        std::copy(norm, norm+3, std::back_inserter(_pt_data));
        _pt_data.emplace_back(colorTexCoord);
        _pt_data.emplace_back(0);
    }

    void addVertices(std::vector<float>&& in_move) {
        _pt_data.insert(_pt_data.end(),
                        std::make_move_iterator(in_move.begin()),
                        std::make_move_iterator(in_move.end()));
    }

    /**
     * Buffers the vertex data onto the GPU.
     */
    void bufferData();

    /**
     * Draws the vertex data.
     */
    void drawObject(GLenum renderAs);
};

class TextBuffer {
private:
    struct _entry {
        float rx, ry, rz;
        std::string text;
        _entry() = default;
        _entry(float x, float y, float z, std::string txt)
            : rx(x), ry(y), rz(z), text(txt) { }
    };

    std::unique_ptr<GLuint> _handle;
    std::vector<_entry> _data;
    size_t _size;
public:
    TextBuffer() : _handle(new GLuint(0)) { };
    ~TextBuffer() {
        if (_handle)
            glDeleteBuffers(1, _handle.get());
    }

    void addText(float x, float y, float z, std::string& text) {
        _data.emplace_back(x, y, z, text);
    }

    void bufferData();

    /**
     * Draws the text.
     */
    void drawObject();

    void clear() {
        _data.clear();
        _size = 0;
    }
};

class GlDrawable {
private:
    std::unordered_map<GLenum, VertexBuffer> buffers[6];
    TextBuffer text_buffer;
public:

    /**
     * Adds a string at the given position in object coordinates.
     */
    void addText(float x, float y, float z, std::string&& text) {
        text_buffer.addText(x, y, z, text);
    }

    void addText(float x, float y, float z, std::string& text) {
        text_buffer.addText(x, y, z, text);
    }

    void addLines(std::vector<float>&& in_move) {
        getBuffer(VertexBuffer::LAYOUT_VTX, GL_LINES)
            .addVertices(std::move(in_move));
    }

    /**
     * Adds a triangle to the drawable object, with the specified face normal
     * and vertex coloring.
     */
    void addTriangle(const double vtx[][3], double (&norm)[3], float (&rgba)[3][4]) {
        float fnorm[3] = { (float) norm[0], (float) norm[1], (float) norm[2] };
        for (int i = 0; i < 3; i++) {
            float fvert[3] = { (float) vtx[i][0], (float) vtx[i][1], (float) vtx[i][2] };
            getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL_COLOR,
                  GL_TRIANGLES).addVertex(fvert, fnorm, rgba[i]);
        }
    }

    /**
     * Adds a triangle to the drawable object, with the specified face normal
     * and color texture coordinates.
     */
    void addTriangle(const double vtx[][3], double (&norm)[3], float (&texcoord)[3]) {
        float fnorm[3] = { (float) norm[0], (float) norm[1], (float) norm[2] };
        for (int i = 0; i < 3; i++) {
            float fvert[3] = { (float) vtx[i][0], (float) vtx[i][1], (float) vtx[i][2] };
            getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL_TEXTURE0,
                  GL_TRIANGLES).addVertex(fvert, fnorm, texcoord[i]);
        }
    }

    /**
     * Adds a quadrilateral to the drawable object, with the specified face normal
     * and vertex coloring.
     */
    void addQuad(const double (&vtx)[4][3], double (&norm)[3], float (&rgba)[4][4]) {
        float fnorm[3] = { (float) norm[0], (float) norm[1], (float) norm[2] };
        int indices[] = {0, 1, 2, 0, 2, 3};
        for (int i : indices) {
            float fvert[3] = { (float) vtx[i][0], (float) vtx[i][1], (float) vtx[i][2] };
            getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL_COLOR,
                  GL_TRIANGLES).addVertex(fvert, fnorm, rgba[i]);
        }
    }

    /**
     * Adds a quadrilateral to the drawable object, with the specified face normal
     * and color.
     */
    void addQuadFace(const double (&vtx)[4][3], double (&norm)[3], float (&rgba)[4]) {
        float fnorm[3] = { (float) norm[0], (float) norm[1], (float) norm[2] };
        int indices[] = {0, 1, 2, 0, 2, 3};
        for (int i : indices) {
            float fvert[3] = { (float) vtx[i][0], (float) vtx[i][1], (float) vtx[i][2] };
            getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL_COLOR,
                  GL_TRIANGLES).addVertex(fvert, fnorm, rgba);
        }
    }


    /**
     * Adds a quadrilateral to the drawable object, with the specified face normal
     * and color texture coordinates.
     */
    void addQuad(const double (&vtx)[4][3], double (&norm)[3], float (&texcoord)[4]) {
        float fnorm[3] = { (float) norm[0], (float) norm[1], (float) norm[2] };
        int indices[] = {0, 1, 2, 0, 2, 3};
        for (int i : indices) {
            float fvert[3] = { (float) vtx[i][0], (float) vtx[i][1], (float) vtx[i][2] };
            getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL_TEXTURE0,
                  GL_TRIANGLES).addVertex(fvert, fnorm, texcoord[i]);
        }
    }

    void addArrow(double px, double py, double pz,
                  double vx, double vy, double vz,
                  double length, double cone_scale);

    void addArrow(double px, double py, double pz,
                  double vx, double vy, double vz,
                  double length, double cone_scale,
                  GlColor color);

    void addShape(GLenum primitive,
                  VertexBuffer::array_layout layout,
                  std::vector<float>&& points) {
            getBuffer(layout, primitive)
                .addVertices(std::move(points));
    }

    VertexBuffer& getBuffer(VertexBuffer::array_layout layout, GLenum shape) {
        auto loc = buffers[layout].find(shape);
        if (loc != buffers[layout].end()) {
            return loc->second;
        }
        buffers[layout].emplace(std::make_pair(shape, VertexBuffer(layout)));
        return buffers[layout].at(shape);
    }

    GlBuilder createBuilder() {
        return GlBuilder(this);
    }

    /**
     * Clears the drawable object.
     */
    void clear() {
        for (int i = 0; i < 6; i++) {
            for (auto& pair : buffers[i]) {
                pair.second.clear();
            }
        }
        text_buffer.clear();
    }
    
    /**
     * Buffers the drawable object onto the GPU.
     */
    void buffer() {
        for (int i = 0; i < 6; i++) {
            for (auto& pair : buffers[i]) {
                pair.second.bufferData();
            }
        }
        text_buffer.bufferData();
    }

    /**
     * Draws the object.
     */
    void draw() {
        for (int i = 0; i < 6; i++) {
            for (auto& pair : buffers[i]) {
                pair.second.drawObject(pair.first);
            }
        }
        text_buffer.drawObject();
    }
};

}
#endif
