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
#include <array>
#include <iostream>
#include <memory>

#include "platform_gl.hpp"

using namespace std;

namespace gl3 {

enum array_layout {
    LAYOUT_VTX = 0,
    LAYOUT_VTX_NORMAL,
    LAYOUT_VTX_COLOR,
    LAYOUT_VTX_TEXTURE0,
    LAYOUT_VTX_NORMAL_COLOR,
    LAYOUT_VTX_NORMAL_TEXTURE0,
    NUM_LAYOUTS
};

struct alignas(16) Vertex
{
    std::array<float, 3> coord;

    static void setupAttribLayout(); 
    static void clearAttribLayout() { }
    static const int layout = LAYOUT_VTX;
};

struct alignas(16) VertexColor
{
    std::array<float, 3> coord;
    std::array<uint8_t, 4> color;

    static void setupAttribLayout();
    static void clearAttribLayout();
    static const int layout = LAYOUT_VTX_COLOR;
};

struct alignas(16) VertexTex
{
    std::array<float, 3> coord;
    std::array<float, 2> texCoord;

    static void setupAttribLayout();
    static void clearAttribLayout();
    static const int layout = LAYOUT_VTX_TEXTURE0;
};

struct alignas(16) VertexNorm
{
    std::array<float, 3> coord;
    std::array<float, 3> norm;

    static void setupAttribLayout();
    static void clearAttribLayout();
    static const int layout = LAYOUT_VTX_NORMAL;
};

struct alignas(16) VertexNormColor
{
    std::array<float, 3> coord;
    std::array<float, 3> norm;
    std::array<uint8_t, 4> color;

    static void setupAttribLayout();
    static void clearAttribLayout();
    static const int layout = LAYOUT_VTX_NORMAL_COLOR;
};

struct alignas(16) VertexNormTex
{
    std::array<float, 3> coord;
    std::array<float, 3> norm;
    std::array<float, 2> texCoord;

    static void setupAttribLayout();
    static void clearAttribLayout();
    static const int layout = LAYOUT_VTX_NORMAL_TEXTURE0;
};

inline std::array<uint8_t, 4> ColorU8(float r, float g, float b, float a) {
    return {
        (r >= 1.0) ? (uint8_t) 255 : (uint8_t)(r * 256.),
        (g >= 1.0) ? (uint8_t) 255 : (uint8_t)(g * 256.),
        (b >= 1.0) ? (uint8_t) 255 : (uint8_t)(b * 256.),
        (a >= 1.0) ? (uint8_t) 255 : (uint8_t)(a * 256.),
    };
}

inline std::array<uint8_t, 4> ColorU8(float rgba[]) {
    return {
        (rgba[0] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[0] * 256.),
        (rgba[1] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[1] * 256.),
        (rgba[2] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[2] * 256.),
        (rgba[3] >= 1.0) ? (uint8_t) 255 : (uint8_t)(rgba[3] * 256.),
    };
}

class GlDrawable;

/**
 * Crude fixed-function OpenGL emulation helper.
 */
class GlBuilder
{
    GlDrawable * parent_buf;
    GLenum render_as;
    int count;

    bool is_line;

    bool use_norm;
    bool use_color;
    bool use_tex;

    struct _vertex {
        std::array<float, 3> coords;
        std::array<float, 3> norm;
        std::array<uint8_t, 4> color;
        std::array<float, 2> texcoord;
    };

    _vertex saved[3];
    _vertex curr;

    void saveVertex(const _vertex& v);

public:
    GlBuilder(GlDrawable * buf)
        : parent_buf(buf)
        , count(0) 
        , is_line(false)
        , use_norm(false)
        , use_color(false)
        , use_tex(false) { }

    void glBegin(GLenum e) {
        if (e == GL_LINES || e == GL_LINE_STRIP || e == GL_LINE_LOOP) {
            is_line = true;
        } else {
            is_line = false;
        }
        render_as = e;
        count = 0;
    }
    
    void glEnd() {
        //create degenerate primitives if necessary
        if (render_as == GL_LINES && count % 2 != 0) {
            saveVertex(curr);
        }
        if (render_as == GL_TRIANGLES) {
            for (int i = 0; i < count % 3; i++) {
                saveVertex(curr);
            }
        }
        if (render_as == GL_LINE_LOOP && count > 2) {
            //link first and last vertex
            saveVertex(saved[0]);
            saveVertex(saved[1]);
        }
        count = 0;
    }

    void glVertex3d(double x, double y, double z) {
        curr.coords = { (float) x, (float) y, (float) z };
        if (render_as == GL_LINES || render_as == GL_TRIANGLES) {
            //Lines and triangles are stored as-is
            saveVertex(curr);
        } else if (is_line) {
            //LINE_LOOP and LINE_STRIP: each vertex creates a new line with the
            //previous vertex
            if (count == 0) {
                saved[0] = curr;
            } else {
                saveVertex(saved[1]);
                saveVertex(curr);
            }
            saved[1] = curr;
        } else if (render_as == GL_QUADS) {
            if (count % 4 == 3) {
                //split on 0-2 diagonal
                saveVertex(saved[0]);
                saveVertex(saved[1]);
                saveVertex(saved[2]);
                saveVertex(saved[0]);
                saveVertex(saved[2]);
                saveVertex(curr);
                count -= 4;
            } else {
                saved[count] = curr;
            }
        } else {
            //TriangleStrip, TriangleFan, Polygon
            if (count >= 2) {
                saveVertex(saved[0]);
                saveVertex(saved[1]);
                saveVertex(curr);
                if (render_as == GL_TRIANGLE_STRIP) {
                    //pop earliest vertex
                    saved[0] = saved[1];
                }
                saved[1] = curr;
            } else {
                saved[count] = curr;
            }
        }
        count++;
    }

    void glVertex3dv(const double * d) {
        glVertex3d(d[0], d[1], d[2]);
    }

    void glNormal3d(double nx, double ny, double nz) {
        use_norm = true;
        curr.norm = { (float) nx, (float) ny, (float) nz };
    }

    void glNormal3dv(const double * d) { glNormal3d(d[0], d[1], d[2]); }

    void glColor4f(float r, float g, float b, float a) {
        if (count == 0) {
            use_color = true;
            use_tex = false;
        }
        curr.color = {
            (r >= 1.0) ? (uint8_t) 255 : (uint8_t)(r * 256.),
            (g >= 1.0) ? (uint8_t) 255 : (uint8_t)(g * 256.),
            (b >= 1.0) ? (uint8_t) 255 : (uint8_t)(b * 256.),
            (a >= 1.0) ? (uint8_t) 255 : (uint8_t)(a * 256.),
        };
    }

    void glColor3f(float r, float g, float b) { glColor4f(r, g, b, 1.f); }

    void glColor4fv(float * cv) { glColor4f(cv[0], cv[1], cv[2], cv[3]); }
    
    void glTexCoord2f(float coord_u, float coord_v) {
        if (count == 0) {
            use_tex = true;
            use_color = false;
        }
        curr.texcoord = { coord_u, coord_v };
    }
};

class IVertexBuffer
{
public:
    virtual ~IVertexBuffer() { }
    virtual void clear() = 0;
    virtual void buffer() = 0;
    virtual void draw() = 0;

    virtual size_t count() const = 0;
    virtual GLenum get_shape() const = 0;
};

template<typename T>
class VertexBuffer : public IVertexBuffer
{
private:
    GLenum _shape;
    std::vector<T> _data;
    std::unique_ptr<GLuint> _handle;
    size_t _buffered_size;
    size_t _allocated_size;

public:
    VertexBuffer(GLenum shape)
        : _shape(shape)
        , _handle(new GLuint)
        , _buffered_size(0)
        , _allocated_size(0) {
        glGenBuffers(1, _handle.get());
    }

    ~VertexBuffer() {
        if (_handle)
            glDeleteBuffers(1, _handle.get());
    }

    VertexBuffer(VertexBuffer&&) = default;
    VertexBuffer& operator = (VertexBuffer&&) = default;

    /**
     * Returns the number of vertices buffered on the GPU.
     */
    virtual size_t count() const { return _buffered_size; }

    /**
     * Gets the primitive type contained by the vertex buffer.
     */
    virtual GLenum get_shape() const { return _shape; }

    /**
     * Clears the buffer of all data.
     */
    virtual void clear() {
        _data.clear();
        _buffered_size = 0;
    }

    /**
     * Buffers vertex data onto the GPU.
     */
    virtual void buffer() {
        if (_data.empty()) {
            return;
        }
        glBindBuffer(GL_ARRAY_BUFFER, *_handle);
        if (_allocated_size >= _data.size()) {
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(T) * _data.size(), _data.data());
        } else {
            glBufferData(GL_ARRAY_BUFFER, sizeof(T) * _data.size(), _data.data(), GL_DYNAMIC_DRAW);
            _allocated_size = _data.size();
        }
        _buffered_size = _data.size();
    }

    /**
     * Draws the vertex data buffered on the GPU.
     */
    virtual void draw() {
        if (_buffered_size == 0) {
            return;
        }
        glBindBuffer(GL_ARRAY_BUFFER, *_handle);
        T::setupAttribLayout();
        glDrawArrays(_shape, 0, _buffered_size);
        T::clearAttribLayout();
    }

    /**
     * Add a vertex to the buffer.
     */
    void addVertex(const T& vertex) {
        _data.emplace_back(vertex);
    }
};

class TextBuffer
{
public:
    struct Entry {
        float rx, ry, rz;
        std::string text;
        Entry() = default;
        Entry(float x, float y, float z, const std::string& txt)
            : rx(x), ry(y), rz(z), text(txt) { }
    };
    typedef std::vector<Entry>::iterator Iterator;
    typedef std::vector<Entry>::const_iterator ConstIterator;
private:
    std::unique_ptr<GLuint> _handle;
    std::vector<Entry> _data;
    size_t _size;
public:
    TextBuffer() : _handle(new GLuint(0)) { };
    ~TextBuffer() {
        if (_handle)
            glDeleteBuffers(1, _handle.get());
    }

    /**
     * Adds a text element at the specified local space (pre-transform) coordinates.
     */
    void addText(float x, float y, float z, const std::string& text) {
        _data.emplace_back(x, y, z, text);
    }

    /**
     * Gets an iterator referencing the first text entry.
     */
    Iterator begin() { return _data.begin(); }
    ConstIterator begin() const { return _data.begin(); }

    /**
     * Gets an iterator pointing to the last text entry.
     */
    Iterator end() { return _data.end(); }
    ConstIterator end() const { return _data.end(); }

    /**
     * Buffers the text onto the GPU.
     */
    void buffer();

    /**
     * Gets the width and height of the bounding box containing the rendered
     * text.
     */
    static void getObjectSize(const std::string& text, int& w, int& h);

    /**
     * Draws the text buffered onto the GPU.
     */
    void draw();

    /**
     * Clears the text buffer.
     */
    void clear() {
        _data.clear();
        _size = 0;
    }
};

class IDrawHook {
public:
    virtual void preDraw(const IVertexBuffer *) = 0;
    virtual void postDraw(const IVertexBuffer *) = 0;

    virtual void preDraw(const TextBuffer&) = 0;
    virtual void postDraw(const TextBuffer&) = 0;
};

class GlDrawable
{
private:
    static IDrawHook * buf_hook;
    const static size_t NUM_SHAPES = 2;
    std::unique_ptr<IVertexBuffer> buffers[NUM_LAYOUTS][NUM_SHAPES];
    TextBuffer text_buffer;

    friend class GlBuilder;

    template<typename Vert>
    VertexBuffer<Vert> * getBuffer(GLenum shape) {
        int idx = -1;
        if (shape == GL_LINES) {
            idx = 0;
        } else if (shape == GL_TRIANGLES) {
            idx = 1;
        } else {
            return nullptr;
        }
        if (!buffers[Vert::layout][idx]) {
            buffers[Vert::layout][idx].reset(new VertexBuffer<Vert>(shape));
        }
        VertexBuffer<Vert> * buf = static_cast<VertexBuffer<Vert>*>(buffers[Vert::layout][idx].get());
        return buf;
    }
public:
    /**
     * Sets a global draw hook to be called before and after each vertex buffer
     * draw.
     */
    static void setDrawHook(IDrawHook * h) { buf_hook = h; }

    /**
     * Adds a string at the given position in object coordinates.
     */
    void addText(float x, float y, float z, const std::string& text) {
        text_buffer.addText(x, y, z, text);
    }

    template<typename Vert>
    void addLine(const Vert& v1, const Vert& v2) {
        getBuffer<Vert>(GL_LINES)->addVertex(v1);
        getBuffer<Vert>(GL_LINES)->addVertex(v2);
    }

    template<typename Vert>
    void addTriangle(const Vert& v1, const Vert& v2, const Vert& v3) {
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v2);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
    }

    template<typename Vert>
    void addQuad(const Vert& v1, const Vert& v2, const Vert& v3, const Vert& v4) {
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v2);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v1);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v3);
        getBuffer<Vert>(GL_TRIANGLES)->addVertex(v4);
    }

    void addCone(float x, float y, float z,
                 float vx, float vy, float vz,
                 float cone_scale = 0.075);

    GlBuilder createBuilder() {
        return GlBuilder(this);
    }

    /**
     * Clears the drawable object.
     */
    void clear() {
        for (int i = 0; i < NUM_LAYOUTS; i++) {
            for (int j = 0; j < NUM_SHAPES; j++) {
                if (buffers[i][j]) {
                    buffers[i][j]->clear();
                }
            }
        }
        text_buffer.clear();
    }
    
    /**
     * Buffers the drawable object onto the GPU.
     */
    void buffer() {
        for (int i = 0; i < NUM_LAYOUTS; i++) {
            for (int j = 0; j < NUM_SHAPES; j++) {
                if (buffers[i][j]) {
                    buffers[i][j]->buffer();
                }
            }
        }
        text_buffer.buffer();
    }

    /**
     * Draws the object.
     */
    void draw() {
        for (int i = 0; i < NUM_LAYOUTS; i++) {
            for (int j = 0; j < NUM_SHAPES; j++) {
                if (!buffers[i][j])
                    continue;
                if (GlDrawable::buf_hook) {
                    GlDrawable::buf_hook->preDraw(buffers[i][j].get());
                }
                buffers[i][j]->draw();
                if (GlDrawable::buf_hook) {
                    GlDrawable::buf_hook->postDraw(buffers[i][j].get());
                }
            }
        }
        if (GlDrawable::buf_hook) {
            GlDrawable::buf_hook->preDraw(text_buffer);
            text_buffer.draw();
            GlDrawable::buf_hook->postDraw(text_buffer);
        } else {
            text_buffer.draw();
        }
    }
};

}
#endif

