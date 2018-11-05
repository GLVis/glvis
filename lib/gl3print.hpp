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

#ifndef __GLVIS_GL3PRINT__
#define __GLVIS_GL3PRINT__

#include "aux_vis.hpp"
#include "glstate.hpp"
#include "aux_gl3.hpp"
#include "gl2ps.h"

#include <vector>
#include <iostream>

namespace gl3
{

struct FeedbackVertex
{
    float pos[4];
    float color[4];
    float clipCoord;
};

void processTriangleTransformFeedback(FeedbackVertex * buf, size_t numVerts);
void processLineTransformFeedback(FeedbackVertex * buf, size_t numVerts);

class GL2PSFeedbackHook : public IDrawHook
{
private:
    GLuint _feedback_buf;
public:
    GL2PSFeedbackHook() {
        glGenBuffers(1, &_feedback_buf);
    }

    ~GL2PSFeedbackHook() {
        if (_feedback_buf != 0) {
            glDeleteBuffers(1, &_feedback_buf);
        }
    }

    GL2PSFeedbackHook(GL2PSFeedbackHook&& other)
        : _feedback_buf(other._feedback_buf) {
        other._feedback_buf = 0;
    }

    GL2PSFeedbackHook& operator = (GL2PSFeedbackHook&& other) {
        if (this != &other) {
            _feedback_buf = other._feedback_buf;
            other._feedback_buf = 0;
        }
        return *this;
    }

    void preDraw(const IVertexBuffer * d) {
        // Allocate buffer and setup feedback
        int buf_size = d->count() * sizeof(FeedbackVertex);
        glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, _feedback_buf);
        glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER,
                     buf_size, nullptr, GL_STATIC_READ);
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, _feedback_buf);
        // Draw objects while capturing vertices
        glEnable(GL_RASTERIZER_DISCARD);
        glBeginTransformFeedback(d->get_shape());
    }

    void postDraw(const IVertexBuffer * d) {
        int buf_size = d->count() * sizeof(FeedbackVertex);
        glEndTransformFeedback();
        glDisable(GL_RASTERIZER_DISCARD);
        // Read buffer
        FeedbackVertex * fb_buf = nullptr;
        fb_buf = new FeedbackVertex[d->count()];
        glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER,
                           0, buf_size, fb_buf);
        if (d->get_shape() == GL_TRIANGLES) {
            processTriangleTransformFeedback(fb_buf, d->count());
        } else if (d->get_shape() == GL_LINES) {
            processLineTransformFeedback(fb_buf, d->count());
        } else { //shape == GL_POINTS/other?
            std::cerr << "Warning: Unhandled primitive type during transform feedback parsing.";
        }
        delete [] fb_buf;
    }

    void preDraw(const TextBuffer& t) {
        glEnable(GL_RASTERIZER_DISCARD);
        GLint vp[4];
        GetGlState()->getViewport(vp);
        for (const auto& entry : t) {
            glm::vec3 raster = glm::project(glm::vec3(entry.rx, entry.ry, entry.rz),
                                            GetGlState()->modelView.mtx,
                                            GetGlState()->projection.mtx,
                                            glm::vec4(0, 0, vp[2], vp[3]));
            GL2PSvertex v = { raster.x, raster.y, raster.z,
                              0, 0, 0, 1 };
            gl2psForceRasterPos(&v);
            gl2psText(entry.text.c_str(), "Times", 8);
        }
    }

    void postDraw(const TextBuffer& t) {
        glDisable(GL_RASTERIZER_DISCARD);
    }

};

}

#endif /* __GLVIS_GLPRINT__ */
