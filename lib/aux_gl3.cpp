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

#include "aux_gl3.hpp"
#include "aux_vis.hpp"
#include <iostream>
#include <utility>
using namespace gl3;

void LineBuilder::glVertex3d(double x, double y, double z) {
#ifdef GLVIS_OGL3
    int offset = (has_color && !has_stipple) ? 7 : 3;
    if (count >= 2 && (render_as == GL_LINE_STRIP || render_as == GL_LINE_LOOP)) {
        pts.reserve(offset * 2 + pts.size());
        //append last-inserted point
        std::copy_n(pts.end() - offset, offset, std::back_inserter(pts));
    } else {
        pts.reserve(offset + pts.size());
    }
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(z);
    if (has_color && !has_stipple) {
        pts.insert(pts.end(), color, color + 4);
    }
    count++;
#else
    ::glVertex3d(x, y, z);
#endif
}

void LineBuilder::glColor3f(float r, float g, float b) {
#ifdef GLVIS_OGL3
    this->color[0] = r;
    this->color[1] = g;
    this->color[2] = b;
    this->color[3] = 1.0;
#else
    ::glColor3f(r,g,b);
#endif
}

void LineBuilder::glColor4fv(float * color) {
#ifdef GLVIS_OGL3
    this->color[0] = color[0];
    this->color[1] = color[1];
    this->color[2] = color[2];
    this->color[3] = color[3];
#else
    ::glColor4fv(color);
#endif
}

void LineBuilder::glEnd() {
#ifdef GLVIS_OGL3
    if (pts.size() < 2
        || (pts.size() == 2 && render_as == GL_LINE_LOOP)) {
        pts.clear();
        return;
    }
    if (render_as == GL_LINE_LOOP) {
        int offset = (has_color && !has_stipple) ? 7 : 3;
        //connect first and last points
        pts.reserve(offset * 2 + pts.size());
        std::copy_n(pts.end() - offset, offset, std::back_inserter(pts));
        std::copy_n(pts.begin(), offset, std::back_inserter(pts));
    }
    if (has_stipple) {
        VertexBuffer& toInsert = parent_buf->getBuffer(VertexBuffer::LAYOUT_VTX, GL_LINES);
        toInsert._pt_data.insert(toInsert._pt_data.end(),
                                 std::make_move_iterator(pts.begin()),
                                 std::make_move_iterator(pts.end()));
    } else if (has_color) {
        VertexBuffer& toInsert = parent_buf->getBuffer(VertexBuffer::LAYOUT_VTX_COLOR, GL_LINES);
        toInsert._pt_data.insert(toInsert._pt_data.end(),
                                 std::make_move_iterator(pts.begin()),
                                 std::make_move_iterator(pts.end()));
    } else {
        VertexBuffer& toInsert = parent_buf->getBuffer(VertexBuffer::LAYOUT_VTX, GL_LINES);
        toInsert._pt_data.insert(toInsert._pt_data.end(),
                                 std::make_move_iterator(pts.begin()),
                                 std::make_move_iterator(pts.end()));
    }
    //if we've std::moved the data, pts is junked
    pts.clear();
    count = 0;
#else
    ::glEnd();
#endif
}

void VertexBuffer::bufferData() {
    if (_pt_data.empty()) {
        return;
    }
    glBindBuffer(GL_ARRAY_BUFFER, *_handle);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _pt_data.size(), _pt_data.data(), GL_STATIC_DRAW);
    _buffered_size = _pt_data.size();
}

void VertexBuffer::drawObject(GLenum renderAs) {
    if (_buffered_size == 0) {
        return;
    }
    glBindBuffer(GL_ARRAY_BUFFER, *_handle);
    glEnableClientState(GL_VERTEX_ARRAY);
    switch (_layout) {
        case LAYOUT_NONE:
            cerr << "ERROR: Unable to buffer vertex data." << endl;
            return;
        case LAYOUT_VTX:
            glVertexPointer(3, GL_FLOAT, 0, 0);
            glDrawArrays(renderAs, 0, _buffered_size / 3);
            break;
        case LAYOUT_VTX_COLOR:
            glVertexPointer(3, GL_FLOAT, sizeof(float) * 7, 0);
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 7);
            glDisableClientState(GL_COLOR_ARRAY);
            break;
        case LAYOUT_VTX_TEXTURE0:
            glVertexPointer(3, GL_FLOAT, sizeof(float) * 8, 0);
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 8);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            break;
        case LAYOUT_VTX_NORMAL_COLOR:
            glVertexPointer(3, GL_FLOAT, sizeof(float) * 10, 0);
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 3));
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 10);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            break;
        case LAYOUT_VTX_NORMAL_TEXTURE0:
            glVertexPointer(3, GL_FLOAT, sizeof(float) * 8, 0);
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 3));
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 8);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            break;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

