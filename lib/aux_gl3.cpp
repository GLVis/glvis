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
    if (count > 2 && (render_as == GL_LINE_STRIP || render_as == GL_LINE_LOOP)) {
        pts.reserve(offset * 2 + pts.size());
        //append last-inserted point
        std::copy_n(pts.end() - offset - 1, offset, std::back_inserter(pts));
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
        std::copy_n(pts.end() - offset - 1, offset, std::back_inserter(pts));
        std::copy_n(pts.begin(), offset, std::back_inserter(pts));
    }
    if (has_stipple) {
        std::move(pts.begin(), pts.end(), std::back_inserter(parent_buf->texcoord_data));
    } else if (has_color) {
        std::move(pts.begin(), pts.end(), std::back_inserter(parent_buf->color_data));
    } else {
        std::move(pts.begin(), pts.end(), std::back_inserter(parent_buf->pt_data));
    }
    //if we've std::moved the data, pts is junked
    pts.clear();
    count = 0;
#else
    ::glEnd();
#endif
}

void VertexBuffer::BufferData() {
    if (!handles_created) { init(); }
    if (!pt_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        pt_cnt = pt_data.size() / 3;
    }
    if (!color_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 10, 0);
        glNormalPointer(GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 3));
        glColorPointer(4, GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 6));
        color_cnt = color_data.size() / 10;
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[2]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 7, 0);
        glNormalPointer(GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
        glTexCoordPointer(1, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 6));
        texcoord_cnt = texcoord_data.size() / 7;
    }
}

void VertexBuffer::DrawObject(GLenum renderAs) {
    if (!handles_created) { init(); }
    glEnableClientState(GL_VERTEX_ARRAY);
    if (!pt_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
        glDrawArrays(renderAs, 0, pt_cnt);
    }
    if (!color_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glDrawArrays(renderAs, 0, color_cnt);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[2]);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glDrawArrays(renderAs, 0, texcoord_cnt);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
}

void TextBuffer::BufferData() {
    // Stub since we're just drawing directly
    LineBuffer::DrawObject();
}

void TextBuffer::DrawObject(GLenum renderAs) {
    DrawObject();
}

void TextBuffer::DrawObject() {
    LineBuffer::DrawObject();
    if (entries.size() == 0) { return; }
#ifndef GLVIS_USE_FREETYPE
    cerr << "Can't use text buffer object without Freetype" << endl;
#else
    for (auto& str_obj : entries) {
        glRasterPos3f(str_obj.x, str_obj.y, str_obj.z);
        DrawBitmapText(str_obj.text.c_str());
    }
#endif
}

void LineBuffer::BufferData() {
    if (!handles_created) { init(); }
    if (!pt_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        pt_cnt = pt_data.size() / 3;
    }
    if (!color_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 7, 0);
        glColorPointer(4, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
        color_cnt = color_data.size() / 7;
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[2]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        texcoord_cnt = texcoord_data.size() / 3;
    }
}

void LineBuffer::DrawObject(GLenum renderAs) {
    DrawObject();
}

void LineBuffer::DrawObject() {
    if (!handles_created) { init(); }
    glEnableClientState(GL_VERTEX_ARRAY);
    if (!pt_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[0]);
        glDrawArrays(GL_LINES, 0, pt_cnt);
    }
    if (!color_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[1]);
        glEnableClientState(GL_COLOR_ARRAY);
        glDrawArrays(GL_LINES, 0, color_cnt);
        glDisableClientState(GL_COLOR_ARRAY);
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_handles[2]);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glDrawArrays(GL_LINES, 0, texcoord_cnt);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
}
