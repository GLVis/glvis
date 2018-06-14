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
        parent_buf->texcoord_data.insert(parent_buf->texcoord_data.end(),
                                         std::make_move_iterator(pts.begin()),
                                         std::make_move_iterator(pts.end()));
    } else if (has_color) {
        parent_buf->color_data.insert(parent_buf->color_data.end(),
                                      std::make_move_iterator(pts.begin()),
                                      std::make_move_iterator(pts.end()));
    } else {
        parent_buf->pt_data.insert(parent_buf->pt_data.end(),
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

void VertexBuffer::BufferData() {
    if (!pt_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(0));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        pt_cnt = pt_data.size() / 3;
    }
    if (!color_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(1));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * color_data.size(), color_data.data(), GL_STATIC_DRAW);
        color_cnt = color_data.size() / 10;
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(2));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        texcoord_cnt = texcoord_data.size() / 8;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void VertexBuffer::DrawObject(GLenum renderAs) {
    if (!pt_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(0));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glDrawArrays(renderAs, 0, pt_cnt);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (!color_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(1));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 10, 0);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 3));
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 6));
        glDrawArrays(renderAs, 0, color_cnt);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (!texcoord_data.empty()) {
        GetGlState()->setModeColorTexture();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(2));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 8, 0);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 3));
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 6));
        glDrawArrays(renderAs, 0, texcoord_cnt);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void TextBuffer::SetText(float x, float y, float z, std::string text) {
    entries.emplace_back(x, y, z, std::move(text));
}

void TextBuffer::BufferData() {
    LineBuffer::BufferData();
}

void TextBuffer::DrawObject(GLenum renderAs) {
    this->DrawObject();
}

void TextBuffer::DrawObject() {
    LineBuffer::DrawObject();
    if (entries.size() == 0) { return; }
#ifndef GLVIS_USE_FREETYPE
    cerr << "Can't use text buffer object without Freetype" << endl;
#else
    for (auto& str_obj : entries) {
        DrawBitmapText(str_obj.str.c_str(), str_obj.x, str_obj.y, str_obj.z);
    }
#endif
}

void LineBuffer::BufferData() {
    if (!pt_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(0));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pt_data.size(), pt_data.data(), GL_STATIC_DRAW);
        pt_cnt = pt_data.size() / 3;
    }
    if (!color_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(1));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * color_data.size(), color_data.data(), GL_STATIC_DRAW);
        color_cnt = color_data.size() / 7;
    }
    if (!texcoord_data.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(2));
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * texcoord_data.size(), texcoord_data.data(), GL_STATIC_DRAW);
        texcoord_cnt = texcoord_data.size() / 3;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void LineBuffer::DrawObject(GLenum renderAs) {
    this->DrawObject();
}

void LineBuffer::DrawObject() {
    if (!pt_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(0));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glDrawArrays(GL_LINES, 0, pt_cnt);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (!color_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(1));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(float) * 7, 0);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
        glDrawArrays(GL_LINES, 0, color_cnt);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (!texcoord_data.empty()) {
        GetGlState()->setModeColor();
        glBindBuffer(GL_ARRAY_BUFFER, vbo->get(2));
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        //glLineStipple(1, 255);
        //glEnable(GL_LINE_STIPPLE);
        glDrawArrays(GL_LINES, 0, texcoord_cnt);
        //glDisable(GL_LINE_STIPPLE);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
