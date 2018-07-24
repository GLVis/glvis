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
#include "openglvis.hpp"
#include <algorithm>
#include <iostream>
#include <utility>
using namespace gl3;

void GlBuilder::glEnd() {
#ifdef GLVIS_OGL3
    int pts_stride = 4;
    VertexBuffer::array_layout dst_layout = VertexBuffer::LAYOUT_VTX;
    if (is_line) {
        if (pts.size() < 2
            || (pts.size() == 2 && render_as == GL_LINE_LOOP)) {
            pts.clear();
            return;
        }
        if (use_color) {
            dst_layout = VertexBuffer::LAYOUT_VTX_COLOR;
        } else if (use_color_tex) {
            dst_layout = VertexBuffer::LAYOUT_VTX_TEXTURE0;
        }
#ifdef GLVIS_DEBUG
        if (pts.size() % pts_stride > 0) {
            cerr << "WARNING: GlBuilder stride does not cleanly divide into pts." << endl;
        }
#endif
        if (render_as == GL_LINE_LOOP) {
            //connect first and last points
            pts.reserve(pts_stride * 2 + pts.size());
            std::copy_n(pts.end() - pts_stride, pts_stride, std::back_inserter(pts));
            std::copy_n(pts.begin(), pts_stride, std::back_inserter(pts));
        }
        VertexBuffer& toInsert = parent_buf->getBuffer(dst_layout, GL_LINES);
        toInsert._pt_data.insert(toInsert._pt_data.end(),
                                 std::make_move_iterator(pts.begin()),
                                 std::make_move_iterator(pts.end()));
    } else {
        if (pts.size() < 3) {
            pts.clear();
            return;
        }
        dst_layout = VertexBuffer::LAYOUT_VTX_NORMAL;
        pts_stride = 6;
        if (use_color) {
            dst_layout = VertexBuffer::LAYOUT_VTX_NORMAL_COLOR;
            pts_stride = 8;
        } else if (use_color_tex) {
            dst_layout = VertexBuffer::LAYOUT_VTX_NORMAL_TEXTURE0;
            pts_stride = 8;
        }
#ifdef GLVIS_DEBUG
        if (pts.size() % pts_stride > 0) {
            cerr << "WARNING: GlBuilder stride does not cleanly divide into pts." << endl;
        }
#endif
        VertexBuffer& toInsert = parent_buf->getBuffer(dst_layout, GL_TRIANGLES);
        if (render_as == GL_POLYGON || render_as == GL_TRIANGLE_FAN) {
            //fan out from point 0
            int fan = pts_stride;
            while (fan + pts_stride * 2 <= pts.size()) {
                toInsert._pt_data.insert(toInsert._pt_data.end(),
                                         pts.begin(),
                                         pts.begin() + pts_stride);
                toInsert._pt_data.insert(toInsert._pt_data.end(),
                                         pts.begin() + fan,
                                         pts.begin() + fan + pts_stride * 2);
                fan += pts_stride;
            }
        } else if (render_as == GL_QUADS) {
            while (!pts.empty()) {
                toInsert._pt_data.insert(toInsert._pt_data.end(), pts.begin(), pts.begin() + pts_stride * 3);
                toInsert._pt_data.insert(toInsert._pt_data.end(), pts.begin(), pts.begin() + pts_stride);
                toInsert._pt_data.insert(toInsert._pt_data.end(), pts.begin() + pts_stride * 2, pts.begin() + pts_stride * 4);
                pts.erase(pts.begin(), pts.begin() + pts_stride * 4);
            }
        } else {
            toInsert._pt_data.insert(toInsert._pt_data.end(),
                                     std::make_move_iterator(pts.begin()),
                                     std::make_move_iterator(pts.end()));
        }
    }
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
    if (_allocated_size >= _pt_data.size()) {
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * _pt_data.size(), _pt_data.data());
    } else {
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _pt_data.size(), _pt_data.data(), GL_DYNAMIC_DRAW);
        _allocated_size = _pt_data.size();
    }
    _buffered_size = _pt_data.size();
}

//TODO: switch between vertexpointer and vertexattribpointer depending on available GL version
void VertexBuffer::drawObject(GLenum renderAs) {
    if (_buffered_size == 0) {
        return;
    }
    if (_layout == LAYOUT_VTX_TEXTURE0 || _layout == LAYOUT_VTX_NORMAL_TEXTURE0) {
        GetGlState()->setModeColorTexture();
    } else {
        GetGlState()->setModeColor();
    }

    glBindBuffer(GL_ARRAY_BUFFER, *_handle);
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    int loc_vtx = GetGlState()->getAttribLoc(GlState::ATTR_VERTEX);
    int loc_nor = GetGlState()->getAttribLoc(GlState::ATTR_NORMAL);
    int loc_color = GetGlState()->getAttribLoc(GlState::ATTR_COLOR);
    int loc_tex = GetGlState()->getAttribLoc(GlState::ATTR_TEXCOORD0);
    switch (_layout) {
        case LAYOUT_VTX:
            //glVertexPointer(3, GL_FLOAT, 0, 0);
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 4, 0);
            glDrawArrays(renderAs, 0, _buffered_size / 4);
            break;
        case LAYOUT_VTX_NORMAL:
            //glVertexPointer(3, GL_FLOAT, sizeof(float) * 6, 0);
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
            //glNormalPointer(GL_FLOAT, sizeof(float) * 6, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 6);
            GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
            break;
        case LAYOUT_VTX_COLOR:
            //glVertexPointer(3, GL_FLOAT, sizeof(float) * 7, 0);
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 4, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
            //glColorPointer(4, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float) * 4, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 4);
            GetGlState()->disableAttribArray(GlState::ATTR_COLOR);
            break;
        case LAYOUT_VTX_TEXTURE0:
            //glVertexPointer(3, GL_FLOAT, sizeof(float) * 5, 0);
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 4, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD0);
            //glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 5, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_tex, 1, GL_FLOAT, GL_FALSE, sizeof(float) * 4, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 4);
            GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD0);
            break;
        case LAYOUT_VTX_NORMAL_COLOR:
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
            //glNormalPointer(GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 3));
            GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
            //glColorPointer(4, GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 6));
            glVertexAttribPointer(loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(float) * 8, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 8);
            GetGlState()->disableAttribArray(GlState::ATTR_COLOR);
            GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
            break;
        case LAYOUT_VTX_NORMAL_TEXTURE0:
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
            //glNormalPointer(GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 3));
            GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD0);
            //glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 8, (void*)(sizeof(float) * 6));
            glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 8);
            GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD0);
            GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
            break;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void TextBuffer::bufferData() {
    std::vector<float> buf_data;
    int offset = 0;
    float tex_w = GetFont()->getAtlasWidth();
    float tex_h = GetFont()->getAtlasHeight();
    for (auto& e : _data) {
        float x = 0.f, y = 0.f;
        e.w = 0;
        e.h = 0;
        for (char c : e.text) {
            GlVisFont::glyph g = GetFont()->GetTexChar(c);
            float cur_x = x + g.bear_x;
            float cur_y = -y - g.bear_y;
            x += g.adv_x;
            y += g.adv_y;
            if (!g.w || !g.h) {
                continue;
            }
            float tris[] = {
                e.rx, e.ry, e.rz, cur_x,       -cur_y,       g.tex_x,               0,           0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y,       g.tex_x + g.w / tex_w, 0,           0,
                e.rx, e.ry, e.rz, cur_x,       -cur_y - g.h, g.tex_x,               g.h / tex_h, 0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y,       g.tex_x + g.w / tex_w, 0,           0,
                e.rx, e.ry, e.rz, cur_x,       -cur_y - g.h, g.tex_x,               g.h / tex_h, 0,
                e.rx, e.ry, e.rz, cur_x + g.w, -cur_y - g.h, g.tex_x + g.w / tex_w, g.h / tex_h, 0
            };
            buf_data.insert(buf_data.end(), tris, tris + 8 * 6 * sizeof(float));
            e.w = (int)(cur_x + g.w);
            e.h = std::max(e.h, (int)g.h);
        }
    }
    _size = buf_data.size() / 8;
    if (*_handle == 0) {
        glGenBuffers(1, _handle.get());
    }
    glBindBuffer(GL_ARRAY_BUFFER, *_handle);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * buf_data.size(), buf_data.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void TextBuffer::drawObject() {
    if (_size == 0 || *_handle == 0) {
        return;
    }
    GetGlState()->setModeRenderText();
    
    int loc_rast_vtx = GetGlState()->getAttribLoc(GlState::ATTR_VERTEX);
    int loc_txt_vtx = GetGlState()->getAttribLoc(GlState::ATTR_TEXT_VERTEX);
    int loc_tex = GetGlState()->getAttribLoc(GlState::ATTR_TEXCOORD1);

    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXT_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD1);
    
    glBindBuffer(GL_ARRAY_BUFFER, *_handle);

    glVertexAttribPointer(loc_rast_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, 0);
    glVertexAttribPointer(loc_txt_vtx, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 3));
    //glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 4, (void*)(sizeof(float) * 2));
    glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 5));
    glDrawArrays(GL_TRIANGLES, 0, _size);
    
    GetGlState()->disableAttribArray(GlState::ATTR_TEXT_VERTEX);
    GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GetGlState()->setModeColor();
}

void GlDrawable::addCone(float x, float y, float z,
                         float vx, float vy, float vz,
                         float cone_scale) {
    VertexBuffer& buf = getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL, GL_TRIANGLES);
    double rhos  = sqrt (vx*vx+vy*vy+vz*vz);
    float phi   = acos(vz/rhos);
    float theta = atan2 (vy, vx);

    glm::mat4 mtx(1.0);
    mtx = glm::translate(mtx, glm::vec3(x, y, z));
    mtx = glm::scale(mtx, glm::vec3(cone_scale));
    mtx = glm::translate(mtx, glm::vec3(x, y, z));
    mtx = glm::rotate(mtx, theta, glm::vec3(0.f, 0.f, 1.f));
    mtx = glm::rotate(mtx, phi, glm::vec3(0.f, 1.f, 0.f));
    glm::mat3 norm(mtx);
    norm = glm::inverseTranspose(norm);

    glm::vec3 start_vtx = glm::vec3(mtx * glm::vec4(0.f, 0.f, 0.f, 1.f));
    glm::vec3 start_norm = glm::vec3(norm * glm::vec3(0.f, 0.f, 1.f));

    glm::vec3 base_pts[] = {
        glm::vec3(mtx * glm::vec4(1, 0, -4, 1)),
        glm::vec3(mtx * glm::vec4(cos(2*M_PI/4), sin(2*M_PI/4), -4, 1)),
        glm::vec3(mtx * glm::vec4(cos(4*M_PI/4), sin(4*M_PI/4), -4, 1)),
        glm::vec3(mtx * glm::vec4(cos(6*M_PI/4), sin(6*M_PI/4), -4, 1)),
    };

    float nz = (1.0/4.0);
    glm::vec3 base_norms[] = {
        glm::vec3(norm * glm::vec3(1, 0, nz)),
        glm::vec3(norm * glm::vec3(cos(2*M_PI/4), sin(2*M_PI/4), nz)),
        glm::vec3(norm * glm::vec3(cos(4*M_PI/4), sin(4*M_PI/4), nz)),
        glm::vec3(norm * glm::vec3(cos(6*M_PI/4), sin(6*M_PI/4), nz)),
    };
    
    float* orig = glm::value_ptr(start_vtx);
    float* orig_n = glm::value_ptr(start_norm);
    float* base[4] = {
        glm::value_ptr(base_pts[0]),
        glm::value_ptr(base_pts[1]),
        glm::value_ptr(base_pts[2]),
        glm::value_ptr(base_pts[3])
    };
    float* base_n[4] = {
        glm::value_ptr(base_norms[0]),
        glm::value_ptr(base_norms[1]),
        glm::value_ptr(base_norms[2]),
        glm::value_ptr(base_norms[3])
    };

    std::vector<float> cone_pts;
    for (int i = 0; i < 4; i++) {
        buf.addVertexNorm({orig[0],   orig[1],   orig[2]},
                          {orig_n[0], orig_n[1], orig_n[2]});
        buf.addVertexNorm({base[i][0],   base[i][1],   base[i][2]},
                          {base_n[i][0], base_n[i][1], base_n[i][2]});
        buf.addVertexNorm({base[(i+1)%4][0],   base[(i+1)%4][1],   base[(i+1)%4][2]},
                          {base_n[(i+1)%4][0], base_n[(i+1)%4][1], base_n[(i+1)%4][2]});
    }
}

