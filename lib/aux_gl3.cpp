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
#include <iostream>
#include <cstddef>

using namespace gl3;

void Vertex::setupAttribLayout() {
    GetGlState()->setModeColor();
    int loc_vtx = GlState::ATTR_VERTEX;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, coord));
}

void VertexColor::setupAttribLayout() {
    GetGlState()->setModeColor();
    int loc_vtx = GlState::ATTR_VERTEX;
    int loc_color = GlState::ATTR_COLOR;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(VertexColor), (void*)(void*)offsetof(VertexColor, coord));
    glVertexAttribPointer(loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(VertexColor), (void*)offsetof(VertexColor, color));
}

void VertexColor::clearAttribLayout() {
    GetGlState()->disableAttribArray(GlState::ATTR_COLOR);
}

void VertexTex::setupAttribLayout() {
    GetGlState()->setModeColorTexture();
    int loc_vtx = GlState::ATTR_VERTEX;
    int loc_tex = GlState::ATTR_TEXCOORD0;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD0);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(VertexTex), (void*)offsetof(VertexTex, coord));
    glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(VertexTex), (void*)offsetof(VertexTex, texCoord));
}

void VertexTex::clearAttribLayout() {
    GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD0);
}

void VertexNorm::setupAttribLayout() {
    GetGlState()->setModeColor();
    int loc_vtx = GlState::ATTR_VERTEX;
    int loc_nor = GlState::ATTR_NORMAL;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNorm), (void*)offsetof(VertexNorm, coord));
    glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNorm), (void*)offsetof(VertexNorm, norm));
}

void VertexNorm::clearAttribLayout() {
    GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
}

void VertexNormColor::setupAttribLayout() {
    GetGlState()->setModeColor();
    int loc_vtx = GlState::ATTR_VERTEX;
    int loc_nor = GlState::ATTR_NORMAL;
    int loc_color = GlState::ATTR_COLOR;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
    GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormColor), (void*)offsetof(VertexNormColor, coord));
    glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormColor), (void*)offsetof(VertexNormColor, norm));
    glVertexAttribPointer(loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(VertexNormColor), (void*)offsetof(VertexNormColor, color));
}

void VertexNormColor::clearAttribLayout() {
    GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
    GetGlState()->disableAttribArray(GlState::ATTR_COLOR);
}

void VertexNormTex::setupAttribLayout() {
    GetGlState()->setModeColorTexture();
    int loc_vtx = GlState::ATTR_VERTEX;
    int loc_nor = GlState::ATTR_NORMAL;
    int loc_tex = GlState::ATTR_TEXCOORD0;
    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD0);
    glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormTex), (void*)offsetof(VertexNormTex, coord));
    glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormTex), (void*)offsetof(VertexNormTex, norm));
    glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(VertexNormTex), (void*)offsetof(VertexNormTex, texCoord));
}

void VertexNormTex::clearAttribLayout() {
    GetGlState()->disableAttribArray(GlState::ATTR_NORMAL);
    GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD0);
}

void TextBuffer::buffer() {
    std::vector<float> buf_data;
    float tex_w = GetFont()->getAtlasWidth();
    float tex_h = GetFont()->getAtlasHeight();
    for (auto& e : _data) {
        float x = 0.f, y = 0.f;
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
            buf_data.insert(buf_data.end(), tris, tris + 8 * 6);
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

void TextBuffer::getObjectSize(const std::string& text, int& w, int& h) {
    float x = 0.f;
    w = 0.f, h = 0.f;
    for (char c : text) {
        GlVisFont::glyph g = GetFont()->GetTexChar(c);
        float cur_x = x + g.bear_x;
        x += g.adv_x;
        if (!g.w || !g.h) {
            continue;
        }
        w = (int)(cur_x + g.w);
        h = std::max(h, (int)g.h);
    }
}

void TextBuffer::draw() {
    if (_size == 0 || *_handle == 0) {
        return;
    }
    GetGlState()->setModeRenderText();
    
    int loc_rast_vtx = GlState::ATTR_VERTEX;
    int loc_txt_vtx = GlState::ATTR_TEXT_VERTEX;
    int loc_tex = GlState::ATTR_TEXCOORD1;

    GetGlState()->enableAttribArray(GlState::ATTR_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXT_VERTEX);
    GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD1);
    
    glBindBuffer(GL_ARRAY_BUFFER, *_handle);

    glVertexAttribPointer(loc_rast_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, 0);
    glVertexAttribPointer(loc_txt_vtx, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 3));
    glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 5));
    glDrawArrays(GL_TRIANGLES, 0, _size);
    
    GetGlState()->disableAttribArray(GlState::ATTR_TEXT_VERTEX);
    GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GetGlState()->setModeColor();
}

IDrawHook * GlDrawable::buf_hook = nullptr;

void GlDrawable::addCone(float x, float y, float z,
                         float vx, float vy, float vz,
                         float cone_scale) {
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
        addTriangle(
            VertexNorm {
                {orig[0],   orig[1],   orig[2]},
                {orig_n[0], orig_n[1], orig_n[2]}
            },
            VertexNorm {
                {base[i][0],   base[i][1],   base[i][2]},
                {base_n[i][0], base_n[i][1], base_n[i][2]}
            },
            VertexNorm {
                {base[(i+1)%4][0],   base[(i+1)%4][1],   base[(i+1)%4][2]},
                {base_n[(i+1)%4][0], base_n[(i+1)%4][1], base_n[(i+1)%4][2]}
            }
        );
    }
}

void GlBuilder::saveVertex(const GlBuilder::_vertex& v) {
    GLenum dst_buf = is_line ? GL_LINES : GL_TRIANGLES;
    if (!use_norm) {
        if (use_color) {
            parent_buf->getBuffer<VertexColor>(dst_buf)
                      ->addVertex(VertexColor{v.coords, v.color});
        } else if (use_tex) {
            parent_buf->getBuffer<VertexTex>(dst_buf)
                      ->addVertex(VertexTex{v.coords, v.texcoord});
        } else {
            parent_buf->getBuffer<Vertex>(dst_buf)
                      ->addVertex(Vertex{v.coords});
        }
    } else {
        if (use_color) {
            parent_buf->getBuffer<VertexNormColor>(dst_buf)
                      ->addVertex(VertexNormColor{v.coords, v.norm, v.color});
        } else if (use_tex) {
            parent_buf->getBuffer<VertexNormTex>(dst_buf)
                      ->addVertex(VertexNormTex{v.coords, v.norm, v.texcoord});
        } else {
            parent_buf->getBuffer<VertexNorm>(dst_buf)
                      ->addVertex(VertexNorm{v.coords, v.norm});
        }
    }
}

