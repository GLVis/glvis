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
            pts_stride = 8;
        } else if (use_color_tex) {
            dst_layout = VertexBuffer::LAYOUT_VTX_TEXTURE0;
            pts_stride = 5;
        }
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
            pts_stride = 10;
        } else if (use_color_tex) {
            dst_layout = VertexBuffer::LAYOUT_VTX_NORMAL_TEXTURE0;
            pts_stride = 8;
        }
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
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * _pt_data.size(), _pt_data.data(), GL_STATIC_DRAW);
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
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 8, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
            //glColorPointer(4, GL_FLOAT, sizeof(float) * 7, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_color, 4, GL_FLOAT, GL_FALSE, sizeof(float) * 8, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 8);
            GetGlState()->disableAttribArray(GlState::ATTR_COLOR);
            break;
        case LAYOUT_VTX_TEXTURE0:
            //glVertexPointer(3, GL_FLOAT, sizeof(float) * 5, 0);
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 5, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_TEXCOORD0);
            //glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 5, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_tex, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 5, (void*)(sizeof(float) * 3));
            glDrawArrays(renderAs, 0, _buffered_size / 5);
            GetGlState()->disableAttribArray(GlState::ATTR_TEXCOORD0);
            break;
        case LAYOUT_VTX_NORMAL_COLOR:
            glVertexAttribPointer(loc_vtx, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 10, 0);
            GetGlState()->enableAttribArray(GlState::ATTR_NORMAL);
            //glNormalPointer(GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 3));
            glVertexAttribPointer(loc_nor, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 10, (void*)(sizeof(float) * 3));
            GetGlState()->enableAttribArray(GlState::ATTR_COLOR);
            //glColorPointer(4, GL_FLOAT, sizeof(float) * 10, (void*)(sizeof(float) * 6));
            glVertexAttribPointer(loc_color, 4, GL_FLOAT, GL_FALSE, sizeof(float) * 10, (void*)(sizeof(float) * 6));
            glDrawArrays(renderAs, 0, _buffered_size / 10);
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

/*
void ArrowHelper(double px, double py, double pz,
                 double vx, double vy, double vz,
                 double length, double cone_scale,
                 double cone[][3], double normal[][3]) {
   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   double phi = acos(vz/rhos), theta = atan2(vy, vx);
   const int n = 8;
   const double step = 2*M_PI/n, nz = (1.0/4.0);
   double point = step;
   int i, j, k;

   cone[0][0] = 0;          cone[0][1] = 0; cone[0][2] = 1;
   cone[1][0] = cone_scale; cone[1][1] = 0; cone[1][2] = -4*cone_scale + 1;
   normal[0][0] = 0.0/cone_scale;
   normal[0][1] = 0.0/cone_scale;
   normal[0][2] = 1.0/cone_scale;
   normal[1][0] = 1.0/cone_scale;
   normal[1][1] = 0.0/cone_scale;
   normal[1][2] = nz/cone_scale;

   for (i=2; i<n+1; i++)
   {
      normal[i][0] = cos(point)/cone_scale;
      normal[i][1] = sin(point)/cone_scale;
      normal[i][2] = nz/cone_scale;

      cone[i][0] = cos(point)*cone_scale;
      cone[i][1] = sin(point)*cone_scale;
      cone[i][2] = -4*cone_scale + 1;
      point += step;
   }
   cone[n+1][0] = cone_scale; cone[n+1][1] = 0; cone[n+1][2] =-4*cone_scale + 1;
   normal[n+1][0] = 1.0/cone_scale;
   normal[n+1][1] = 0.0/cone_scale;
   normal[n+1][2] = nz/cone_scale;

   cone[n+2][0] = 0; cone[n+2][1] = 0; cone[n+2][2] = 0;
   cone[n+3][0] = 0; cone[n+3][1] = 0; cone[n+3][2] = 1;

   if (arrow_scaling_type == 0)
   {
      length = rhos;
   }

   // double xc = 0.5*(x[0]+x[1]), yc = 0.5*(y[0]+y[1]), zc = 0.5*(z[0]+z[1]);
   double coord[3];
   // double rlen = length/rhos;

   // px = (px-xc)*xscale;  py = (py-yc)*yscale;  pz = (pz-zc)*zscale;
   // vx *= rlen*xscale;    vy *= rlen*yscale;    vz *= rlen*zscale;

   if (arrow_type == 1)
      for (i=0; i<n+4; i++)
      {
         cone[i][2] -= 0.5;
      }

   double M[3][3]= {{cos(theta)*cos(phi), -sin(theta),  cos(theta)*sin(phi)},
      {sin(theta)*cos(phi),  cos(theta),  sin(theta)*sin(phi)},
      {          -sin(phi),          0.,             cos(phi)}
   };
   double v[3] = { M[0][2]/xscale, M[1][2]/yscale, M[2][2]/zscale };
   length /= sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

   for (i=0; i<n+4; i++)
   {
      for (j=0; j<3; j++)
      {
         coord[j] = cone[i][j] * length;
      }

      for (k=0; k<3; k++)
      {
         cone[i][k] = 0.;
         for (j=0; j<3; j++)
         {
            cone[i][k] += M[k][j] * coord[j];
         }
      }
      // cone[i][0] = (cone[i][0] + px)/xscale + xc;
      // cone[i][1] = (cone[i][1] + py)/yscale + yc;
      // cone[i][2] = (cone[i][2] + pz)/zscale + zc;
      cone[i][0] = cone[i][0]/xscale + px;
      cone[i][1] = cone[i][1]/yscale + py;
      cone[i][2] = cone[i][2]/zscale + pz;
   }

   for (i=0; i<=n+1; i++)
   {
      for (j=0; j<3; j++)
      {
         coord[j] = normal[i][j];
      }

      for (k=0; k<3; k++)
      {
         normal[i][k] = 0.;
         for (j=0; j<3; j++)
         {
            normal[i][k] += M[k][j] * coord[j];
         }
      }
      normal[i][0] *= xscale;
      normal[i][1] *= yscale;
      normal[i][2] *= zscale;
   }

}

void GlDrawable::addArrow(double px, double py, double pz,
                          double vx, double vy, double vz,
                          double length, double cone_scale) {
    double rhos = sqrt (vx*vx+vy*vy+vz*vz);
    if (rhos == 0.0)
    {
        return;
    }
    const int n = 8;
    double cone[n+4][3], normal[n+2][3];
    ArrowHelper(px, py, pz, vx, vy, vz,
                length, cone_scale, cone, normal);
    VertexBuffer& buf_tris = this->getBuffer(VertexBuffer::LAYOUT_VTX_NORMAL, GL_TRIANGLES);
    VertexBuffer& buf_lines = this->getBuffer(VertexBuffer::LAYOUT_VTX, GL_LINES);
    for (i = 2; i <= n+1; i++) {
        buf_tris.addVertex(cone[0], normal[0]);
        buf_tris.addVertex(cone[i-1], normal[i-1]);
        buf_tris.addVertex(cone[i], normal[i]);
    }
    buf_lines.addVertex(cone[n+2]);
    buf_lines.addVertex(cone[n+3]);
}

void GlDrawable::addArrow(double px, double py, double pz,
                          double vx, double vy, double vz,
                          double length, double cone_scale,
                          GlColor color) {
    double rhos = sqrt (vx*vx+vy*vy+vz*vz);
    if (rhos == 0.0)
    {
        return;
    }

    const int n = 8;
    double cone[n+4][3], normal[n+2][3];
    ArrowHelper(px, py, pz, vx, vy, vz,
            length, cone_scale, cone, normal);
    VertexBuffer::array_layout buf_layout;
    if (color.use_texture) {
        buf_layout = VertexBuffer::LAYOUT_VTX_NORMAL_TEXTURE0;
    } else {
        buf_layout = VertexBuffer::LAYOUT_VTX_NORMAL_COLOR;
    }

    VertexBuffer& buf_tris = this->getBuffer(buf_layout, GL_TRIANGLES);
    VertexBuffer& buf_lines = this->getBuffer(buf_layout, GL_LINES);
    for (i = 2; i <= n+1; i++) {
        if (color.use_texture) {
            buf_tris.addVertex(cone[0], normal[0], color.texcoord);
            buf_tris.addVertex(cone[i-1], normal[i-1], color.texcoord);
            buf_tris.addVertex(cone[i], normal[i], color.texcoord);
        } else {
            buf_tris.addVertex(cone[0], normal[0], color.rgba);
            buf_tris.addVertex(cone[i-1], normal[i-1], color.rgba);
            buf_tris.addVertex(cone[i], normal[i], color.rgba);
        }
    }
    if (color.use_texture) {
        buf_lines.addVertex(cone[n+2], normal[n+1], color.texcoord);
        buf_lines.addVertex(cone[n+3], normal[n+1], color.texcoord);
    } else {
        buf_lines.addVertex(cone[n+2], normal[n+1], color.rgba);
        buf_lines.addVertex(cone[n+3], normal[n+1], color.rgba);
    }
}

*/
