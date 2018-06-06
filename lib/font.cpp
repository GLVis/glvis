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

#include <vector>

#include "font.hpp"
#include "aux_vis.hpp"
#include "platform_gl.hpp"

struct vert_tex2d {
    float x, y;
    float u, v;
    vert_tex2d() = default;
    vert_tex2d(float x, float y, float u, float v)
        : x(x), y(y), u(u), v(v) { }
};


bool GlVisFont::LoadFont(const char* path, int font_size) {
    if (!init) {
        return false;
    }
    if (font_init) {
        glDeleteTextures(1, &font_tex);
        FT_Done_Face(face);
        font_init = false;
    }
    if (FT_New_Face(library, path, 0, &face)) {
        cout << "GLVis: Cannot open font file: " << path << endl;
        return false;
    }
    int ppi_w, ppi_h;
    GetAppWindow()->getDpi(ppi_w, ppi_h);
    if (FT_Set_Char_Size(face, 0, font_size*64, ppi_w, ppi_h)) {
        cout << "GLVis: Cannot set font height: " << font_size << " pts"
             << endl;
        FT_Done_Face(face);
        return false;
    }

    //generate atlas
    int w = 0, h = 0;
    for (int c = 32; c < 128; c++) {
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
            cout << "GLVis: Cannot load glyph: " << (char) c << endl;
            continue;
        }
        w += face->glyph->bitmap.width;
        if (h < face->glyph->bitmap.rows) {
            h = face->glyph->bitmap.rows;
        }
    }
    tex_w = w;
    tex_h = h;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, w, h, 0, GL_ALPHA, GL_UNSIGNED_BYTE, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    int x = 0;
    for (int c = 32; c < 128; c++) {
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
            cout << "GLVis: Cannot load glyph: " << (char) c << endl;
            continue;
        }
        glTexSubImage2D(GL_TEXTURE_2D,
                        0,
                        x, 0,
                        face->glyph->bitmap.width,
                        face->glyph->bitmap.rows,
                        GL_ALPHA,
                        GL_UNSIGNED_BYTE,
                        face->glyph->bitmap.buffer);
        font_chars[c] = {
            face->glyph->bitmap.width,
            face->glyph->bitmap.rows,
            face->glyph->bitmap_left,
            face->glyph->bitmap_top,
            (int)(face->glyph->advance.x >> 6),
            (int)(face->glyph->advance.y >> 6),
            (float) x / w
        };
        x += face->glyph->bitmap.width;
    }
    font_init = true;
    return true;
}

uint32_t GlVisFont::BufferText(std::string& str) {
    std::vector<vert_tex2d> coordData(str.size() * 6);
    float x = 0.0, y = 0.0;
    for (char& c : str) {
        glyph g = font_chars[(uint8_t) c];
        float cur_x = x + g.bear_x;
        float cur_y = -y - g.bear_y;
        x += g.adv_x;
        y += g.adv_y;
        if (!g.w || !g.h) {
            continue;
        }
        coordData.emplace_back(cur_x,       -cur_y,       g.tex_x,               0);
        coordData.emplace_back(cur_x + g.w, -cur_y,       g.tex_x + g.w / tex_w, 0);
        coordData.emplace_back(cur_x,       -cur_y - g.h, g.tex_x,               g.h / tex_h);
        coordData.emplace_back(cur_x + g.w, -cur_y,       g.tex_x + g.w / tex_w, 0);
        coordData.emplace_back(cur_x,       -cur_y - g.h, g.tex_x,               g.h / tex_h);  
        coordData.emplace_back(cur_x + g.w, -cur_y - g.h, g.tex_x + g.w / tex_w, g.h / tex_h);
    }
    GLuint buf;
    glGenBuffers(1, &buf);
    glBindBuffer(GL_ARRAY_BUFFER, buf);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vert_tex2d) * coordData.size(), coordData.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return buf;
}

void GlVisFont::RenderBuffer(uint32_t buf) {
    
    GLint viewport[4];
    GLdouble raster[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_CURRENT_RASTER_POSITION, raster);

    glPushAttrib(GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, viewport[2], 0, viewport[3], -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glTranslated(raster[0], raster[1], raster[2]);

    glEnable(GL_TEXTURE_2D);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, buf);
    glVertexPointer(2, GL_FLOAT, sizeof(float) * 4, 0);
    glTexCoordPointer(2, GL_FLOAT, sizeof(float) * 4, (void*)(sizeof(float) * 2));
    GLint size = 0;
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
    glDrawArrays(GL_TRIANGLES, 0, size / (sizeof(float) * 4)); 
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glPopAttrib();

}

