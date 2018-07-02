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

#ifndef FONT_HPP
#define FONT_HPP

#ifdef GLVIS_USE_FREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H

#include <string>
#include <iostream>

using namespace std;


class GlVisFont
{
public:
    struct glyph {
        uint32_t w, h;
        int32_t bear_x, bear_y;
        int32_t adv_x, adv_y;
        float tex_x;
    };
private:
   bool init;
   bool font_init;
   glyph font_chars[256];
   float tex_w;
   float tex_h;
   uint32_t font_tex;

   FT_Library  library;
   FT_Face     face;
public:

    bool LoadFont(const char* path, int font_size);

    glyph GetTexChar(char c) {
        return font_chars[(uint8_t) c];
    }

    uint32_t BufferText(std::string& str);

    void RenderBuffer(uint32_t buf, double x, double y, double z);

    GlVisFont()
        : init(false)
        , font_init(false) {
        if (FT_Init_FreeType(&library)) {
            cout << "GLVis: Can not initialize FreeType library!" << endl;
        }
        init = true;
    }

    ~GlVisFont() {
        if (init)
            FT_Done_FreeType(library);
    }

    bool isFontLoaded() { return font_init; }

    float getAtlasWidth() { return tex_w; }

    float getAtlasHeight() { return tex_h; }
};
#endif /* GLVIS_USE_FREETYPE */
#endif /* FONT_HPP */
