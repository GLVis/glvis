// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_FONT_HPP
#define GLVIS_FONT_HPP

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H

#include <string>
#include <iostream>

#include "gl/platform_gl.hpp"

using namespace std;

class GlVisFont
{
public:
   struct glyph
   {
      uint32_t w, h;
      int32_t bear_x, bear_y;
      float adv_x, adv_y;
      float tex_x;
   };
private:
   bool init;
   bool font_init;

   GLenum alpha_channel;

   glyph font_chars[256];
   float tex_w;
   float tex_h;
   uint32_t font_tex;

   FT_Library  library;
   FT_Face     face;
   bool        face_has_kerning;
public:

   bool LoadFont(const std::string& path, int font_index, int font_size);

   const glyph &GetTexChar(char c) const
   {
      return font_chars[(uint8_t) c];
   }

   /// Get the width and height of the bounding box containing the rendered text
   void getObjectSize(const std::string& text, int& w, int& h);

   GlVisFont()
      : init(false),
        font_init(false),
        face_has_kerning(false)
   {
      if (FT_Init_FreeType(&library))
      {
         cout << "GLVis: Can not initialize FreeType library!" << endl;
      }
      init = true;
   }

   ~GlVisFont()
   {
      if (init)
      {
         FT_Done_FreeType(library);
      }
   }

   float GetKerning(char cprev, char c)
   {
      if (!face_has_kerning || cprev == '\0')
      {
         return 0;
      }
      FT_UInt glyph_prev = FT_Get_Char_Index(face, cprev);
      FT_UInt glyph_curr = FT_Get_Char_Index(face, c);

      FT_Vector delta;
      FT_Get_Kerning(face, glyph_prev, glyph_curr,
                     FT_KERNING_DEFAULT, &delta);
      return delta.x / 64.f;
   }

   bool isFontLoaded() { return font_init; }

   float getFontLineSpacing() const { return face->size->metrics.height/64.f; }

   float getFontDescender() const { return face->size->metrics.descender/64.f; }

   float getAtlasWidth() { return tex_w; }

   float getAtlasHeight() { return tex_h; }

   unsigned getFontTex() { return font_tex; }

   void setAlphaChannel(GLenum alpha) { alpha_channel = alpha; }
};

#endif /* GLVIS_FONT_HPP */
