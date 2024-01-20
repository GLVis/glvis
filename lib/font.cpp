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

#include <vector>

#include "font.hpp"
#include "aux_vis.hpp"

struct vert_tex2d
{
   float x, y;
   float u, v;

   vert_tex2d(float x, float y, float u, float v)
      : x(x), y(y), u(u), v(v) { }
};


const int font_scale = 64;
bool GlVisFont::LoadFont(const std::string& path, int font_index, int font_size)
{
   if (!init)
   {
      return false;
   }
   if (font_init)
   {
      glDeleteTextures(1, &font_tex);
      FT_Done_Face(face);
      font_init = false;
   }
   if (FT_New_Face(library, path.c_str(), font_index, &face))
   {
      cout << "GlVisFont::LoadFont(): Cannot open font file: " << path << endl;
      return false;
   }
   face_has_kerning = FT_HAS_KERNING(face);
   int ppi_w, ppi_h;
   GetAppWindow()->getDpi(ppi_w, ppi_h);
   const bool use_fixed_ppi_h = true;
   if (use_fixed_ppi_h)
   {
      double ratio = double(ppi_w)/ppi_h;
      ppi_h = GetAppWindow()->isHighDpi() ? 192 : 96;
      ppi_w = ratio*ppi_h + 0.5;
#ifdef GLVIS_DEBUG
      cout << "Fonts use fixed ppi: " << ppi_w << " x " << ppi_h << endl;
#endif
   }
   if (FT_Set_Char_Size(face, 0, font_size*font_scale, ppi_w, ppi_h))
   {
      cout << "GlVisFont::LoadFont(): Cannot set font height: " << font_size
           << " pts" << endl;
      FT_Done_Face(face);
      return false;
   }
#ifdef GLVIS_DEBUG
   cout << "Loaded font: " << path << ", index: " << font_index
        << ", height: " << font_size << endl;
#endif

   // generate atlas
   size_t w = 0, h = 0;
   for (int c = 32; c < 128; c++)
   {
      if (FT_Load_Char(face, c, FT_LOAD_RENDER))
      {
         cout << "GlVisFont::LoadFont(): Cannot load glyph: " << (char) c << endl;
         continue;
      }
      w += face->glyph->bitmap.width + 2;
      if (h < size_t(face->glyph->bitmap.rows))
      {
         h = face->glyph->bitmap.rows;
      }
   }
   tex_w = w;
   tex_h = h + 2;

#ifdef GLVIS_DEBUG
   cout << "GlVisFont::LoadFont(): font texture dimensions are ("
        << tex_w << ", " << tex_h << ")" << endl;
#endif

   glGenTextures(1, &font_tex);

   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, font_tex);
   std::vector<uint8_t> zeros(tex_w * tex_h, 0);
   GLenum alpha_internal_fmt = alpha_channel == GL_RED ? GL_R8 : GL_ALPHA;
   glTexImage2D(GL_TEXTURE_2D, 0, alpha_internal_fmt, tex_w, tex_h, 0,
                alpha_channel,
                GL_UNSIGNED_BYTE, zeros.data());

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

   int x = 0;
   for (int c = 32; c < 128; c++)
   {
      if (FT_Load_Char(face, c, FT_LOAD_RENDER))
      {
         cout << "GlVisFont::LoadFont(): Cannot load glyph: " << (char) c << endl;
         continue;
      }
      FT_GlyphSlot glyph_slot = face->glyph;
      glTexSubImage2D(GL_TEXTURE_2D,
                      0,
                      x + 1, 1,
                      glyph_slot->bitmap.width,
                      glyph_slot->bitmap.rows,
                      alpha_channel,
                      GL_UNSIGNED_BYTE,
                      glyph_slot->bitmap.buffer);
      // To improve the rendering quality in some cases, we can utilize the
      // 'lsb_delta' and 'rsb_delta' fields in 'glyph_slot'. These will need to
      // be stored in the GlVisFont::glyph struct and used in the method
      // bufferToDevice() of class CoreGLDevice and class FFGLDevice.
      font_chars[c] =
      {
         (unsigned)(glyph_slot->bitmap.width + 2),
         (unsigned)(glyph_slot->bitmap.rows + 2),
         glyph_slot->bitmap_left,
         glyph_slot->bitmap_top,
         glyph_slot->advance.x / 64.f,
         glyph_slot->advance.y / 64.f,
         x / (float) w
      };
      x += glyph_slot->bitmap.width + 2;
   }
   font_init = true;
   // glEnable(GL_TEXTURE_2D);
   glActiveTexture(GL_TEXTURE0);
   return true;
}

void GlVisFont::getObjectSize(const std::string& text, int& w, int& h)
{
   float pen_x = 0.f, pen_y = 0.f;
   char prev_c = '\0';
   int min_x = INT_MAX, max_x = INT_MIN;
   int min_y = INT_MAX, max_y = INT_MIN;
   for (char c : text)
   {
      const glyph &g = GetTexChar(c);
      if (!g.w || !g.h) { continue; }
      pen_x += GetKerning(prev_c, c);
      // note: both g.w and g.h include padding of 2
      // note: the x-direction bounds are not pixel-tight (include bear_x/adv_x)
      // note: the y-direction bounds are pixel-tight
      min_x = std::min(min_x, (int)floorf(pen_x));
      max_x = std::max(max_x, (int)ceilf (pen_x + g.adv_x));
      min_y = std::min(min_y, (int)floorf(pen_y + g.bear_y - (g.h - 2)));
      max_y = std::max(max_y, (int)ceilf (pen_y + g.bear_y));
      pen_x += g.adv_x;
      pen_y += g.adv_y;
      prev_c = c;
   }
   w = max_x - min_x;
   h = max_y - min_y;
}
