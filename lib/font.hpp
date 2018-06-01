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
#include <fontconfig/fontconfig.h>

#include <string>
#include <iostream>

using namespace std;

#include "platform_gl.hpp"

class GLVisFont
{
private:
   int         init;

   FT_Library  library;
   FT_Face     face;

   FT_Bool     use_kerning;
   FT_Glyph   *glyph;
   FT_UInt     alloc_glyphs, num_glyphs;
   FT_BBox     bbox;

   int            image_width, image_height;
   unsigned char *image;

#ifdef GLVIS_OGL3
   //TODO: if you want unicode characters, don't do this
   
#endif

   void LoadSequence(const char *text)
   {
      int          err;
      FT_UInt      len, glyph_idx, prev_glyph, i;
      FT_Vector    pen, delta;

      FreeSeq();

      len      = strlen(text);
      glyph    = new FT_Glyph[len];
      alloc_glyphs = len;

      pen.x = pen.y = 0;
      prev_glyph = i = 0;
      for (FT_UInt ic = 0; ic < len; ic++)
      {
         glyph_idx = FT_Get_Char_Index(face, text[ic]);

         if (use_kerning && prev_glyph && glyph_idx)
         {
            FT_Get_Kerning(face, prev_glyph, glyph_idx,
                           FT_KERNING_DEFAULT, &delta);
            pen.x += delta.x;
         }

         err = FT_Load_Glyph(face, glyph_idx, FT_LOAD_DEFAULT);
         if (err) { continue; }

         err = FT_Get_Glyph(face->glyph, &glyph[i]);
         if (err) { continue; }

         FT_Glyph_Transform(glyph[i], 0, &pen);

         pen.x += face->glyph->advance.x;
         prev_glyph = glyph_idx;

         i++;
      }
      num_glyphs = i;
   }

   void GetBBox()
   {
      FT_BBox glyph_bbox;

      if (num_glyphs <= 0)
      {
         bbox.xMin = bbox.yMin = 0;
         bbox.xMax = bbox.yMax = 0;
         return;
      }

      FT_Glyph_Get_CBox(glyph[0], FT_GLYPH_BBOX_PIXELS, &bbox);
      for (FT_UInt i = 1; i < num_glyphs; i++)
      {
         FT_Glyph_Get_CBox(glyph[i], FT_GLYPH_BBOX_PIXELS, &glyph_bbox);

         if (glyph_bbox.xMin < bbox.xMin)
         {
            bbox.xMin = glyph_bbox.xMin;
         }
         if (glyph_bbox.xMax > bbox.xMax)
         {
            bbox.xMax = glyph_bbox.xMax;
         }
         if (glyph_bbox.yMin < bbox.yMin)
         {
            bbox.yMin = glyph_bbox.yMin;
         }
         if (glyph_bbox.yMax > bbox.yMax)
         {
            bbox.yMax = glyph_bbox.yMax;
         }
      }
   }

   void FreeSeq()
   {
      if (alloc_glyphs <= 0)
      {
         return;
      }

      for (FT_UInt i = 0; i < num_glyphs; i++)
      {
         FT_Done_Glyph(glyph[i]);
      }

      delete [] glyph;

      glyph = NULL;

      alloc_glyphs = num_glyphs = 0;
   }

   static int FindFontFile(const char *font_patterns[], int num_patterns,
                           string &font_file)
   {
      FcObjectSet *os  = NULL;
      FcFontSet   *fs  = NULL;
      FcPattern   *pat = NULL;
      FcChar8     *s   = NULL;
#ifdef GLVIS_DEBUG
      FcChar8     *fnt;
#endif
      FcResult     res;

      if (num_patterns <= 0)
      {
         return -2;
      }

      if (!FcInit())
      {
         return -1;
      }

      os = FcObjectSetBuild(FC_FAMILY, FC_STYLE, FC_FILE, (void *)(NULL));

      for (int i = 0; i < num_patterns && !s; i++)
      {
         pat = FcNameParse((FcChar8 *)font_patterns[i]);
         if (!pat)
         {
            continue;
         }

         fs = FcFontList(0, pat, os);
         FcPatternDestroy(pat);
         if (!fs)
         {
            continue;
         }

#ifdef GLVIS_DEBUG
         if (fs->nfont > 1)
         {
            cout <<
                 "Font pattern '" << font_patterns[i] << "' matched"
                 " multiple fonts:\n";
            for (int j = 0; j < fs->nfont; j++)
            {
               fnt = FcNameUnparse(fs->fonts[j]);
               cout << fnt << endl;
               free(fnt);
            }
            cout << "-----" << endl;
         }
#endif

         for (int j = 0; j < fs->nfont; j++)
         {
            res = FcPatternGetString(fs->fonts[j], FC_FILE, 0, &s);
            if (res == FcResultMatch && s)
            {
               font_file = (char *)s;
#ifdef GLVIS_DEBUG
               fnt = FcNameUnparse(fs->fonts[j]);
               cout << "Using font: " << fnt << endl;
               free(fnt);
#endif
               break;
            }
         }

         FcFontSetDestroy(fs);
      }

      if (os)
      {
         FcObjectSetDestroy(os);
      }

      FcFini();

      return s ? 0 : -2;
   }

public:
   GLVisFont() { init = 0; }

   int Initialized() const { return init; }

   void Init(const char *font_patterns[], int num_patterns, int height)
   {
      SetFont(font_patterns, num_patterns, height);
   }

   int SetFont(const char *font_patterns[], int num_patterns, int height)
   {
      string font_file;

      if (FindFontFile(font_patterns, num_patterns, font_file))
      {
         if (!init)
         {
            init = -4;
         }
         return -4;
      }

#ifdef GLVIS_DEBUG
      cout << "Using font file: " << font_file << endl;
#endif

      return SetFontFile(font_file.c_str(), height);
   }

   int SetFontFile(const char *font_file, int height);

   int Render(const char *text)
   {
      int            err;
      FT_Glyph       a_glyph;
      FT_BitmapGlyph bit_glyph;
      FT_Bitmap     *bitmap;

      if (init <= 0)
      {
         return 1;
      }

      delete [] image;
      image = NULL;

      LoadSequence(text);

      if (num_glyphs <= 0)
      {
         return 2;
      }

      GetBBox();

      const int pad = 2;

      image_width  = (bbox.xMax - bbox.xMin) + 2*pad;
      image_height = (bbox.yMax - bbox.yMin) + 2*pad;

      if (image_width <= 0 || image_height <= 0)
      {
         return 3;
      }

      GLfloat col[4];
      glGetFloatv(GL_CURRENT_COLOR, col);

      image = new unsigned char[4*image_width*image_height];
      for (int i = 0; i < 4*image_width*image_height; i++)
      {
         image[i] = 0;
      }

      for (FT_UInt g = 0; g < num_glyphs; g++)
      {
         a_glyph = glyph[g];

         err = FT_Glyph_To_Bitmap(&a_glyph, FT_RENDER_MODE_NORMAL, 0, 0);
         if (err) { continue; }

         bit_glyph = (FT_BitmapGlyph)a_glyph;
         bitmap = &bit_glyph->bitmap;

         int off_i = bit_glyph->left - bbox.xMin + pad;
         int off_j = image_height - bit_glyph->top + bbox.yMin - pad;

         for (int j = 0; j < (int) bitmap->rows; j++)
         {
            int im_j = image_height - 1 - (j + off_j);
            if (im_j < 0 || im_j >= image_height)
            {
#ifdef GLVIS_DEBUG
               cout <<
                    "GLVisFont::Render : outside 'y' range!\n"
                    "   text            = " << text            << "\n"
                    "   j               = " << j               << "\n"
                    "   bitmap->width   = " << bitmap->width   << "\n"
                    "   bitmap->rows    = " << bitmap->rows    << "\n"
                    "   off_i           = " << off_i           << "\n"
                    "   off_j           = " << off_j           << "\n"
                    "   bit_glyph->left = " << bit_glyph->left << "\n"
                    "   bit_glyph->top  = " << bit_glyph->top  << endl;
#endif
               continue;
            }
            for (int i = 0; i < (int) bitmap->width; i++)
            {
               int im_i = i + off_i;
               if (im_i < 0 || im_i >= image_width)
               {
#ifdef GLVIS_DEBUG
                  cout <<
                       "GLVisFont::Render : outside 'x' range!\n"
                       "   text            = " << text            << "\n"
                       "   i               = " << i               << "\n"
                       "   j               = " << j               << "\n"
                       "   bitmap->width   = " << bitmap->width   << "\n"
                       "   bitmap->rows    = " << bitmap->rows    << "\n"
                       "   off_i           = " << off_i           << "\n"
                       "   off_j           = " << off_j           << "\n"
                       "   bit_glyph->left = " << bit_glyph->left << "\n"
                       "   bit_glyph->top  = " << bit_glyph->top  << endl;
#endif
                  continue;
               }
               int off = 4*(im_i + im_j*image_width);
               int val = bitmap->buffer[i + j*bitmap->width];
               image[off + 0] = (unsigned char)(col[0]*255);
               image[off + 1] = (unsigned char)(col[1]*255);
               image[off + 2] = (unsigned char)(col[2]*255);
               image[off + 3] = val;
            }
         }

         FT_Done_Glyph(a_glyph);
      }

      FreeSeq();

      return 0;
   }

   const unsigned char *GetImage() const { return image; }
   int GetImageWidth() const { return image_width; }
   int GetImageHeight() const { return image_height; }

   ~GLVisFont()
   {
      if (init <= 0)
      {
         return;
      }

      FreeSeq();

      delete [] image;

      FT_Done_Face(face);
      FT_Done_FreeType(library);
   }
};
#endif /* GLVIS_USE_FREETYPE */
#endif /* FONT_HPP */
