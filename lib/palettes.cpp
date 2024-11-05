// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "palettes.hpp"
#include "gl/renderer.hpp"

#include <cmath>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mutex>

using namespace std;


int PaletteState::ChoosePalette()
{
   const int buflen = 256;
   char buffer[buflen];
   int pal;
   cout << "Choose a palette:\n";
   Palettes->printSummary();
   cout << "\n ---> [" << curr_palette+1 << "] " << flush;

   cin.getline (buffer, buflen);
   cin.getline (buffer, buflen);

   if (buffer[0])
   {
      sscanf(buffer, "%i", &pal);
   }
   else
   {
      pal = curr_palette+1;
   }

   if (pal < 1)
   {
      pal = 1;
   }
   else if (pal > (int)Palettes->NumPalettes())
   {
      pal = (int)Palettes->NumPalettes();
   }

   return pal-1;
}


// Generates a discrete texture from the given palette.
void PaletteState::ToTextureDiscrete(double * palette, size_t plt_size,
                                     GLuint tex)
{
   vector<array<float,4>> texture_buf(plt_size);

   if (RepeatPaletteTimes > 0)
   {
      for (size_t i = 0; i < plt_size; i++)
      {
         texture_buf[i] =
         {
            (float) palette[3*i],
            (float) palette[3*i+1],
            (float) palette[3*i+2],
            1.0
         };
      }
   }
   else
   {
      for (size_t i = 0; i < plt_size; i++)
      {
         texture_buf[i] =
         {
            (float) palette[3*(plt_size-1-i)+0],
            (float) palette[3*(plt_size-1-i)+1],
            (float) palette[3*(plt_size-1-i)+2],
            1.0
         };
      }
   }
   if (PaletteNumColors > 1 && (plt_size > (size_t)PaletteNumColors))
   {
      texture_buf.resize(PaletteNumColors);
      for (int i = 0; i < PaletteNumColors; i++)
      {
         int plt_i = i * plt_size / (PaletteNumColors-1);
         if (i >= PaletteNumColors - 1)
         {
            plt_i = plt_size - 1;
         }
         if (RepeatPaletteTimes < 0)
         {
            plt_i = plt_size-1-plt_i;
         }
         texture_buf[i] =
         {
            (float) palette[3*plt_i],
            (float) palette[3*plt_i+1],
            (float) palette[3*plt_i+2],
            1.0
         };
      }
      plt_size = PaletteNumColors;
   }
   glBindTexture(GL_TEXTURE_2D, tex);
   glTexImage2D(GL_TEXTURE_2D,
                0,
                rgba_internal,
                plt_size,
                1,
                0,
                GL_RGBA,
                GL_FLOAT,
                texture_buf.data());

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

// Generates a smooth texture from the given palette.
void PaletteState::ToTextureSmooth(double * palette, size_t plt_size,
                                   GLuint tex)
{
   vector<array<float,4>> texture_buf(MaxTextureSize);
   glBindTexture(GL_TEXTURE_2D, tex);

   size_t textureSize = MaxTextureSize;
   if (plt_size * abs(RepeatPaletteTimes) <= textureSize)
   {
      int flip_start = RepeatPaletteTimes < 0;
      for (int rpt = 0; rpt < abs(RepeatPaletteTimes); rpt++)
      {
         for (size_t i = 0; i < plt_size; i++)
         {
            // flip = 0: p_i = i
            // flip = 1: p_i = plt_size-1-i
            int p_i = (flip_start + rpt) % 2 == 0 ? i : plt_size - 1 - i;
            texture_buf[i + plt_size * rpt] =
            {
               (float) palette[3*p_i],
               (float) palette[3*p_i + 1],
               (float) palette[3*p_i + 2],
               1.0
            };
         }
      }
      glTexImage2D(GL_TEXTURE_2D, 0, rgba_internal,
                   plt_size * abs(RepeatPaletteTimes), 1,
                   0, GL_RGBA, GL_FLOAT, texture_buf.data());
   }
   else
   {
      for (size_t i = 0; i < textureSize; i++)
      {
         double t = double(i) / textureSize - 1;
         t *= 0.999999999 * (plt_size - 1) * abs(RepeatPaletteTimes);
         int j = floor(t);
         t -= j;
         int p_i;
         if (((j / (plt_size-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
             ((j / (plt_size-1)) % 2 == 1 && RepeatPaletteTimes < 0))
         {
            p_i = j % (plt_size - 1);
         }
         else
         {
            p_i = plt_size - 2 - j % (plt_size - 1);
            t = 1.0 - t;
         }
         texture_buf[i] =
         {
            (float)((1.0-t) * palette[3*p_i] + t * palette[3*(p_i+1)]),
            (float)((1.0-t) * palette[3*p_i+1] + t * palette[3*(p_i+1)+1]),
            (float)((1.0-t) * palette[3*p_i+2] + t * palette[3*(p_i+1)+2]),
            1.0
         };
      }
      glTexImage2D(GL_TEXTURE_2D, 0, rgba_internal,
                   textureSize, 1,
                   0, GL_RGBA, GL_FLOAT, texture_buf.data());
   }

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

PaletteState::PaletteState()
   : first_init(false)
   ,Palettes(&BasePalettes)
   , palette_tex(BasePalettes.NumPalettes()
                )
{}

static std::mutex init_mtx;

void PaletteState::Init()
{
   if (!first_init)
   {

      glGetIntegerv(GL_MAX_TEXTURE_SIZE, &MaxTextureSize);
      if (MaxTextureSize < 4096)
      {
         cerr << "Warning: GL_MAX_TEXTURE_SIZE is less than 4096." << endl;
      }
      MaxTextureSize = std::min(MaxTextureSize, 4096);
      {
         std::lock_guard<std::mutex> lk{init_mtx};
      }

      vector<array<GLuint, 2>> paletteTexIds(Palettes->MAX_PALETTES);
      GLuint alphaTexId;

      glGenTextures(Palettes->NumPalettes() * 2, &(paletteTexIds[0][0]));
      glGenTextures(1, &alphaTexId);

      for (size_t ipal = 0; ipal < Palettes->NumPalettes(); ipal++)
      {
         palette_tex[ipal][0] = paletteTexIds[ipal][0];
         palette_tex[ipal][1] = paletteTexIds[ipal][1];
      }
      alpha_tex = alphaTexId;

      GLenum alpha_internal;
      if (gl3::GLDevice::useLegacyTextureFmts())
      {
         alpha_internal = GL_ALPHA;
         alpha_channel = GL_ALPHA;
         rgba_internal = GL_RGBA;
      }
      else
      {
         // WebGL 2 requires sized internal format for float texture
         alpha_internal = GL_R32F;
         alpha_channel = GL_RED;
         rgba_internal = GL_RGBA32F;
      }
      // set alpha texture to 1.0
      std::vector<float> alphaTexData(MaxTextureSize);
      std::fill(alphaTexData.begin(), alphaTexData.end(), 1.0f);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, alpha_tex);
      glTexImage2D(GL_TEXTURE_2D, 0, alpha_internal, MaxTextureSize, 1, 0,
                   alpha_channel, GL_FLOAT, alphaTexData.data());
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      glActiveTexture(GL_TEXTURE0);
      first_init = true;
   }

   for (int i = 0; i < (int)Palettes->NumPalettes(); i++)
   {
      ToTextureDiscrete(GetData(i),
                        GetSize(i),
                        palette_tex[i][0]);
      ToTextureSmooth(GetData(i),
                      GetSize(i),
                      palette_tex[i][1]);
   }
}



double PaletteState::GetColorCoord(double val, double min, double max)
{
   // static double eps = 1e-24;
   static const double eps = 0.0;
   if (use_logscale)
   {
      if (val < min)
      {
         val = min;
      }
      if (val > max)
      {
         val = max;
      }
      return log(fabs(val/(min+eps))) / (log(fabs(max/(min+eps)))+eps);
   }
   else
   {
      return ((val-min)/(max-min));
   }
}

void PaletteState::GetColorFromVal(double val, float * rgba)
{
   int palSize = GetSize();
   const double* palData = GetData();
   val *= 0.999999999 * ( palSize - 1 ) * abs(RepeatPaletteTimes);
   int i = (int) floor( val );
   double t = val - i;

   const double* pal;
   if (((i / (palSize-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
       ((i / (palSize-1)) % 2 == 1 && RepeatPaletteTimes < 0))
   {
      pal = palData + 3 * ( i % (palSize-1) );
   }
   else
   {
      pal = palData + 3 * ( (palSize-2) - i % (palSize-1) );
      t = 1.0 - t;
   }
   rgba[0] = (1.0 - t) * pal[0] + t * pal[3];
   rgba[1] = (1.0 - t) * pal[1] + t * pal[4];
   rgba[2] = (1.0 - t) * pal[2] + t * pal[5];
   rgba[3] = 1.f;
}

double * PaletteState::GetData(int pidx)
{
   // return RGB_Palettes[curr_palette];
   if (pidx == -1)
   {
      return Palettes->get(curr_palette)->as_rgb_array();
   }
   return Palettes->get(pidx)->as_rgb_array();
}

void PaletteState::GenerateAlphaTexture(float matAlpha, float matAlphaCenter)
{
   std::vector<float> alphaTexData(MaxTextureSize);
   if (matAlpha >= 1.0)
   {
      // transparency off
      std::fill(alphaTexData.begin(), alphaTexData.end(), 1.0f);
   }
   else
   {
      for (int i = 0; i < MaxTextureSize; i++)
      {
         double val = double(2*i + 1)/(2*MaxTextureSize); // midpoint of texel
         if (matAlphaCenter > 1.0)
         {
            alphaTexData[i] = matAlpha * std::exp(-(matAlphaCenter)*std::abs(val - 1.0));
         }
         else if (matAlphaCenter < 0.0)
         {
            alphaTexData[i] = matAlpha * std::exp((matAlphaCenter - 1.0)*std::abs(val));
         }
         else
         {
            alphaTexData[i] = matAlpha * std::exp(-std::abs(val - matAlphaCenter));
         }
      }
   }
   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, alpha_tex);
   glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, MaxTextureSize, 1, alpha_channel,
                   GL_FLOAT, alphaTexData.data());
   glActiveTexture(GL_TEXTURE0);
}

void PaletteState::NextIndex()
{
   SetIndex((curr_palette + 1) % Palettes->NumPalettes());
}

void PaletteState::PrevIndex()
{
   SetIndex((curr_palette == 0) ? Palettes->NumPalettes() - 1 :
            curr_palette - 1);
}

int PaletteState::SelectNewRGBPalette()
{
   int pal = ChoosePalette();

   SetIndex(pal);

   return pal;
}

int PaletteState::NumPalettes()
{
   return Palettes->NumPalettes();
}

int PaletteState::GetSize(int pidx) const
{
   if (pidx == -1)
   {
      return Palettes->get(curr_palette)->size();
   }
   return Palettes->get(pidx)->size();
}
