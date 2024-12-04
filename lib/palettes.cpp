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

#include <cmath>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


int PaletteState::ChoosePalette()
{
   const int buflen = 256;
   char buffer[buflen];
   int pal;
   cout << "Choose a palette:\n";
   Palettes->PrintSummary();
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
   else if (pal > Palettes->NumPalettes())
   {
      pal = Palettes->NumPalettes();
   }

   return pal-1;
}

PaletteState::PaletteState()
   : Palettes(&BasePalettes)
   , textures(Palettes->NumPalettes())
{
   // Init the palette textures (will creates new texture ids)
   InitTextures();
   // Generate the textures
   GenerateTextures();
}

void PaletteState::InitTextures()
{
   int N = Palettes->NumPalettes();

   // Initialize both discrete [0] and smooth [1] textures
   // Texture constructor will assign texture ids
   for (int i = 0; i < N; i++)
   {
      const Palette* pal = static_cast<const Palette*>(Palettes->Get(i));
      textures[i][0] = Texture(pal, Texture::TextureType::DISCRETE,
                               RepeatPaletteTimes, PaletteNumColors);
      textures[i][1] = Texture(pal, Texture::TextureType::SMOOTH,
                               RepeatPaletteTimes, PaletteNumColors);
   }

   // Init the alpha texture (TODO: simplify this using Texture class)
   // Generate the texture id
   GLuint alphaTexId;
   glGenTextures(1, &alphaTexId);
   alpha_tex = alphaTexId;

   // set alpha texture to 1.0
   std::vector<float> alphaTexData(Texture::max_texture_size);
   std::fill(alphaTexData.begin(), alphaTexData.end(), 1.0f);
   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, alpha_tex);
   glTexImage2D(GL_TEXTURE_2D, 0, Texture::alpha_internal,
                Texture::max_texture_size, 1, 0,
                Texture::alpha_channel, GL_FLOAT, alphaTexData.data());
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

   glActiveTexture(GL_TEXTURE0);

}

void PaletteState::GenerateTextures()
{
   for (int i = 0; i < Palettes->NumPalettes(); i++)
   {
      textures[i][0].GenerateGLTexture(RepeatPaletteTimes, PaletteNumColors);
      textures[i][1].GenerateGLTexture(RepeatPaletteTimes, PaletteNumColors);
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
   Palette* pal = GetPalette();
   val *= 0.999999999 * ( palSize - 1 ) * abs(RepeatPaletteTimes);
   int i = (int) floor( val );
   float t = float(val) - i;
   int idx = 0;

   // const float* pal;
   if (((i / (palSize-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
       ((i / (palSize-1)) % 2 == 1 && RepeatPaletteTimes < 0))
   {
      idx = i % (palSize-1);
   }
   else
   {
      idx = (palSize-2) - i % (palSize-1);
      t = 1.0 - t;
   }
   RGBAf color1 = pal->Color(idx);
   RGBAf color2 = pal->Color(idx+1);
   rgba[0] = (1.0 - t) * color1.r + t * color2.r;
   rgba[1] = (1.0 - t) * color1.g + t * color2.g;
   rgba[2] = (1.0 - t) * color1.b + t * color2.b;
   rgba[3] = (1.0 - t) * color1.a + t * color2.a;
}

Palette* PaletteState::GetPalette(int pidx) const
{
   if (pidx == -1)
   {
      return Palettes->Get(curr_palette);
   }
   return Palettes->Get(pidx);
}

void PaletteState::GenerateAlphaTexture(float matAlpha, float matAlphaCenter)
{
   std::vector<float> alphaTexData(Texture::max_texture_size);
   if (matAlpha >= 1.0)
   {
      // transparency off
      std::fill(alphaTexData.begin(), alphaTexData.end(), 1.0f);
   }
   else
   {
      for (int i = 0; i < Texture::max_texture_size; i++)
      {
         double val = double(2*i + 1)/(2*Texture::max_texture_size); // midpoint of texel
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
   glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Texture::max_texture_size, 1,
                   Texture::alpha_channel,
                   GL_FLOAT, alphaTexData.data());
   glActiveTexture(GL_TEXTURE0);
}

void PaletteState::SetIndex(int num)
{
   curr_palette = num;
   cout << "Palette: " << num << ") " << Palettes->Get(curr_palette)->name << endl;
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
      return Palettes->Get(curr_palette)->Size();
   }
   return Palettes->Get(pidx)->Size();
}
