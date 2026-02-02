// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
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
   : initialized(false)
   , Palettes(&BasePalettes)
   , textures(Palettes->NumPalettes())
{
}

void PaletteState::InitTextures()
{
   cout << "Initializing textures." << endl;
   int N = Palettes->NumPalettes();
   textures.resize(N);

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

   // Initialize the global alpha texture
   alpha_tex = Texture(1.0, 0.5);

   // Set flag
   initialized = true;
}

void PaletteState::GenerateTextures(bool reinitialize)
{
   if (!initialized || reinitialize)
   {
      InitTextures();
   }

   for (int i = 0; i < Palettes->NumPalettes(); i++)
   {
      textures[i][0].UpdateParameters(RepeatPaletteTimes, PaletteNumColors);
      textures[i][0].Generate();
      textures[i][1].UpdateParameters(RepeatPaletteTimes, PaletteNumColors);
      textures[i][1].Generate();
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
   alpha_tex.UpdateAlphaParameters(matAlpha, matAlphaCenter);
   alpha_tex.Generate();
}

void PaletteState::SetIndex(int num)
{
   if ((num >= 0) && (num < Palettes->NumPalettes()))
   {
      curr_palette = num;
      cout << "Palette: " << num+1 << ") " << Palettes->Get(curr_palette)->name
           << endl;
   }
   else
   {
      cout << "Palette index " << num+1 << " is out of bounds." << endl;
   }
}

void PaletteState::SetByName(const std::string& palette_name)
{
   int num = Palettes->GetIndexByName(palette_name);
   if ((num >= 0) && (num < Palettes->NumPalettes()))
   {
      curr_palette = num;
      cout << "Palette: " << num+1 << ") " << Palettes->Get(curr_palette)->name
           << endl;
   }
   else
   {
      cout << "Palette " << palette_name << " is not defined." << endl;
   }
}

bool PaletteState::UseDefaultIndex()
{
   int num = Palettes->GetDefault();
   if ((num >= 0) && (num < Palettes->NumPalettes()))
   {
      curr_palette = num;
      cout << "Palette: " << num+1 << ") " << Palettes->Get(curr_palette)->name
           << endl;
      return true;
   }
   return false;
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
