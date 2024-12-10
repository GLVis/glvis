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

#ifndef GLVIS_PALETTES_HPP
#define GLVIS_PALETTES_HPP
#include "palettes_base.hpp"
#include <vector>
#include <array>

class PaletteState
{
public:
   PaletteState();
   /// Initializes the palette textures (and defines texture ids).
   void InitTextures();
   /// Binds the discrete version of the current palette texture.
   void UseDiscrete() { use_smooth = 0; }
   /// Binds the smooth version of the current palette texture.
   void UseSmooth() { use_smooth = 1; }
   /// Gets whether the smooth texture is being used (1 = true)
   int GetSmoothSetting() { return use_smooth; }
   /// Sets the palette texture to bind.
   void SetIndex(int num);
   int GetCurrIndex() const { return curr_palette; }
   void NextIndex();
   void PrevIndex();
   int ChoosePalette();
   int SelectNewRGBPalette();
   /// Gets a pointer to a palette (default = current palette).
   Palette* GetPalette(int pidx = -1) const;
   /// Gets the total number of colors in the current palette color array.
   int GetSize(int pidx = -1) const;
   /// Gets the number of colors used in the current palette color array.
   int GetNumColors(int pal = -1) const
   { return PaletteNumColors ? PaletteNumColors : GetSize(pal); }
   /// Sets the number of colors to use in the current palette color array.
   void SetNumColors(int numColors) { PaletteNumColors = numColors; }
   int GetRepeatTimes() const { return RepeatPaletteTimes; }
   void SetRepeatTimes(int rpt) { RepeatPaletteTimes = rpt; }

   void SetUseLogscale(bool logscale) { use_logscale = logscale; }
   bool GetUseLogscale() { return use_logscale; }
   double GetColorCoord(double val, double min, double max);
   void GetColorFromVal(double val, float* rgba);

   GLuint GetColorTexture() const
   { return textures[curr_palette][use_smooth].Get(); }
   GLuint GetAlphaTexture() const { return alpha_tex.Get(); }
   /// Generates new textures with the same ids, using current settings
   void GenerateTextures(bool reinitialize = false);
   void GenerateAlphaTexture(float matAlpha, float matAlphaCenter);

   int NumPalettes();

private:
   bool initialized;
   PaletteRegistry* Palettes;

   // Regular (rgba) textures
   std::vector<std::array<Texture,2>> textures;
   // Global alpha texture - blended with other textures
   Texture alpha_tex;

   int curr_palette = 2;
   int use_smooth = 0;
   int RepeatPaletteTimes = 1;
   int PaletteNumColors = 0;

   bool use_logscale = false;

};

#endif
