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

#ifndef GLVIS_PALETTES_HPP
#define GLVIS_PALETTES_HPP
#include "gl/types.hpp"
#include <vector>
#include <array>

class PaletteState
{
public:
   PaletteState();
   /// Initializes the palette textures.
   void Init();
   /// Binds the discrete version of the current palette texture.
   void UseDiscrete() { use_smooth = 0; }
   /// Binds the smooth version of the current palette texture.
   void UseSmooth() { use_smooth = 1; }
   /// Gets whether the smooth texture is being used (1 = true)
   int GetSmoothSetting() { return use_smooth; }
   /// Sets the palette texture to bind.
   void SetIndex(int num) { curr_palette = num; }
   int GetCurrIndex() const { return curr_palette; }
   void NextIndex();
   void PrevIndex();
   int ChoosePalette();
   int SelectNewRGBPalette();
   /// Gets the data in the palette color array.
   const double* GetData() const;
   /// Gets the total number of colors in the current palette color array.
   int GetSize(int pal = -1) const;
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
   { return palette_tex[curr_palette][use_smooth]; }
   GLuint GetAlphaTexture() const { return alpha_tex; }
   void GenerateAlphaTexture(float matAlpha, float matAlphaCenter);
private:
   void ToTextureDiscrete(double * palette, size_t plt_size, GLuint tex);
   void ToTextureSmooth(double * palette, size_t plt_size, GLuint tex);
   using TexHandle = gl3::resource::TextureHandle;

   std::vector<std::array<TexHandle,2>> palette_tex;
   TexHandle alpha_tex;

   int curr_palette = 2;
   int use_smooth = 0;
   int RepeatPaletteTimes = 1;
   int PaletteNumColors = 0;

   bool use_logscale = false;

   bool first_init;
   int MaxTextureSize;
   GLenum alpha_channel;
   GLenum rgba_internal;
};

#endif
