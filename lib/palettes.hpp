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
#include "gl/types.hpp"
#include <vector>
#include <array>

// Cast a double in range [0,1] to a uint8_t
uint8_t as_uint8(double x) {
   if (x >= 0 && x <= 1.0)
   {
      return static_cast<uint8_t>(x * 255.0f);
   }
   else
   {
      throw std::out_of_range("Value out of range [0, 1]");
   }
}
uint8_t as_uint8(int x) {
   if (x >= 0 && x <= 255)
   {
      return static_cast<uint8_t>(x);
   }
   else
   {
      throw std::out_of_range("Value out of range [0, 255]");
   }
}
struct RGB
{
   uint8_t r,g,b;

   RGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
   RGB(int r, int g, int b) :
      r(as_uint8(r)),
      g(as_uint8(g)),
      b(as_uint8(b)) {};
   RGB(double r, double g, double b) :
      r(as_uint8(r)),
      g(as_uint8(g)),
      b(as_uint8(b)) {};

   void print();
   array<uint8_t,3> as_array();
};

struct Palette
{
   string name;
   std::vector<RGB> colors;
   Palette(const string& name) : name(name) {};

   void addColor(double r, double g, double b);
   void addColor(int r, int g, int b);
   void print();
   void reverse();
   int size();
   double* as_double_array();
};

class PaletteManager
{
public:
   std::vector<Palette> palettes;
   void addPalette(Palette palette);
   void addPalette(string name);
   int get_index_by_name(string name);
   Palette* get_name(string name);
   Palette* get(int idx);
   int size();
   void print();
   void load(istream &pal);
private:
};

class PaletteState
{
public:
   PaletteState(PaletteManager* palettes);
   /// Palettes
   PaletteManager* palettes;
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
