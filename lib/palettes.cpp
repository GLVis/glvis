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

struct RGB {
   uint8_t r,g,b;

   RGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
   RGB(int r, int g, int b) :
      r(as_uint8(r)),
      g(as_uint8(g)),
      b(as_uint8(b)) {}
   RGB(double r, double g, double b) :
      r(as_uint8(r)),
      g(as_uint8(g)),
      b(as_uint8(b)) {}

   void print() {
      cout << +r << " " << +g << " " << +b << endl;
   }
   array<uint8_t,3> as_array() {
      return array<uint8_t,3>{r,g,b};
   }
};


struct Palette {
   string name;
   std::vector<RGB> colors;

   Palette(const string& name) : name(name) {}

   void addColor(double r, double g, double b) {
      colors.push_back(RGB(r,g,b));
   }
   void addColor(int r, int g, int b) {
      colors.push_back(RGB(r,g,b));
   }

   void print() {
      cout << name << " (size=" << size() << ")" << endl;;
      for (auto color : colors) {
         color.print();
      }
      cout << endl;
   }

   void reverse() {
      std::reverse(colors.begin(), colors.end());
   }

   int size() {
      return colors.size();
   }
};

class PaletteManager {
private:

public:
   std::vector<Palette> palettes;

   void addPalette(Palette palette) {
      // palette name is unique || container is empty
      if (get_index_by_name(palette.name) == -1 || palettes.empty()) {
         palettes.push_back(palette);
      }
   }

   void addPalette(string name) {
      if (get_index_by_name(name) == -1 || palettes.empty()) {
         palettes.push_back(Palette(name));
      }
   }

   int get_index_by_name(string name) {
      for (int i = 0; i < size(); i++) {
         if (palettes[i].name == name) {
            return i;
         }
      }
      return -1;
   }

   int size() {
      return palettes.size();
   }

   void print() {
      for (Palette cmap : palettes) {
         cmap.print();
         cout << endl;
      }
   }

   void load(istream &pal) {
      string word, palname, channeltype;
      int idx = -1;

      // read initializing commands
      while (1) {
         pal >> ws;
         if (!pal.good()) {
            // cout << "Error in palette" << endl;
            break;
         }
         if (pal.peek() == '#') {
            getline(pal, word);
            continue;
         }
         pal >> word;
         if (word == "palette") {
            pal >> palname >> channeltype;
            addPalette(palname);
            idx = get_index_by_name(palname);
            cout << "Reading palette: (" << idx << ") " << palname << endl;
         }
         else if (channeltype == "float" && idx != -1) {
            float r, g, b;
            r = stof(word);
            pal >> g >> b;
            palettes[idx].addColor(r,g,b);
         }
         else if (channeltype == "int" && idx != -1) {
            int r, g, b;
            r = stoi(word);
            pal >> g >> b;
            palettes[idx].addColor(r,g,b);
         }
      }
      cout << "Finished loading palettes from file" << endl;
   }
};


// const int Num_RGB_Palettes = 71;
// const int RGB_Palettes_Sizes[Num_RGB_Palettes] =
// double *RGB_Palettes[Num_RGB_Palettes] =
// const char *RGB_Palettes_Names[Num_RGB_Palettes] =


int PaletteState::ChoosePalette()
{
   const int buflen = 256;
   char buffer[buflen];
   int pal;
   cout << "Choose a palette:\n";
   for (pal = 0; pal < palettes->size(); pal++)
   {
      cout << setw(4) << pal+1 << ") " << RGB_Palettes_Names[pal];
      if ((pal+1)%5 == 0)
      {
         cout << '\n';
      }
   }
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
   else if (pal > palettes->size())
   {
      pal = palettes->size();
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

PaletteState::PaletteState(PaletteManager* palettes)
   : palettes(palettes)
   , palette_tex(palettes->size())
   , first_init(false)
{
}

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
         // Init_Palettes();
      }

      GLuint paletteTexIds[palettes->size()][2];
      GLuint alphaTexId;

      glGenTextures(palettes->size() * 2, &(paletteTexIds[0][0]));
      glGenTextures(1, &alphaTexId);

      for (int ipal = 0; ipal < palettes->size(); ipal++)
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

   for (int i = 0; i < palettes->size(); i++)
   {
      ToTextureDiscrete(RGB_Palettes[i], RGB_Palettes_Sizes[i],
                        palette_tex[i][0]);
      ToTextureSmooth(RGB_Palettes[i], RGB_Palettes_Sizes[i],
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
   int palSize = RGB_Palettes_Sizes[curr_palette];
   const double* palData = RGB_Palettes[curr_palette];
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

const double * PaletteState::GetData() const
{
   return RGB_Palettes[curr_palette];
}

int PaletteState::GetSize(int pal) const
{
   if (pal == -1)
   {
      return RGB_Palettes_Sizes[curr_palette];
   }
   return RGB_Palettes_Sizes[pal];
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
   SetIndex((curr_palette + 1) % palettes->size());
}

void PaletteState::PrevIndex()
{
   SetIndex((curr_palette == 0) ? palettes->size() - 1 :
            curr_palette - 1);
}

int PaletteState::SelectNewRGBPalette()
{
   int pal = ChoosePalette();

   SetIndex(pal);

   return pal;
}
