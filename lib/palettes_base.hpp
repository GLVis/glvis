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

#ifndef GLVIS_PALETTESBASE_HPP
#define GLVIS_PALETTESBASE_HPP

#include "gl/types.hpp"

#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <iomanip>

using namespace std;

struct RGBAf
{
   float r, g, b, a;

   constexpr RGBAf(float r = 0.0, float g = 0.0, float b = 0.0, float a = 1.0)
      : r(r), g(g), b(b), a(a) {}

   void print(bool printalpha = false, ostream& os = cout) const
   {
      os << fixed << setprecision(6)
         << setw(10) << r << " "
         << setw(10) << g << " "
         << setw(10) << b;
      if (printalpha)
      {
         os << " " << setw(10) << a;
      }
   }

   array<float, 4> as_array() const
   {
      return {r, g, b, a};
   }

};


class Palette
{
public:
   const string name;

   /// Constructors
   Palette(const string& name) : name(name) {}

   /// Constructor from vector of RGBAf
   Palette(const string& name, const vector<RGBAf>& colors)
      : name(name), colors(colors) {}

   /// Constructor from Nx3 array
   template <size_t N>
   Palette(const string& name, const array<array<float,3>,N>& arr) : name(name)
   {
      colors.reserve(N);
      for (size_t i = 0; i < N; ++i)
      {
         colors[i] = RGBAf(arr[i][0], arr[i][1], arr[i][2]);
      }
   }

   /// Constructor from Nx4 array
   template <size_t N>
   Palette(const string& name, const array<array<float,4>,N>& arr) : name(name)
   {
      colors.reserve(N);
      for (size_t i = 0; i < N; ++i)
      {
         colors[i] = RGBAf(arr[i][0], arr[i][1], arr[i][2], arr[i][3]);
      }
   }

   int size() const { return colors.size(); }

   /// Add color to palette
   void addColor(float r, float g, float b, float a = 1.0)
   {
      colors.push_back(RGBAf(r, g, b, a));
   }

   /// Print this palette
   void print(ostream& os = cout) const
   {
      os << "palette " << name << " RGBf" << endl;
      for (const auto& color : colors)
      {
         color.print(false, os);
         os << endl;
      }
      os << endl;
   }

   /// Get color at index i (optionally, use reversed order)
   RGBAf color(int i, bool reversed = false) const
   {
      int j = reversed ? size() - 1 - i : i;
      return colors[j];
   }

   vector<array<float,4>> data(bool reversed = false) const
   {
      vector<array<float,4>> rgba_data(size());
      for (int i = 0; i < size(); ++i)
      {
         rgba_data[i] = color(i, reversed).as_array();
      }
      return rgba_data;
   }

private:
   vector<RGBAf> colors;
};

struct Texture
{
   /// The palette to create a texture of
   Palette* const palette;
   /// Repeat the palette multiple times (negative for reverse); cannot be 0
   int Nrepeat;
   /// Reverse the palette
   bool reversed;
   /// Number of colors to discretize with (0 uses the original number of colors)
   int Ncolors;
   /// Is texture smooth or discrete?
   bool smooth;
   /// Texture size
   int size;
   /// Max texture size
   int MAX_TEXTURE_SIZE;
   /// Texture data
   vector<array<float,4>> texture;

   Texture(Palette* palette, int Nrepeat_ = 1, int Ncolors_ = 0,
           bool smooth = false)
      : palette(palette)
   {
      // Get the maximum texture size
      glGetIntegerv(GL_MAX_TEXTURE_SIZE, &MAX_TEXTURE_SIZE);
      if (MAX_TEXTURE_SIZE < 4096)
      {
         cerr << "Warning: GL_MAX_TEXTURE_SIZE is less than 4096." << endl;
      }
      // Is limiting to 4096 necessary?
      MAX_TEXTURE_SIZE = min(MAX_TEXTURE_SIZE, 4096);
      // Nrepeat cannot be 0; we also extract the sign
      reversed = Nrepeat_ < 0;
      Nrepeat = Nrepeat_ == 0 ? 1 : abs(Nrepeat_);
      // Ncolors must be positive
      Ncolors = Ncolors_ <= 0 ? palette->size() : Ncolors_;

      generate();
   }

   /// Generate the texture
   void generate()
   {
      // original palette size
      int plt_size = palette->size();
      // Set the texture size
      size = Nrepeat * Ncolors;
      if (size > MAX_TEXTURE_SIZE)
      {
         cerr << "Warning: Texture size "
              << "(" << size << ")" << " exceeds maximum "
              << "(" << MAX_TEXTURE_SIZE << ")" << endl;
         if (Ncolors >= MAX_TEXTURE_SIZE)
         {
            Ncolors = MAX_TEXTURE_SIZE;
            Nrepeat = 1;
            size = Nrepeat * Ncolors;
         }
         else
         {
            Nrepeat = MAX_TEXTURE_SIZE / Ncolors;
            size = Nrepeat * Ncolors;
         }
      }
      texture.clear();
      texture.resize(size);

      // generate the discrete texture
      // indices: plt_size x Nrepeat -> size
      if (!smooth)
      {
         for (int rpt = 0; rpt < Nrepeat; rpt++)
         {
            bool reverse = (reversed + rpt) % 2 != 0;
            for (int i = 0; i < Ncolors; i++)
            {
               int j = 0.999999 * i * plt_size / (Ncolors - 1);
               texture[rpt*Ncolors + i] = palette->color(j, reverse).as_array();
            }
         }
      }
      // smooth texture interpolate colors
      else
      {
         for (int rpt = 0; rpt < Nrepeat; rpt++)
         {
            bool reverse = (reversed + rpt) % 2 != 0;
            for (int i = 0; i < Ncolors; i++)
            {
               float t = 0.999999 * i * (plt_size - 1) / (Ncolors - 1);
               int j = floor(t);
               t -= j;
               array<float,4> col1 = palette->color(j, reverse).as_array();
               array<float,4> col2 = palette->color(j+1, reverse).as_array();
               texture[rpt*Ncolors + i] =
               {
                  (1-t) * col1[0] + t * col2[0],
                  (1-t) * col1[1] + t * col2[1],
                  (1-t) * col1[2] + t * col2[2],
                  (1-t) * col1[3] + t * col2[3]
               };
            }
         }

      }
   }
};


// Behaves like make_unique (only available in >= c++14)
template<typename T, typename... Args>
std::unique_ptr<T> as_unique(Args&&... args)
{
   return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// PaletteRegistry with a vector of unique_ptr<Palette>. Besides holding
// the palettes, this should be stateless.
class PaletteRegistry
{
private:
   vector<unique_ptr<Palette>> palettes;

   int get_index_by_name(const string& name) const
   {
      for (int i = 0; i < NumPalettes(); i++)
      {
         if (get(i)->name == name)
         {
            return i;
         }
      }
      return -1;
   }

public:
   // empty constructor
   PaletteRegistry() {}

   PaletteRegistry(const vector<Palette>& paletteRefs)
   {
      for (const Palette& palette : paletteRefs)
      {
         if (check_name(palette.name))
         {
            palettes.push_back(as_unique<Palette>(palette));
         }
      }
   }

   void addPalette(Palette& palette)
   {
      if (check_name(palette.name))
      {
         palettes.push_back(as_unique<Palette>(palette));
      }
   }

   void addPalette(const string& name)
   {
      if (check_name(name))
      {
         palettes.push_back(as_unique<Palette>(name));
      }
   }

   bool check_name(const string& name) const
   {
      // palette name is unique || container is empty
      if (get_index_by_name(name) == -1 || palettes.empty())
      {
         return true;
      }
      else
      {
         cout << "Palette with name: '" << name << "' already exists in registry.";
         return false;
      }
   }

   // get by index
   Palette* get(int index) const
   {
      if (0 <= index && index <= NumPalettes()-1)
      {
         return palettes[index].get();
      }
      cout << "Palette (index = " << index+1 << ") out of range. Available palettes:"
           << endl;
      this->printSummary();
      return palettes.back().get();
   }

   // get by name
   Palette* get(const string& name) const
   {
      int idx = get_index_by_name(name);
      if (idx != -1)
      {
         return palettes[idx].get();
      }
      cout << "Palette (name = " << name << ") not found. Available palettes:" <<
           endl;
      this->printSummary();
      return palettes.back().get();
   }

   void printSummary(ostream& os = cout) const
   {
      for (int i = 0; i < NumPalettes(); i++)
      {
         os << setw(3) << i+1 << ") "
            << left << setw(12) << get(i)->name << right;
         if ((i+1)%5 == 0)
         {
            os << endl;
         }
      }
      os << endl;
   }

   void printAll(ostream& os = cout) const
   {
      for (int i = 0; i < NumPalettes(); i++)
      {
         get(i)->print(os);
      }
   }

   int NumPalettes() const
   {
      return palettes.size();
   }

   void load(const string& palette_filename)
   {

      ifstream pfile(palette_filename);
      if (!pfile)
      {
         cout << "Could not open palette file: " << palette_filename << endl;
         return;
      }
      string word, palname, channeltype;
      int idx = -1;

      // read initializing commands
      while (1)
      {
         pfile >> ws;
         if (!pfile.good())
         {
            break;
         }
         if (pfile.peek() == '#')
         {
            getline(pfile, word);
            continue;
         }
         pfile >> word;
         if (word == "palette")
         {
            pfile >> palname >> channeltype;
            idx = get_index_by_name(palname);
            if (idx == -1)
            {
               addPalette(palname);
               idx = get_index_by_name(palname);
               cout << "Reading palette: (" << idx+1 << ") " << palname << endl;
            }
            else
            {
               cout << "Error reading palette: " << palname
                    << ". Palette with same name already exists." << endl;
               break;
            }
         }
         else if (channeltype == "RGBf" && idx != -1)
         {
            float r, g, b;
            r = stof(word);
            pfile >> g >> b;
            get(idx)->addColor(r,g,b);
         }
         else if (channeltype == "RGBAf" && idx != -1)
         {
            float r, g, b, a;
            r = stof(word);
            pfile >> g >> b >> a;
            get(idx)->addColor(r,g,b,a);
         }
         else
         {
            cout << "Error reading palette file: " << palette_filename << endl;
            break;
         }
      }
      cout << "Finished loading palettes from file: " << palette_filename << endl;
   }

};


extern PaletteRegistry BasePalettes;

#endif
