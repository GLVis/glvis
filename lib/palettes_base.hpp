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

   void print(ostream& os = cout) const;

   array<float, 4> as_array() const { return {r, g, b, a}; }

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
   Palette(const string& name, const array<array<float,3>,N>& arr);

   /// Constructor from Nx4 array
   template <size_t N>
   Palette(const string& name, const array<array<float,4>,N>& arr);

   /// Get size
   int size() const { return colors.size(); }

   /// Add color to palette
   void addColor(float r, float g, float b, float a = 1.0);

   /// Print each color of this palette to a stream
   void print(ostream& os = cout) const;

   /// Get color at index i (optionally, use reversed order)
   RGBAf color(int i, bool reversed = false) const;

   /// Get all colors as a vector of float arrays
   vector<array<float,4>> data(bool reversed = false) const;

private:
   vector<RGBAf> colors;
};

/// Generates the texture data for a given palette, to be used in OpenGL
class Texture
{
public:
   /// The palette to create a texture of
   Palette* const palette;
   /// Repeat the palette multiple times (negative for reverse); cannot be 0
   int Nrepeat_;
   /// Number of colors to discretize with (0 uses the original number of colors)
   int Ncolors_;
   /// Is texture smooth or discrete?
   bool smooth;
   /// Constructor - generates texture
   Texture(Palette* palette, int Nrepeat_ = 1, int Ncolors_ = 0,
           bool smooth = false);
   /// Texture size
   int size() const { return texture_data.size(); }
   /// Get texture data
   const vector<array<float,4>>& texture() const { return texture_data; }
   /// If true, all colors in palette are read in reverse
   bool isReversed() const { return Nrepeat_ < 0; }
   /// Generates the texture data
   void generate();
private:
   int MAX_TEXTURE_SIZE;
   vector<array<float,4>> texture_data;
};


/// Behaves like make_unique (only available in >= c++14)
template<typename T, typename... Args>
std::unique_ptr<T> as_unique(Args&&... args)
{
   return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/// Holds a vector of unique_ptr<Palette>. Besides holding the
/// palettes this should be stateless (put state in PaletteState)
class PaletteRegistry
{
private:
   vector<unique_ptr<Palette>> palettes;

   /// Find the index of a palette by name
   int get_index_by_name(const string& name) const;

public:
   /// Empty constructor
   PaletteRegistry() {}

   /// Constructor via a const vector of Palettes; if name already exists, skip
   /// Used for loading compiled palettes (i.e. `palette_definitions.cpp`)
   PaletteRegistry(const vector<Palette>& paletteRefs);

   /// Adds an existing palette to the registry
   void addPalette(Palette& palette);

   /// Create a new palette with the given name and add it to the registry
   void addPalette(const string& name);

   /// Returns true if name is unique
   bool check_name(const string& name) const;

   /// Get a palette pointer by index; if not found, returns last palette
   Palette* get(int index) const;

   /// Get a palette pointer by name; if not found, returns last palette
   Palette* get(const string& name) const;

   /// Prints a summary (index + name) of all palettes
   void printSummary(ostream& os = cout) const;

   /// Prints all colors for all palettes
   void printAll(ostream& os = cout) const;

   /// Number of palettes in the registry
   int NumPalettes() const { return palettes.size(); }

   /* Loads palette(s) from a file. Format is:

      palette <palette_name> <RGBf/RGBAf>
      <r1> <g1> <b1> [<a1>]
      <r2> <g2> <b2> [<a2>]
      ...

      see `share/palettes-crameri.txt` for an example */
   void load(const string& palette_filename);
};


extern PaletteRegistry BasePalettes;

#endif
