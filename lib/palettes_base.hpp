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

   void Print(ostream& os = cout) const;

   array<float, 4> AsArray() const { return {r, g, b, a}; }
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
   int Size() const { return colors.size(); }

   /// Add color to palette
   void AddColor(float r, float g, float b, float a = 1.0);

   /// Print each color of this palette to a stream
   void Print(ostream& os = cout) const;

   /// Get color at index i (optionally, use reversed order)
   RGBAf Color(int i, bool reversed = false) const;

   /// Get all colors as a vector of float arrays
   vector<array<float,4>> GetData(bool reversed = false) const;

   /// Are any alpha != 1.0?
   bool IsTranslucent() const;

private:
   vector<RGBAf> colors;
};

using TexHandle = gl3::resource::TextureHandle;
/// Generates the texture data for a given palette, to be used in OpenGL
class Texture
{
public:
   // What type of texture is this (discrete, smooth, alphamap)
   enum class TextureType
   {
      DISCRETE = 0,
      SMOOTH = 1,
      ALPHAMAP = 2,
   };

   /// Empty constructor
   Texture() {}

   /// Constructor - generates texture
   Texture(const Palette* palette,
           TextureType textype = TextureType::DISCRETE,
           int cycles = 1, int colors = 0);

   /// Constructor for alphamap
   Texture(float matAlpha, float matAlphaCenter);

   /// Get texture size.
   int Size() { return tsize; }

   /// Get the GL texture
   GLuint Get() const { return texture; }

   /// Generate the GL texture and binds it to `texture`
   void Generate();

   /// Update alpha/regular texture parameters
   void UpdateAlphaParameters(float matAlpha, float matAlphaCenter);
   void UpdateParameters(int cycles, int colors);

private:
   /// The palette to create a texture of
   const Palette* palette;
   // What type of texture is this (discrete, smooth, alphamap)
   TextureType textype;

   /// GL static parameters
   static GLenum alpha_internal;
   static GLenum alpha_channel;
   static GLenum rgba_internal;
   static GLenum rgba_channel;
   static int max_texture_size;

   /// Repeat the palette multiple times (negative for reverse); cannot be 0
   int nrepeat;
   /// Number of colors to discretize with (0 uses the original number of colors)
   int ncolors;
   /// Is the texture reversed?
   bool reversed;
   /// Texture size
   int tsize;
   /// The GL texture
   TexHandle texture;
   /// Only used for alphamap
   float alpha;
   float alpha_center;

   /// Initialize Static GL parameters
   static void InitStaticGL();

   /// Generate alpha/regular texture data
   vector<float> GenerateAlphaTextureData();
   vector<array<float,4>> GenerateTextureData();

   /// Set the number of cycles
   void SetCycles(int cycles);

   /// Set the number of colors
   void SetColors(int colors);

   /// Update the texture size (may change ncolors and/or cycles if too large)
   void UpdateTextureSize();
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
   int GetIndexByName(const string& name) const;

public:
   /// Empty constructor
   PaletteRegistry() {}

   /// Constructor via a const vector of Palettes; if name already exists, skip
   /// Used for loading compiled palettes (i.e. `palettes_default.cpp`)
   PaletteRegistry(const vector<Palette>& paletteRefs);

   /// Adds an existing palette to the registry
   void AddPalette(Palette& palette);

   /// Create a new palette with the given name and add it to the registry
   void AddPalette(const string& name);

   /// Returns true if name is unique
   bool IsNameUnique(const string& name) const;

   /// Get a palette pointer by index; if not found, returns last palette
   Palette* Get(int index) const;

   /// Get a palette pointer by name; if not found, returns last palette
   Palette* Get(const string& name) const;

   /// Prints a summary (index + name) of all palettes
   void PrintSummary(ostream& os = cout) const;

   /// Prints all colors for all palettes
   void PrintAll(ostream& os = cout) const;

   /// Number of palettes in the registry
   int NumPalettes() const { return palettes.size(); }

   /// Loads palette(s) from a file.
   /**  Format is:
      palette <palette_name> <RGBf/RGBAf>
      <r1> <g1> <b1> [<a1>]
      <r2> <g2> <b2> [<a2>]
      ...

      see `share/palettes-crameri.txt` for an example */
   void Load(const string& palette_filename);
};


extern PaletteRegistry BasePalettes;

#endif
