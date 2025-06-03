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

#include "palettes_base.hpp"
#include "gl/renderer.hpp"

#include <iomanip>
#include <fstream>

void RGBAf::Print(ostream& os) const
{
   os << std::fixed << std::setprecision(6)
      << std::setw(10) << r << " "
      << std::setw(10) << g << " "
      << std::setw(10) << b << " "
      << std::setw(10) << a;
}

template <size_t N>
Palette::Palette(const string& name,
                 const array<array<float,3>,N>& arr) : name(name)
{
   colors.resize(N);
   for (size_t i = 0; i < N; ++i)
   {
      colors[i] = RGBAf(arr[i][0], arr[i][1], arr[i][2]);
   }
}

template <size_t N>
Palette::Palette(const string& name,
                 const array<array<float,4>,N>& arr) : name(name)
{
   colors.resize(N);
   for (size_t i = 0; i < N; ++i)
   {
      colors[i] = RGBAf(arr[i][0], arr[i][1], arr[i][2], arr[i][3]);
   }
}


void Palette::AddColor(float r, float g, float b, float a)
{
   colors.push_back(RGBAf(r, g, b, a));
}


void Palette::Print(ostream& os) const
{
   os << "palette " << name << " RGBf" << endl;
   for (const auto& color : colors)
   {
      color.Print(os);
      os << endl;
   }
   os << endl;
}

RGBAf Palette::Color(int i, bool reversed) const
{
   int j = reversed ? Size() - 1 - i : i;
   return colors[j];
}

vector<array<float,4>> Palette::GetData(bool reversed) const
{
   vector<array<float,4>> rgba_data(Size());
   for (int i = 0; i < Size(); ++i)
   {
      rgba_data[i] = Color(i, reversed).AsArray();
   }
   return rgba_data;
}

bool Palette::IsTranslucent() const
{
   for (const auto& color : colors)
   {
      if (color.a != 1.0) { return true; }
   }
   return false;
}

// Initialize GL parameters
int Texture::max_texture_size = -1;
// WebGL 2 requires sized internal format for float texture
GLenum Texture::alpha_internal = GL_R32F;
GLenum Texture::alpha_channel = GL_RED;
GLenum Texture::rgba_internal = GL_RGBA32F;
GLenum Texture::rgba_channel = GL_RGBA;

Texture::Texture(const Palette* palette,
                 TextureType textype,
                 int cycles, int colors)
   : palette(palette), textype(textype)
{
   // Initialize static GL parameters
   InitStaticGL();

   // Input sanitization/init
   UpdateParameters(cycles, colors);

   // Generate the texture id
   GLuint texid;
   glGenTextures(1, &texid);
   texture = texid;
}

Texture::Texture(float matAlpha, float matAlphaCenter)
   : palette(nullptr), textype(TextureType::ALPHAMAP)
{
   // Initialize static GL parameters
   InitStaticGL();

   // Input sanitization/init
   UpdateAlphaParameters(matAlpha, matAlphaCenter);

   // Generate the texture id
   GLuint texid;
   glGenTextures(1, &texid);
   texture = texid;

   // set alpha texture to 1.0
   std::vector<float> texture_data(Texture::max_texture_size);
   std::fill(texture_data.begin(), texture_data.end(), 1.0f);
   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, texture);
   glTexImage2D(GL_TEXTURE_2D, 0, Texture::alpha_internal,
                Texture::max_texture_size, 1, 0,
                Texture::alpha_channel, GL_FLOAT, texture_data.data());
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

   glActiveTexture(GL_TEXTURE0);
}

void Texture::InitStaticGL()
{
   if (Texture::max_texture_size < 0)
   {
      glGetIntegerv(GL_MAX_TEXTURE_SIZE, &Texture::max_texture_size);

      if (gl3::GLDevice::useLegacyTextureFmts())
      {
         Texture::alpha_internal = GL_ALPHA;
         Texture::alpha_channel = GL_ALPHA;
         Texture::rgba_internal = GL_RGBA;
      }
   }
}

vector<float> Texture::GenerateAlphaTextureData()
{
   vector<float> texture_data(Texture::max_texture_size);
   if (alpha >= 1.0)
   {
      // transparency off
      std::fill(texture_data.begin(), texture_data.end(), 1.0f);
   }
   else
   {
      for (int i = 0; i < Texture::max_texture_size; i++)
      {
         float val = float(2*i + 1)/(2*Texture::max_texture_size); // midpoint of texel
         texture_data[i] = alpha * std::exp(-std::abs(val - alpha_center));
      }
   }
   return texture_data;
}

vector<array<float,4>> Texture::GenerateTextureData()
{
   // Original palette size
   int plt_size = palette->Size();

   // Initialize the texture data
   vector<array<float,4>> texture_data(tsize);

   // Discrete texture
   if ( textype == TextureType::DISCRETE )
   {
      // Generate the texture data
      for (int rpt = 0; rpt < nrepeat; rpt++)
      {
         bool reverse = (reversed + rpt) % 2 != 0;
         for (int i = 0; i < ncolors; i++)
         {
            int j = std::min(i * plt_size / (ncolors - 1), plt_size - 1);
            texture_data[rpt*ncolors + i] = palette->Color(j, reverse).AsArray();
         }
      }
   }
   // Smooth texture (interpolates colors)
   else if ( textype == TextureType::SMOOTH )
   {
      // Generate the texture data
      for (int rpt = 0; rpt < nrepeat; rpt++)
      {
         bool reverse = (reversed + rpt) % 2 != 0;
         for (int i = 0; i < ncolors; i++)
         {
            float t = i * (plt_size - 1) / (ncolors - 1);
            int j = std::min((int)t, plt_size - 2);
            t -= j;
            array<float,4> col1 = palette->Color(j, reverse).AsArray();
            array<float,4> col2 = palette->Color(j+1, reverse).AsArray();
            texture_data[rpt*ncolors + i] =
            {
               (1-t) * col1[0] + t * col2[0],
               (1-t) * col1[1] + t * col2[1],
               (1-t) * col1[2] + t * col2[2],
               (1-t) * col1[3] + t * col2[3]
            };
         }
      }
   }
   return texture_data;
}

void Texture::SetCycles(int cycles)
{
   if (cycles == 0)
   {
      cycles = 1;
   }
   reversed = cycles < 0;
   nrepeat = abs(cycles);
}

void Texture::SetColors(int colors)
{
   ncolors = colors <= 0 ? palette->Size() : colors;
}

void Texture::UpdateTextureSize()
{
   tsize = nrepeat * ncolors;
   if (tsize > Texture::max_texture_size)
   {
      cerr << "Warning: Texture size "
           << "(" << tsize << ")" << " exceeds maximum "
           << "(" << Texture::max_texture_size << ")" << endl;
      if (ncolors >= Texture::max_texture_size)
      {
         ncolors = Texture::max_texture_size;
         nrepeat = 1;
      }
      else
      {
         nrepeat = Texture::max_texture_size / ncolors;
      }
      tsize = nrepeat * ncolors;
   }
}

void Texture::UpdateAlphaParameters(float matAlpha, float matAlphaCenter)
{
   alpha = std::max(0.0f, std::min(1.0f, matAlpha));
   alpha_center = std::max(0.0f, std::min(1.0f, matAlphaCenter));
   tsize = Texture::max_texture_size;
}

void Texture::UpdateParameters(int cycles, int colors)
{
   SetCycles(cycles);
   SetColors(colors);
   UpdateTextureSize();
}

void Texture::Generate()
{
   // Get texture data and formats (different for alpha textures)
   // Define the texture image
   if ( textype == TextureType::ALPHAMAP )
   {
      vector<float> texture_data = GenerateAlphaTextureData();
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, texture);
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
                      tsize, 1,
                      Texture::alpha_channel, GL_FLOAT, texture_data.data());
      glActiveTexture(GL_TEXTURE0);
   }
   else
   {
      vector<array<float,4>> texture_data = GenerateTextureData();
      glBindTexture(GL_TEXTURE_2D, texture);
      glTexImage2D(GL_TEXTURE_2D, 0, Texture::rgba_internal,
                   tsize, 1, 0,
                   Texture::rgba_channel, GL_FLOAT, texture_data.data());
   }

   // Discrete
   if ( textype == TextureType::DISCRETE)
   {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   }
   // Smooth
   else if ( textype == TextureType::SMOOTH )
   {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   }
}

int PaletteRegistry::GetIndexByName(const string& name) const
{
   for (int i = 0; i < NumPalettes(); i++)
   {
      if (Get(i)->name == name) { return i; }
   }
   return -1;
}

PaletteRegistry::PaletteRegistry(const vector<Palette>& paletteRefs)
{
   for (const Palette& palette : paletteRefs)
   {
      if (IsNameUnique(palette.name))
      {
         palettes.push_back(as_unique<Palette>(palette));
      }
   }
}

void PaletteRegistry::AddPalette(Palette& palette)
{
   if (IsNameUnique(palette.name))
   {
      palettes.push_back(as_unique<Palette>(palette));
   }
}

void PaletteRegistry::AddPalette(const string& name)
{
   if (IsNameUnique(name))
   {
      palettes.push_back(as_unique<Palette>(name));
   }
}

bool PaletteRegistry::IsNameUnique(const string& name) const
{
   // palette name is unique || container is empty
   if (GetIndexByName(name) == -1 || palettes.empty())
   {
      return true;
   }
   else
   {
      cout << "Palette with name: '" << name << "' already exists in registry.";
      return false;
   }
}

Palette* PaletteRegistry::Get(int index) const
{
   if (0 <= index && index <= NumPalettes()-1)
   {
      return palettes[index].get();
   }
   cout << "Palette (index = " << index+1 << ") out of range. Available palettes:"
        << endl;
   this->PrintSummary();
   return palettes.back().get();
}

Palette* PaletteRegistry::Get(const string& name) const
{
   int idx = GetIndexByName(name);
   if (idx != -1)
   {
      return palettes[idx].get();
   }
   cout << "Palette (name = " << name << ") not found. Available palettes:" <<
        endl;
   this->PrintSummary();
   return palettes.back().get();
}

void PaletteRegistry::PrintSummary(ostream& os) const
{
   for (int i = 0; i < NumPalettes(); i++)
   {
      os << setw(3) << i+1 << ") "
         << left << setw(12) << Get(i)->name << right;
      if ((i+1)%5 == 0)
      {
         os << endl;
      }
   }
   os << endl;
}

void PaletteRegistry::PrintAll(ostream& os) const
{
   for (int i = 0; i < NumPalettes(); i++)
   {
      Get(i)->Print(os);
   }
}

void PaletteRegistry::Load(const string& palette_filename)
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
         idx = GetIndexByName(palname);
         if (idx != -1)
         {
            cout << "Warning: palette name <" << palname
                 << "> already exists. Overriding..." << endl;
            palettes.erase(palettes.begin() + idx);
         }
         // Add new palette
         AddPalette(palname);
         idx = GetIndexByName(palname);
         cout << "Reading palette: (" << idx+1 << ") " << palname << endl;
      }
      else if (channeltype == "RGBf" && idx != -1)
      {
         float r, g, b;
         r = stof(word);
         pfile >> g >> b;
         Get(idx)->AddColor(r,g,b);
      }
      else if (channeltype == "RGBAf" && idx != -1)
      {
         float r, g, b, a;
         r = stof(word);
         pfile >> g >> b >> a;
         Get(idx)->AddColor(r,g,b,a);
      }
      else
      {
         cout << "Error reading palette file: " << palette_filename << endl;
         break;
      }
   }
   cout << "Finished loading palettes from file: " << palette_filename << endl;
}
