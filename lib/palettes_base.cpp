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

#include "palettes_base.hpp"
#include "palettes_default.cpp"


void RGBAf::Print(ostream& os) const
{
   os << fixed << setprecision(6)
      << setw(10) << r << " "
      << setw(10) << g << " "
      << setw(10) << b << " "
      << setw(10) << a;
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
   int j = reversed ? size() - 1 - i : i;
   return colors[j];
}

vector<array<float,4>> Palette::GetData(bool reversed) const
{
   vector<array<float,4>> rgba_data(size());
   for (int i = 0; i < size(); ++i)
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

// Initialize
int Texture::max_texture_size = 4096;

Texture::Texture(const Palette* palette, int Nrepeat_, int Ncolors_,
                 bool smooth) : palette(palette), Nrepeat_(Nrepeat_),
   Ncolors_(Ncolors_), smooth(smooth)
{
   if (Texture::max_texture_size == 0)
   {
      glGetIntegerv(GL_max_texture_size, &Texture::max_texture_size);
   }
   // Generate the texture data
   Generate();
}

void Texture::Generate()
{
   // Nrepeat cannot be 0; we also extract the sign
   bool reversed = IsReversed();
   int Nrepeat = Nrepeat_ == 0 ? 1 : abs(Nrepeat_);
   // Ncolors must be positive
   int Ncolors = Ncolors_ <= 0 ? palette->size() : Ncolors_;

   // Original palette size
   int plt_size = palette->size();
   // Set the texture size
   int tsize = Nrepeat * Ncolors;
   if (tsize > Texture::max_texture_size)
   {
      cerr << "Warning: Texture size "
           << "(" << tsize << ")" << " exceeds maximum "
           << "(" << Texture::max_texture_size << ")" << endl;
      if (Ncolors >= Texture::max_texture_size)
      {
         Ncolors = Texture::max_texture_size;
         Nrepeat = 1;
         tsize = Nrepeat * Ncolors;
      }
      else
      {
         Nrepeat = Texture::max_texture_size / Ncolors;
         tsize = Nrepeat * Ncolors;
      }
   }
   texture_data.clear();
   texture_data.resize(tsize);

   // Generate the discrete texture data
   if (!smooth)
   {
      for (int rpt = 0; rpt < Nrepeat; rpt++)
      {
         bool reverse = (reversed + rpt) % 2 != 0;
         for (int i = 0; i < Ncolors; i++)
         {
            int j = std::min(i * plt_size / (Ncolors - 1), plt_size - 1);
            texture_data[rpt*Ncolors + i] = palette->Color(j, reverse).AsArray();
         }
      }
   }
   // Generate the smooth texture data (interpolates colors)
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
            array<float,4> col1 = palette->Color(j, reverse).AsArray();
            array<float,4> col2 = palette->Color(j+1, reverse).AsArray();
            texture_data[rpt*Ncolors + i] =
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
         if (idx == -1)
         {
            AddPalette(palname);
            idx = GetIndexByName(palname);
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
