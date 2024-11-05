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

#ifndef GLVIS_BASEPALETTES_HPP
#define GLVIS_BASEPALETTES_HPP

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


struct Palette
{
   const string name;
   vector<RGBAf> colors;

   // --- Constructors ---
   Palette(const string& name) : name(name) {}
   Palette(const string& name, const vector<RGBAf>& colors)
      : name(name), colors(colors) {}

   // from Nx3 array
   template <size_t N>
   // Palette(const string& name, const float (&array)[N][3]) : name(name) {
   Palette(const string& name, const array<array<float,3>,N>& arr) : name(name)
   {
      for (size_t i = 0; i < N; ++i)
      {
         colors.push_back(RGBAf(arr[i][0], arr[i][1], arr[i][2]));
      }
   }

   // from Nx4 array
   template <size_t N>
   Palette(const string& name, const array<array<float,4>,N>& arr) : name(name)
   {
      for (size_t i = 0; i < N; ++i)
      {
         colors.push_back(RGBAf(arr[i][0], arr[i][1], arr[i][2], arr[i][3]));
      }
   }

   // --- Getters ---
   size_t size() const
   {
      return colors.size();
   }

   // --- Add color ---
   void addColor(float r, float g, float b, float a = 1.0)
   {
      colors.push_back(RGBAf(r, g, b, a));
   }

   // --- Print ---
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

   // helper function - not used in glvis.cpp
   void printAsCPP(const string& filename, const string& varname) const
   {
      ofstream outfile(filename, ios::app); // Open file in append mode
      if (!outfile)
      {
         cerr << "Could not open file: " << filename << endl;
         return;
      }
      outfile << "const Palette " << varname << " = Palette(\"" << name << "\", {" <<
              endl;
      for (const auto& color : colors)
      {
         outfile << "   {" << color.r << ", " << color.g << ", " << color.b << ", " <<
                 color.a << "},\n";
      }
      outfile << "});" << endl;
      outfile << endl;
      outfile.close();
   }

   shared_ptr<Palette> shared() const
   {
      return make_shared<Palette>(*this);
   }

   double* as_rgb_array() const
   {
      size_t N = colors.size();
      double* arr = new double[N * 3];
      for (size_t i = 0; i < N; ++i)
      {
         arr[i * 3 + 0] = colors[i].r;
         arr[i * 3 + 1] = colors[i].g;
         arr[i * 3 + 2] = colors[i].b;
         // arr[i * 4 + 3] = colors[i].a;
      }
      return arr;
   }
};


// PaletteRegistry with a vector of shared pointers to Palette
// Besides holding the palettes, this should be stateless
class PaletteRegistry
{
private:
   vector<shared_ptr<Palette>> palettes;

public:
   // empty constructor
   PaletteRegistry() {}

   PaletteRegistry(const vector<Palette>& paletteRefs)
   {
      for (const auto& palette : paletteRefs)
      {
         palettes.push_back(palette.shared());
      }
   }

   void addPalette(const Palette& palette)
   {
      // palette name is unique || container is empty
      if (get_index_by_name(palette.name) == -1 || palettes.empty())
      {
         palettes.push_back(palette.shared());
      }
   }

   void addPalette(const string& name)
   {
      // palette name is unique || container is empty
      if (get_index_by_name(name) == -1 || palettes.empty())
      {
         addPalette(Palette(name));
      }
   }

   // get by index
   shared_ptr<Palette> get(size_t index) const
   {
      if (0 <= index && index <= NumPalettes()-1)
      {
         return palettes[index];
      }
      cout << "Palette (index = " << index+1 << ") out of range. Available palettes:"
           << endl;
      this->printSummary();
      return palettes[NumPalettes()-1];
   }

   // get by name
   shared_ptr<Palette> get(const string& name) const
   {
      int idx = get_index_by_name(name);
      if (idx != -1)
      {
         return palettes[idx];
      }
      cout << "Palette (name = " << name << ") not found. Available palettes:" <<
           endl;
      this->printSummary();
      return palettes[NumPalettes()-1];
   }

   void printSummary(ostream& os = cout) const
   {
      size_t idx = 1;
      for (const auto& palette : palettes)
      {
         os << setw(3) << idx << ") "
            << left << setw(12) << palette->name << right;
         if (idx%5 == 0)
         {
            os << endl;
         }
         idx++;
      }
      os << endl;
   }

   void printAll(ostream& os = cout) const
   {
      for (const auto& palette : palettes)
      {
         palette->print(os);
      }
   }

   size_t NumPalettes() const
   {
      return palettes.size();
   }

   int get_index_by_name(const string& name) const
   {
      for (size_t i = 0; i < NumPalettes(); i++)
      {
         if (palettes[i]->name == name)
         {
            return i;
         }
      }
      return -1;
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
            palettes[idx]->addColor(r,g,b);
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