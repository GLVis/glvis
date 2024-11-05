// base_palettes.hpp
#ifndef GLVIS_BASEPALETTES_HPP
#define GLVIS_BASEPALETTES_HPP

#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <vector>
#include <cctype>
#include <iomanip> // Include iomanip for std::setw and std::setfill

using namespace std;

struct RGBAf {
   float r, g, b, a;

   constexpr RGBAf(float r = 0.0, float g = 0.0, float b = 0.0, float a = 1.0)
      : r(r), g(g), b(b), a(a) {}

   void print() const {
      cout << fixed << setprecision(6)
           << setw(10) << r << " "
           << setw(10) << g << " "
           << setw(10) << b << " "
           << setw(10) << a;// << endl;
   }

   array<float, 4> as_array() const {
      return {r, g, b, a};
   }

};


struct Palette {
   vector<RGBAf> colors;
   const string name;

   // --- Constructors ---
   Palette(const string& name) : name(name) {}
   Palette(const string& name, const vector<RGBAf>& colors)
      : name(name), colors(colors) {}

   // from Nx3 array
   template <size_t N>
   // Palette(const string& name, const float (&array)[N][3]) : name(name) {
   Palette(const string& name, const array<array<float,3>,N>& arr) : name(name) {
      for (size_t i = 0; i < N; ++i) {
         colors.push_back(RGBAf(arr[i][0], arr[i][1], arr[i][2]));
      }
   }

   // from Nx4 array
   template <size_t N>
   Palette(const string& name, const array<array<float,4>,N>& arr) : name(name) {
      for (size_t i = 0; i < N; ++i) {
         colors.push_back(RGBAf(arr[i][0], arr[i][1], arr[i][2], arr[i][3]));
      }
   }

   // --- Getters ---
   size_t size() const {
      return colors.size();
   }

   // --- Add color ---
   void addColor(float r, float g, float b, float a = 1.0) {
      colors.push_back(RGBAf(r, g, b, a));
   }

   // --- Print ---
   void print() const {
      cout << "palette " << name << "\n";
      for (const auto& color : colors) {
         color.print();
         cout << endl;
      }
      cout << endl;
   }

   // helper function
   void printAsCPP(const string& filename, const string& varname) const {
      ofstream outfile(filename, ios::app); // Open file in append mode
      if (!outfile) {
         cerr << "Could not open file: " << filename << endl;
         return;
      }
      outfile << "const Palette " << varname << " = Palette(\"" << name << "\", {" << endl;
      for (const auto& color : colors) {
         outfile << "   {" << color.r << ", " << color.g << ", " << color.b << ", " << color.a << "},\n";
      }
      outfile << "});" << endl;
      outfile << endl;
      outfile.close();
   }

   shared_ptr<Palette> shared() const {
      return make_shared<Palette>(*this);
   }

   double* as_array() const {
      size_t N = colors.size();
      double* arr = new double[N * 4];
      for (size_t i = 0; i < N; ++i) {
         arr[i * 4 + 0] = colors[i].r;
         arr[i * 4 + 1] = colors[i].g;
         arr[i * 4 + 2] = colors[i].b;
         arr[i * 4 + 3] = colors[i].a;
      }
      return arr;
   }
};


// PaletteRegistry with a vector of shared pointers to Palette
// Besides holding the palettes, this should be stateless
class PaletteRegistry {
private:
   vector<shared_ptr<Palette>> palettes;

public:
   // empty constructor
   PaletteRegistry() {}

   PaletteRegistry(const vector<Palette>& paletteRefs) {
      for (const auto& palette : paletteRefs) {
         palettes.push_back(palette.shared());
      }
   }

   // PaletteRegistry(const vector<shared_ptr<Palette>>& paletteRefs)
   //    : palettes(paletteRefs) {}

   void addPalette(const Palette& palette) {
      // palette name is unique || container is empty
      if (get_index_by_name(palette.name) == -1 || palettes.empty()) {
         palettes.push_back(palette.shared());
      }
   }

   void addPalette(const string& name) {
      // palette name is unique || container is empty
      if (get_index_by_name(name) == -1 || palettes.empty()) {
         addPalette(Palette(name));
      }
   }

   // get by index
   shared_ptr<Palette> get(size_t index) const {
      if (1 <= index && index <= NumPalettes()) {
         return palettes[index-1];
      }
      cout << "Palette (index = " << index << ") out of range. Available palettes:" << endl;
      this->printSummary();
      return palettes[NumPalettes()-1];
   }

   // get by name
   shared_ptr<Palette> get(const string& name) const {
      int idx = get_index_by_name(name);
      if (idx != -1) {
         return palettes[idx];
      }
      cout << "Palette (name = " << name << ") not found. Available palettes:" << endl;
      this->printSummary();
      return palettes[NumPalettes()-1];
   }

   void printSummary() const {
      size_t idx = 1;
      for (const auto& palette : palettes) {
         cout << idx << ") "
              << left << setw(12) << palette->name << right;
         if (idx%5 == 0)
         {
            cout << '\n';
         }
         idx++;
      }
      cout << endl;
   }

   void printAll() const {
      for (const auto& palette : palettes) {
         palette->print();
      }
   }

   size_t NumPalettes() const {
      return palettes.size();
   }

   int get_index_by_name(const string& name) const {
      for (int i = 0; i < NumPalettes(); i++) {
         if (palettes[i]->name == name) {
            return i;
         }
      }
      return -1;
   }

   void load(const string& palette_filename) {

      ifstream pfile(palette_filename);
      if (!pfile)
      {
         cout << "Could not open palette file: " << palette_filename << endl;
         return;
      }
      string word, palname, channeltype;
      int idx = -1;

      // read initializing commands
      while (1) {
         pfile >> ws;
         if (!pfile.good()) {
            // cout << "Error in palette" << endl;
            break;
         }
         if (pfile.peek() == '#') {
            getline(pfile, word);
            continue;
         }
         pfile >> word;
         if (word == "palette") {
            pfile >> palname >> channeltype;
            idx = get_index_by_name(palname);
            if (idx == -1) {
               addPalette(palname);
               idx = get_index_by_name(palname);
               cout << "Reading palette: (" << idx+1 << ") " << palname << endl;
            } else {
               cout << "Error reading palette: " << palname
                    << ". Palette with same name already exists." << endl;
               break;
            }
         }
         else if (channeltype == "RGBf" && idx != -1) {
            float r, g, b;
            r = stof(word);
            pfile >> g >> b;
            palettes[idx]->addColor(r,g,b);
         }
         else {
            cout << "Error reading palette file: " << palette_filename << endl;
            break;
         }
      }
      cout << "Finished loading palettes from file." << endl;
   }

};


extern PaletteRegistry BasePalettes;

#endif