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

#ifndef GLVIS_GLTF_HPP
#define GLVIS_GLTF_HPP

#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <fstream>
#include <limits>
#include <memory>


class glTF_Builder
{
public:
   typedef std::array<float,2> vec2f;
   typedef std::array<float,3> vec3f;
   typedef std::array<float,4> color4f;
   typedef std::vector<float> vecnf;

protected:
   template <typename T>
   struct node_type
   {
      bool valid;
      std::string key;
      T value;
   };

   typedef node_type<bool>        node_bool;
   typedef node_type<unsigned>    node_unsigned;
   typedef node_type<float>       node_float;
   typedef node_type<std::string> node_string;
   typedef node_type<vec3f>       node_vec3f;
   typedef node_type<color4f>     node_color4f;
   typedef node_type<vecnf>       node_vecnf;

   struct struct_buffer
   {
      node_string uri = { false, "uri", "" };
      node_unsigned byteLength = { false, "byteLength", 0 };

      std::unique_ptr<std::ofstream> file;
   };
   struct struct_buffer_view
   {
      node_unsigned buffer = { false, "buffer", 0 };
      node_unsigned byteOffset = { false, "byteOffset", 0 };
      node_unsigned byteLength = { false, "byteLength", 0 };
      node_unsigned byteStride = { false, "byteStride", 0};
      node_unsigned target = { false, "target", 0 };
   };
   struct struct_accessor
   {
      node_unsigned bufferView = { false, "bufferView", 0 };
      node_unsigned byteOffset = { false, "byteOffset", 0 };
      node_unsigned componentType = { false, "componentType", 0 };
      node_unsigned count = { false, "count", 0 };
      node_string type = { false, "type", "" };
      node_vecnf min = { false, "min", {} };
      node_vecnf max = { false, "max", {} };
      node_bool normalized = { false, "normalized", false };
   };
   struct struct_image
   {
      node_string uri = { false, "uri", "" };
      node_string name = { false, "name", "" };
   };
   struct struct_sampler
   {
      node_unsigned magFilter = { false, "magFilter", 0 };
      node_unsigned minFilter = { false, "minFilter", 0 };
      node_unsigned wrapS = { false, "wrapS", 0 };
      node_unsigned wrapT = { false, "wrapT", 0 };
   };
   struct struct_texture
   {
      node_unsigned sampler = { false, "sampler", 0 };
      node_unsigned source = { false, "source", 0 };
   };
   struct struct_texture_info
   {
      node_unsigned index = { false, "index", 0 };
      node_unsigned texCoord = { false, "texCoord", 0 };
   };
   struct struct_pbrMetallicRoughness
   {
      node_color4f baseColorFactor =
      { false, "baseColorFactor", {{0.f, 0.f, 0.f, 0.f}} };
      node_type<struct_texture_info> baseColorTexture =
      { false, "baseColorTexture", {} };
      node_float metallicFactor = { false, "metallicFactor", 0.f };
      node_float roughnessFactor = { false, "roughnessFactor", 0.f };
   };
   struct struct_material
   {
      node_type<struct_pbrMetallicRoughness> pbrMetallicRoughness =
      { false, "pbrMetallicRoughness", {} };
      node_bool doubleSided = { false, "doubleSided", false };
      node_string name = { false, "name", "" };
   };
   struct struct_attributes
   {
      node_unsigned POSITION = { false, "POSITION", 0 };
      node_unsigned NORMAL = { false, "NORMAL", 0 };
      node_unsigned TEXCOORD_0 = { false, "TEXCOORD_0", 0 };
      node_unsigned COLOR_0 = { false, "COLOR_0", 0 };

   };
   struct struct_primitive
   {
      node_type<struct_attributes> attributes = { false, "attributes", {} };
      node_unsigned indices = { false, "indices", 0 };
      node_unsigned material = { false, "material", 0 };
      node_unsigned mode = { false, "mode", 0 };
   };
   struct struct_mesh
   {
      std::vector<struct_primitive> primitives;
      node_string name = { false, "name", "" };
   };
   struct struct_node
   {
      node_unsigned mesh = { false, "mesh", 0 };
      node_vec3f scale = { false, "scale", {{0.f, 0.f, 0.f}} };
      node_vec3f translation = { false, "translation", {{0.f, 0.f, 0.f}} };
      node_string name = { false, "name", "" };
   };

   // begin: printing functions
   static const char *sep(size_t i) { return (i==0) ? "" : ","; }
   template <typename T>
   static void print_node(std::ostream &out,
                          int &pfx_counter,
                          const std::string &pfx,
                          const node_type<T> &n)
   {
      if (n.valid)
      {
         out << sep(pfx_counter++) << pfx << '"' << n.key << "\" : ";
         print(out, n.value);
      }
   }

   static void print(std::ostream &out, const bool &v) { out << v; }
   static void print(std::ostream &out, const unsigned &v) { out << v; }
   static void print(std::ostream &out, const float &v) { out << v; }
   static void print(std::ostream &out, const std::string &v)
   { out << '"' << v << '"'; }
   template <typename T, size_t s>
   static void print(std::ostream &out, const std::array<T,s> &v)
   {
      out << '[';
      for (size_t i = 0; i != s; ++i) { out << sep(i) << ' ' << v[i]; }
      out << " ]";
   }
   template <typename T>
   static void print(std::ostream &out, const std::vector<T> &v)
   {
      out << '[';
      for (size_t i = 0; i != v.size(); ++i) { out << sep(i) << ' ' << v[i]; }
      out << " ]";
   }
   static void print(std::ostream &out, const struct_attributes &a)
   {
      out << '{';
      int pos = 0;
      print_node(out, pos, "\n        ", a.POSITION);
      print_node(out, pos, "\n        ", a.NORMAL);
      print_node(out, pos, "\n        ", a.TEXCOORD_0);
      print_node(out, pos, "\n        ", a.COLOR_0);
      out << "\n      }";
   }
   static void print(std::ostream &out, const struct_texture_info &ti)
   {
      out << '{';
      int pos = 0;
      print_node(out, pos, "\n        ", ti.index);
      print_node(out, pos, "\n        ", ti.texCoord);
      out << "\n      }";
   }
   static void print(std::ostream &out, const struct_pbrMetallicRoughness &pbr)
   {
      out << '{';
      int pos = 0;
      print_node(out, pos, "\n      ", pbr.baseColorFactor);
      print_node(out, pos, "\n      ", pbr.baseColorTexture);
      print_node(out, pos, "\n      ", pbr.metallicFactor);
      print_node(out, pos, "\n      ", pbr.roughnessFactor);
      out << "\n    }";
   }
   // end: printing functions

   const std::string file_prefix;

   std::vector<struct_buffer> buffers;
   std::vector<struct_buffer_view> buffer_views;
   std::vector<struct_accessor> accessors;
   std::vector<struct_image> images;
   std::vector<struct_sampler> samplers;
   std::vector<struct_texture> textures;
   std::vector<struct_material> materials;
   std::vector<struct_mesh> meshes;
   std::vector<struct_node> nodes;

public:
   static constexpr unsigned INVALID_ID = std::numeric_limits<unsigned>::max();
   typedef struct { unsigned id; } buffer_id;
   typedef struct { unsigned id; } buffer_view_id;
   typedef struct { unsigned id; } accessor_id;
   typedef struct { unsigned id; } image_id;
   typedef struct { unsigned id; } sampler_id;
   typedef struct { unsigned id; } texture_id;
   typedef struct { unsigned id; } material_id;
   typedef struct { unsigned id; } mesh_id;
   typedef struct { unsigned id; } node_id;

   // buffer view option
   enum struct target_type
   {
      ARRAY_BUFFER = 34962,
      ELEMENT_ARRAY_BUFFER = 34963
   };

   // accessor option
   enum struct component_type
   {
      BYTE           = 5120,
      UNSIGNED_BYTE  = 5121,
      SHORT          = 5122,
      UNSIGNED_SHORT = 5123,
      UNSIGNED_INT   = 5125,
      FLOAT          = 5126
   };
   // accessor option
   enum struct tensor_type
   {
      SCALAR = 0, VEC2, VEC3, VEC4, MAT2, MAT3, MAT4
   };
   // string constants corresponding to the tensor_type constants
   static const char *tensorTypes[];

   // sampler option: magnification filter
   enum struct mag_filter { NEAREST = 9728, LINEAR = 9729 };
   // sampler option: minification filter
   enum struct min_filter
   {
      NEAREST = 9728, LINEAR = 9729, NEAREST_MIPMAP_NEAREST = 9984,
      LINEAR_MIPMAP_NEAREST = 9985, NEAREST_MIPMAP_LINEAR = 9986,
      LINEAR_MIPMAP_LINEAR = 9987
   };
   // sampler option: S/T (or U/V) wrapping mode
   enum struct wrap_type
   {
      CLAMP_TO_EDGE = 33071, MIRRORED_REPEAT = 33648, REPEAT = 10497
   };

   struct pbr_matallic_roughness
   {
      bool haveTexture;  // if true, baseColorTexture must be defined
      color4f baseColorFactor;
      texture_id baseColorTexture;
      float metallicFactor;
      float roughnessFactor;
   };


   glTF_Builder(const std::string &filePrefix)
      : file_prefix(filePrefix)
   { }

   buffer_id addBuffer(const std::string &bufferName);

   buffer_view_id addBufferView(buffer_id buffer,
                                const void *data,
                                size_t byteLength,
                                size_t byteStride,
                                size_t byteAlign,
                                target_type target);

   // can be called after the call to addBufferView() that created the given
   // bufferView but only before the next call to addBufferView()
   void appendToBufferView(buffer_view_id bufferView,
                           const void *data,
                           size_t byteLength);

   // count must be >= 1
   accessor_id addAccessor(buffer_view_id bufferView,
                           size_t byteOffset,
                           component_type componentType,
                           size_t count,
                           tensor_type tensorType);

   // count must be >= 1
   accessor_id addAccessorVec2f(buffer_view_id bufferView,
                                size_t byteOffset,
                                size_t count,
                                vec2f min,
                                vec2f max);

   // count must be >= 1
   accessor_id addAccessorVec3f(buffer_view_id bufferView,
                                size_t byteOffset,
                                size_t count,
                                vec3f min,
                                vec3f max);

   image_id addImage(const std::string &imageName,
                     int width,
                     int height,
                     const color4f *pixels);

   sampler_id addSampler(mag_filter magFilter = mag_filter::NEAREST,
                         min_filter minFilter = min_filter::NEAREST,
                         wrap_type wrapS = wrap_type::CLAMP_TO_EDGE,
                         wrap_type wrapT = wrap_type::CLAMP_TO_EDGE);

   texture_id addTexture(sampler_id sampler, image_id source);

   material_id addMaterial(const std::string &materialName,
                           const pbr_matallic_roughness &pbrMetallicRoughness,
                           bool doubleSided = false);

   mesh_id addMesh(const std::string &meshName);

   void addMeshTriangles(mesh_id mesh,
                         accessor_id vertexPositions,
                         accessor_id vertexNormals,
                         accessor_id vertexTexCoords0,
                         accessor_id vertexIndices,
                         material_id material);

   void addMeshLines(mesh_id mesh,
                     accessor_id vertexPositions,
                     accessor_id vertexTexcoords0,
                     accessor_id vertexColors0,
                     material_id material);

   node_id addNode(const std::string &nodeName);

   void addNodeMesh(node_id node, mesh_id mesh);

   void addNodeScale(node_id node, vec3f scale);

   void addNodeTranslation(node_id node, vec3f translation);

   void getMaterialPBRMR(material_id material,
                         pbr_matallic_roughness &pbr_mr_copy);

   int writeFile();
};

#endif // GLVIS_GLTF_HPP
