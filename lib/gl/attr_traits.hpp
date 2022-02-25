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

#ifndef GLVIS_ATTR_TRAITS_HPP
#define GLVIS_ATTR_TRAITS_HPP

#include "types.hpp"
#include "renderer_core.hpp"
#include <type_traits>

namespace gl3
{

struct AttrNone
{
   static void setup() { }
   static void setupLegacy(void*) { }
   static void clear() { }
   static void clearLegacy() { }
   enum { exists = false };
};

// Base class to generate vertex attribute setup functions.
template<
   typename TV, typename TAttr, TAttr TV::*Attrib, typename TAttrInfo>
struct AttrBase
{
   constexpr static bool NormalizeAttr =
      std::is_integral<typename TAttr::value_type>::value;

   constexpr static TAttr* getAttrOffset()
   {
      return &(((TV*)0)->*Attrib);
   }

   // Sets up vertex attributes in the currently-bound buffer object.
   static void setup()
   {
      glEnableVertexAttribArray(TAttrInfo::ShaderIdx);
      glVertexAttribPointer(TAttrInfo::ShaderIdx,
                            std::tuple_size<TAttr>::value,
                            TAttrInfo::AttrGLType,
                            NormalizeAttr,
                            sizeof(TV),
                            (void*) getAttrOffset());
   }

   // Sets up client-side vertex pointers for the given buffer.
   static void setupLegacy(TV* buffer)
   {
      glEnableClientState(TAttrInfo::FFArrayIdx);
      TAttrInfo::FFSetupFunc(std::tuple_size<TAttr>::value,
                             TAttrInfo::AttrGLType,
                             sizeof(TV),
                             (char*) buffer + (size_t) getAttrOffset());
   }

   // Disables the attribute array.
   static void clear()
   {
      glDisableVertexAttribArray(TAttrInfo::ShaderIdx);
   }

   // Disables the client-side vertex array.
   static void clearLegacy()
   {
      glDisableClientState(TAttrInfo::FFArrayIdx);
   }

   enum { exists = true };
};

// Default attribute traits for vertex types. Provides no-op setup/clear
// functions if an attribute doesn't exist.
template<typename TV, typename = int>
struct AttrCoord : AttrNone { };

template<typename TV, typename = int>
struct AttrNormal : AttrNone { };

template<typename TV, typename = int>
struct AttrColor : AttrNone { };

template<typename TV, typename = int>
struct AttrTexcoord : AttrNone { };

// Template specializations for attribute traits. If an attribute exists in a
// vertex, generates setup and clear functions which setup OpenGL vertex
// attributes.
//
// Each attribute specialization defines static parameters which are passed to
// the AttrBase base class via CRTP to generate the attribute setup functions:
//  - AttrGLType: the OpenGL type of the attribute data
//  - ShaderIdx: the index of the generic vertex attribute
//  - FFArrayIdx: the index of the client-side vertex attribute array
//  - FFSetupFunc: the function to use when setting up the FF array pointer;
//      this can either be a direct pointer to a gl*Pointer function or a custom
//      function
template<typename TV>
struct AttrCoord<TV, decltype((void)TV::coord, 0)>
: AttrBase<TV, decltype(TV::coord), &TV::coord,
AttrCoord<TV, decltype((void)TV::coord, 0)>>
{
   const static GLenum AttrGLType = GL_FLOAT;
   const static int ShaderIdx = CoreGLDevice::ATTR_VERTEX;
   const static GLenum FFArrayIdx = GL_VERTEX_ARRAY;
   constexpr static auto FFSetupFunc = glVertexPointer;
};

template<typename TV>
struct AttrNormal<TV, decltype((void)TV::norm, 0)>
: AttrBase<TV, decltype(TV::norm), &TV::norm,
AttrNormal<TV, decltype((void)TV::norm, 0)>>
{
   const static GLenum AttrGLType = GL_FLOAT;
   const static int ShaderIdx = CoreGLDevice::ATTR_NORMAL;
   const static GLenum FFArrayIdx = GL_NORMAL_ARRAY;
   static void FFSetupFunc(GLint /*size*/, GLenum type, GLsizei stride,
                           const GLvoid* ptr)
   {
      glNormalPointer(type, stride, ptr);
   }
};

template<typename TV>
struct AttrColor<TV, decltype((void)TV::color, 0)>
: AttrBase<TV, decltype(TV::color), &TV::color,
AttrColor<TV, decltype((void)TV::color, 0)>>
{
   const static GLenum AttrGLType = GL_UNSIGNED_BYTE;
   const static int ShaderIdx = CoreGLDevice::ATTR_COLOR;
   const static GLenum FFArrayIdx = GL_COLOR_ARRAY;
   constexpr static auto FFSetupFunc = glColorPointer;
};

template<typename TV>
struct AttrTexcoord<TV, decltype((void)TV::texCoord, 0)>
: AttrBase<TV, decltype(TV::texCoord), &TV::texCoord,
AttrTexcoord<TV, decltype((void)TV::texCoord, 0)>>
{
   const static GLenum AttrGLType = GL_FLOAT;
   const static int ShaderIdx = CoreGLDevice::ATTR_TEXCOORD0;
   const static GLenum FFArrayIdx = GL_TEXTURE_COORD_ARRAY;
   constexpr static auto FFSetupFunc = glTexCoordPointer;
};

}

#endif
