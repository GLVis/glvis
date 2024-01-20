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

#ifndef GLVIS_MATERIAL_HPP
#define GLVIS_MATERIAL_HPP
#include <array>

struct Material
{
   std::array<float, 4> ambient;
   std::array<float, 4> diffuse;
   std::array<float, 4> specular;
   float shininess;
};

struct Light
{
   std::array<float, 4> position;
   std::array<float, 4> diffuse;
   std::array<float, 4> specular;
};

#endif
