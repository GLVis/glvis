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

#include <array>

#include "material.hpp"

Material materials[5] =
{
   {
      { 0.8f, 0.8f, 0.8f, 1.0f },
      { 0.8f, 0.8f, 0.8f, 1.0f },
      { 1.0f, 1.0f, 1.0f, 1.0f },
      100.0f
   },
   {
      { 0.3f, 0.3f, 0.3f, 1.0f },
      { 0.7f, 0.7f, 0.7f, 1.0f },
      { 0.8f, 0.8f, 0.8f, 1.0f },
      20.0f
   },
   {
      { 0.3f, 0.3f, 0.3f, 1.0f },
      { 1.0f, 1.0f, 1.0f, 1.0f },
      { 0.0f, 0.0f, 0.0f, 1.0f },
      0.0f
   },
   {
      { 0.24725f, 0.1995f, 0.0745f, 1.0f },
      { 0.75164f, 0.60648f, 0.22648f, 1.0f },
      { 0.628281f, 0.555802f, 0.366065f, 1.0f },
      51.2f
   },
   {
      { 0.0f, 0.0f, 0.0f, 1.0f },
      { 0.8f, 0.8f, 0.8f, 1.0f },
      { 0.1f, 0.1f, 0.1f, 1.0f },
      1.0f
   }
};

Light lights[] =
{
   { { 1.0f, 1.0f, 1.0f, 0.0f }, { 0.9f, 0.9f, 0.9f, 1.0f }, { 0.8f, 0.8f, 0.8f, 1.0f } },
   { { 0.5f, 0.5f, 1.0f, 0.0f }, { 0.5f, 0.5f, 0.5f, 1.0f }, { 1.0f, 1.0f, 1.0f, 1.0f } },
   { { 0.0f, 0.0f, 1.0f, 0.0f }, { 0.5f, 0.5f, 0.5f, 1.0f }, { 0.0f, 0.0f, 0.0f, 1.0f } },
   { { 0.0f, 0.0f, 1.0f, 0.0f }, { 0.7f, 0.7f, 0.7f, 1.0f }, { 0.6f, 0.6f, 0.6f, 1.0f } }
};

std::array<float,4> amb_setting[] =
{
   { 0.3f, 0.3f, 0.3f, 1.0f },
   { 0.5f, 0.5f, 0.5f, 1.0f },
   { 0.5f, 0.5f, 0.5f, 1.0f },
   { 0.5f, 0.5f, 0.5f, 1.0f },
   { 0.5f, 0.5f, 0.5f, 1.0f }
};

Light lights_4[] =
{
   { { 1.0f, 0.0f, 1.0f, 0.0f }, { 0.4f, 0.0f, 0.0f, 1.0f }, { 0.3f, 0.3f, 0.3f, 1.0f } },
   { { 1.0f, 1.0f, 1.0f, 0.0f }, { 0.0f, 0.4f, 0.0f, 1.0f }, { 0.3f, 0.3f, 0.3f, 1.0f } },
   { { 0.0f, 1.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 0.4f, 1.0f }, { 0.3f, 0.3f, 0.3f, 1.0f } }
};
