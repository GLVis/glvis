// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "material.hpp"
#include <array>
Material materials[5] =
{
    {
        { 0.8, 0.8, 0.8, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        { 1.0, 1.0, 1.0, 1.0 },
        100
    },
    {
        { 0.3, 0.3, 0.3, 1.0 },
        { 0.7, 0.7, 0.7, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        20
    },
    {
        { 0.3, 0.3, 0.3, 1.0 },
        { 1.0, 1.0, 1.0, 1.0 },
        { 0.0, 0.0, 0.0, 1.0 },
        0
    },
    {
        { 0.24725, 0.1995, 0.0745, 1.0 },
        { 0.75164, 0.60648, 0.22648, 1.0 },
        { 0.628281, 0.555802, 0.366065, 1.0 }, 
        51.2
    },
    {
        { 0.0, 0.0, 0.0, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        { 0.1, 0.1, 0.1, 1.0 },
        1.0
    }
};

Light lights[] =
{
    { { 1.0, 1.0, 1.0, 0.0 }, { 0.9, 0.9, 0.9, 1.0 }, { 0.8, 0.8, 0.8, 1.0 } },
    { { 0.5, 0.5, 1.0, 0.0 }, { 0.5, 0.5, 0.5, 1.0 }, { 1.0, 1.0, 1.0, 1.0 } },
    { { 0.0, 0.0, 1.0, 0.0 }, { 0.5, 0.5, 0.5, 1.0 }, { 0.0, 0.0, 0.0, 1.0 } },
    { { 0.0, 0.0, 1.0, 0.0 }, { 0.7, 0.7, 0.7, 1.0 }, { 0.6, 0.6, 0.6, 1.0 } }
};

std::array<float,4> amb_setting[] = 
{
    { 0.3, 0.3, 0.3, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 }
};

Light lights_4[] = 
{
    { { 1.0, 0.0, 1.0, 0.0 }, { 0.4, 0.0, 0.0, 1.0 }, { 0.3, 0.3, 0.3, 1.0 } },
    { { 1.0, 1.0, 1.0, 0.0 }, { 0.0, 0.4, 0.0, 1.0 }, { 0.3, 0.3, 0.3, 1.0 } },
    { { 0.0, 1.0, 1.0, 0.0 }, { 0.0, 0.0, 0.4, 1.0 }, { 0.3, 0.3, 0.3, 1.0 } }
};
