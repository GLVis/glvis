R"(
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

uniform bool containsText;
uniform bool useColorTex;

uniform sampler2D fontTex;
uniform sampler2D alphaTex;
uniform sampler2D colorTex;

varying vec4 fColor;
varying float fClipCoord;

uniform bool useClipPlane;

void main()
{
   gl_FragColor = fColor;
}
)"
