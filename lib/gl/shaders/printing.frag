R"(
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
