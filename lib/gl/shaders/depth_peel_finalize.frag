R"(
// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

uniform sampler2D lastFrontColor;
uniform sampler2D lastBackColor;

// adapted from "Order Independent Transparency with Dual Depth Peeling":
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.193.3485&rep=rep1&type=pdf
// and https://github.com/tsherif/webgl2examples/blob/master/oit-dual-depth-peeling.html

void main()
{
    ivec2 iQuadCoord = ivec2(gl_FragCoord.xy);
    vec4 frontColor = texelFetch(lastFrontColor, iQuadCoord, 0);
    vec4 backColor  = texelFetch(lastBackColor, iQuadCoord, 0);
    float alphaMultiplier = 1.0 - frontColor.a;

    gl_FragColor.rgb = frontColor.rgb + alphaMultiplier * backColor.rgb;
    gl_FragColor.a = frontColor.a + backColor.a;
}
)"
