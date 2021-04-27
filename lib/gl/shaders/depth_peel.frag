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

uniform sampler2D alphaTex;
uniform sampler2D colorTex;
uniform sampler2D lastDepthTex;
uniform sampler2D lastFrontColorTex;

varying vec3 fNormal;
varying vec3 fPosition;
varying vec4 fColor;
varying vec2 fTexCoord;

uniform bool useClipPlane;
varying float fClipVal;

void fragmentClipPlane()
{
   if (useClipPlane && fClipVal < 0.0)
   {
      discard;
   }
}

#define MAX_DEPTH 1.0

// location 0: depth
// location 1: front color
// location 2: back color

// adapted from "Order Independent Transparency with Dual Depth Peeling":
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.193.3485&rep=rep1&type=pdf
// and https://github.com/tsherif/webgl2examples/blob/master/oit-dual-depth-peeling.html

void main()
{
   fragmentClipPlane();

   ivec2 iFragCoord = ivec2(gl_FragCoord.xy);
   vec2 lastDepths = texelFetch(lastDepthTex, iFragCoord, 0).xy;
   vec4 lastFrontColor = texelFetch(lastFrontColorTex, iFragCoord, 0);
   float nearestDepth = -lastDepths.x;
   float farthestDepth = lastDepths.y;
   float thisDepth = gl_FragCoord.z;
   float alphaMultiplier = 1.0 - lastFrontColor.a;

   gl_FragData[1] = lastFrontColor;
   gl_FragData[2] = vec4(0.0);

   if (thisDepth < nearestDepth || thisDepth > farthestDepth)
   {
       gl_FragData[0].xy = vec2(-MAX_DEPTH);
       return;
   }

   if (thisDepth > nearestDepth && thisDepth < farthestDepth)
   {
       gl_FragData[0].xy = vec2(-thisDepth, thisDepth);
       return;
   }

   gl_FragData[0].xy = vec2(-MAX_DEPTH);
   vec4 color = fColor * texture2D(colorTex, vec2(fTexCoord));
   color = blinnPhong(fPosition, fNormal, color);
#ifdef USE_ALPHA
   color.a *= texture2D(alphaTex, vec2(fTexCoord)).a;
#else
   color.a *= texture2D(alphaTex, vec2(fTexCoord)).r;
#endif
    
    if (thisDepth == nearestDepth)
    {
        gl_FragData[1].rgb += color.rgb * color.a * alphaMultiplier;
        gl_FragData[1].a = 1.0 - alphaMultiplier * (1.0 - color.a);
    }
    else
    {
        gl_FragData[2] += color;
    }
}
)"
