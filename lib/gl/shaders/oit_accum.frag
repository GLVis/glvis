R"(
// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
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

void main()
{
   fragmentClipPlane();
   vec4 color = fColor * texture2D(colorTex, vec2(fTexCoord));
   color = blinnPhong(fPosition, fNormal, color);
#ifdef USE_ALPHA
   color.a *= texture2D(alphaTex, vec2(fTexCoord)).a;
#else
   color.a *= texture2D(alphaTex, vec2(fTexCoord)).r;
#endif

   float alpha = clamp(color.a, 0.0, 1.0);
   float weight = clamp(pow(alpha, 4.0) * 1000.0 + 0.01, 0.01, 3000.0);
   gl_FragColor = vec4(color.rgb * alpha * weight, alpha * weight);
}
)"

