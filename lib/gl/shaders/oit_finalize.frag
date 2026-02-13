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

uniform sampler2D sceneTex;
uniform sampler2D accumTex;
uniform sampler2D revealTex;

varying vec2 vUv;

void main()
{
   vec3 scene = texture2D(sceneTex, vUv).rgb;
   vec4 accum = texture2D(accumTex, vUv);
   float reveal = texture2D(revealTex, vUv).r;

   float transAlpha = 1.0 - reveal;
   vec3 transColor = (accum.a > 1e-5) ? (accum.rgb / accum.a) : vec3(0.0);
   vec3 outColor = transColor * transAlpha + scene * reveal;
   gl_FragColor = vec4(outColor, 1.0);
}
)"

