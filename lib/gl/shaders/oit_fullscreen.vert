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

varying vec2 vUv;

void main()
{
   vec2 pos;
   if (gl_VertexID == 0)
   {
      pos = vec2(-1.0, -1.0);
   }
   else if (gl_VertexID == 1)
   {
      pos = vec2(3.0, -1.0);
   }
   else
   {
      pos = vec2(-1.0, 3.0);
   }
   vUv = 0.5 * pos + 0.5;
   gl_Position = vec4(pos, 0.0, 1.0);
}
)"

