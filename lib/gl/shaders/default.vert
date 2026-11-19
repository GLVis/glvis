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

attribute vec3 vertex;
attribute vec2 textVertex;
attribute vec4 color;
attribute vec3 normal;
attribute vec2 texCoord0;

attribute float line_orientation;
attribute vec3 line_prev_vtx;
attribute vec3 line_next_vtx;

uniform bool containsText;

uniform bool expandLines;
uniform float lineWidth;
uniform float aspectRatio;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform mat4 textProjMatrix;
uniform mat3 normalMatrix;

uniform vec4 clipPlane;

varying vec3 fNormal;
varying vec3 fPosition;
varying vec4 fColor;
varying vec2 fTexCoord;

varying float fClipVal;
varying float fLineEdgeDist;

void setupClipPlane(in float dist)
{
   fClipVal = dist;
}

void main()
{
   vec4 pos = modelViewMatrix * vec4(vertex, 1.0);
   fPosition = pos.xyz;
   fNormal = normalize(normalMatrix * normal);
   fColor = color;
   fTexCoord = texCoord0.xy;
   fLineEdgeDist = 0.0; // Default: not a line
   setupClipPlane(dot(vec4(pos.xyz, 1.0), clipPlane));
   pos = projectionMatrix * pos;
   gl_Position = pos;
   if (containsText)
   {
      vec4 textOffset = textProjMatrix * vec4(textVertex, 0.0, 0.0);
      gl_Position += vec4((textOffset.xy * pos.w), -0.005, 0.0);
   }
   // Line expansion (adapted from: https://github.com/mattdesl/webgl-lines)
   if (expandLines)
   {
      fLineEdgeDist = line_orientation; // ±1 at edges, 0 at center

      // Transform to screen space
      mat4 mvp = projectionMatrix * modelViewMatrix;
      vec4 prev_clip = mvp * vec4(line_prev_vtx, 1.0);
      vec4 next_clip = mvp * vec4(line_next_vtx, 1.0);

      vec2 curr_scrn = pos.xy / pos.w;
      vec2 prev_scrn = prev_clip.xy / prev_clip.w;
      vec2 next_scrn = next_clip.xy / next_clip.w;

      curr_scrn.x *= aspectRatio;
      prev_scrn.x *= aspectRatio;
      next_scrn.x *= aspectRatio;

      float width = lineWidth;

      // Compute line segment direction
      float dist_to_prev = distance(vertex, line_prev_vtx);
      float dist_to_next = distance(vertex, line_next_vtx);

      vec2 dir;
      vec2 tangent_offset = vec2(0.0);

      if (dist_to_prev < 0.0001)
      {
         // Start cap (rounded)
         dir = normalize(next_scrn - curr_scrn);
         tangent_offset = -dir * lineWidth * 0.5;
      }
      else if (dist_to_next < 0.0001)
      {
         // End cap (rounded)
         dir = normalize(curr_scrn - prev_scrn);
         tangent_offset = dir * lineWidth * 0.5;
      }
      else
      {
         // Line join
         vec2 dirA = normalize(curr_scrn - prev_scrn);
         vec2 dirB = normalize(next_scrn - curr_scrn);

         float dot_dirs = dot(dirA, dirB);
         if (dot_dirs > 0.99)
         {
            // Nearly straight line
            dir = normalize(dirA + dirB);
         }
         else
         {
            // Miter join with bevel fallback for sharp angles
            vec2 tangent = normalize(dirA + dirB);
            vec2 miter = vec2(-tangent.y, tangent.x);
            vec2 normal = vec2(-dirA.y, dirA.x);
            float miter_length = 1.0 / max(dot(miter, normal), 0.001);

            if (miter_length > 2.5)
            {
               // Bevel join for sharp angles (< ~75°)
               dir = dirA;
               width = lineWidth;
            }
            else
            {
               dir = tangent;
               width = lineWidth * miter_length;
            }
         }
      }

      vec2 line_normal = vec2(-dir.y, dir.x);
      line_normal.x /= aspectRatio;
      tangent_offset.x /= aspectRatio;
      vec4 offset = vec4(line_normal * line_orientation * width * pos.w / 2.0 + tangent_offset * pos.w, 0.0, 0.0);
      gl_Position += offset;
   }
})"
