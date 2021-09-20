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
   setupClipPlane(dot(vec4(pos.xyz, 1.0), clipPlane));
   pos = projectionMatrix * pos;
   gl_Position = pos;
   if (containsText)
   {
      vec4 textOffset = textProjMatrix * vec4(textVertex, 0.0, 0.0);
      gl_Position += vec4((textOffset.xy * pos.w), -0.005, 0.0);
   }
   if (expandLines)
   {
      // Compute screen-space coordinates for current and previous segments
      mat4 mvp = projectionMatrix * modelViewMatrix;
      vec4 prev_clip = mvp * vec4(line_prev_vtx, 1.0);
      vec4 next_clip = mvp * vec4(line_next_vtx, 1.0);

      vec2 curr_scrn = pos.xy / pos.w;
      vec2 prev_scrn = prev_clip.xy / prev_clip.w;
      vec2 next_scrn = next_clip.xy / next_clip.w;

      // Correct for current aspect ratio
      curr_scrn.x *= aspectRatio;
      prev_scrn.x *= aspectRatio;
      next_scrn.x *= aspectRatio;

      float width = lineWidth;

      // Get direction and normal of line segment
      vec2 dir;
      if (vertex == line_prev_vtx)
      {
         dir = normalize(next_scrn - curr_scrn);
      }
      else if (vertex == line_next_vtx)
      {
         dir = normalize(curr_scrn - prev_scrn);
      }
      else
      {
         dir = normalize(curr_scrn - prev_scrn);
         vec2 perp = vec2(-dir.y, dir.x);

         dir += normalize(next_scrn - curr_scrn);
         dir = normalize(dir);
         vec2 miter = vec2(-dir.y, dir.x);
         width = lineWidth / dot(miter, perp);
      }
      vec2 line_normal = vec2(-dir.y, dir.x);
      vec4 offset = vec4(line_normal * line_orientation * width / 2.0, 0.0, 0.0);
      gl_Position += offset;
   }
})"
