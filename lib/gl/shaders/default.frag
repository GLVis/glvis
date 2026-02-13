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
varying float fLineEdgeDist;
uniform bool useLineAA;
uniform float lineWidth;

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

   // Fragment shader antialiasing for lines (combined with MSAA when enabled)
   if (useLineAA && abs(fLineEdgeDist) > 0.001)
   {
      float dist = abs(fLineEdgeDist);
      float delta = fwidth(dist);

      // Edge transition with asymmetric smoothing
      float edge0 = 1.0 - delta * 1.6;
      float edge1 = 1.0 + delta * 0.6;
      float t = clamp((dist - edge0) / (edge1 - edge0), 0.0, 1.0);

      // Smootherstep interpolation (C2 continuous)
      float smootherT = t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
      float alpha = 1.0 - smootherT;
      alpha = alpha * (1.0 + alpha * 0.03);

      // Dual-ring multi-sampling (8 inner + 4 outer samples)
      float dx = dFdx(fLineEdgeDist);
      float dy = dFdy(fLineEdgeDist);

      const float r1 = 0.3535533905932738; // Inner ring: 1/sqrt(8)
      float s1 = abs(fLineEdgeDist + dx * r1);
      float s2 = abs(fLineEdgeDist - dx * r1);
      float s3 = abs(fLineEdgeDist + dy * r1);
      float s4 = abs(fLineEdgeDist - dy * r1);
      float s5 = abs(fLineEdgeDist + (dx + dy) * r1 * 0.707);
      float s6 = abs(fLineEdgeDist + (dx - dy) * r1 * 0.707);
      float s7 = abs(fLineEdgeDist - (dx + dy) * r1 * 0.707);
      float s8 = abs(fLineEdgeDist - (dx - dy) * r1 * 0.707);

      const float r2 = 0.5; // Outer ring
      float s9 = abs(fLineEdgeDist + (dx + dy) * r2);
      float s10 = abs(fLineEdgeDist + (dx - dy) * r2);
      float s11 = abs(fLineEdgeDist - (dx + dy) * r2);
      float s12 = abs(fLineEdgeDist - (dx - dy) * r2);

      // Apply smootherstep to all 12 samples
      float t1 = clamp((s1 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t2 = clamp((s2 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t3 = clamp((s3 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t4 = clamp((s4 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t5 = clamp((s5 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t6 = clamp((s6 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t7 = clamp((s7 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t8 = clamp((s8 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t9 = clamp((s9 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t10 = clamp((s10 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t11 = clamp((s11 - edge0) / (edge1 - edge0), 0.0, 1.0);
      float t12 = clamp((s12 - edge0) / (edge1 - edge0), 0.0, 1.0);

      float a1 = 1.0 - t1 * t1 * t1 * (t1 * (t1 * 6.0 - 15.0) + 10.0);
      float a2 = 1.0 - t2 * t2 * t2 * (t2 * (t2 * 6.0 - 15.0) + 10.0);
      float a3 = 1.0 - t3 * t3 * t3 * (t3 * (t3 * 6.0 - 15.0) + 10.0);
      float a4 = 1.0 - t4 * t4 * t4 * (t4 * (t4 * 6.0 - 15.0) + 10.0);
      float a5 = 1.0 - t5 * t5 * t5 * (t5 * (t5 * 6.0 - 15.0) + 10.0);
      float a6 = 1.0 - t6 * t6 * t6 * (t6 * (t6 * 6.0 - 15.0) + 10.0);
      float a7 = 1.0 - t7 * t7 * t7 * (t7 * (t7 * 6.0 - 15.0) + 10.0);
      float a8 = 1.0 - t8 * t8 * t8 * (t8 * (t8 * 6.0 - 15.0) + 10.0);
      float a9 = 1.0 - t9 * t9 * t9 * (t9 * (t9 * 6.0 - 15.0) + 10.0);
      float a10 = 1.0 - t10 * t10 * t10 * (t10 * (t10 * 6.0 - 15.0) + 10.0);
      float a11 = 1.0 - t11 * t11 * t11 * (t11 * (t11 * 6.0 - 15.0) + 10.0);
      float a12 = 1.0 - t12 * t12 * t12 * (t12 * (t12 * 6.0 - 15.0) + 10.0);

      // Weighted average: center(3x) + inner_ring(8×1x) + outer_ring(4×0.6x)
      alpha = (alpha * 3.0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + (a9 + a10 + a11 + a12) * 0.6) / 13.4;

      // Contrast-adaptive sharpening based on sample variance
      float variance = 0.0;
      variance += (a1 - alpha) * (a1 - alpha);
      variance += (a2 - alpha) * (a2 - alpha);
      variance += (a3 - alpha) * (a3 - alpha);
      variance += (a4 - alpha) * (a4 - alpha);
      variance += (a5 - alpha) * (a5 - alpha);
      variance += (a6 - alpha) * (a6 - alpha);
      variance += (a7 - alpha) * (a7 - alpha);
      variance += (a8 - alpha) * (a8 - alpha);
      variance += (a9 - alpha) * (a9 - alpha);
      variance += (a10 - alpha) * (a10 - alpha);
      variance += (a11 - alpha) * (a11 - alpha);
      variance += (a12 - alpha) * (a12 - alpha);
      variance /= 12.0;

      float sharpness = smoothstep(0.0, 0.10, variance);
      float sharpenAmount = sharpness * 0.25;
      float deviation = alpha - 0.5;
      alpha = clamp(alpha + deviation * sharpenAmount, 0.0, 1.0);

      // Subtle outer feather for smooth edges
      if (dist > 0.88 && dist < 1.25)
      {
         float featherDist = (dist - 0.88) / 0.37;
         float featherAlpha = exp(-featherDist * featherDist * 3.5);
         alpha = alpha * 0.88 + featherAlpha * 0.12;
      }

      // Sub-pixel coverage correction
      float pixelWidth = lineWidth * 1000.0;
      if (pixelWidth < 1.0)
      {
         alpha *= mix(pixelWidth, 1.0, pixelWidth * pixelWidth);
      }

      // Gamma correction for sRGB
      alpha = pow(clamp(alpha, 0.0, 1.0), 1.0 / 2.2);

      // Dual-frequency dithering to reduce banding
      vec2 screenPos = gl_FragCoord.xy;
      float noise1 = fract(52.9829189 * fract(0.06711056 * screenPos.x + 0.00583715 * screenPos.y));
      float noise2 = fract(31.5491234 * fract(0.09127341 * screenPos.x + 0.04283951 * screenPos.y));
      float noise = noise1 * 0.6 + noise2 * 0.4;
      alpha = clamp(alpha + (noise - 0.5) * 0.012, 0.0, 1.0);

      color.a *= alpha;
   }

   gl_FragColor = color;
}
)"
