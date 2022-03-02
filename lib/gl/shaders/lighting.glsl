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

struct Material {
   vec4 specular;
   float shininess;
};

struct PointLight {
   vec4 position;
   vec4 diffuse;
   vec4 specular;
};

uniform int num_lights;
uniform PointLight lights[3];
uniform vec4 g_ambient;

uniform Material material;

vec4 blinnPhong(in vec3 pos, in vec3 norm, in vec4 color)
{
   if (num_lights == 0)
   {
      return color;
   }
   // approximation of gl_FrontFacing from:
   // https://stackoverflow.com/questions/24375171
   if (dot(norm, pos) >= 0.0)
   {
      // invert normal direction for back-facing polygons
      norm *= -1.0;
   }
   // calculate ambient
   vec4 lit_color = g_ambient * color;
   for (int i = 0; i < 3; i++)
   {
      if (i >= num_lights)
      {
         break;
      }
      vec3 light_dir;
      if (lights[i].position.w == 0.0)
      {
         // directional light - no attenuation
         light_dir = normalize(lights[i].position.xyz);
      }
      else
      {
         light_dir = normalize(lights[i].position.xyz - pos);
      }

      // calculate diffuse
      lit_color += lights[i].diffuse * color * max(dot(norm, light_dir), 0.0);

      // calculate specular
      vec3 half_v = normalize(light_dir - pos);
      float specular_factor = max(dot(half_v, norm), 0.0);
      lit_color += lights[i].specular * material.specular * pow(specular_factor, material.shininess);
   }
   return lit_color;
}
)"
