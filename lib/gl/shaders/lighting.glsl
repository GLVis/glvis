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

vec4 blinnPhong(in vec3 pos, in vec3 norm, in vec4 color) {
    if (num_lights == 0) {
	return color;
    }
    // approximation of gl_FrontFacing from:
    // https://stackoverflow.com/questions/24375171
    if (dot(norm, pos) >= 0.0) {
        // invert normal direction for back-facing polygons
        norm *= -1.0;
    }
    // calculate ambient
    vec4 lit_color = g_ambient * color;
    for (int i = 0; i < 3; i++)
    {
	if (i >= num_lights)
	    break;
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
