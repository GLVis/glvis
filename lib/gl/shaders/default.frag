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

uniform sampler2D alphaTex;
uniform sampler2D colorTex;

varying vec3 fNormal;
varying vec3 fPosition;
varying vec4 fColor;
varying vec2 fTexCoord;

uniform bool useClipPlane;
varying float fClipVal;

void fragmentClipPlane() {
    if (useClipPlane && fClipVal < 0.0) {
	discard;
    }
}

void main()
{
    fragmentClipPlane();
    vec4 color = fColor * texture2D(colorTex, vec2(fTexCoord));
    color = blinnPhong(fPosition, fNormal, color);
#ifdef GL_ES
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).a;
#else
    color.a *= texture2D(alphaTex, vec2(fTexCoord)).r;
#endif
    gl_FragColor = color;
}
)"
