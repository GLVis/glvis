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

attribute vec3 vertex;
attribute vec2 textVertex;
attribute vec4 color;
attribute vec3 normal;
attribute vec2 texCoord0;

uniform bool containsText;

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

void setupClipPlane(in float dist) {
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
    if (containsText) {
	vec4 textOffset = textProjMatrix * vec4(textVertex, 0.0, 0.0);
	gl_Position += vec4((textOffset.xy * pos.w), -0.005, 0.0);
    }
})"
