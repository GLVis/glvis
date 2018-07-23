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
#ifndef MATERIAL_HPP
#define MATERIAL_HPP

struct Material
{
    float ambient[4];
    float diffuse[4];
    float specular[4];
    float shininess;
};

struct Light
{
    float position[4];
    float diffuse[4];
    float specular[4];
};

void Set_Material();

void Set_Light();

int Next_Material_And_Light();

void Set_Material_And_Light(int,int);

float Set_Black_Material();

void Set_Background();

void Toggle_Background();

void Set_Transparency();

void Remove_Transparency();

int  Get_AntiAliasing();
void Set_AntiAliasing();
void Remove_AntiAliasing();

double Get_LineWidth();
void Set_LineWidth(double);
double Get_MS_LineWidth();
void Set_MS_LineWidth(double);

#endif
