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

#include "platform_gl.hpp"
#include "material.hpp"
#include <string>
#ifdef GLVIS_DEBUG
#include <iostream>
#endif
#include "aux_vis.hpp"

#include "gl2ps.h"

extern int visualize;
extern int GetMultisample();

int Num_Materials = 5;
int Current_Material = 3;
int AntiAliasing = 0;
double NM_LineWidth = 1.0;
#ifdef GLVIS_MS_LINEWIDTH
double MS_LineWidth = GLVIS_MS_LINEWIDTH;
#else
double MS_LineWidth = 1.4;
#endif

Material materials[5] =
{
    {
        { 0.8, 0.8, 0.8, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        { 1.0, 1.0, 1.0, 1.0 },
        100
    },
    {
        { 0.3, 0.3, 0.3, 1.0 },
        { 0.7, 0.7, 0.7, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        20
    },
    {
        { 0.3, 0.3, 0.3, 1.0 },
        { 1.0, 1.0, 1.0, 1.0 },
        { 0.0, 0.0, 0.0, 1.0 },
        0
    },
    {
        { 0.24725, 0.1995, 0.0745, 1.0 },
        { 0.75164, 0.60648, 0.22648, 1.0 },
        { 0.628281, 0.555802, 0.366065, 1.0 }, 
        51.2
    },
    {
        { 0.0, 0.0, 0.0, 1.0 },
        { 0.8, 0.8, 0.8, 1.0 },
        { 0.1, 0.1, 0.1, 1.0 },
        1.0
    }
};

Light lights[] =
{
    { { 1.0, 1.0, 1.0, 0.0 }, { 0.9, 0.9, 0.9, 1.0 }, { 0.8, 0.8, 0.8, 1.0 } },
    { { 0.5, 0.5, 1.0, 0.0 }, { 0.5, 0.5, 0.5, 1.0 }, { 1.0, 1.0, 1.0, 1.0 } },
    { { 0.0, 0.0, 1.0, 0.0 }, { 0.5, 0.5, 0.5, 1.0 }, { 0.0, 0.0, 0.0, 1.0 } },
    { { 0.0, 0.0, 1.0, 0.0 }, { 0.7, 0.7, 0.7, 1.0 }, { 0.6, 0.6, 0.6, 1.0 } }
};

float amb_setting[][4] = 
{
    { 0.3, 0.3, 0.3, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.5, 0.5, 0.5, 1.0 },
    { 0.2, 0.2, 0.2, 1.0 }
};

Light lights_4[] = 
{
    { { 1.0, 0.0, 1.0, 0.0 }, { 0.3, 0.3, 0.3, 1.0 }, { 0.4, 0.0, 0.0, 0.0 } },
    { { 1.0, 1.0, 1.0, 0.0 }, { 0.3, 0.3, 0.3, 1.0 }, { 0.0, 0.4, 0.0, 0.0 } },
    { { 0.0, 1.0, 1.0, 0.0 }, { 0.3, 0.3, 0.3, 1.0 }, { 0.0, 0.0, 0.4, 0.0 } }
};

void Set_Material()
{
    GetGlState()->setMaterial(materials[Current_Material]);
}

void Set_Light()
{
    GetGlState()->setGlobalAmbLight(amb_setting[Current_Material]);
    if (Current_Material == 4) {
        for (int i = 0; i < 3; i++) {
            GetGlState()->setLight(i, lights_4[i]);
        }
        GetGlState()->setNumLights(3);
    } else {
        GetGlState()->setLight(0, lights[Current_Material]);
        GetGlState()->setNumLights(1);
    }
}

int Next_Material_And_Light()
{
   Current_Material = (Current_Material + 1) % Num_Materials;
   Set_Material();
   Set_Light();

   return Current_Material+1;
}

void Set_Material_And_Light (int m, int l)
{
   Current_Material = l;
   Set_Light();
   Current_Material = m;
   Set_Material();
}

int Background = 1;

float Set_Black_Material()
{
    Material m = {
        { 0.0, 0.0, 0.0, 1.0 },
        { 0.0, 0.0, 0.0, 1.0 },
        { 0.0, 0.0, 0.0, 1.0 },
        0.0
    };
   switch (Background)
   {
      case 0:
         // glColor3f (0.75, 0.75, 0.75);
         glColor3f (1., 1., 1.);
         GetGlState()->setMaterial(m);
         return 1;

      case 1:
         glColor3f (0, 0, 0);
         GetGlState()->setMaterial(m);
         return 0;
   }
}

void Set_Background()
{
   switch (Background)
   {
      case 0:
         glClearColor (0.0, 0.0, 0.0, 1.0);
         break;
      case 1:
         glClearColor (1.0, 1.0, 1.0, 1.0);
         break;
   }
}

void Toggle_Background()
{
   Background = 1 - Background;
}

void Set_Transparency()
{
   if (AntiAliasing == 0)
   {
#ifndef GLVIS_OGL3
      gl2psEnable(GL2PS_BLEND);
#endif
      GetGlState()->enableBlend();
   }
   glDepthMask(GL_FALSE);
}

void Remove_Transparency()
{
   if (AntiAliasing == 0)
   {
      GetGlState()->disableBlend(); 
#ifndef GLVIS_OGL3
      gl2psDisable(GL2PS_BLEND);
#endif
   }
   glDepthMask(GL_TRUE);
}

int Get_AntiAliasing()
{
   return AntiAliasing;
}

void Set_AntiAliasing()
{
   if (GetMultisample() > 0)
   {
      glEnable(GL_MULTISAMPLE);

#ifdef GL_MULTISAMPLE_FILTER_HINT_NV
      std::string s = (const char *)glGetString(GL_EXTENSIONS);
      if (s.find("GL_NV_multisample_filter_hint") != std::string::npos)
      {
         glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
      }
#endif
   }

#ifndef GLVIS_OGL3
   gl2psEnable(GL2PS_BLEND);
#endif
   //glEnable(GL_BLEND);
   // glDisable(GL_DEPTH_TEST);
   //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   // glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
   GetGlState()->enableBlend();

   glLineWidth(MS_LineWidth);

   // In order for polygon smoothing to work nicely, the polygons need to be
   // sorted by depth. Therefore we leave polygon smoothing to multisampling.
   if (0)
   {
      glEnable(GL_POLYGON_SMOOTH);
      glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
   }

   std::string vendor = (const char *)glGetString(GL_VENDOR);
   if (vendor.find("ATI") == std::string::npos || GetMultisample() <= 0)
   {
      // In order for line smoothing to blend nicely with polygons, we draw
      // the lines after the polygons.
      glEnable(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   }
   else
   {
#ifdef GLVIS_DEBUG
      std::cout <<
                "Found 'ATI' in the GL vendor string:"
                " using line smoothing via multisampling." << std::endl;
#endif
   }

   glEnable(GL_POINT_SMOOTH);
   glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

   /* add fog
      glEnable(GL_AUTO_NORMAL);
      glEnable(GL_NORMALIZE);
      glEnable(GL_FOG);
      GLfloat fogColor[4] = {1.0, 1.0, 1.0, 1.0};
      GLint fogMode;
      fogMode = GL_EXP2;
      glFogi (GL_FOG_MODE, fogMode);
      glFogfv (GL_FOG_COLOR, fogColor);
      glFogf (GL_FOG_DENSITY, 2.0);
      glHint (GL_FOG_HINT, GL_NICEST);
   */

   AntiAliasing = 1;
}

void Remove_AntiAliasing()
{
   glDisable(GL_POLYGON_SMOOTH);
   glDisable(GL_LINE_SMOOTH);
   glDisable(GL_POINT_SMOOTH);

   glLineWidth(NM_LineWidth);

   // glEnable(GL_DEPTH_TEST);
   GetGlState()->disableBlend();
#ifndef GLVIS_OGL3
   gl2psDisable(GL2PS_BLEND);
#endif

   if (GetMultisample() > 0)
   {
#ifdef GL_MULTISAMPLE_FILTER_HINT_NV
      std::string s = (const char *)glGetString(GL_EXTENSIONS);
      if (s.find("GL_NV_multisample_filter_hint") != std::string::npos)
      {
         glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_FASTEST);
      }
#endif

      glDisable(GL_MULTISAMPLE);
   }

   // glDisable(GL_FOG);

   AntiAliasing = 0;
}

double Get_LineWidth()
{
   return NM_LineWidth;
}

void Set_LineWidth(double lw)
{
   NM_LineWidth = lw;
#ifdef GLVIS_DEBUG
   std::cout << "Normal LineWidth set to " << NM_LineWidth << std::endl;
#endif
   if (AntiAliasing == 0 && visualize)
   {
      glLineWidth(NM_LineWidth);
   }
}

double Get_MS_LineWidth()
{
   return MS_LineWidth;
}

void Set_MS_LineWidth(double lw)
{
   MS_LineWidth = lw;
#ifdef GLVIS_DEBUG
   std::cout << "Multisample LineWidth set to " << MS_LineWidth << std::endl;
#endif
   if (AntiAliasing == 1 && visualize)
   {
      glLineWidth(MS_LineWidth);
   }
}
