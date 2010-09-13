// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.googlecode.com.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <string>

int Num_Materials = 5;
int Current_Material = 3;

void Set_Material()
{
   switch(Current_Material)
   {
   case 0:
   {
      GLfloat mdiff[] = { 0.8, 0.8, 0.8, 1.0 };
      GLfloat mspec[] = { 1.0, 1.0, 1.0, 1.0 };
      GLfloat mshin[] = { 100.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);

      // if (glIsEnabled (GL_BLEND))
      //    glDisable(GL_BLEND);
   }
   break;

   case 1:
   {
      GLfloat mdiff[] = { 0.7, 0.7, 0.7, 1.0 };
      GLfloat mambi[] = { 0.3, 0.3, 0.3, 1.0 };
      GLfloat mspec[] = { 0.8, 0.8, 0.8, 1.0 };
      GLfloat mshin[] = { 20.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);

      // if (glIsEnabled (GL_BLEND))
      //    glDisable(GL_BLEND);
   }
   break;

   case 2:
   {
      GLfloat mdiff[] = { 1.0, 1.0, 1.0, 1.0 };
      GLfloat mambi[] = { 0.3, 0.3, 0.3, 1.0 };
      GLfloat mspec[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mshin[] = { 0.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);

      // if (glIsEnabled (GL_BLEND))
      //    glDisable(GL_BLEND);
   }
   break;

   case 3:
   {
      GLfloat mdiff[] = { 0.75164, 0.60648, 0.22648, 1.0 };
      GLfloat mambi[] = { 0.24725, 0.1995, 0.0745, 1.0 };
      GLfloat mspec[] = { 0.628281, 0.555802, 0.366065, 1.0 };
      GLfloat mshin[] = { 51.2 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);

      // if (glIsEnabled (GL_BLEND))
      //    glDisable(GL_BLEND);
   }
   break;

   case 4:
   {
      GLfloat mdiff[] = { 0.8, 0.8, 0.8, 1.0 };
      GLfloat mambi[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mspec[] = { 0.1, 0.1, 0.1, 1.0 };
      GLfloat mshin[] = { 1.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);

      // glEnable(GL_BLEND);
      // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }
   break;

   }

   GLfloat memis[] = { 0.0, 0.0, 0.0, 1.0 };
   glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, memis);
}

void Set_Light()
{
   glLoadIdentity();

   switch(Current_Material)
   {
   case 0:
   {
      GLfloat light[] = { 1.0, 1.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);

      GLfloat ambvals[]  = { 0.1f, 0.1f, 0.1f, 1.0f };
      GLfloat diffvals[] = { 0.9f, 0.9f, 0.9f, 1.0f };
      GLfloat specvals[] = { 0.8f, 0.8f, 0.8f, 1.0f };
      glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
      GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

      if (glIsEnabled (GL_LIGHT1))
         glDisable(GL_LIGHT1);
      if (glIsEnabled (GL_LIGHT2))
         glDisable(GL_LIGHT2);
   }
   break;

   case 1:
   {
      GLfloat light[] = { 0.5, 0.5, 1.0, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);

      GLfloat ambvals[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
      GLfloat diffvals[] = { 0.5f, 0.5f, 0.5f, 1.0f };
      GLfloat specvals[] = { 1.0f, 1.0f, 1.0f, 1.0f };
      glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
      GLfloat lmodel_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

      if (glIsEnabled (GL_LIGHT1))
         glDisable(GL_LIGHT1);
      if (glIsEnabled (GL_LIGHT2))
         glDisable(GL_LIGHT2);
   }
   break;

   case 2:
   {
      GLfloat light[] = { 0.0, 0.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);

      GLfloat ambvals[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
      GLfloat diffvals[] = { 0.5f, 0.5f, 0.5f, 1.0f };
      GLfloat specvals[] = { 0.0f, 0.0f, 0.0f, 1.0f };
      glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
      GLfloat lmodel_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

      if (glIsEnabled (GL_LIGHT1))
         glDisable(GL_LIGHT1);
      if (glIsEnabled (GL_LIGHT2))
         glDisable(GL_LIGHT2);
   }
   break;

   case 3:
   {
      GLfloat light[] = { 0.0, 0.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);

      GLfloat ambvals[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
      GLfloat diffvals[] = { 0.7f, 0.7f, 0.7f, 1.0f };
      GLfloat specvals[] = { 0.6f, 0.6f, 0.6f, 1.0f };
      glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
      GLfloat lmodel_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

      if (glIsEnabled (GL_LIGHT1))
         glDisable(GL_LIGHT1);
      if (glIsEnabled (GL_LIGHT2))
         glDisable(GL_LIGHT2);
   }
   break;

   case 4:
   {
      GLfloat light[] = { 1.0, 0.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);
      GLfloat light1[] = { 1.0, 1.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT1, GL_POSITION, light1);
      GLfloat light2[] = { 0.0, 1.0, 1.0, 0.0 };
      glLightfv(GL_LIGHT2, GL_POSITION, light2);

      GLfloat specvals[] = { 0.3f, 0.3f, 0.3f, 1.0f };

      GLfloat diffvals[] = { 0.4f, 0.0f, 0.0f, 0.0f };
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);

      GLfloat diffvals1[] = { 0.0f, 0.4f, 0.0f, 0.0f };
      glLightfv(GL_LIGHT1, GL_DIFFUSE, diffvals1);
      glLightfv(GL_LIGHT1, GL_SPECULAR, specvals);

      GLfloat diffvals2[] = { 0.0f, 0.0f, 0.4f, 0.0f };
      glLightfv(GL_LIGHT2, GL_DIFFUSE, diffvals2);
      glLightfv(GL_LIGHT2, GL_SPECULAR, specvals);

      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      glEnable(GL_LIGHT2);
   }
   break;

   }

   // Enable textures with specular color
   glLightModeli (GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
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

void Set_Black_Material()
{
   switch (Background)
   {
   case 0:
   {
      // glColor3f (0.75, 0.75, 0.75);
      glColor3f (1., 1., 1.);
      GLfloat mdiff[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mambi[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mspec[] = { 0.0, 0.0, 0.0, 1.0 };
      // GLfloat memis[] = { 0.75, 0.75, 0.75, 1.0 };
      GLfloat memis[] = { 0., 0., 0., 1.0 };
      GLfloat mshin[] = { 0.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, memis);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);
   }
   break;

   case 1:
   {
      glColor3f (0, 0, 0);
      GLfloat mdiff[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mambi[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mspec[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat memis[] = { 0.0, 0.0, 0.0, 1.0 };
      GLfloat mshin[] = { 0.0 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mdiff);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mambi);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mspec);
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, memis);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mshin);
      break;
   }
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
   glEnable (GL_BLEND);
   glDepthMask (GL_FALSE);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Remove_Transparency()
{
   glDisable (GL_BLEND);
   glDepthMask (GL_TRUE);
}

void Set_AntiAliasing()
{
#ifdef GLVIS_MULTISAMPLE
   glEnable(GL_MULTISAMPLE);

#ifdef GL_MULTISAMPLE_FILTER_HINT_NV
   std::string s = (const char *)glGetString(GL_EXTENSIONS);
   if (s.find("GL_NV_multisample_filter_hint") != std::string::npos)
      glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
#endif
#endif

   glEnable(GL_BLEND);
   // glDisable(GL_DEPTH_TEST);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   // glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);

   glLineWidth(1.4);

   glEnable(GL_POLYGON_SMOOTH);
   glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

   glEnable(GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

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
}

void Remove_AntiAliasing()
{
   glDisable(GL_POLYGON_SMOOTH);
   glDisable(GL_LINE_SMOOTH);
   glDisable(GL_POINT_SMOOTH);

   glLineWidth(1.);

   // glEnable(GL_DEPTH_TEST);
   glDisable(GL_BLEND);

#ifdef GLVIS_MULTISAMPLE
#ifdef GL_MULTISAMPLE_FILTER_HINT_NV
   std::string s = (const char *)glGetString(GL_EXTENSIONS);
   if (s.find("GL_NV_multisample_filter_hint") != std::string::npos)
      glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_FASTEST);
#endif

   glDisable(GL_MULTISAMPLE);
#endif

   // glDisable(GL_FOG);
}
