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
#include <math.h>
#include "openglvis.hpp"

VisualizationScene :: VisualizationScene ()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity ();
   glGetDoublev (GL_MODELVIEW_MATRIX, translmat);
   glRotatef(-60.0, 1.0f, 0.0f, 0.0f);
   glRotatef(-40.0, 0.0f, 0.0f, 1.0f);
   glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
   xscale = yscale = zscale = 1;
   spinning = print = movie = 0;
   OrthogonalProjection = 0;
   ViewAngle = 45;
   ViewScale = 1;
   ViewCenterX = 0.0;
   ViewCenterY = 0.0;
}

VisualizationScene :: ~VisualizationScene () {;}

void VisualizationScene::Rotate(double anglex, double angley)
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity ();
   glRotatef(anglex, 0.0f, 1.0f, 0.0f);
   glRotatef(angley, 1.0f, 0.0f, 0.0f);
   glMultMatrixd (rotmat);
   glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
}

void VisualizationScene::Translate(double _x, double _y, double _z)
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity ();
   glTranslatef (_x, -_y, _z);
   glMultMatrixd (translmat);
   glGetDoublev (GL_MODELVIEW_MATRIX, translmat);
}

void VisualizationScene::Scale(double s)
{
   Scale (s, s, s);
}

void VisualizationScene::Scale(double s1, double s2, double s3){
   /*
     glMatrixMode (GL_MODELVIEW);
     glLoadIdentity ();
     glScaled (s1, s2, s3);
     glMultMatrixd (rotmat);
     glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
   */
   xscale *= s1;
   yscale *= s2;
   zscale *= s3;
}

void VisualizationScene::Prepare(){}

void VisualizationScene::CenterObject(){}

void VisualizationScene::CenterObject2D(){}

void VisualizationScene::PrepareAxes(){}

void VisualizationScene::Draw(){}

void VisualizationScene::SetView(double theta, double phi)
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity ();
   glGetDoublev (GL_MODELVIEW_MATRIX, translmat);

   glRotatef(-theta, 1.0f, 0.0f, 0.0f);
   glRotatef(-phi, 0.0f, 0.0f, 1.0f);
   glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
}

void VisualizationScene::Zoom(double factor)
{
   if (OrthogonalProjection)
   {
      ViewScale *= factor;
   }
   else
   {
      double va = ViewAngle * ( M_PI / 360.0 );
      ViewAngle = atan( tan( va ) / factor ) * (360.0 / M_PI);
   }
}
