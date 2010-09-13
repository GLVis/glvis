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

#ifndef GLVIS_OPENGLVIS
#define GLVIS_OPENGLVIS

// Visualization header file

class VisualizationScene
{

protected:

   // How to scale the visualized object(s)
   double xscale, yscale, zscale;

public:
   VisualizationScene();
   virtual ~VisualizationScene();

   int spinning, OrthogonalProjection, print, movie;
   double ViewAngle, ViewScale;
   double ViewCenterX, ViewCenterY;

   /// Bounding box.
   double x[2], y[2], z[2];

   double rotmat[16];
   double translmat[16];

   virtual void Prepare();
   virtual void CenterObject();
   virtual void CenterObject2D();
   virtual void PrepareAxes();
   virtual void Draw();

   void SetView(double theta, double phi);
   void Zoom(double factor);

   void Rotate (double anglex, double angley);
   void Translate (double x, double y, double z=0);
   void Scale (double s);
   void Scale (double s1, double s2, double s3);

   /// This is set by SetVisualizationScene
   int view;
};

#endif
