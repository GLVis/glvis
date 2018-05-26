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

#ifndef GLVIS_OPENGLVIS
#define GLVIS_OPENGLVIS

#include <cmath>

#include "sdl.hpp"

// Visualization header file

// Some inline functions

inline void LinearCombination(const double a, const double x[],
                              const double b, const double y[], double z[])
{
   z[0] = a*x[0] + b*y[0];
   z[1] = a*x[1] + b*y[1];
   z[2] = a*x[2] + b*y[2];
}

inline double InnerProd(const double a[], const double b[])
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void CrossProd(const double a[], const double b[], double cp[])
{
   cp[0] = a[1] * b[2] - a[2] * b[1];
   cp[1] = a[2] * b[0] - a[0] * b[2];
   cp[2] = a[0] * b[1] - a[1] * b[0];
}

inline int Normalize(double v[])
{
   double len = sqrt(InnerProd(v, v));
   if (len > 0.0)
   {
      len = 1.0 / len;
   }
   else
   {
      return 1;
   }
   for (int i = 0; i < 3; i++)
   {
      v[i] *= len;
   }
   return 0;
}

inline int ProjectVector(double v[], const double n[])
{
   // project 'v' on the plane with normal given by 'n' and  then normalize 'v'
   LinearCombination(InnerProd(n, n), v, -InnerProd(v, n), n, v);
   return Normalize(v);
}


class Camera
{
private:
   double eye[3], dir[3], up[3];
   double left[3];

   void MoveEye(double dist, const double dir_[])
   { LinearCombination(1.0, eye, dist, dir_, eye); }

public:
   Camera() { Reset(); }

   void Reset();
   void Set(const double cam[]);

   const double *GetEye() { return eye; }
   const double *GetDir() { return dir; }
   const double *GetUp() { return up; }
   const double *GetLeft() { CrossProd(up, dir, left); return left; }

   void MoveForwardBackward(double dist) { MoveEye(dist, dir); }
   void MoveLeftRight(double dist) { MoveEye(dist, GetLeft()); }
   void MoveUpDown(double dist) { MoveEye(dist, up); }

   void TiltLeftRight(double angle);

   void TurnLeftRight(double angle);
   void TurnUpDown(double angle);

   void GLMultRotMatrix();
   void GLMultTransposeRotMatrix();
   void GLMultMatrix();
   void Print();
};

class VisualizationScene
{
protected:
   // How to scale the visualized object(s)
   double xscale, yscale, zscale;

   SdlWindow * wnd;
public:
   VisualizationScene();
   virtual ~VisualizationScene();

   int spinning, OrthogonalProjection, print, movie;
   double ViewAngle, ViewScale;
   double ViewCenterX, ViewCenterY;

   Camera cam;

   /// Bounding box.
   double x[2], y[2], z[2];

   double rotmat[16];
   double translmat[16];

   virtual void Draw() = 0;

   void SetView(double theta, double phi);
   void Zoom(double factor);

   void Rotate(double angle, double x, double y, double z);
   void PreRotate(double angle, double x, double y, double z);

   void Rotate(double angley, double anglex);
   void Translate(double x, double y, double z = 0.0);
   void Scale(double s);
   void Scale(double s1, double s2, double s3);

   void CenterObject();
   void CenterObject2D();

   void ModelView();

   /// This is set by SetVisualizationScene
   int view;
};

#endif
