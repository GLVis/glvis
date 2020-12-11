// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_OPENGLVIS
#define GLVIS_OPENGLVIS

#include <cmath>
#include "sdl.hpp"
#include "gl/types.hpp"
#include "material.hpp"

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

   void Print();

   glm::mat4 RotMatrix();
   glm::mat4 TransposeRotMatrix();
   glm::mat4 TranslateMatrix();
};

class VisualizationScene
{
protected:
   // How to scale the visualized object(s)
   double xscale, yscale, zscale;

   SdlWindow * wnd;

   glm::mat4 proj_mtx;

   enum
   {
      BG_BLK = 0,
      BG_WHITE = 1
   } background;

   const Material BLK_MAT =
   {
      {{ 0.0, 0.0, 0.0, 1.0 }},
      {{ 0.0, 0.0, 0.0, 1.0 }},
      {{ 0.0, 0.0, 0.0, 1.0 }},
      0.0
   };

   std::array<float, 4> _l0_pos;
   bool _use_cust_l0_pos;
   int light_mat_idx;
   bool use_light;

   gl3::RenderParams GetMeshDrawParams();
   glm::mat4 GetModelViewMtx();

   std::array<float, 4> GetLineColor()
   {
      if (background == BG_BLK)
      {
         return { 1.f, 1.f, 1.f, 1.f };
      }
      else
      {
         return { 0.f, 0.f, 0.f, 1.f };
      }
   }

public:
   VisualizationScene();
   virtual ~VisualizationScene();

   int spinning, OrthogonalProjection, print, movie;
   double ViewAngle, ViewScale;
   double ViewCenterX, ViewCenterY;

   Camera cam;

   /// Bounding box.
   double x[2], y[2], z[2];

   glm::mat4 rotmat;
   glm::mat4 translmat;

   virtual gl3::SceneInfo GetSceneObjs() = 0;

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

   void SetProjectionMtx(glm::mat4 projection) { proj_mtx = projection; }
   void SetLightMatIdx(unsigned i);
   int GetLightMatIdx() { return light_mat_idx; }

   void SetLight0CustomPos(std::array<float, 4> pos);
   void ToggleBackground();

   /// This is set by SetVisualizationScene
   int view;
};

#endif
