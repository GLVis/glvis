// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_OPENGLVIS_HPP
#define GLVIS_OPENGLVIS_HPP

#include <cmath>
#include "sdl.hpp"
#include "font.hpp"
#include "gl/types.hpp"
#include "material.hpp"
#include "palettes.hpp"
#include "mfem.hpp"
#include "geom_utils.hpp"

// Visualization header file

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

   GlVisFont* font = nullptr;

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

   void MySetColor (gl3::GlBuilder& builder, double val, double min, double max)
   {
      MySetColor(builder, palette.GetColorCoord(val, min, max));
   }

   void MySetColor (gl3::GlBuilder& builder, double val)
   {
      if (val < 0.0) { val = 0.0; }
      if (val > 1.0) { val = 1.0; }

      builder.glTexCoord2f(val, 1.0);
   }

   // We only need 3 points, but the array is 4x3
   void DrawTriangle(gl3::GlDrawable& buff,
                     const double (&pts)[4][3], const double (&cv)[4],
                     const double minv, const double maxv);

   void DrawQuad(gl3::GlDrawable& buff,
                 const double (&pts)[4][3], const double (&cv)[4],
                 const double minv, const double maxv);

   void DrawPatch(gl3::GlDrawable& buff, const mfem::DenseMatrix &pts,
                  mfem::Vector &vals, mfem::DenseMatrix &normals,
                  const int n, const mfem::Array<int> &ind, const double minv,
                  const double maxv, const int normals_opt = 0);

public:
   VisualizationScene();
   virtual ~VisualizationScene();

   int spinning, OrthogonalProjection, print, movie;
   double ViewAngle, ViewScale;
   double ViewCenterX, ViewCenterY;

   Camera cam;
   PaletteState palette;

   /// Bounding box.
   double x[2], y[2], z[2];

   glm::mat4 rotmat;
   glm::mat4 translmat;

   float matAlpha = 1.0;
   float matAlphaCenter = 0.5;

   /// Needs to be called before GetSceneObjs()
   virtual void UpdateWindowSize(int w, int h, int gl_w, int gl_h) = 0;
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

   // Toggles between orthogonal and perspective projections.
   void ToggleProjectionMode()
   { OrthogonalProjection = !OrthogonalProjection; }

   void SetProjectionMtx(glm::mat4 projection) { proj_mtx = projection; }
   void SetLightMatIdx(unsigned i);
   int GetLightMatIdx() { return light_mat_idx; }

   void SetLight0CustomPos(std::array<float, 4> pos);
   void ToggleBackground();
   std::array<float, 4> GetBackgroundColor()
   {
      if (background == BG_BLK)
      {
         return { 0.f, 0.f, 0.f, 1.f };
      }
      else
      {
         return { 1.f, 1.f, 1.f, 1.f };
      }
   }

   void GenerateAlphaTexture()
   { palette.GenerateAlphaTexture(matAlpha, matAlphaCenter); }

   void DecrementAlpha();
   void IncrementAlpha();
   void DecrementAlphaCenter();
   void IncrementAlphaCenter();

   /// This is set by SetVisualizationScene
   int view;
};

#endif
