// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
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
#include "gl/types.hpp"
#include "material.hpp"
#include "palettes.hpp"
#include "mfem.hpp"
#include "geom_utils.hpp"
#include "sdl.hpp"
#include "gltf.hpp"

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

   /// Draw a 3D triangle in physical space with a central triangle removed. The
   /// cut is controlled by value of cut_lambda. See keys Ctrl+F3/F4. Similar to
   /// CutReferenceTriangle in lib/vssolution3d.cpp.
   void DrawCutTriangle(gl3::GlDrawable& buff,
                        const double (&pts)[4][3], const double (&cv)[4],
                        const double minv, const double maxv);

   /// Draw a 3D quad in physical space with a central square removed. The cut
   /// is controlled by the value of cut_lambda. See keys Ctrl+F3/F4. Similar to
   /// CutReferenceSquare in lib/vssolution3d.cpp.
   void DrawCutQuad(gl3::GlDrawable& buff,
                    const double (&pts)[4][3], const double (&cv)[4],
                    const double minv, const double maxv);

   void DrawPatch(gl3::GlDrawable& buff, const mfem::DenseMatrix &pts,
                  mfem::Vector &vals, mfem::DenseMatrix &normals,
                  const int n, const mfem::Array<int> &ind, const double minv,
                  const double maxv, const int normals_opt = 0);

   glTF_Builder::material_id AddPaletteMaterial(glTF_Builder &bld);
   glTF_Builder::material_id AddBlackMaterial(glTF_Builder &bld);
   glTF_Builder::material_id AddPaletteLinesMaterial(
      glTF_Builder &bld, glTF_Builder::material_id palette_mat);

   glTF_Builder::node_id AddModelNode(glTF_Builder &bld,
                                      const std::string &nodeName);

   // returns number of triangles
   int AddTriangles(glTF_Builder &bld,
                    glTF_Builder::mesh_id mesh,
                    glTF_Builder::buffer_id buffer,
                    glTF_Builder::material_id material,
                    const gl3::GlDrawable &gl_drawable);

   // returns number of lines
   int AddLines(glTF_Builder &bld,
                glTF_Builder::mesh_id mesh,
                glTF_Builder::buffer_id buffer,
                glTF_Builder::material_id material,
                const gl3::GlDrawable &gl_drawable);

public:
   VisualizationScene();
   virtual ~VisualizationScene();

   int spinning, OrthogonalProjection, print, movie;
   double ViewAngle, ViewScale;
   double ViewCenterX, ViewCenterY;

   Camera cam;
   PaletteState palette;

   /// Bounding box.
   struct
   {
      double x[2], y[2], z[2];
   } bb;

   /// Amount of face cutting with keys Ctrl-F3/F4 (0: no cut, 1: cut to edges)
   double cut_lambda;
   /// Have the reference geometries been updated for the cut?
   bool cut_updated;

   glm::mat4 rotmat;
   glm::mat4 translmat;

   float matAlpha = 1.0;
   float matAlphaCenter = 0.5;

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

   void GenerateAlphaTexture()
   { palette.GenerateAlphaTexture(matAlpha, matAlphaCenter); }

   /// This is set by SetVisualizationScene
   int view;
};

#endif
