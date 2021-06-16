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

#include <iostream>
#include "openglvis.hpp"
#include "material.hpp"
#include "aux_vis.hpp"

using namespace mfem;

const int NUM_MATERIALS = 5;
extern Material materials[5];
extern Light lights[];
extern std::array<float, 4> amb_setting[];
extern Light lights_4[];

void Camera::Reset()
{
   static const double cam[9] =
   {
      0.0, 0.0,  2.5,   // eye
      0.0, 0.0, -1.0,   // dir
      0.0, 1.0,  0.0
   }; // up

   Set(cam);
}

void Camera::Set(const double cam[])
{
   eye[0] = cam[0]; eye[1] = cam[1]; eye[2] = cam[2];
   dir[0] = cam[3]; dir[1] = cam[4]; dir[2] = cam[5];
   up [0] = cam[6]; up [1] = cam[7]; up [2] = cam[8];
   Normalize(dir);
   ProjectVector(up, dir);
}

void Camera::TiltLeftRight(double angle)
{
   LinearCombination(cos(angle), up, sin(angle), GetLeft(), up);
   ProjectVector(up, dir);
}

void Camera::TurnLeftRight(double angle)
{
   LinearCombination(cos(angle), dir, sin(angle), GetLeft(), dir);
   Normalize(dir);
   ProjectVector(up, dir);
}

void Camera::TurnUpDown(double angle)
{
   double old_dir[3] = { dir[0], dir[1], dir[2] };
   double c = cos(angle), s = sin(angle);
   LinearCombination( c, old_dir, s, up, dir);
   LinearCombination(-s, old_dir, c, up, up);
   Normalize(dir);
   ProjectVector(up, dir);
}

glm::mat4 Camera::RotMatrix()
{
   GetLeft();

   double mat[16] =
   {
      -left[0], up[0], -dir[0], 0.0,
      -left[1], up[1], -dir[1], 0.0,
      -left[2], up[2], -dir[2], 0.0,
      0.0, 0.0, 0.0, 1.0
   };
   return glm::make_mat4(mat);
}

glm::mat4 Camera::TransposeRotMatrix()
{
   GetLeft();
   double mat_t[16] =
   {
      -left[0], -left[1], -left[2], 0.0,
      up[0],    up[1],    up[2],   0.0,
      -dir[0],  -dir[1],  -dir[2],  0.0,
      0.0, 0.0, 0.0, 1.0
   };
   return glm::make_mat4(mat_t);
}

glm::mat4 Camera::TranslateMatrix()
{
   glm::mat4 rotmtx = RotMatrix();
   return glm::translate(rotmtx, glm::vec3(-eye[0], -eye[1], -eye[2]));
}

void Camera::Print()
{
   std::cout <<
             "camera " << eye[0] << ' ' << eye[1] << ' ' << eye[2] << "\n"
             "       " << dir[0] << ' ' << dir[1] << ' ' << dir[2] << "\n"
             "       " <<  up[0] << ' ' <<  up[1] << ' ' <<  up[2] << '\n'
             << std::endl;
}

VisualizationScene::VisualizationScene()
{
   translmat = glm::mat4(1.0);
   rotmat = glm::mat4(1.0);
   rotmat = glm::rotate(rotmat, glm::radians(-60.f), glm::vec3(1.f, 0.f, 0.f));
   rotmat = glm::rotate(rotmat, glm::radians(-40.f), glm::vec3(0.f, 0.f, 1.f));
   xscale = yscale = zscale = 1;
   spinning = print = movie = 0;
   OrthogonalProjection = 0;
   ViewAngle = 45;
   ViewScale = 1;
   ViewCenterX = 0.0;
   ViewCenterY = 0.0;

   background = BG_WHITE;
   _use_cust_l0_pos = false;
   light_mat_idx = 3;
   use_light = true;

   palette.Init();
}

VisualizationScene::~VisualizationScene() {}

void VisualizationScene::Init(GLVisWindow* wnd)
{
   font = wnd->getFont();
}

void VisualizationScene
::DrawTriangle(gl3::GlDrawable& buff,
               const double (&pts)[4][3], const double (&cv)[4],
               const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], nor))
   {
      return;
   }

   std::array<float, 2> texcoord[3];
   std::array<float, 3> fpts[3];
   std::array<float, 3> fnorm = {(float) nor[0], (float) nor[1], (float) nor[2]};

   for (int i = 0; i < 3; i++)
   {
      float pal_coord = palette.GetColorCoord(cv[i], minv, maxv);
      texcoord[i] = { pal_coord, 1.0 };
      fpts[i] = {(float) pts[i][0], (float) pts[i][1], (float) pts[i][2]};
   }
   buff.addTriangle<gl3::VertexNormTex>(
   {fpts[0], fnorm, texcoord[0]},
   {fpts[1], fnorm, texcoord[1]},
   {fpts[2], fnorm, texcoord[2]});
}

void VisualizationScene
::DrawQuad(gl3::GlDrawable& buff,
           const double (&pts)[4][3], const double (&cv)[4],
           const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], nor))
   {
      return;
   }

   std::array<float, 2> texcoord[4];
   std::array<float, 3> fpts[4];
   std::array<float, 3> fnorm = {(float) nor[0], (float) nor[1], (float) nor[2]};

   for (int i = 0; i < 4; i++)
   {
      float pal_coord = palette.GetColorCoord(cv[i], minv, maxv);
      texcoord[i] = { pal_coord, 1.0 };
      fpts[i] = {(float) pts[i][0], (float) pts[i][1], (float) pts[i][2]};
   }
   buff.addQuad<gl3::VertexNormTex>(
   {fpts[0], fnorm, texcoord[0]},
   {fpts[1], fnorm, texcoord[1]},
   {fpts[2], fnorm, texcoord[2]},
   {fpts[3], fnorm, texcoord[3]});
}

void VisualizationScene
::DrawPatch(gl3::GlDrawable& drawable, const DenseMatrix &pts, Vector &vals,
            DenseMatrix &normals,
            const int n, const Array<int> &ind, const double minv,
            const double maxv, const int normals_opt)
{
   gl3::GlBuilder poly = drawable.createBuilder();
   double na[3];

   if (normals_opt == 1 || normals_opt == -2)
   {
      normals.SetSize(3, pts.Width());
      normals = 0.;
      for (int i = 0; i < ind.Size(); i += n)
      {
         int j;
         if (n == 3)
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), na);
         else
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), &pts(0, ind[i+3]), na);
         if (j == 0)
            for ( ; j < n; j++)
               for (int k = 0; k < 3; k++)
               {
                  normals(k, ind[i+j]) += na[k];
               }
      }
   }

   if (normals_opt != 0 && normals_opt != -1)
   {
      std::vector<gl3::VertexNormTex> vertices;
      std::vector<int> indices;
      vertices.reserve(pts.Size());
      indices.reserve(ind.Size());
      for (int i = 0; i < pts.Width(); i++)
      {
         vertices.emplace_back(
            gl3::VertexNormTex
         {
            {(float) pts(0, i), (float) pts(1, i), (float) pts(2, i)},
            {(float) normals(0, i), (float) normals(1, i), (float) normals(2, i)},
            {(float) palette.GetColorCoord(vals(i), minv, maxv), 1.0 }
         });
      }
      if (normals_opt > 0)
      {
         for (int i = 0; i < ind.Size(); i++)
         {
            indices.emplace_back(ind[i]);
         }
      }
      else
      {
         for (int i = ind.Size()-1; i >= 0; i--)
         {
            indices.emplace_back(ind[i]);
         }
      }
      if (n == 3)
      {
         drawable.addTriangleIndexed(vertices, indices);
      }
      else
      {
         drawable.addQuadIndexed(vertices, indices);
      }
   }
   else
   {
      if (n == 3)
      {
         poly.glBegin(GL_TRIANGLES);
      }
      else
      {
         poly.glBegin(GL_QUADS);
      }
      for (int i = 0; i < ind.Size(); i += n)
      {
         int j;
         if (n == 3)
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), na);
         else
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), &pts(0, ind[i+3]), na);
         if (j == 0)
         {
            if (normals_opt == 0)
            {
               poly.glNormal3dv(na);
               for ( ; j < n; j++)
               {
                  MySetColor(poly, vals(ind[i+j]), minv, maxv);
                  poly.glVertex3dv(&pts(0, ind[i+j]));
               }
            }
            else
            {
               poly.glNormal3d(-na[0], -na[1], -na[2]);
               for (j = n-1; j >= 0; j--)
               {
                  MySetColor(poly, vals(ind[i+j]), minv, maxv);
                  poly.glVertex3dv(&pts(0, ind[i+j]));
               }
            }
         }
      }
      poly.glEnd();
   }
}



gl3::RenderParams VisualizationScene::GetMeshDrawParams()
{
   gl3::RenderParams params = {};
   params.model_view.mtx = GetModelViewMtx();
   params.projection.mtx = proj_mtx;
   params.mesh_material = materials[light_mat_idx];
   params.num_pt_lights = 0;
   if (use_light)
   {
      if (light_mat_idx == 4)
      {
         for (int i = 0; i < 3; i++)
         {
            params.lights[i] = lights_4[i];
         }
         params.num_pt_lights = 3;
      }
      else
      {
         params.lights[0] = lights[light_mat_idx];
         params.num_pt_lights = 1;
      }
   }
   if (_use_cust_l0_pos)
   {
      params.lights[0].position = _l0_pos;
   }
   params.light_amb_scene = amb_setting[light_mat_idx];
   params.static_color = GetLineColor();
   return params;
}

void VisualizationScene::SetLightMatIdx(unsigned i)
{
   if (i < NUM_MATERIALS)
   {
      light_mat_idx = i;
      _use_cust_l0_pos = false;
   }
}

void VisualizationScene::SetLight0CustomPos(std::array<float, 4> pos)
{
   _l0_pos = pos;
   _use_cust_l0_pos = true;
}

void VisualizationScene::ToggleBackground()
{
   if (background == BG_BLK)
   {
      background = BG_WHITE;
   }
   else
   {
      background = BG_BLK;
   }
}

void VisualizationScene::Rotate(double angle, double x, double y, double z)
{
   gl3::GlMatrix rot_tmp;
   rot_tmp.identity();
   rot_tmp.mult(cam.TransposeRotMatrix());
   rot_tmp.rotate(angle, x, y, z);
   rot_tmp.mult(cam.RotMatrix());
   rot_tmp.mult(rotmat);
   rotmat = rot_tmp.mtx;
}

void VisualizationScene::PreRotate(double angle, double x, double y, double z)
{
   rotmat = glm::rotate<float>(rotmat, glm::radians(angle), glm::vec3(x,y,z));
}

void VisualizationScene::Rotate(double angley, double anglex)
{
   gl3::GlMatrix rot_tmp;
   rot_tmp.identity();
   rot_tmp.mult(cam.TransposeRotMatrix());
   rot_tmp.rotate(angley, 0.0, 1.0, 0.0);
   rot_tmp.rotate(anglex, 1.0, 0.0, 0.0);
   rot_tmp.mult(cam.RotMatrix());
   rot_tmp.mult(rotmat);
   rotmat = rot_tmp.mtx;

}

void VisualizationScene::Translate(double _x, double _y, double _z)
{
   gl3::GlMatrix trans_tmp;
   trans_tmp.identity();
   trans_tmp.translate(_x, -_y, _z);
   trans_tmp.mult(translmat);
   translmat = trans_tmp.mtx;
}

void VisualizationScene::Scale(double s)
{
   Scale (s, s, s);
}

void VisualizationScene::Scale(double s1, double s2, double s3)
{
   xscale *= s1;
   yscale *= s2;
   zscale *= s3;
}

void VisualizationScene::CenterObject()
{
   gl3::GlMatrix tmp_mtx;
   tmp_mtx.identity();
   translmat = tmp_mtx.mtx;
   tmp_mtx.rotate(-60.f, 1.f, 0.f, 0.f);
   tmp_mtx.rotate(-40.f, 0.f, 0.f, 1.f);
   rotmat = tmp_mtx.mtx;
}

void VisualizationScene::CenterObject2D()
{
   translmat = glm::mat4(1.0);
   rotmat = glm::mat4(1.0);
}

void VisualizationScene::SetView(double theta, double phi)
{
   gl3::GlMatrix tmp_mtx;
   tmp_mtx.identity();
   translmat = tmp_mtx.mtx;
   tmp_mtx.rotate(-theta, 1.f, 0.f, 0.f);
   tmp_mtx.rotate(-phi, 0.f, 0.f, 1.f);
   rotmat = tmp_mtx.mtx;
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

glm::mat4 VisualizationScene::GetModelViewMtx()
{
   gl3::GlMatrix modelView;
   modelView.identity();
   modelView.mult(cam.TranslateMatrix());
   modelView.mult(translmat);
   modelView.mult(rotmat);
   modelView.scale(xscale, yscale, zscale);
   modelView.translate(-(x[0]+x[1])/2, -(y[0]+y[1])/2, -(z[0]+z[1])/2);
   return modelView.mtx;
}

void VisualizationScene::DecrementAlpha()
{
   matAlpha -= 0.05;
   if (matAlpha < 0.0)
   {
      matAlpha = 0.0;
   }
   GenerateAlphaTexture();
}

void VisualizationScene::IncrementAlpha()
{
   matAlpha += 0.05;
   if (matAlpha > 1.0)
   {
      matAlpha = 1.0;
   }
   GenerateAlphaTexture();
}

void VisualizationScene::DecrementAlphaCenter()
{
   matAlphaCenter -= 0.25;
   GenerateAlphaTexture();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << matAlphaCenter << endl;
#endif
}

void VisualizationScene::IncrementAlphaCenter()
{
   matAlphaCenter += 0.25;
   GenerateAlphaTexture();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << matAlphaCenter << endl;
#endif
}

