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

#include <iostream>
#include "platform_gl.hpp"
#include "openglvis.hpp"
#include "material.hpp"
#include "aux_vis.hpp"

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

glm::mat4 Camera::TranslateMatrix() {
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
    gl = GetGlState();
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
}

VisualizationScene::~VisualizationScene() {}

void VisualizationScene::Rotate(double angle, double x, double y, double z)
{
    GlMatrix rot_tmp;
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
    GlMatrix rot_tmp;
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
    GlMatrix trans_tmp;
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
    GlMatrix tmp_mtx;
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
    GlMatrix tmp_mtx;
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

void VisualizationScene::ModelView()
{
    gl->modelView.identity();
    gl->modelView.mult(cam.TranslateMatrix());
    gl->modelView.mult(translmat);
    gl->modelView.mult(rotmat);
    gl->modelView.scale(xscale, yscale, zscale);
    gl->modelView.translate(-(x[0]+x[1])/2, -(y[0]+y[1])/2, -(z[0]+z[1])/2);
    gl->projection.mtx = _projmat;
    Set_Light();
    gl->loadMatrixUniforms();
}
