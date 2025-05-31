// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include <cstdlib>
#include <cmath>

#include <unordered_set>

#include <iomanip>
#include <sstream>
#include <limits>
using namespace std;

#include "vsdata.hpp"
#include "aux_vis.hpp"
#include "material.hpp"
#include "palettes.hpp"
#include "threads.hpp"

const char *strings_off_on[] = { "off", "on" };

void VisualizationSceneScalarData::FixValueRange()
{
   double am = fabs(minv);
   if (am < fabs(maxv)) { am = fabs(maxv); }
   if (float(am) < 100*numeric_limits<float>::min()) { am = 1e-3; }
   if ((maxv-minv) < 1e-5*am)
   {
      // Shading quality may be bad since OpenGL uses single precision. We
      // should probably pre-scale the solution before feeding it to OpenGL
      int old_prec = cout.precision(12);
      cout << "[minv,maxv] = " << "[" << minv << "," << maxv
           << "] (maxv-minv = " << maxv-minv << ")\n --> ";
      minv -= 0.49999e-5*am;
      maxv += 0.50001e-5*am;
      cout << "[" << minv << "," << maxv << "]" << endl;
      cout.precision(old_prec);
   }
}

int VisualizationSceneScalarData::GetFunctionAutoRefineFactor(GridFunction &gf)
{
   Mesh *mesh = gf.FESpace()->GetMesh();
   const int order = gf.FESpace()->GetMaxElementOrder();

   // check for integral elements
   const int dim = mesh->Dimension();
   const FiniteElementCollection *fec = gf.FESpace()->FEColl();
   if (fec && fec->GetMapType(dim) == FiniteElement::INTEGRAL)
   {
      cout << "Warning: integral elements are non-polynomial in the physical space,\n"
           << "         consider increasing the refinement by the key 'o'."
           << endl;
   }

   return std::max(order, 1);
}

int VisualizationSceneScalarData::GetAutoRefineFactor()
{
   const int dim = mesh->Dimension();
   const int ne = (dim == 3)?(mesh->GetNBE()):(mesh->GetNE());

   // determine the refinement based on the order of the mesh and grid function
   int order_ref = GetFunctionAutoRefineFactor();

   // mesh
   const FiniteElementSpace *nfes = mesh->GetNodalFESpace();
   if (nfes)
   {
      const int order = nfes->GetMaxElementOrder();
      order_ref = std::max(order_ref, order);
   }

   // limit the total number of vertices
   int auto_ref_surf_vert = ne * (order_ref+1) * (order_ref+1);
   auto_ref_surf_vert = std::min(std::max(auto_ref_surf_vert,
                                          auto_ref_min_surf_vert), auto_ref_max_surf_vert);

   // approach the given number of vertices
   int ref = 1;
   while (ref < auto_ref_max && ne*(ref+2)*(ref+2) <= auto_ref_surf_vert)
   { ref++; }

   if (ref < order_ref)
   {
      cout << "Warning: the automatic refinement does not resolve the data fully,\n"
           << "         consider increasing the refinement by the key 'o'."
           << endl;
   }

   return ref;
}

void VisualizationSceneScalarData::DoAutoscale(bool prepare)
{
   if (autoscale == 1)
   {
      FindNewBoxAndValueRange(prepare);
   }
   else if (autoscale == 2)
   {
      FindNewValueRange(prepare);
   }
   else if (autoscale == 3)
   {
      FindMeshBox(prepare);
   }
}

void VisualizationSceneScalarData::DoAutoscaleValue(bool prepare)
{
   if (autoscale == 1 || autoscale == 3)
   {
      FindNewBoxAndValueRange(prepare);
   }
   else
   {
      FindNewValueRange(prepare);
   }
}

template<typename T>
static std::array<float, 3> ToVec3(const T* vec)
{
   return { (float) vec[0], (float) vec[1], (float) vec[2] };
}

void VisualizationSceneScalarData::Cone(gl3::GlDrawable& buf,
                                        glm::mat4 xfrm,
                                        double cval)
{
   const int n = 8;
   const double step = 2*M_PI/n;
   const double nz = (1.0/4.0);
   double point = step;

   glm::mat3 normXfrm = glm::inverseTranspose(glm::mat3(xfrm));

   glm::vec3 start1vtx = glm::vec3(xfrm * glm::vec4(0, 0, 0, 1));
   glm::vec3 start1norm = glm::vec3(normXfrm * glm::vec3(0, 0, 1));
   glm::vec3 start2vtx = glm::vec3(xfrm * glm::vec4(1, 0, -4, 1));
   glm::vec3 start2norm = glm::vec3(normXfrm * glm::vec3(1, 0, nz));

   int indices[n*3];
   for (int i = 0; i < n; i++)
   {
      indices[3*i] = 0;
      indices[3*i+1] = i+1;
      indices[3*i+2] = i+2;
   }
   if (cval == HUGE_VAL)
   {
      gl3::VertexNorm verts[n+2];

      verts[0].coord = ToVec3(glm::value_ptr(start1vtx));
      verts[0].norm = ToVec3(glm::value_ptr(start1norm));
      verts[1].coord = ToVec3(glm::value_ptr(start2vtx));
      verts[1].norm = ToVec3(glm::value_ptr(start2norm));
      for (int i = 2; i <= n+1; i++)
      {
         glm::vec3 baseVtx = glm::vec3(xfrm * glm::vec4(cos(point), sin(point), -4, 1));
         glm::vec3 baseNorm = glm::vec3(normXfrm * glm::vec3(cos(point), sin(point),
                                                             nz));
         verts[i].coord = ToVec3(glm::value_ptr(baseVtx));
         verts[i].norm = ToVec3(glm::value_ptr(baseNorm));
         point += step;
      }
      buf.addTriangleIndexed(n+2, verts, n*3, indices);
   }
   else
   {
      float colortex = palette.GetColorCoord(cval, minv, maxv);
      gl3::VertexNormTex verts[n+2];
      verts[0].coord = ToVec3(glm::value_ptr(start1vtx));
      verts[0].norm = ToVec3(glm::value_ptr(start1norm));
      verts[0].texCoord = colortex;
      verts[1].coord = ToVec3(glm::value_ptr(start2vtx));
      verts[1].norm = ToVec3(glm::value_ptr(start2norm));
      verts[1].texCoord = colortex;
      for (int i = 2; i <= n+1; i++)
      {
         glm::vec3 baseVtx = glm::vec3(xfrm * glm::vec4(cos(point), sin(point), -4, 1));
         glm::vec3 baseNorm = glm::vec3(normXfrm * glm::vec3(cos(point), sin(point),
                                                             nz));
         verts[i].coord = ToVec3(glm::value_ptr(baseVtx));
         verts[i].norm = ToVec3(glm::value_ptr(baseNorm));
         verts[i].texCoord = colortex;
         point += step;
      }
      buf.addTriangleIndexed(n+2, verts, n*3, indices);
   }
}

// Draw an arrow starting at point (px, py, pz) with orientation (vx, vy, vz)
// and length "length".
void VisualizationSceneScalarData::Arrow3(gl3::GlDrawable& buf,
                                          double px, double py, double pz,
                                          double vx, double vy, double vz,
                                          double length,
                                          double cone_scale,
                                          double cval)
{
   double xc = 0.5*(bb.x[0]+bb.x[1]);
   double yc = 0.5*(bb.y[0]+bb.y[1]);
   double zc = 0.5*(bb.z[0]+bb.z[1]);
   glm::mat4 xfrm(1.0);
   xfrm = glm::translate(xfrm, glm::vec3(xc, yc, zc));
   xfrm = glm::scale(xfrm, glm::vec3(1.0/xscale, 1.0/yscale, 1.0/zscale));

   double rlen = length/sqrt(vx*vx+vy*vy+vz*vz);
   px = (px-xc)*xscale;
   py = (py-yc)*yscale;
   pz = (pz-zc)*zscale;
   vx *= rlen*xscale;
   vy *= rlen*yscale;
   vz *= rlen*zscale;

   if (arrow_scaling_type == 0)
   {
      length = sqrt(vx*vx+vy*vy+vz*vz);
   }

   xfrm = glm::translate(xfrm, glm::vec3(px, py, pz));

   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   float phi   = acos(vz/rhos);
   float theta;
   theta = atan2 (vy, vx);

   xfrm = glm::rotate(xfrm, theta, glm::vec3(0, 0, 1));
   xfrm = glm::rotate(xfrm, phi, glm::vec3(0, 1, 0));

   xfrm = glm::scale(xfrm, glm::vec3(length));


   if (arrow_type == 1)
   {
      xfrm = glm::translate(xfrm, glm::vec3(0, 0, -0.5));
   }

   glm::vec4 pt1 = xfrm * glm::vec4(0, 0, 0, 1);
   glm::vec4 pt2 = xfrm * glm::vec4(0, 0, 1, 1);

   if (cval == HUGE_VAL)
   {
      buf.addLine<gl3::Vertex>(
         gl3::Vertex{ToVec3(glm::value_ptr(pt1))},
         gl3::Vertex{ToVec3(glm::value_ptr(pt2))});
   }
   else
   {
      double colortex = palette.GetColorCoord(cval, minv, maxv);
      buf.addLine<gl3::VertexTex>(
         gl3::VertexTex{ToVec3(glm::value_ptr(pt1)), (float)colortex},
         gl3::VertexTex{ToVec3(glm::value_ptr(pt2)), (float)colortex});
   }

   xfrm = glm::translate(xfrm, glm::vec3(0, 0, 1));
   xfrm = glm::scale(xfrm, glm::vec3(cone_scale));

   if (cone_scale > 0.0)
   {
      Cone(buf, xfrm, cval);
   }
}

void VisualizationSceneScalarData::Arrow2(gl3::GlDrawable& buf,
                                          double px, double py, double pz,
                                          double vx, double vy, double vz,
                                          double length,
                                          double cone_scale,
                                          double cval)
{
   glm::mat4 xfrm(1.0);
   xfrm = glm::translate(xfrm, glm::vec3(px, py, pz));

   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   float phi   = acos(vz/rhos);
   float theta;
   theta = atan2 (vy, vx);

   xfrm = glm::rotate(xfrm, theta, glm::vec3(0, 0, 1));
   xfrm = glm::rotate(xfrm, phi, glm::vec3(0, 1, 0));

   xfrm = glm::scale(xfrm, glm::vec3(length));

   glm::vec4 pt1 = xfrm * glm::vec4(0, 0, 0, 1);
   glm::vec4 pt2 = xfrm * glm::vec4(0, 0, 1, 1);

   if (cval == HUGE_VAL)
   {
      buf.addLine<gl3::Vertex>(
         gl3::Vertex{ToVec3(glm::value_ptr(pt1))},
         gl3::Vertex{ToVec3(glm::value_ptr(pt2))});
   }
   else
   {
      double colortex = palette.GetColorCoord(cval, minv, maxv);
      buf.addLine<gl3::VertexTex>(
         gl3::VertexTex{ToVec3(glm::value_ptr(pt1)), (float)colortex},
         gl3::VertexTex{ToVec3(glm::value_ptr(pt2)), (float)colortex});
   }

   xfrm = glm::translate(xfrm, glm::vec3(0, 0, 1));
   xfrm = glm::scale(xfrm, glm::vec3(cone_scale));

   Cone(buf, xfrm, cval);
}

void VisualizationSceneScalarData::Arrow(gl3::GlDrawable& buf,
                                         double px, double py, double pz,
                                         double vx, double vy, double vz,
                                         double length,
                                         double cone_scale,
                                         double cval)
{
   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   if (rhos == 0.0)
   {
      return;
   }
   double phi = acos(vz/rhos), theta = atan2(vy, vx);
   constexpr int n = 8;
   const double step = 2*M_PI/n, nz = (1.0/4.0);
   double point = step, cone[n+4][3], normal[n+2][3];

   cone[0][0] = 0;          cone[0][1] = 0; cone[0][2] = 1;
   cone[1][0] = cone_scale; cone[1][1] = 0; cone[1][2] = -4*cone_scale + 1;
   normal[0][0] = 0.0/cone_scale;
   normal[0][1] = 0.0/cone_scale;
   normal[0][2] = 1.0/cone_scale;
   normal[1][0] = 1.0/cone_scale;
   normal[1][1] = 0.0/cone_scale;
   normal[1][2] = nz/cone_scale;

   for (int i=2; i<n+1; i++)
   {
      normal[i][0] = cos(point)/cone_scale;
      normal[i][1] = sin(point)/cone_scale;
      normal[i][2] = nz/cone_scale;

      cone[i][0] = cos(point)*cone_scale;
      cone[i][1] = sin(point)*cone_scale;
      cone[i][2] = -4*cone_scale + 1;
      point += step;
   }
   cone[n+1][0] = cone_scale; cone[n+1][1] = 0; cone[n+1][2] =-4*cone_scale + 1;
   normal[n+1][0] = 1.0/cone_scale;
   normal[n+1][1] = 0.0/cone_scale;
   normal[n+1][2] = nz/cone_scale;

   cone[n+2][0] = 0; cone[n+2][1] = 0; cone[n+2][2] = 0;
   cone[n+3][0] = 0; cone[n+3][1] = 0; cone[n+3][2] = 1;

   if (arrow_scaling_type == 0)
   {
      length = rhos;
   }

   // double xc = 0.5*(x[0]+x[1]), yc = 0.5*(y[0]+y[1]), zc = 0.5*(z[0]+z[1]);
   double coord[3];
   // double rlen = length/rhos;

   // px = (px-xc)*xscale;  py = (py-yc)*yscale;  pz = (pz-zc)*zscale;
   // vx *= rlen*xscale;    vy *= rlen*yscale;    vz *= rlen*zscale;

   if (arrow_type == 1)
      for (int i=0; i<n+4; i++)
      {
         cone[i][2] -= 0.5;
      }

   double M[3][3]= {{cos(theta)*cos(phi), -sin(theta),  cos(theta)*sin(phi)},
      {sin(theta)*cos(phi),  cos(theta),  sin(theta)*sin(phi)},
      {          -sin(phi),          0.,             cos(phi)}
   };
   double v[3] = { M[0][2]/xscale, M[1][2]/yscale, M[2][2]/zscale };
   length /= sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

   for (int i=0; i<n+4; i++)
   {
      for (int j=0; j<3; j++)
      {
         coord[j] = cone[i][j] * length;
      }

      for (int k=0; k<3; k++)
      {
         cone[i][k] = 0.;
         for (int j=0; j<3; j++)
         {
            cone[i][k] += M[k][j] * coord[j];
         }
      }
      // cone[i][0] = (cone[i][0] + px)/xscale + xc;
      // cone[i][1] = (cone[i][1] + py)/yscale + yc;
      // cone[i][2] = (cone[i][2] + pz)/zscale + zc;
      cone[i][0] = cone[i][0]/xscale + px;
      cone[i][1] = cone[i][1]/yscale + py;
      cone[i][2] = cone[i][2]/zscale + pz;
   }

   for (int i=0; i<=n+1; i++)
   {
      for (int j=0; j<3; j++)
      {
         coord[j] = normal[i][j];
      }

      for (int k=0; k<3; k++)
      {
         normal[i][k] = 0.;
         for (int j=0; j<3; j++)
         {
            normal[i][k] += M[k][j] * coord[j];
         }
      }
      normal[i][0] *= xscale;
      normal[i][1] *= yscale;
      normal[i][2] *= zscale;
   }

   int indices[n*3];
   for (int i = 0; i < n; i++)
   {
      indices[3*i] = 0;
      indices[3*i+1] = i+1;
      indices[3*i+2] = i+2;
   }
   if (cval == HUGE_VAL)
   {
      gl3::VertexNorm verts[n+2];
      for (int i = 0; i <= n+1; i++)
      {
         verts[i].coord = ToVec3(cone[i]);
         verts[i].norm = ToVec3(normal[i]);
      }
      buf.addTriangleIndexed(n+2, verts, n*3, indices);
      buf.addLine<gl3::Vertex>(
         gl3::Vertex{ToVec3(cone[n+2])},
         gl3::Vertex{ToVec3(cone[n+3])});

   }
   else
   {
      float colortex = palette.GetColorCoord(cval, minv, maxv);
      gl3::VertexNormTex verts[n+2];
      for (int i = 0; i <= n+1; i++)
      {
         verts[i].coord = ToVec3(cone[i]);
         verts[i].norm = ToVec3(normal[i]);
         verts[i].texCoord = colortex;
      }
      buf.addTriangleIndexed(n+2, verts, n*3, indices);
      buf.addLine<gl3::VertexTex>(
         gl3::VertexTex{ToVec3(cone[n+2]), colortex},
         gl3::VertexTex{ToVec3(cone[n+3]), colortex});
   }
}

void VisualizationSceneScalarData::SetColorbarNumberFormat(int precision,
                                                           char format,
                                                           bool showsign)
{
   colorbar_formatter = NumberFormatter(precision, format, showsign);
   // The first two arguments are required but I don't think they are used?
   PrepareColorBar(0,0);
}

void VisualizationSceneScalarData::SetColorbarNumberFormat(string formatting)
{
   colorbar_formatter = NumberFormatter(formatting);
   // The first two arguments are required but I don't think they are used?
   PrepareColorBar(0,0);
}

void VisualizationSceneScalarData::PrepareColorBar (double minval,
                                                    double maxval,
                                                    Array<double> *mesh_level,
                                                    Array<double> *lsurf_levels)
{

   int i;

   float miny;
   float maxy;
   float minx;
   float maxx;
   float posz = -4.0;

   if (OrthogonalProjection)
   {
      miny = -.65;
      maxy =  .65;
      minx = 0.73;
      maxx = 0.80;
   }
   else
   {
      miny = -1.;
      maxy =  1.;
      minx =  1.2;
      maxx =  1.3;
   }
   color_bar.clear();
   color_bar.addQuad<gl3::VertexTex>(
   {{minx, miny, posz}, 0.f},
   {{maxx, miny, posz}, 0.f},
   {{maxx, maxy, posz}, 1.f},
   {{minx, maxy, posz}, 1.f}
   );

   static const int border = 2;

   if (border == 1)
   {
      color_bar.addLines<gl3::Vertex>(
      {
         {minx, miny, posz}, {maxx, miny, posz},
         {maxx, miny, posz}, {maxx, maxy, posz},
         {maxx, maxy, posz}, {minx, maxy, posz},
         {minx, maxy, posz}, {minx, miny, posz}
      });
   }
   else if (border == 2)
   {
      color_bar.addLine<gl3::Vertex>({minx, miny, posz}, {minx, maxy, posz});
      color_bar.addLine<gl3::Vertex>({maxx, miny, posz}, {maxx, maxy, posz});
   }

   if (lsurf_levels)
   {
      for (i = 0; i < lsurf_levels->Size(); i++)
      {
         float Y = miny + (maxy - miny) * LogUVal((*lsurf_levels)[i]);
         color_bar.addLine<gl3::Vertex>({minx, Y, posz}, {maxx, Y, posz});
      }
   }
   if (mesh_level)
   {
      for (i = 0; i < mesh_level->Size(); i++)
      {
         float Y = miny + (maxy - miny) * LogUVal((*mesh_level)[i]);
         color_bar.addLine<gl3::Vertex>({minx, Y, posz}, {maxx, Y, posz});
      }
   }

   const double text_x = maxx + 0.4*(maxx-minx);
   double val;
   double Y;
   if (!mesh_level)
   {
      for (i = 0; i <= 4; i++)
      {
         Y = miny + i * (maxy-miny) / 4;

         val = ULogVal(i / 4.0);

         color_bar.addText(text_x,Y,posz,colorbar_formatter(val));
      }
   }
   else
   {
      for (i = 0; i < mesh_level->Size(); i++)
      {
         val = (*mesh_level)[i];
         Y = miny + (maxy - miny) * LogUVal(val);

         color_bar.addText(text_x,Y,posz,colorbar_formatter(val));
      }
   }

   if (lsurf_levels)
   {
      for (i = 0; i < lsurf_levels->Size(); i++)
      {
         val = (*lsurf_levels)[i];
         Y = miny + (maxy - miny) * LogUVal(val);

         color_bar.addText(text_x,Y,posz,colorbar_formatter(val));
      }
   }
   updated_bufs.emplace_back(&color_bar);
}

// Draw a centered caption at the top (visible with the colorbar)
void VisualizationSceneScalarData::PrepareCaption()
{
   bool empty = win.plot_caption.empty();
   colorbar = (colorbar ? empty+1 : !empty);

   string caption(win.plot_caption);
   if (!win.extra_caption.empty())
   {
      caption += " (" + win.extra_caption + ")";
   }

   caption_buf.clear();
   caption_buf.addText(0, 0, 0, caption);
   updated_bufs.emplace_back(&caption_buf);
   GetFont()->getObjectSize(caption, caption_w, caption_h);
}

static thread_local VisualizationSceneScalarData *vsdata;
static thread_local Window *window;

void KeycPressed(GLenum state)
{
   if (state & KMOD_ALT)
   {
      cout << "Setting colorbar number formatting..." << endl;
      int default_precision = 4;
      char default_format = 'd';
      bool default_showsign = false;

      int precision = prompt<int>("Enter precision (4): ",
      &default_precision, [](int p) { return p>=0; });
      char format =
         prompt<char>("Enter format [(d)efault, (f)ixed, (s)cientific] (d): ",
      &default_format, [](char c) { return c=='d' || c=='f' || c=='s'; });
      bool showsign = prompt<bool>("Show sign? [(1)true, (0)false] (0): ",
                                   &default_showsign);
      vsdata->SetColorbarNumberFormat(precision, format, showsign);
      SendExposeEvent();
   }
   else
   {
      vsdata->ToggleDrawColorbar();
      SendExposeEvent();
   }
}

void KeyCPressed()
{
   cout << "Enter new caption: " << flush;
   std::getline(cin, window->plot_caption);
   vsdata->PrepareCaption(); // turn on or off the caption
   SendExposeEvent();
}

void KeySPressed()
{
   vsdata -> ToggleScaling();
   SendExposeEvent();
}

void KeyaPressed()
{
   vsdata -> ToggleDrawAxes();
   SendExposeEvent();
}

void Key_Mod_a_Pressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      static const char *autoscale_modes[] = { "off", "on", "value", "mesh" };
      int autoscale = vsdata->GetAutoscale();
      autoscale = (autoscale + 1)%4;
      cout << "Autoscale: " << flush;
      vsdata->SetAutoscale(autoscale);
      cout << autoscale_modes[autoscale] << endl;
      SendExposeEvent();
   }
   else if (state & KMOD_ALT)
   {
      cout << "Setting axes number formatting..." << endl;
      int default_precision = 4;
      char default_format = 'd';
      bool default_showsign = false;

      int precision = prompt<int>("Enter precision (4): ",
      &default_precision, [](int p) { return p>=0; });
      char format =
         prompt<char>("Enter format [(d)efault, (f)ixed, (s)cientific] (d): ",
      &default_format, [](char c) { return c=='d' || c=='f' || c=='s'; });
      bool showsign = prompt<bool>("Show sign? [(1)true, (0)false] (0): ",
                                   &default_showsign);
      vsdata->SetAxisNumberFormat(precision, format, showsign);
      SendExposeEvent();
   }
   else
   {
      vsdata->ToggleDrawAxes();
      SendExposeEvent();
   }
}

void KeyHPressed()
{
   cout << vsdata->GetHelpString() << flush;
}

void KeylPressed()
{
   vsdata -> ToggleLight();
   SendExposeEvent();
}

void KeyLPressed()
{
   vsdata->ToggleLogscale(true);
   SendExposeEvent();
}

void KeyrPressed()
{
   window->vs -> spinning = 0;
   RemoveIdleFunc(MainLoop);
   vsdata -> CenterObject();

   window->vs -> ViewAngle = 45.0;
   window->vs -> ViewScale = 1.0;
   window->vs -> ViewCenterX = 0.0;
   window->vs -> ViewCenterY = 0.0;
   window->vs->cam.Reset();
   vsdata -> key_r_state = 0;
   SendExposeEvent();
}

void KeyRPressed()
{
   window->vs->spinning = 0;
   RemoveIdleFunc(MainLoop);
   vsdata->Toggle2DView();
   SendExposeEvent();
}

void KeypPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      KeyCtrlP();
   }
   else
   {
      window->vs->palette.NextIndex();
      SendExposeEvent();
   }
}

void KeyPPressed()
{
   window->vs->palette.PrevIndex();
   SendExposeEvent();
}

static void KeyF5Pressed()
{
   int n;
   double min, max;

   cout << "Enter min : " << flush;
   cin >> min;
   cout << "Enter max : " << flush;
   cin >> max;
   cout << "Enter n : " << flush;
   cin >> n;

   vsdata -> SetLevelLines (min, max, n, 0);

   vsdata -> UpdateLevelLines();
   SendExposeEvent();
}

void KeyF6Pressed()
{
   int RepeatPaletteTimes = vsdata->palette.GetRepeatTimes();
   cout << "Palette is repeated " << RepeatPaletteTimes << " times.\n"
        << "(Negative value means the palette is flipped.)\n"
        << "Enter new value: " << flush;
   cin >> RepeatPaletteTimes;
   if (RepeatPaletteTimes == 0)
   {
      RepeatPaletteTimes = 1;
   }
   cout << "Palette will be repeated " << RepeatPaletteTimes
        << " times now.\n\n";
   vsdata->palette.SetRepeatTimes(RepeatPaletteTimes);

   int pal = vsdata->palette.ChoosePalette();

   int colors_used = vsdata->palette.GetNumColors(pal);
   int palette_size = vsdata->palette.GetSize(pal);
   cout << "\nPalette is using " << colors_used << " colors.\n"
        << "Enter new value (0 = use original " << palette_size
        << " colors): " << flush;
   cin >> colors_used;
   if (colors_used == 1) { colors_used = 0; }
   vsdata->palette.SetNumColors(colors_used);

   vsdata->palette.GenerateTextures();
   vsdata->palette.SetIndex(pal);

   colors_used = vsdata->palette.GetNumColors();
   cout << "Palette will be using " << colors_used << " colors now.\n";

   vsdata->EventUpdateColors();
   SendExposeEvent();
}

void KeyF7Pressed(GLenum state)
{
   if (state & KMOD_SHIFT)
   {
      cout << "Current bounding box:\n"
           << "   min: (" << vsdata->bb.x[0] << ',' << vsdata->bb.y[0] << ','
           << vsdata->bb.z[0] << ")\n"
           << "   max: (" << vsdata->bb.x[1] << ',' << vsdata->bb.y[1] << ','
           << vsdata->bb.z[1] << ")\n"
           << "Enter new bounding box:\n"
           << "x_min = " << flush;
      cin >> vsdata->bb.x[0];
      cout << "y_min = " << flush;
      cin >> vsdata->bb.y[0];
      cout << "z_min = " << flush;
      cin >> vsdata->bb.z[0];
      cout << "x_max = " << flush;
      cin >> vsdata->bb.x[1];
      cout << "y_max = " << flush;
      cin >> vsdata->bb.y[1];
      cout << "z_max = " << flush;
      cin >> vsdata->bb.z[1];
      cout << "New bounding box:\n"
           << "   min: (" << vsdata->bb.x[0] << ',' << vsdata->bb.y[0] << ','
           << vsdata->bb.z[0] << ")\n"
           << "   max: (" << vsdata->bb.x[1] << ',' << vsdata->bb.y[1] << ','
           << vsdata->bb.z[1] << ")\n" << flush;
      vsdata->UpdateBoundingBox();
      SendExposeEvent();
   }
   else
   {
      cout << "[minv,maxv] = [" << vsdata->GetMinV() << "," << vsdata->GetMaxV()
           << "]  maxv-minv = " << vsdata->GetMaxV()-vsdata->GetMinV() << "\n"
           << "New value for minv: " << flush;
      cin >> vsdata->GetMinV();
      cout << "New value for maxv: " << flush;
      cin >> vsdata->GetMaxV();
      vsdata->UpdateValueRange(true);
      SendExposeEvent();
   }
}

void KeyBackslashPressed()
{
   float x, y, z, w;

   cout << "Enter light source position\n(0,0,1,w) - from camera\n"
        "(0,1,0,w) - from above\n(1,0,0,w) - from the right\n"
        "w = 0/1  defines directional/spot light\n";
   cout << "x = " << flush;
   cin >> x;
   cout << "y = " << flush;
   cin >> y;
   cout << "z = " << flush;
   cin >> z;
   cout << "w = " << flush;
   cin >> w;

   vsdata->SetLight0CustomPos({x, y, z, w});
   SendExposeEvent();
}

void KeyTPressed()
{
   int ml;

   ml = (vsdata->GetLightMatIdx() + 1) % 5;
   vsdata->SetLightMatIdx(ml);
   SendExposeEvent();
   cout << "New material/light : " << ml << endl;
}

void KeygPressed()
{
   vsdata->ToggleBackground();
   vsdata->PrepareAxes();
   vsdata->EventUpdateBackground();
   SendExposeEvent();
}

void KeyGPressed()
{
   vsdata->glTF_Export();
}

void KeyF1Pressed()
{
   vsdata->PrintState();
}

void KeyF2Pressed()
{
   vsdata -> EventUpdateColors();
   vsdata -> PrepareLines();
   // vsdata->CPPrepare();
   SendExposeEvent();
}

void KeykPressed()
{
   window->vs->matAlpha -= 0.05;
   if (window->vs->matAlpha < 0.0)
   {
      window->vs->matAlpha = 0.0;
   }
   window->vs->GenerateAlphaTexture();
   SendExposeEvent();
}

void KeyKPressed()
{
   window->vs->matAlpha += 0.05;
   if (window->vs->matAlpha > 1.0)
   {
      window->vs->matAlpha = 1.0;
   }
   window->vs->GenerateAlphaTexture();
   SendExposeEvent();
}

void KeyAPressed()
{
   bool curr_aa = window->wnd->getRenderer().getAntialiasing();
   window->wnd->getRenderer().setAntialiasing(!curr_aa);

   cout << "Multisampling/Antialiasing: "
        << strings_off_on[!curr_aa ? 1 : 0] << endl;

   // vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyCommaPressed()
{
   window->vs->matAlphaCenter -= 0.25;
   // vsdata -> EventUpdateColors();
   window->vs->GenerateAlphaTexture();
   SendExposeEvent();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << window->vs->matAlphaCenter << endl;
#endif
}

void KeyLessPressed()
{
   window->vs->matAlphaCenter += 0.25;
   // vsdata -> EventUpdateColors();
   window->vs->GenerateAlphaTexture();
   SendExposeEvent();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << window->vs->matAlphaCenter << endl;
#endif
}

void KeyGravePressed()
{
   vsdata->ToggleRuler();
   SendExposeEvent();
}

void KeyTildePressed()
{
   vsdata->RulerPosition();
   SendExposeEvent();
}

void KeyToggleTexture()
{
   vsdata->ToggleTexture();
   SendExposeEvent();
}

void VisualizationSceneScalarData::PrintLogscale(bool warn)
{
   if (warn)
   {
      cout << "The range [" << minv << ',' << maxv
           << "] is not appropriate for logarithmic scale!" << endl;
   }
   cout << "Logarithmic scale: " << strings_off_on[logscale ? 1 : 0]
        << endl;
}

void VisualizationSceneScalarData::ToggleLogscale(bool print)
{
   if (logscale || LogscaleRange())
   {
      logscale = !logscale;
      palette.SetUseLogscale(logscale);
      SetLogA();
      SetLevelLines(minv, maxv, nl);
      UpdateLevelLines();
      EventUpdateColors();
      if (print)
      {
         PrintLogscale(false);
      }
   }
   else if (print)
   {
      PrintLogscale(true);
   }
   PrepareRuler();
}

void VisualizationSceneScalarData::ToggleRuler()
{
   ruler_on = (ruler_on + 1) % 3;
   PrepareRuler();
}

void VisualizationSceneScalarData::RulerPosition()
{
   cout << "Current ruler position: (" << ruler_x << ','
        << ruler_y << ',' << ruler_z << ")\n";
   cout << "x = " << flush; cin >> ruler_x;
   cout << "y = " << flush; cin >> ruler_y;
   cout << "z = " << flush; cin >> ruler_z;
   if (ruler_x < bb.x[0])
   {
      ruler_x = bb.x[0];
   }
   else if (ruler_x > bb.x[1])
   {
      ruler_x = bb.x[1];
   }
   if (ruler_y < bb.y[0])
   {
      ruler_y = bb.y[0];
   }
   else if (ruler_y > bb.y[1])
   {
      ruler_y = bb.y[1];
   }
   if (ruler_z < bb.z[0])
   {
      ruler_z = bb.z[0];
   }
   else if (ruler_z > bb.z[1])
   {
      ruler_z = bb.z[1];
   }
   cout << "New ruler position: (" << ruler_x << ','
        << ruler_y << ',' << ruler_z << ")" << endl;
   PrepareRuler();
}

void VisualizationSceneScalarData::PrepareRuler(bool log_z)
{
   float pos_z = LogVal(ruler_z, log_z);
   float x_f[2] = {(float) bb.x[0], (float) bb.x[1]};
   float y_f[2] = {(float) bb.y[0], (float) bb.y[1]};
   float z_f[2] = {(float) bb.z[0], (float) bb.z[1]};
   float ruler_x_f = (float) ruler_x,
         ruler_y_f = (float) ruler_y;
   ruler_buf.clear();
   if (ruler_on == 2)
   {
      std::array<uint8_t, 4> color = gl3::ColorU8(0.8, 0.8, 0.8, 1.0);
      std::array<float, 3> norm = { 0, 0, 1 };
      ruler_buf.addQuad<gl3::VertexNormColor>(
      {{x_f[0], y_f[0], pos_z},norm,color},
      {{x_f[1], y_f[0], pos_z},norm,color},
      {{x_f[1], y_f[1], pos_z},norm,color},
      {{x_f[0], y_f[1], pos_z},norm,color}
      );

      std::array<float, 3> norm_2 = { 0, 1, 0 };
      ruler_buf.addQuad<gl3::VertexNormColor>(
      {{x_f[0], ruler_y_f, z_f[0]},norm_2,color},
      {{x_f[0], ruler_y_f, z_f[1]},norm_2,color},
      {{x_f[1], ruler_y_f, z_f[1]},norm_2,color},
      {{x_f[1], ruler_y_f, z_f[0]},norm_2,color}
      );

      std::array<float, 3> norm_3 = { 1, 0, 0 };
      ruler_buf.addQuad<gl3::VertexNormColor>(
      {{ruler_x_f, y_f[0], z_f[0]},norm_3,color},
      {{ruler_x_f, y_f[1], z_f[0]},norm_3,color},
      {{ruler_x_f, y_f[1], z_f[1]},norm_3,color},
      {{ruler_x_f, y_f[0], z_f[1]},norm_3,color}
      );

      ruler_buf.addLines<gl3::Vertex>(
      {
         {x_f[0], y_f[0], pos_z}, {x_f[1], y_f[0], pos_z},
         {x_f[1], y_f[0], pos_z}, {x_f[1], y_f[1], pos_z},
         {x_f[1], y_f[1], pos_z}, {x_f[0], y_f[1], pos_z},
         {x_f[0], y_f[1], pos_z}, {x_f[0], y_f[0], pos_z},

         {x_f[0], ruler_y_f, z_f[0]}, {x_f[1], ruler_y_f, z_f[0]},
         {x_f[1], ruler_y_f, z_f[0]}, {x_f[1], ruler_y_f, z_f[1]},
         {x_f[1], ruler_y_f, z_f[1]}, {x_f[0], ruler_y_f, z_f[1]},
         {x_f[0], ruler_y_f, z_f[1]}, {x_f[0], ruler_y_f, z_f[0]},

         {ruler_x_f, y_f[0], z_f[0]}, {ruler_x_f, y_f[1], z_f[0]},
         {ruler_x_f, y_f[1], z_f[0]}, {ruler_x_f, y_f[1], z_f[1]},
         {ruler_x_f, y_f[1], z_f[1]}, {ruler_x_f, y_f[0], z_f[1]},
         {ruler_x_f, y_f[0], z_f[1]}, {ruler_x_f, y_f[0], z_f[0]}
      });
   }

   ruler_buf.addLines<gl3::Vertex>(
   {
      {x_f[0], ruler_y_f, pos_z},
      {x_f[1], ruler_y_f, pos_z},
      {ruler_x_f, y_f[0], pos_z},
      {ruler_x_f, y_f[1], pos_z},
      {ruler_x_f, ruler_y_f, z_f[0]},
      {ruler_x_f, ruler_y_f, z_f[1]}
   });

   updated_bufs.emplace_back(&ruler_buf);
}

void VisualizationSceneScalarData::Toggle2DView()
{
   gl3::GlMatrix newrot;
   newrot.identity();
   translmat = newrot.mtx;

   switch (key_r_state)
   {
      case 0:
         break;

      case 1:
         newrot.rotate(-90.0, 1.0f, 0.0f, 0.0f);
         break;

      case 2:
         newrot.rotate(-90.0, 1.0f, 0.0f, 0.0f);
         newrot.rotate(-90.0, 0.0f, 0.0f, 1.0f);
         break;

      case 3:
         newrot.rotate(-90.0, 1.0f, 0.0f, 0.0f);
         newrot.rotate(-180.0, 0.0f, 0.0f, 1.0f);
         break;

      case 4:
         newrot.rotate(-90.0, 1.0f, 0.0f, 0.0f);
         newrot.rotate(-270.0, 0.0f, 0.0f, 1.0f);
         break;

      case 5:
         newrot.rotate(180.0, 1.0f, 0.0f, 0.0f);
         break;
   }

   // if (window->vs -> view != 2) // make 'R' work the same in 2D and 3D
   key_r_state = (key_r_state+1)%6;

   rotmat = newrot.mtx;
}

gl3::SceneInfo VisualizationSceneScalarData::GetSceneObjs()
{
   int w, h;
   wnd->getWindowSize(w, h);
   gl3::SceneInfo scene {};

   gl3::RenderParams params {};
   params.model_view.identity();
   params.projection.identity();
   params.mesh_material = BLK_MAT;
   params.num_pt_lights = 0;
   params.static_color = this->GetLineColor();
   params.use_clip_plane = false;
   params.contains_translucent = true;

   if (colorbar)
   {
      // add color bar to draw list
      params.projection.mtx = proj_mtx;
      scene.queue.emplace_back(params, &color_bar);
      params.projection.identity();
   }
   if (colorbar == 1)
   {
      // caption size is in screen pixels and needs to be centered with
      // GL pixel size
      int gl_w, gl_h;
      wnd->getGLDrawSize(gl_w, gl_h);
      // add caption to draw list
      double v_pos = 2.;
      double line_h = GetFont()->getFontLineSpacing();
      params.model_view.translate(-(double)caption_w / gl_w,
                                  1.0 - 2 * v_pos * line_h / gl_h, 0.0);
      scene.queue.emplace_back(params, &caption_buf);
   }
   params.contains_translucent = true;
   if (drawaxes && drawaxes != 3)
   {
      // add coordinate cross to draw list
      params.projection.ortho(-1.,1.,-1.,1.,-2.,2.);
      params.model_view.identity();
      params.model_view.translate(-1, -1, 0.0);
      params.model_view.scale(40.0 / w, 40.0 / h, 1);
      params.model_view.translate(2.0, 2.0, 0.0);
      params.model_view.mult(cam.RotMatrix());
      params.model_view.mult(rotmat);
      scene.queue.emplace_back(params, &coord_cross_buf);
   }
   params.projection.mtx = proj_mtx;
   params.model_view.mtx = GetModelViewMtx();
   if (drawaxes)
   {
      // add axes to draw list
      scene.queue.emplace_back(params, &axes_buf);
   }
   params.contains_translucent = false;
   if (ruler_on)
   {
      // add ruler to draw list
      params = GetMeshDrawParams();
      params.use_clip_plane = false;
      params.contains_translucent = false;
      scene.queue.emplace_back(params, &ruler_buf);
   }
   return scene;
}

void VisualizationSceneScalarData::ProcessUpdatedBufs(gl3::SceneInfo& scene)
{
   std::unordered_set<gl3::GlDrawable*> bufs_in_scene;
   for (const auto& qelem : scene.queue)
   {
      bufs_in_scene.insert(qelem.second);
   }

   auto it = updated_bufs.begin();
   while (it != updated_bufs.end())
   {
      if (bufs_in_scene.find(*it) != bufs_in_scene.end())
      {
         scene.needs_buffering.emplace_back(*it);
         it = updated_bufs.erase(it);
      }
      else
      {
         ++it;
      }
   }
}

void VisualizationSceneScalarData::glTF_ExportBox(
   glTF_Builder &bld,
   glTF_Builder::buffer_id buffer,
   glTF_Builder::material_id black_mat)
{
   auto white_mat = glTF_Builder::material_id{glTF_Builder::INVALID_ID};

   auto box_node = AddModelNode(bld, "Box");
   auto box_mesh = bld.addMesh("Box Mesh");
   bld.addNodeMesh(box_node, box_mesh);

   int nlines = AddLines(
                   bld,
                   box_mesh,
                   buffer,
                   (drawaxes != 3) ? black_mat : white_mat,
                   axes_buf);
   if (nlines == 0)
   {
      cout << "glTF export: no box found to export!" << endl;
   }
}

void VisualizationSceneScalarData::glTF_ExportElements(
   glTF_Builder &bld,
   glTF_Builder::buffer_id buffer,
   glTF_Builder::material_id palette_mat,
   const gl3::GlDrawable &gl_drawable)
{
   auto elements_node = AddModelNode(bld, "Elements");
   auto elements_mesh = bld.addMesh("Elements Mesh");
   bld.addNodeMesh(elements_node, elements_mesh);

   int ntria = AddTriangles(
                  bld,
                  elements_mesh,
                  buffer,
                  palette_mat,
                  gl_drawable);
   if (ntria == 0)
   {
      cout << "glTF export: no elements found to export!" << endl;
   }
}

void VisualizationSceneScalarData::glTF_ExportMesh(
   glTF_Builder &bld,
   glTF_Builder::buffer_id buffer,
   glTF_Builder::material_id black_mat,
   const gl3::GlDrawable &gl_drawable)
{
   auto lines_node = AddModelNode(bld, "Lines");
   auto lines_mesh = bld.addMesh("Lines Mesh");
   bld.addNodeMesh(lines_node, lines_mesh);

   int nlines = AddLines(
                   bld,
                   lines_mesh,
                   buffer,
                   black_mat,
                   gl_drawable);
   if (nlines == 0)
   {
      cout << "glTF export: no mesh/level lines found to export!" << endl;
   }
}

void VisualizationSceneScalarData::glTF_Export()
{
   cout << "glTF export is not yet implemented for this visualization mode."
        << endl;
}

void VisualizationSceneScalarData::ToggleTexture()
{
   int isSmooth = palette.GetSmoothSetting();
   if (isSmooth)
   {
      palette.UseDiscrete();
      cout << "Texture type : discrete" << endl;
   }
   else
   {
      palette.UseSmooth();
      cout << "Texture type : smooth" << endl;
   }
}

void VisualizationSceneScalarData::SetAutoscale(int _autoscale)
{
   if (autoscale != _autoscale)
   {
      autoscale = _autoscale;
      DoAutoscale(true);
   }
}

VisualizationSceneScalarData::VisualizationSceneScalarData(
   Window &win_, bool init) : VisualizationScene(*win_.wnd), win(win_)
{
   mesh = win.data_state.mesh.get();
   mesh_coarse = win.data_state.mesh_quad.get();
   sol  = &win.data_state.sol;

   if (init)
   {
      Init();
   }
}

void VisualizationSceneScalarData::Init()
{
   vsdata = this;
   window = &win;

   arrow_type = arrow_scaling_type = 0;
   scaling = 0;
   drawaxes = colorbar = 0;
   auto_ref_max = 16;
   auto_ref_min_surf_vert = 100000;
   auto_ref_max_surf_vert = 2000000;
   minv = 0.0;
   maxv = 1.0;
   logscale = false;
   SetLogA();
   PrepareCaption(); // turn on or off the caption

   CuttingPlane = NULL;

   key_r_state = 0;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      wnd->setOnKeyDown('l', KeylPressed);
      wnd->setOnKeyDown('L', KeyLPressed);

      wnd->setOnKeyDown('s', KeySPressed);

      // wnd->setOnKeyDown('a', KeyaPressed);
      wnd->setOnKeyDown('a', Key_Mod_a_Pressed);
      wnd->setOnKeyDown('A', KeyAPressed);

      wnd->setOnKeyDown('r', KeyrPressed);
      wnd->setOnKeyDown('R', KeyRPressed);

      wnd->setOnKeyDown('p', KeypPressed);
      wnd->setOnKeyDown('P', KeyPPressed);

      wnd->setOnKeyDown('h', KeyHPressed);
      wnd->setOnKeyDown('H', KeyHPressed);

      wnd->setOnKeyDown(SDLK_F5, KeyF5Pressed);
      wnd->setOnKeyDown(SDLK_F6, KeyF6Pressed);
      wnd->setOnKeyDown(SDLK_F7, KeyF7Pressed);

      wnd->setOnKeyDown(SDLK_BACKSLASH, KeyBackslashPressed);
      wnd->setOnKeyDown('t', KeyTPressed);
      wnd->setOnKeyDown('T', KeyTPressed);

      wnd->setOnKeyDown('g', KeygPressed);
      wnd->setOnKeyDown('G', KeyGPressed);

      wnd->setOnKeyDown('c', KeycPressed);
      wnd->setOnKeyDown('C', KeyCPressed);

      wnd->setOnKeyDown('k', KeykPressed);
      wnd->setOnKeyDown('K', KeyKPressed);

      wnd->setOnKeyDown(SDLK_F1, KeyF1Pressed);
      wnd->setOnKeyDown(SDLK_F2, KeyF2Pressed);

      wnd->setOnKeyDown(SDLK_COMMA, KeyCommaPressed);
      wnd->setOnKeyDown(SDLK_LESS, KeyLessPressed);
      wnd->setOnKeyDown('~', KeyTildePressed);
      wnd->setOnKeyDown('`', KeyGravePressed);

      wnd->setOnKeyDown(SDLK_EXCLAIM, KeyToggleTexture);
   }

   // Set_Light();

   // glEnable (GL_COLOR_MATERIAL);
   // glShadeModel (GL_SMOOTH);

   // gl->enableLight();
   // gl->enableDepthTest();
   // glEnable(GL_AUTO_NORMAL);
   // glEnable(GL_NORMALIZE);

   if (GetMultisample() > 0)
   {
      glDisable(GL_MULTISAMPLE);
   }
   // add black fog
   // glEnable(GL_FOG);
   // GLfloat fogcol[4] = {0,0,0,1};
   // glFogfv(GL_FOG_COLOR, fogcol);
   // glFogf(GL_FOG_DENSITY,1.0f);

   SetLevelLines(minv, maxv, 15);

   FindNewBox(false);
   ruler_on = 0;
   ruler_x = 0.5 * (bb.x[0] + bb.x[1]);
   ruler_y = 0.5 * (bb.y[0] + bb.y[1]);
   ruler_z = 0.5 * (bb.z[0] + bb.z[1]);

   PrepareRuler();

   autoscale = 1;
}

VisualizationSceneScalarData::~VisualizationSceneScalarData()
{
   delete CuttingPlane;
}

void VisualizationSceneScalarData::SetNewScalingFromBox()
{
   // double eps = 1e-12;
   double eps = 0.0;

   // Find the new scaling
   if (scaling)
   {
      // Scale all sides of the box to 1.
      xscale = yscale = zscale = 1.;
      if ((bb.x[1]-bb.x[0])>eps) { xscale /= (bb.x[1]-bb.x[0]); }
      if ((bb.y[1]-bb.y[0])>eps) { yscale /= (bb.y[1]-bb.y[0]); }
      if ((bb.z[1]-bb.z[0])>eps) { zscale /= (bb.z[1]-bb.z[0]); }
   }
   else
   {
      // Find the largest side of the box in xscale
      xscale = bb.x[1]-bb.x[0];
      yscale = bb.y[1]-bb.y[0];
      zscale = bb.z[1]-bb.z[0];
      if (xscale < yscale) { xscale = yscale; }
      if (xscale < zscale) { xscale = zscale; }
      // Set proportional scaling so that the largest side of the box is 1.
      if (xscale > eps)
      {
         xscale = ( 1.0 / xscale );
      }
      else
      {
         xscale = 1.0;
      }
      zscale = yscale = xscale;
   }
}


void VisualizationSceneScalarData::SetValueRange(double min, double max)
{
   minv = min;
   maxv = max;

   UpdateValueRange(true);
}

void VisualizationSceneScalarData::SetAxisLabels(const char * a_x,
                                                 const char * a_y,
                                                 const char * a_z)
{
   a_label_x = a_x;
   a_label_y = a_y;
   a_label_z = a_z;
   PrepareAxes();
}

void VisualizationSceneScalarData::SetAxisNumberFormat(int precision,
                                                       char format,
                                                       bool showsign)
{
   axis_formatter = NumberFormatter(precision, format, showsign);
   PrepareAxes();
}

void VisualizationSceneScalarData::SetAxisNumberFormat(string formatting)
{
   axis_formatter = NumberFormatter(formatting);
   PrepareAxes();
}

void VisualizationSceneScalarData::PrepareAxes()
{
   // std::array<float, 4> blk = GetLineColor();
   axes_buf.clear();

   gl3::GlBuilder bld = axes_buf.createBuilder();
   if (drawaxes == 3)
   {
      // glLineStipple(1, 255);
      bld.glBegin(GL_LINES);
      bld.glColor3f(1., 0., 0.);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[0]);
      bld.glColor3f(0.75, 0.75, 0.75); // bld.glColor4fv(blk.data());
      bld.glVertex3d(bb.x[1], bb.y[0], bb.z[0]);
      bld.glVertex3d(bb.x[0], bb.y[1], bb.z[0]);
      bld.glColor3f(0., 1., 0.);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[0]);
      bld.glEnd();
      // bld.glColor4fv(blk.data());
      // bld.setUseColor(false);
      // bld.glEnable(GL_LINE_STIPPLE);
      bld.glBegin(GL_LINE_STRIP);
      bld.glColor3f(0.75, 0.75, 0.75);
      bld.glVertex3d(bb.x[1], bb.y[0], bb.z[0]);
      bld.glVertex3d(bb.x[1], bb.y[1], bb.z[0]);
      bld.glVertex3d(bb.x[0], bb.y[1], bb.z[0]);
      bld.glEnd();
   }
   else
   {
      // bld.setUseColor(false);
      bld.glBegin(GL_LINE_LOOP);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[0]);
      bld.glVertex3d(bb.x[1], bb.y[0], bb.z[0]);
      bld.glVertex3d(bb.x[1], bb.y[1], bb.z[0]);
      bld.glVertex3d(bb.x[0], bb.y[1], bb.z[0]);
      bld.glEnd();
   }
   bld.glBegin(GL_LINE_LOOP);
   bld.glVertex3d(bb.x[0], bb.y[0], bb.z[1]);
   bld.glVertex3d(bb.x[1], bb.y[0], bb.z[1]);
   bld.glVertex3d(bb.x[1], bb.y[1], bb.z[1]);
   bld.glVertex3d(bb.x[0], bb.y[1], bb.z[1]);
   bld.glEnd();

   if (drawaxes == 3)
   {
      // bld.setUseColor(true);
      // bld.glDisable(GL_LINE_STIPPLE);
      bld.glBegin(GL_LINES);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[1]);
      bld.glColor3f(0., 0., 1.);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[0]);
      bld.glEnd();
      // bld.setUseColor(false);
      // bld.glEnable(GL_LINE_STIPPLE);
      bld.glColor3f(0.75, 0.75, 0.75);
      bld.glBegin(GL_LINES);
   }
   else
   {
      bld.glBegin(GL_LINES);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[0]);
      bld.glVertex3d(bb.x[0], bb.y[0], bb.z[1]);
   }
   bld.glVertex3d(bb.x[1], bb.y[0], bb.z[0]);
   bld.glVertex3d(bb.x[1], bb.y[0], bb.z[1]);
   bld.glVertex3d(bb.x[1], bb.y[1], bb.z[0]);
   bld.glVertex3d(bb.x[1], bb.y[1], bb.z[1]);
   bld.glVertex3d(bb.x[0], bb.y[1], bb.z[0]);
   bld.glVertex3d(bb.x[0], bb.y[1], bb.z[1]);
   bld.glEnd();
   if (drawaxes == 3)
   {
      // bld.glDisable(GL_LINE_STIPPLE);
   }

   // Write the coordinates of the lower left and upper right corner.
   //   glEnable (GL_COLOR_MATERIAL);
   //   GLfloat textcol[3] = { 0, 0, 0 };
   //   glColor3fv (textcol);

   if (drawaxes == 1)
   {
      int desc = GetFont()->getFontDescender();
      int ox = -desc/2;
      int oy = -3*desc/2;
      ostringstream buf;
      buf << "("
          << axis_formatter(bb.x[0]) << ","
          << axis_formatter(bb.y[0]) << ","
          << axis_formatter(bb.z[0]) << ")";
      axes_buf.addText(bb.x[0], bb.y[0], bb.z[0], ox, oy, buf.str());

      ostringstream buf1;
      buf1 << "("
           << axis_formatter(bb.x[1]) << ","
           << axis_formatter(bb.y[1]) << ","
           << axis_formatter(bb.z[1]) << ")";
      axes_buf.addText(bb.x[1], bb.y[1], bb.z[1], ox, oy, buf1.str());
   }
   updated_bufs.emplace_back(&axes_buf);

   constexpr float len = 1.2f;
   constexpr float l = .9f, cl = .27f, cb = l-cl;
   coord_cross_buf.clear();
   coord_cross_buf.addLines<gl3::Vertex>(
   {
      {0, 0, 0}, {0, 0, l},
      {0, 0, 0}, {0, l, 0},
      {0, 0, 0}, {l, 0, 0}
   });
   coord_cross_buf.addCone(0,0,cb,0,0,1, cl);
   coord_cross_buf.addCone(0,cb,0,0,1,0, cl);
   coord_cross_buf.addCone(cb,0,0,1,0,0, cl);
   coord_cross_buf.addText(len, 0.0f, 0.0f, a_label_x);
   coord_cross_buf.addText(0.0f, len, 0.0f, a_label_y);
   coord_cross_buf.addText(0.0f, 0.0f, len, a_label_z);
   updated_bufs.emplace_back(&coord_cross_buf);
}

void VisualizationSceneScalarData::DrawPolygonLevelLines(
   gl3::GlBuilder& builder, double * point, int n, Array<double> &mesh_level,
   bool log_vals)
{
   int l, k, k1;
   double curve, t;
   double p[3];

   for (l = 0; l < mesh_level.Size(); l++)
   {
      // Using GL_LINE_STRIP (explicitly closed for more than 2 points)
      // should produce the same result, however visually the level lines
      // have discontinuities. Using GL_LINE_LOOP does not have that problem.
      builder.glBegin(GL_LINE_LOOP);
      curve = LogVal(mesh_level[l], log_vals);
      for (k = 0; k < n; k++)
      {
         k1 = (k+1)%n;
         if ( (curve <=point[4*k+3] && curve >= point[4*k1+3]) ||
              (curve >=point[4*k+3] && curve <= point[4*k1+3]) )
         {
            if ((curve - point[4*k1+3]) == 0.)
            {
               t = 1.;
            }
            else if ((curve - point[4*k+3]) == 0.)
            {
               t = 0.;
            }
            else
            {
               t = (curve - point[4*k+3]) / (point[4*k1+3]-point[4*k+3]);
            }
            p[0] = (1.0-t)*point[4*k+0]+t*point[4*k1+0];
            p[1] = (1.0-t)*point[4*k+1]+t*point[4*k1+1];
            p[2] = (1.0-t)*point[4*k+2]+t*point[4*k1+2];
            builder.glVertex3d(p[0], p[1], p[2]);
         }
      }
      builder.glEnd();
   }
}

void VisualizationSceneScalarData::SetLevelLines (
   double min, double max, int n, int adj)
{
   int i;
   double t, eps;

   if (min < minv)
   {
      min = minv;
      cout << "min set to minv : " << min << endl;
   }
   if (max > maxv)
   {
      max = maxv;
      cout << "max set to maxv : " << max << endl;
   }

   nl = n;
   level.SetSize(nl+1);
   for (i = 0; i <= nl; i++)
   {
      t = (double) i / nl;
      level[i] = min * (1.0 - t) + t * max;
   }

   if (adj)
   {
      eps = 1.0E-5;
      level[0]  = level[0]  * (1.0 - eps) + level[1]    * eps;
      level[nl] = level[nl] * (1.0 - eps) + level[nl-1] * eps;
   }

   if (logscale)
   {
      for (i = 0; i <= nl; i++)
      {
         level[i] = _ULogVal((level[i] - minv) / (maxv - minv));
      }
   }
}

void VisualizationSceneScalarData::PrintState()
{
   cout << "\nkeys: " << wnd->getSavedKeys() << "\n"
        << "\nlight " << strings_off_on[use_light ? 1 : 0]
        << "\nperspective " << strings_off_on[OrthogonalProjection ? 0 : 1]
        << "\nviewcenter " << ViewCenterX << ' ' << ViewCenterY
        << "\nzoom " << (OrthogonalProjection ? ViewScale :
                         tan(M_PI / 8.) / tan(ViewAngle * (M_PI / 360.0)))
        << "\nvaluerange " << minv << ' ' << maxv;
   const float *r = glm::value_ptr(rotmat);
   ios::fmtflags fmt = cout.flags();
   cout << fixed << showpos
        << "\nrotmat " << r[ 0] << ' ' << r[ 1] << ' ' << r[ 2] << ' ' << r[ 3]
        << "\n       " << r[ 4] << ' ' << r[ 5] << ' ' << r[ 6] << ' ' << r[ 7]
        << "\n       " << r[ 8] << ' ' << r[ 9] << ' ' << r[10] << ' ' << r[11]
        << "\n       " << r[12] << ' ' << r[13] << ' ' << r[14] << ' ' << r[15]
        << '\n' << endl;
   cam.Print();
   cout.flags(fmt);
   mesh->PrintInfo(cout);
}

void VisualizationSceneScalarData::ShrinkPoints(DenseMatrix &pointmat,
                                                int i, int fn, int di)
{
   int dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();

   if (shrink != 1.0)
   {
      if (dim == 2)
      {
         int d, k;

         for (d = 0; d < sdim; d++)
         {
            double cd = 0.0;
            for (k = 0; k < pointmat.Width(); k++)
            {
               cd += pointmat(d,k);
            }
            cd /= pointmat.Width();

            for (k = 0; k < pointmat.Width(); k++)
            {
               pointmat(d,k) = shrink*pointmat(d,k) + (1-shrink)*cd;
            }
         }
      }
      else
      {
         int attr = mesh->GetBdrAttribute(i);
         for (int k = 0; k < pointmat.Width(); k++)
            for (int d = 0; d < sdim; d++)
            {
               pointmat(d,k) = shrink*pointmat(d,k) + (1-shrink)*bdrc(d,attr-1);
            }
      }
   }

   if (shrinkmat != 1.0)
   {
      int attr, elem1, elem2;
      if (dim == 2 || sdim == 2)
      {
         attr = mesh->GetAttribute(i);
      }
      else
      {
         mesh->GetFaceElements(fn, &elem1, &elem2);
         if (di == 0)
         {
            attr = mesh->GetAttribute(elem1);
         }
         else
         {
            attr = mesh->GetAttribute(elem2);
         }
      }

      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < pointmat.Height(); d++)
         {
            pointmat(d,k) = shrinkmat*pointmat(d,k) + (1-shrinkmat)*matc(d,attr-1);
         }
   }
}

void VisualizationSceneScalarData::ComputeBdrAttrCenter()
{
   DenseMatrix pointmat;
   Vector nbdrc(mesh->bdr_attributes.Max());
   int sdim = mesh->SpaceDimension();

   bdrc.SetSize(sdim,mesh->bdr_attributes.Max());
   bdrc = 0.0;
   nbdrc = 0.0;

   for (int i = 0; i < mesh -> GetNBE(); i++)
   {
      mesh->GetBdrPointMatrix(i, pointmat);
      nbdrc(mesh->GetBdrAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < sdim; d++)
         {
            bdrc(d,mesh->GetBdrAttribute(i)-1) += pointmat(d,k);
         }
   }

   for (int i = 0; i < mesh->bdr_attributes.Max(); i++)
      if (nbdrc(i) != 0)
         for (int d = 0; d < sdim; d++)
         {
            bdrc(d,i) /= nbdrc(i);
         }
}

void VisualizationSceneScalarData::ComputeElemAttrCenter()
{
   DenseMatrix pointmat;
   Vector nmatc(mesh->attributes.Max());
   int sdim = mesh->SpaceDimension();

   matc.SetSize(sdim,mesh->attributes.Max());
   matc = 0.0;
   nmatc = 0.0;

   for (int i = 0; i < mesh -> GetNE(); i++)
   {
      mesh->GetPointMatrix(i, pointmat);
      nmatc(mesh->GetAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < sdim; d++)
         {
            matc(d,mesh->GetAttribute(i)-1) += pointmat(d,k);
         }
   }

   for (int i = 0; i < mesh->attributes.Max(); i++)
      if (nmatc(i) != 0)
         for (int d = 0; d < sdim; d++)
         {
            matc(d,i) /= nmatc(i);
         }
}


Plane::Plane(double A,double B,double C,double D)
{
   eqn[0] = A;
   eqn[1] = B;
   eqn[2] = C;
   eqn[3] = D;

   CartesianToSpherical();

   double x[2] = {vsdata -> bb.x[0], vsdata -> bb.x[1]};
   double y[2] = {vsdata -> bb.y[0], vsdata -> bb.y[1]};
   double z[2] = {vsdata -> bb.z[0], vsdata -> bb.z[1]};
   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );

   x0 = (x[0]+x[1])/2.0;
   y0 = (y[0]+y[1])/2.0;
   z0 = (z[0]+z[1])/2.0;

   phi_step = M_PI / 36;
   theta_step = M_PI / 36;
   rho_step = bbox_diam / 200;
}

void Plane::CartesianToSpherical()
{
   rho = sqrt(eqn[0]*eqn[0]+eqn[1]*eqn[1]+eqn[2]*eqn[2]);
   phi = asin(eqn[2]/rho);
   theta = atan2(eqn[1], eqn[0]);
}

void Plane::SphericalToCartesian()
{
   eqn[0] = rho * cos(phi) * cos(theta);
   eqn[1] = rho * cos(phi) * sin(theta);
   eqn[2] = rho * sin(phi);
   eqn[3] = - (eqn[0] * x0 + eqn[1] * y0 + eqn[2] * z0);
}

void Plane::IncreasePhi()
{
   phi += phi_step;
   SphericalToCartesian();
}

void Plane::DecreasePhi()
{
   phi -= phi_step;
   SphericalToCartesian();
}

void Plane::IncreaseTheta()
{
   theta += theta_step;
   SphericalToCartesian();
}

void Plane::DecreaseTheta()
{
   theta -= theta_step;
   SphericalToCartesian();
}

void Plane::IncreaseDistance()
{
   double k = (rho_step) / (rho*rho);
   x0 -= eqn[0] * k;
   y0 -= eqn[1] * k;
   z0 -= eqn[2] * k;
   eqn[3] += rho_step;
   CartesianToSpherical();
}

void Plane::DecreaseDistance()
{
   double k = (rho_step) / (rho*rho);
   x0 += eqn[0] * k;
   y0 += eqn[1] * k;
   z0 += eqn[2] * k;
   eqn[3] -= rho_step;
   CartesianToSpherical();
}
