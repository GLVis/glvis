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

#ifndef GLVIS_VSDATA
#define GLVIS_VSDATA

#include <array>

#include "openglvis.hpp"
#include "mfem.hpp"
#include "aux_gl3.hpp"
using namespace mfem;

extern std::string plot_caption; // defined in glvis.cpp
extern std::string extra_caption; // defined in glvis.cpp

class Plane
{
private:
   double eqn[4];
   double phi, theta, rho;
   double x0,y0,z0;
   void CartesianToSpherical();
   void SphericalToCartesian();

   double bbox_diam;

   double phi_step, theta_step, rho_step;

public:
   Plane(double A,double B,double C,double D);
   inline double * Equation() { return eqn; }
   inline double Transform(double x, double y, double z)
   { return eqn[0]*x+eqn[1]*y+eqn[2]*z+eqn[3]; }
   inline double Transform(double * x)
   { return eqn[0]*x[0]+eqn[1]*x[1]+eqn[2]*x[2]+eqn[3]; }

   void IncreasePhi();
   void DecreasePhi();
   void IncreaseTheta();
   void DecreaseTheta();
   void IncreaseDistance();
   void DecreaseDistance();
};


class VisualizationSceneScalarData : public VisualizationScene
{
protected:
   Mesh   *mesh;
   Vector *sol;

   double minv, maxv;

   std::string a_label_x, a_label_y, a_label_z;

   int scaling, colorbar, drawaxes, axeslist;
   int auto_ref_max, auto_ref_max_surf_elem;

   gl3::GlDrawable axes_buf;
   bool coord_cross_init = false,
        color_bar_init = false;
   gl3::GlDrawable coord_cross, color_bar;

   void Init();

   int arrow_type, arrow_scaling_type;

   int nl;
   Array<double> level;

   int ruler_on;
   double ruler_x, ruler_y, ruler_z;

   // autoscale controls the behavior when the mesh/solution are updated:
   // 0 - do not change the bounding box and the value range
   // 1 - recompute both the bounding box and the value range (default)
   // 2 - recompute only the value range
   // 3 - recompute only the bounding box
   int autoscale;

   bool logscale;

   bool LogscaleRange() { return (minv > 0.0 && maxv > minv); }
   void PrintLogscale(bool warn);

   double log_a, unit_a;
   void SetLogA()
   {
      if (logscale)
      {
         unit_a = 1.0/log(maxv/minv), log_a = (maxv - minv)*unit_a;
      }
      else
      {
         unit_a = 1.0/(maxv - minv), log_a = 1.0;
      }
   }
   double _ULogVal(const double &u) { return minv*pow(maxv/minv, u); }
   double ULogVal(const double &u)
   { return (logscale ? _ULogVal(u) : minv + (maxv - minv)*u); }
   double LogUVal(const double &z)
   {
      return ((logscale && z >= minv && z <= maxv) ?
              (log(z/minv)*unit_a) : (z - minv)*unit_a);
   }
   double _LogVal_(const double &z) { return (log(z/minv)*log_a + minv); }
   double _LogVal(const double &z)
   { return ((z >= minv && z <= maxv) ? _LogVal_(z) : (z)); }
   double LogVal(const double &z, const bool &log_val)
   { return (log_val ? _LogVal(z) : z); }
   double LogVal(const double &z) { return LogVal(z, logscale); }

   void FixValueRange();

public:
   Plane *CuttingPlane;
   int light;
   int key_r_state;
   /** Shrink factor with respect to the center of each element (2D) or the
       center of each boundary attribute (3D) */
   double shrink;
   /// Shrink factor with respect to the element (material) attributes centers
   double shrinkmat;

   VisualizationSceneScalarData()
      : a_label_x("x"), a_label_y("y"), a_label_z("z") {}
   VisualizationSceneScalarData (Mesh & m, Vector & s);

   virtual ~VisualizationSceneScalarData();

   // Determine 'xscale', 'yscale', and 'zscale' using the current bounding
   // box, depending on the value of 'scaling'.
   virtual void SetNewScalingFromBox();

   // Compute the bounding box, call UpdateBoundingBox.
   // In 2D the z range is the value range, so FixValueRange and
   // UpdateValueRange are also called.
   virtual void FindNewBox(bool prepare) = 0;

   // Compute the value range based on the current solution, adjust it by
   // calling FixValueRange and then call UpdateValueRange.
   virtual void FindNewValueRange(bool prepare) = 0;

   // Redefined in 2D to call just FindNewBox
   virtual void FindNewBoxAndValueRange(bool prepare)
   { FindNewBox(prepare); FindNewValueRange(prepare); }

   // Redefined in 2D to update only the x- and y-ranges.
   virtual void FindMeshBox(bool prepare) { FindNewBox(prepare); }

   // Perform autoscaling depending on the value of 'autoscale':
   // 0 - do nothing
   // 1 - call call FindNewBoxAndValueRange
   // 2 - call FindNewValueRange
   // 3 - call FindMeshBox
   void DoAutoscale(bool prepare);
   // Similar to the above but force recomputation of the value range
   void DoAutoscaleValue(bool prepare);

   virtual void Prepare() = 0;
   virtual void PrepareLines() = 0;

   void UpdateBoundingBox() { SetNewScalingFromBox(); PrepareAxes(); }
   virtual void EventUpdateColors() { Prepare(); }
   virtual void UpdateLevelLines() = 0;
   virtual void UpdateValueRange(bool prepare) = 0;
   void SetValueRange(double, double);

   virtual void SetShading(int, bool) = 0;
   virtual void SetRefineFactors(int, int) = 0;
   void SetAutoRefineLimits(int max_ref, int max_surf_elem)
   {
      auto_ref_max = max_ref;
      auto_ref_max_surf_elem = max_surf_elem;
   }
   virtual void AutoRefine() = 0;
   virtual void ToggleAttributes(Array<int> &attr_list) = 0;

   virtual void PrintState();

   Mesh *GetMesh() { return mesh; }

   void DrawColorBar(double minval, double maxval,
                     Array<double> * level = NULL,
                     Array<double> * levels = NULL);
   void DrawCaption();
   void DrawCoordinateCross();

   double &GetMinV() { return minv; }
   double &GetMaxV() { return maxv; }

   void SetLevelLines(double min, double max, int n, int adj = 1);

   void Arrow(gl3::GlBuilder& builder,
              double px, double py, double pz,
              double vx, double vy, double vz, double length,
              double cone_scale = 0.075);
   void Arrow2(double px, double py, double pz,
               double vx, double vy, double vz,
               double length,
               double cone_scale = 0.075);
   void Arrow3(double px, double py, double pz,
               double vx, double vy, double vz,
               double length,
               double cone_scale = 0.075);

   void DrawPolygonLevelLines(double *point, int n, Array<double> &level,
                              bool log_vals, gl3::GlBuilder& builder);

   void DrawPolygonLevelLines(double *point, int n, Array<double> &level,
                              bool log_vals);

   void ToggleLight() { light = !light; }

   void ToggleDrawColorbar()
   {
      // colorbar states are: 0) no colorbar, no caption; 1) colorbar with
      // caption; 2) colorbar without caption.
      static const int next[2][3] = { { 1, 2, 0 }, { 2, 0, 0 } };
      colorbar = next[plot_caption.empty()][colorbar];
   }

   // Turn on or off the caption
   void UpdateCaption()
   {
      bool empty = plot_caption.empty();
      colorbar = (colorbar ? empty+1 : !empty);
   }

   void SetAxisLabels(const char * a_x, const char * a_y, const char * a_z);

   void PrepareAxes();
   void ToggleDrawAxes()
   {
      drawaxes = (drawaxes+1)%4;
      if (drawaxes)
      {
         PrepareAxes();
      }
   }

   void ToggleScaling()
   { scaling = !scaling; SetNewScalingFromBox(); }

   virtual void ToggleLogscale(bool print);

   void ToggleRuler();
   void RulerPosition();
   void DrawRuler(bool log_z = false);

   void ToggleTexture();

   void SetAutoscale(int _autoscale);
   int GetAutoscale() const { return autoscale; }

   /// Shrink the set of points towards attributes centers of gravity
   void ShrinkPoints(DenseMatrix &pointmat, int i, int fn, int di);
   // Centers of gravity based on the boundary/element attributes
   DenseMatrix bdrc, matc;
   /// Compute the center of gravity for each boundary attribute
   void ComputeBdrAttrCenter();
   /// Compute the center of gravity for each element attribute
   void ComputeElemAttrCenter();
};

#endif
