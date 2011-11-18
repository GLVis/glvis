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

#ifndef GLVIS_VSDATA
#define GLVIS_VSDATA

#include "openglvis.hpp"
#include "mfem.hpp"


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
   Mesh * mesh;
   Vector * sol;

   double minv, maxv;

   int scaling, colorbar, drawaxes, axeslist;

   void Init();

   int arrow_type, arrow_scaling_type;

   int nl;
   Array<double> level;

   int ruler_on;
   double ruler_x, ruler_y, ruler_z;

public:
   Plane * CuttingPlane;
   int light;
   int key_r_state;
   /** Shrink factor with respect to the center of each element (2D) or the
       center of each boundary attribute (3D) */
   double shrink;
   /// Shrink factor with respect to the element (material) attributes centers
   double shrinkmat;

   VisualizationSceneScalarData() {}
   VisualizationSceneScalarData (Mesh & m, Vector & s);

   virtual ~VisualizationSceneScalarData();

   virtual void SetNewScalingFromBox();
   virtual void CenterObject();
   virtual void CenterObject2D();
   virtual void ResetScaling();

   virtual void FindNewBox() = 0;

   virtual void PrepareAxes();
   virtual void PrepareLines() = 0;

   virtual void EventUpdateColors () { Prepare(); };
   virtual void UpdateLevelLines() = 0;
   virtual void UpdateValueRange() = 0;
   void SetValueRange(double, double);

   virtual void SetShading(int) = 0;
   virtual void SetRefineFactors(int, int) = 0;
   virtual void ToggleAttributes(Array<int> &attr_list) = 0;

   virtual void PrintState();

   Mesh * GetMesh() { return mesh; }

   void DrawColorBar (double minval, double maxval,
                      Array<double> * level = NULL,
                      Array<double> * levels = NULL);
   void DrawCoordinateCross ();

   double &GetMinV() { return minv; };
   double &GetMaxV() { return maxv; };

   void SetLevelLines (double min, double max, int n, int adj = 1);

   void Arrow(double px, double py, double pz,
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

   void DrawPolygonLevelLines(double * point, int n, Array<double> & level);

   void ToggleLight ()
   { light = !light; }

   void ToggleDrawColorbar ()
   { colorbar = !colorbar; }

   void ToggleDrawAxes ()
   {
      drawaxes = (drawaxes+1)%4;
      if (drawaxes)
         PrepareAxes();
   }

   void ToggleScaling (){
      scaling = !scaling;
      ResetScaling();
   }

   void ToggleRuler();
   void RulerPosition();
   void DrawRuler();

   void ToggleTexture();

   /// Shrink the set of points towards attributes centers of gravity
   void ShrinkPoints(DenseMatrix &pointmat, int i, int fn, int fo);
   // Centers of gravity based on the bounday/element attributes
   DenseMatrix bdrc, matc;
   /// Compute the center of gravity for each boundary attribute
   void ComputeBdrAttrCenter();
   /// Compute the center of gravity for each element attribute
   void ComputeElemAttrCenter();
};

#endif
