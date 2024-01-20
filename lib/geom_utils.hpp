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

#ifndef GLVIS_GEOM_UTILS_HPP
#define GLVIS_GEOM_UTILS_HPP
#include "mfem.hpp"

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

inline int Normalize(mfem::DenseMatrix &normals)
{
   int err = 0;
   for (int i = 0; i < normals.Width(); i++)
   {
      err += Normalize(&normals(0, i));
   }
   return err;
}
// 'v' must define three vectors that add up to the zero vector
// that is v[0], v[1], and v[2] are the sides of a triangle
inline int UnitCrossProd(double v[][3], double nor[])
{
   // normalize the three vectors
   for (int i = 0; i < 3; i++)
      if (Normalize(v[i]))
      {
         return 1;
      }

   // find the pair that forms an angle closest to pi/2, i.e. having
   // the longest cross product:
   double cp[3][3], max_a = 0.;
   int k = 0;
   for (int i = 0; i < 3; i++)
   {
      CrossProd(v[(i+1)%3], v[(i+2)%3], cp[i]);
      double a = sqrt(InnerProd(cp[i], cp[i]));
      if (max_a < a)
      {
         max_a = a, k = i;
      }
   }
   if (max_a == 0.)
   {
      return 1;
   }
   for (int i = 0; i < 3; i++)
   {
      nor[i] = cp[k][i] / max_a;
   }

   return 0;
}

inline int Compute3DUnitNormal(const double p1[], const double p2[],
                               const double p3[], double nor[])
{
   double v[3][3];

   for (int i = 0; i < 3; i++)
   {
      v[0][i] = p2[i] - p1[i];
      v[1][i] = p3[i] - p2[i];
      v[2][i] = p1[i] - p3[i];
   }

   return UnitCrossProd(v, nor);
}

inline int Compute3DUnitNormal (const double p1[], const double p2[],
                                const double p3[], const double p4[], double nor[])
{
   double v[3][3];

   for (int i = 0; i < 3; i++)
   {
      // cross product of the two diagonals:
      /*
        v[0][i] = p3[i] - p1[i];
        v[1][i] = p4[i] - p2[i];
        v[2][i] = (p1[i] + p2[i]) - (p3[i] + p4[i]);
      */

      // cross product of the two vectors connecting the midpoints of the
      // two pairs of opposing sides; this gives the normal vector in the
      // midpoint of the quad:
      v[0][i] = 0.5 * ((p2[i] + p3[i]) - (p1[i] + p4[i]));
      v[1][i] = 0.5 * ((p4[i] + p3[i]) - (p1[i] + p2[i]));
      v[2][i] = p1[i] - p3[i];
   }

   return UnitCrossProd(v, nor);
}

inline int ProjectVector(double v[], const double n[])
{
   // project 'v' on the plane with normal given by 'n' and  then normalize 'v'
   LinearCombination(InnerProd(n, n), v, -InnerProd(v, n), n, v);
   return Normalize(v);
}

#endif // GLVIS_GEOM_UTILS_HPP
