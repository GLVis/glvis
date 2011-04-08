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

#include <stdlib.h>
#include <math.h>
#include <X11/keysym.h>

#include <sstream>
using namespace std;

#include "vsdata.hpp"
#include "aux_vis.hpp"
#include "material.hpp"

#include "gl2ps.h"


// Draw an arrow starting at point (px, py, pz) with orientation (vx, vy, vz)
// and length "length".
void VisualizationSceneScalarData::Arrow3(double px, double py, double pz,
                                          double vx, double vy, double vz,
                                          double length,
                                          double cone_scale)
{
   glPushMatrix();

   double xc = 0.5*(x[0]+x[1]);
   double yc = 0.5*(y[0]+y[1]);
   double zc = 0.5*(z[0]+z[1]);
   glTranslated(xc, yc, zc);
   glScaled(1.0/xscale, 1.0/yscale, 1.0/zscale);

   double rlen = length/sqrt(vx*vx+vy*vy+vz*vz);
   px = (px-xc)*xscale;
   py = (py-yc)*yscale;
   pz = (pz-zc)*zscale;
   vx *= rlen*xscale;
   vy *= rlen*yscale;
   vz *= rlen*zscale;

   if (arrow_scaling_type == 0)
      length = sqrt(vx*vx+vy*vy+vz*vz);

   glTranslated(px, py, pz);

   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   double phi   = acos(vz/rhos);
   double theta;
   theta = atan2 (vy, vx);

   glRotatef(theta*180/M_PI, 0.0f, 0.0f, 1.0f);
   glRotatef(phi*180/M_PI, 0.0f, 1.0f, 0.0f);

   glScaled (length, length, length);

   if (arrow_type == 1)
      glTranslated (0, 0, -0.5);

   glBegin(GL_LINES);
   glVertex3d(0, 0, 0);
   glVertex3d(0, 0, 1);
   glEnd();

   glTranslated (0, 0, 1);
   glScaled (cone_scale, cone_scale, cone_scale);

   Cone();

   glPopMatrix();
}

void VisualizationSceneScalarData::Arrow2(double px, double py, double pz,
                                          double vx, double vy, double vz,
                                          double length,
                                          double cone_scale)
{
   glPushMatrix();
   glTranslated(px, py, pz);

   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   double phi   = acos(vz/rhos);
   double theta;
   theta = atan2 (vy, vx);

   glRotatef(theta*180/M_PI, 0.0f, 0.0f, 1.0f);
   glRotatef(phi*180/M_PI, 0.0f, 1.0f, 0.0f);

   glScaled (length, length, length);

   glBegin(GL_LINES);
   glVertex3d(0, 0, 0);
   glVertex3d(0, 0, 1);
   glEnd();

   glTranslated (0, 0, 1);
   glScaled (cone_scale, cone_scale, cone_scale);

   Cone();

   glPopMatrix();
}

void VisualizationSceneScalarData::Arrow(double px, double py, double pz,
                                         double vx, double vy, double vz,
                                         double length,
                                         double cone_scale)
{
   double rhos = sqrt (vx*vx+vy*vy+vz*vz);
   if (rhos == 0.0)
      return;
   double phi = acos(vz/rhos), theta = atan2(vy, vx);
   const int n = 8;
   const double step = 2*M_PI/n, nz = (1.0/4.0);
   double point = step, cone[n+4][3], normal[n+2][3];
   int i, j, k;

   cone[0][0] = 0;          cone[0][1] = 0; cone[0][2] = 1;
   cone[1][0] = cone_scale; cone[1][1] = 0; cone[1][2] = -4*cone_scale + 1;
   normal[0][0] = 0.0/cone_scale;
   normal[0][1] = 0.0/cone_scale;
   normal[0][2] = 1.0/cone_scale;
   normal[1][0] = 1.0/cone_scale;
   normal[1][1] = 0.0/cone_scale;
   normal[1][2] = nz/cone_scale;

   for(i=2; i<n+1; i++){
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
      length = rhos;

   // double xc = 0.5*(x[0]+x[1]), yc = 0.5*(y[0]+y[1]), zc = 0.5*(z[0]+z[1]);
   double coord[3];
   // double rlen = length/rhos;

   // px = (px-xc)*xscale;  py = (py-yc)*yscale;  pz = (pz-zc)*zscale;
   // vx *= rlen*xscale;    vy *= rlen*yscale;    vz *= rlen*zscale;

   if (arrow_type == 1)
      for(i=0; i<n+4; i++)
         cone[i][2] -= 0.5;

   double M[3][3]={{cos(theta)*cos(phi), -sin(theta),  cos(theta)*sin(phi)},
                   {sin(theta)*cos(phi),  cos(theta),  sin(theta)*sin(phi)},
                   {          -sin(phi),          0.,             cos(phi)}};
   double v[3] = { M[0][2]/xscale, M[1][2]/yscale, M[2][2]/zscale };
   length /= sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

   for(i=0; i<n+4; i++){
      for(j=0; j<3; j++)
         coord[j] = cone[i][j] * length;

      for(k=0; k<3; k++){
         cone[i][k] = 0.;
         for(j=0; j<3; j++)
            cone[i][k] += M[k][j] * coord[j];
      }
      // cone[i][0] = (cone[i][0] + px)/xscale + xc;
      // cone[i][1] = (cone[i][1] + py)/yscale + yc;
      // cone[i][2] = (cone[i][2] + pz)/zscale + zc;
      cone[i][0] = cone[i][0]/xscale + px;
      cone[i][1] = cone[i][1]/yscale + py;
      cone[i][2] = cone[i][2]/zscale + pz;
   }

   for(i=0; i<=n+1; i++){
      for(j=0; j<3; j++)
         coord[j] = normal[i][j];

      for(k=0; k<3; k++){
         normal[i][k] = 0.;
         for(j=0; j<3; j++)
            normal[i][k] += M[k][j] * coord[j];
      }
      normal[i][0] *= xscale;
      normal[i][1] *= yscale;
      normal[i][2] *= zscale;
   }

   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
   glBegin(GL_TRIANGLE_FAN);
   for(i=0; i<=n+1; i++){
      glNormal3dv(normal[i]);
      glVertex3dv(cone[i]);
   }
   glEnd();

   glBegin(GL_LINES);
   glVertex3dv(cone[n+2]);
   glVertex3dv(cone[n+3]);
   glEnd();
}

void VisualizationSceneScalarData::DrawColorBar (double minval, double maxval,
                                                 Array<double> *level,
                                                 Array<double> *levels)
{
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   int i;

   double miny;
   double maxy;
   double minx;
   double maxx;
   double posz = -1.5;

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

   double Y;

   // glEnable (GL_COLOR_MATERIAL);
   glNormal3d (0, 0, 1);

   Set_Material();

   double eps = 1e-24;

   if (GetUseTexture())
   {
      glEnable (GL_TEXTURE_1D);
      glColor4d(1, 1, 1, 1);
      glBegin(GL_QUADS);
      glTexCoord1d(0.);
      glVertex3d(minx, miny, posz);
      glVertex3d(maxx, miny, posz);
      glTexCoord1d(1.);
      glVertex3d(maxx, maxy, posz);
      glVertex3d(minx, maxy, posz);
      glEnd();
      glDisable (GL_TEXTURE_1D);
   }
   else
   {
      const int nquads = 256;
      glBegin(GL_QUAD_STRIP);
      for (int i = 0; i <= nquads; i++)
      {
         const double a = double(i) / nquads;
         Y = (1.0 - a) * miny + a * maxy;
         if (MySetColorLogscale)
         {
            MySetColorLogscale = 0;
            MySetColor (Y,miny,maxy);
            MySetColorLogscale = 1;
         }
         else
            MySetColor (Y,miny,maxy);
         glVertex3d (minx,Y,posz);
         glVertex3d (maxx,Y,posz);
      }
      glEnd();
   }

   Set_Black_Material();
   glBegin(GL_LINES);
   if (levels) {
      for (i = 0; i < levels->Size(); i++) {
         if (MySetColorLogscale)
            Y = miny + (maxy - miny) *
               log(fabs((*levels)[i]/(minv+eps)))/(log(fabs(maxv/(minv+eps)))+eps);
         else
            Y = miny + (maxy - miny)*((*levels)[i]-minv)/(maxv-minv);
         glVertex3d (minx,Y,posz);
         glVertex3d (maxx,Y,posz);
      }
   }
   if (level) {
      for (i = 0; i < level->Size(); i++) {
         if (MySetColorLogscale)
            Y = miny + (maxy - miny)* log(fabs((*level)[i]/(minv+eps))) /
               (log(fabs(maxv/(minv+eps)))+eps);
         else
            Y = miny + (maxy - miny)*((*level)[i]-minv)/(maxv-minv);
         glVertex3d (minx,Y,posz);
         glVertex3d (maxx,Y,posz);
      }
   }
   glEnd();

   // GLfloat textcol[3] = {0,0,0};
   // glColor3fv (textcol);

   glPushAttrib (GL_LIST_BIT);
   extern GLuint fontbase;
   glListBase (fontbase);

   double val;
   ostringstream * buf;

   if (!level)
   {
      for (i = 0; i <= 4; i++)
      {
         Y = miny + i * (maxy-miny) / 4;
         glRasterPos3d (maxx+0.02,Y,posz);

         if (MySetColorLogscale)
            val = minval * pow (maxval / minval, i / 4.0);
         else
            val = minval + i * (maxval-minval) / 4;

         buf = new ostringstream;
         (*buf) << val;
         glCallLists (buf->str().size(), GL_UNSIGNED_BYTE, buf->str().c_str());
         delete buf;
      }
   }
   else
   {
      for (i = 0; i < level->Size(); i++)
      {
         if (MySetColorLogscale)
            Y = miny + (maxy - miny)* log(fabs((*level)[i]/(minv+eps))) /
               (log(fabs(maxv/(minv+eps)))+eps);
         else
            Y = miny + (maxy - miny)*((*level)[i]-minv)/(maxv-minv);
         glRasterPos3d (maxx+0.02,Y,posz);

         val = (*level)[i];

         buf = new ostringstream;
         (*buf) << val;
         glCallLists (buf->str().size(), GL_UNSIGNED_BYTE, buf->str().c_str());
         delete buf;
      }
   }

   if (levels)
   {
      for (i = 0; i < levels->Size(); i++)
      {
         if (MySetColorLogscale)
            Y = miny + (maxy - miny)* log(fabs((*levels)[i]/(minv+eps))) /
               (log(fabs(maxv/(minv+eps)))+eps);
         else
            Y = miny + (maxy - miny)*((*levels)[i]-minv)/(maxv-minv);
         glRasterPos3d (maxx+0.02,Y,posz);

         val = (*levels)[i];

         buf = new ostringstream;
         (*buf) << val;
         glCallLists (buf->str().size(), GL_UNSIGNED_BYTE, buf->str().c_str());
         delete buf;
      }
   }

   glPopAttrib();

   if (print)
   {
      if (!level)
      {
         for (i = 0; i <= 4; i++)
         {
            Y = miny + i * (maxy-miny) / 4;
            glRasterPos3d (maxx+0.02,Y,posz);

            if (MySetColorLogscale)
               val = minval * pow (maxval / minval, i / 4.0);
            else
               val = minval + i * (maxval-minval) / 4;

            buf = new ostringstream;
            (*buf) << val;
            gl2psText(buf->str().c_str(),"Times",8);
            delete buf;
         }
      }
      else
      {
         for (i = 0; i < level->Size(); i++)
         {
            if (MySetColorLogscale)
               Y = miny + (maxy - miny)* log(fabs((*level)[i]/(minv+eps))) /
                  (log(fabs(maxv/(minv+eps)))+eps);
            else
               Y = miny + (maxy - miny)*((*level)[i]-minv)/(maxv-minv);
            glRasterPos3d (maxx+0.02,Y,posz);
            glRasterPos3d (maxx+0.02,Y,posz);

            val = (*level)[i];

            buf = new ostringstream;
            (*buf) << val;
            gl2psText(buf->str().c_str(),"Times",8);
            delete buf;
         }
      }

      if (levels)
      {
         for (i = 0; i < levels->Size(); i++)
         {
            if (MySetColorLogscale)
               Y = miny + (maxy - miny)* log(fabs((*levels)[i]/(minv+eps))) /
                  (log(fabs(maxv/(minv+eps)))+eps);
            else
               Y = miny + (maxy - miny)*((*levels)[i]-minv)/(maxv-minv);
            glRasterPos3d (maxx+0.02,Y,posz);

            val = (*levels)[i];

            buf = new ostringstream;
            (*buf) << val;
            gl2psText(buf->str().c_str(),"Times",8);
            delete buf;
         }
      }
   }

   // glMultMatrixd (rotmat);
   // glScaled(xscale, yscale, zscale);
   // glTranslated(-(x[0]+x[1])/2, -(y[0]+y[1])/2, -(z[0]+z[1])/2);
   glPopMatrix();
}

void VisualizationSceneScalarData::DrawCoordinateCross()
{
   glMatrixMode (GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();

   glMatrixMode (GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   GLint viewport[4];
   glGetIntegerv (GL_VIEWPORT, viewport);

   glTranslatef (-1, -1, 0.0);
   glScaled (40.0 / viewport[2], 40.0 / viewport[3], 1);
   glTranslatef (2.0, 2.0, 0.0);
   glMultMatrixd (rotmat);

   // glEnable (GL_COLOR_MATERIAL);

   // glLineWidth (1.0f);

   float lenx, leny, lenz;

   lenx = leny = lenz = 1;

   Arrow2(0,0,0,1,0,0, 0.9);
   Arrow2(0,0,0,0,1,0, 0.9);
   Arrow2(0,0,0,0,0,1, 0.9);

   glPushAttrib (GL_LIST_BIT);
   extern GLuint fontbase;
   glListBase (fontbase);

   const char *a_labels[] = {"x", "y", "z"};
   glRasterPos3d (lenx, 0.0f, 0.0f);
   if (print) gl2psText(a_labels[0],"Times",8);
   glCallLists (1, GL_UNSIGNED_BYTE, a_labels[0]);

   glRasterPos3d (0.0f, leny, 0.0f);
   if (print) gl2psText(a_labels[1],"Times",8);
   glCallLists (1, GL_UNSIGNED_BYTE, a_labels[1]);

   glRasterPos3d (0.0f, 0.0f, lenz);
   if (print) gl2psText(a_labels[2],"Times",8);
   glCallLists (1, GL_UNSIGNED_BYTE, a_labels[2]);

   glPopAttrib();

   glMatrixMode (GL_PROJECTION);
   glPopMatrix();
   glMatrixMode (GL_MODELVIEW);
   glPopMatrix();

   if (print)
   {
      ostringstream buf4;
      buf4 << "(" << x[0] << "," << y[0] << ","  << z[0] << ")" ;
      glRasterPos3d (x[0], y[0], z[0]);
      if (print) gl2psText(buf4.str().c_str(),"Times",8);

      ostringstream buf5;
      buf5 << "(" << x[1] << "," << y[1] << "," << z[1] << ")" ;
      glRasterPos3d (x[1], y[1], z[1]);
      if (print) gl2psText(buf5.str().c_str(),"Times",8);
   }
}

VisualizationSceneScalarData * vsdata;
extern VisualizationScene  * locscene;

void KeyCPressed()
{
   vsdata -> ToggleDrawColorbar();
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

void KeylPressed()
{
   vsdata -> ToggleLight();
   if (! vsdata -> light)
      glDisable(GL_LIGHTING);
   else
      glEnable(GL_LIGHTING);
   SendExposeEvent();
}

void KeyLPressed()
{
   MySetColorLogscale = !MySetColorLogscale;
   vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyrPressed()
{
   locscene -> spinning = 0;
   RemoveIdleFunc(MainLoop);
   vsdata -> CenterObject();

   locscene -> ViewAngle = 45.0;
   locscene -> ViewScale = 1.0;
   locscene -> ViewCenterX = 0.0;
   locscene -> ViewCenterY = 0.0;
   vsdata -> key_r_state = 0;
   SendExposeEvent();
}

void KeyRPressed()
{
   locscene -> spinning = 0;
   RemoveIdleFunc(MainLoop);
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity();
   glGetDoublev (GL_MODELVIEW_MATRIX,   locscene -> translmat);
   Set_Light();

   switch (vsdata -> key_r_state)
   {
   case 0:
      break;

   case 1:
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);
      break;

   case 2:
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);
      glRotatef(-90.0, 0.0f, 0.0f, 1.0f);
      break;

   case 3:
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);
      glRotatef(-180.0, 0.0f, 0.0f, 1.0f);
      break;

   case 4:
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);
      glRotatef(-270.0, 0.0f, 0.0f, 1.0f);
      break;

   case 5:
      glRotatef(180.0, 1.0f, 0.0f, 0.0f);
      break;
   }

   // if (locscene -> view != 2) // make 'R' work the same in 2D and 3D
   vsdata -> key_r_state = (vsdata -> key_r_state+1)%6;

   glGetDoublev (GL_MODELVIEW_MATRIX,   locscene -> rotmat);
   SendExposeEvent();
}

void Next_RGB_Palette();
int Select_New_RGB_Palette();

void KeyPPressed()
{
   Next_RGB_Palette();
   if (!GetUseTexture())
      vsdata->EventUpdateColors();
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

extern int RepeatPaletteTimes;
void KeyF6Pressed()
{
   cout << "Palette is repeated " << RepeatPaletteTimes << " times.\n"
        << "New value: " << flush;
   cin >> RepeatPaletteTimes;
   if (RepeatPaletteTimes < 1)
      RepeatPaletteTimes = 1;
   else if (RepeatPaletteTimes > 32)
      RepeatPaletteTimes = 32;
   cout << "Palette will be repeated " << RepeatPaletteTimes
        << " times now." << endl;

   Select_New_RGB_Palette();

   if (!GetUseTexture())
      vsdata->EventUpdateColors();
   SendExposeEvent();
}

void KeyF7Pressed()
{
   cout << "[minv,maxv] = [" << vsdata->GetMinV() << "," << vsdata->GetMaxV()
        << "]  maxv-minv = " << vsdata->GetMaxV()-vsdata->GetMinV() << "\n"
        << "New value for minv: " << flush;
   cin >> vsdata->GetMinV();
   cout << "New value for maxv: " << flush;
   cin >> vsdata->GetMaxV();
   vsdata -> UpdateValueRange();
   vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyBackslashPressed()
{
   double x, y, z, w;

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

   glLoadIdentity();
   GLfloat light[] = { x, y, z, w };
   // load modelview matrix before calling glLightfv?
   glLightfv(GL_LIGHT0, GL_POSITION, light);
   SendExposeEvent();
}

void KeyTPressed()
{
   int ml;

   ml = Next_Material_And_Light();
   vsdata -> EventUpdateColors();
   SendExposeEvent();
   cout << "New material/light : " << ml << endl;
}

void KeyGPressed()
{
   Toggle_Background();
   SendExposeEvent();
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
   MatAlpha -= 0.05;
   if (MatAlpha < 0.0)
      MatAlpha = 0.0;
   vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyKPressed()
{
   MatAlpha += 0.05;
   if (MatAlpha > 1.0)
      MatAlpha = 1.0;
   vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyAPressed()
{
   static int a = 0;

   a = 1-a;

   if (a)
      Set_AntiAliasing();
   else
      Remove_AntiAliasing();

   vsdata -> EventUpdateColors();
   SendExposeEvent();
}

void KeyCommaPressed()
{
   MatAlphaCenter -= 0.25;
   vsdata -> EventUpdateColors();
   SendExposeEvent();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << MatAlphaCenter << endl;
#endif
}

void KeyLessPressed()
{
   MatAlphaCenter += 0.25;
   vsdata -> EventUpdateColors();
   SendExposeEvent();
#ifdef GLVIS_DEBUG
   cout << "MatAlphaCenter = " << MatAlphaCenter << endl;
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

void VisualizationSceneScalarData::ToggleRuler()
{
   ruler_on = (ruler_on + 1) % 3;
}

void VisualizationSceneScalarData::RulerPosition()
{
   cout << "Current ruler position: (" << ruler_x << ','
        << ruler_y << ',' << ruler_z << ")\n";
   cout << "x = " << flush; cin >> ruler_x;
   cout << "y = " << flush; cin >> ruler_y;
   cout << "z = " << flush; cin >> ruler_z;
   if (ruler_x < x[0])
      ruler_x = x[0];
   else if (ruler_x > x[1])
      ruler_x = x[1];
   if (ruler_y < y[0])
      ruler_y = y[0];
   else if (ruler_y > y[1])
      ruler_y = y[1];
   if (ruler_z < z[0])
      ruler_z = z[0];
   else if (ruler_z > z[1])
      ruler_z = z[1];
   cout << "New ruler position: (" << ruler_x << ','
        << ruler_y << ',' << ruler_z << ")" << endl;
}

void VisualizationSceneScalarData::DrawRuler()
{
   if (ruler_on)
   {
      glBegin(GL_LINES);
      glVertex3d(x[0], ruler_y, ruler_z);
      glVertex3d(x[1], ruler_y, ruler_z);
      glVertex3d(ruler_x, y[0], ruler_z);
      glVertex3d(ruler_x, y[1], ruler_z);
      glVertex3d(ruler_x, ruler_y, z[0]);
      glVertex3d(ruler_x, ruler_y, z[1]);
      glEnd();

      if (ruler_on == 2)
      {
         glBegin(GL_LINE_LOOP);
         glVertex3d(x[0], y[0], ruler_z);
         glVertex3d(x[1], y[0], ruler_z);
         glVertex3d(x[1], y[1], ruler_z);
         glVertex3d(x[0], y[1], ruler_z);
         glEnd();

         glBegin(GL_LINE_LOOP);
         glVertex3d(x[0], ruler_y, z[0]);
         glVertex3d(x[1], ruler_y, z[0]);
         glVertex3d(x[1], ruler_y, z[1]);
         glVertex3d(x[0], ruler_y, z[1]);
         glEnd();

         glBegin(GL_LINE_LOOP);
         glVertex3d(ruler_x, y[0], z[0]);
         glVertex3d(ruler_x, y[1], z[0]);
         glVertex3d(ruler_x, y[1], z[1]);
         glVertex3d(ruler_x, y[0], z[1]);
         glEnd();

         Set_Material();
         if (light)
            glEnable(GL_LIGHTING);
         glBegin(GL_QUADS);
         glColor4d(0.8, 0.8, 0.8, 1.0);
         if (light)
            glNormal3d(0, 0, 1);
         glVertex3d(x[0], y[0], ruler_z);
         glVertex3d(x[1], y[0], ruler_z);
         glVertex3d(x[1], y[1], ruler_z);
         glVertex3d(x[0], y[1], ruler_z);

         if (light)
            glNormal3d(0, 1, 0);
         glVertex3d(x[0], ruler_y, z[0]);
         glVertex3d(x[0], ruler_y, z[1]);
         glVertex3d(x[1], ruler_y, z[1]);
         glVertex3d(x[1], ruler_y, z[0]);

         if (light)
            glNormal3d(1, 0, 0);
         glVertex3d(ruler_x, y[0], z[0]);
         glVertex3d(ruler_x, y[1], z[0]);
         glVertex3d(ruler_x, y[1], z[1]);
         glVertex3d(ruler_x, y[0], z[1]);
         glEnd();
         glDisable(GL_LIGHTING);
         Set_Black_Material();
      }
   }
}

void VisualizationSceneScalarData::ToggleTexture()
{
   SetUseTexture((GetUseTexture()+1)%3);
   if (GetUseTexture() != 2)
      EventUpdateColors();
}

VisualizationSceneScalarData::VisualizationSceneScalarData(
   Mesh & m, Vector & s)
{
   mesh = &m;
   sol  = &s;

   Init();
}

void VisualizationSceneScalarData::Init()
{
   vsdata = this;

   arrow_type = arrow_scaling_type = 0;
   scaling = 0;
   light   = 1;
   drawaxes = colorbar = 0;

   CuttingPlane = NULL;

   key_r_state = 0;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      auxKeyFunc (AUX_l, KeylPressed);
      auxKeyFunc (AUX_L, KeyLPressed);

      auxKeyFunc (AUX_s, KeySPressed);

      auxKeyFunc (AUX_a, KeyaPressed);
      auxKeyFunc (AUX_A, KeyAPressed);

      auxKeyFunc (AUX_r, KeyrPressed);
      auxKeyFunc (AUX_R, KeyRPressed);

      auxKeyFunc (AUX_p, KeyPPressed);

      auxKeyFunc (XK_F5, KeyF5Pressed);
      auxKeyFunc (XK_F6, KeyF6Pressed);
      auxKeyFunc (XK_F7, KeyF7Pressed);

      auxKeyFunc (XK_backslash, KeyBackslashPressed);
      auxKeyFunc (AUX_t, KeyTPressed);
      auxKeyFunc (AUX_T, KeyTPressed);

      auxKeyFunc (AUX_g, KeyGPressed);
      auxKeyFunc (AUX_G, KeyGPressed);

      auxKeyFunc (AUX_c, KeyCPressed);
      auxKeyFunc (AUX_C, KeyCPressed);

      auxKeyFunc (AUX_k, KeykPressed);
      auxKeyFunc (AUX_K, KeyKPressed);

      auxKeyFunc (XK_F1, KeyF1Pressed);
      auxKeyFunc (XK_F2, KeyF2Pressed);

      auxKeyFunc (XK_comma, KeyCommaPressed);
      auxKeyFunc (XK_less, KeyLessPressed);
      auxKeyFunc (XK_grave, KeyGravePressed);
      auxKeyFunc (XK_asciitilde, KeyTildePressed);

      auxKeyFunc (XK_exclam, KeyToggleTexture);
   }

   Set_Light();

   glEnable (GL_COLOR_MATERIAL);
   glShadeModel (GL_SMOOTH);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glDepthFunc(GL_LEQUAL);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_AUTO_NORMAL);
   glEnable(GL_NORMALIZE);

#ifdef GLVIS_MULTISAMPLE
   glDisable(GL_MULTISAMPLE);
#endif

   // add black fog
   // glEnable(GL_FOG);
   // GLfloat fogcol[4] = {0,0,0,1};
   // glFogfv(GL_FOG_COLOR, fogcol);
   // glFogf(GL_FOG_DENSITY,1.0f);

   axeslist   = glGenLists (1);
   FindNewBox();
   ruler_on = 0;
   ruler_x = 0.5 * (x[0] + x[1]);
   ruler_y = 0.5 * (y[0] + y[1]);
   ruler_z = 0.5 * (z[0] + z[1]);
   PrepareAxes();
}

VisualizationSceneScalarData::~VisualizationSceneScalarData()
{
   glDeleteLists (axeslist, 1);
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
      if ((x[1]-x[0])>eps) xscale /= (x[1]-x[0]);
      if ((y[1]-y[0])>eps) yscale /= (y[1]-y[0]);
      if ((z[1]-z[0])>eps) zscale /= (z[1]-z[0]);
   }
   else
   {
      // Find the largest side of the box in xscale
      xscale = x[1]-x[0];
      yscale = y[1]-y[0];
      zscale = z[1]-z[0];
      if (xscale < yscale) xscale = yscale;
      if (xscale < zscale) xscale = zscale;
      // Set proportional scaling so that the largest side of the box is 1.
      if (xscale > eps)
         xscale = ( 1.0 / xscale );
      else
         xscale = 1.0;
      zscale = yscale = xscale;
   }
}


void VisualizationSceneScalarData::SetValueRange(double min, double max)
{
   // if (max < min)
   //   return;

   minv = min;
   maxv = max;

   UpdateValueRange();
   // EventUpdateColors();
}

void VisualizationSceneScalarData::ResetScaling()
{
   SetNewScalingFromBox();
}

void VisualizationSceneScalarData::CenterObject()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity();
   glGetDoublev (GL_MODELVIEW_MATRIX, translmat);

   Set_Light();

   // scaling = 1;
   // SetNewScalingFromBox();

   glRotatef(-60.0, 1.0f, 0.0f, 0.0f);
   glRotatef(-40.0, 0.0f, 0.0f, 1.0f);
   glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
}

void VisualizationSceneScalarData::CenterObject2D()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity();
   glGetDoublev (GL_MODELVIEW_MATRIX, translmat);
   Set_Light();
   glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
}

void VisualizationSceneScalarData::PrepareAxes()
{
   glNewList (axeslist, GL_COMPILE);

   // Set_Black_Material();

   glBegin (GL_LINE_LOOP);
   glVertex3d(x[0], y[0], z[0]);
   glVertex3d(x[1], y[0], z[0]);
   glVertex3d(x[1], y[1], z[0]);
   glVertex3d(x[0], y[1], z[0]);
   glEnd();

   glBegin (GL_LINE_LOOP);
   glVertex3d(x[0], y[0], z[1]);
   glVertex3d(x[1], y[0], z[1]);
   glVertex3d(x[1], y[1], z[1]);
   glVertex3d(x[0], y[1], z[1]);
   glEnd();

   glBegin (GL_LINES);
   glVertex3d(x[0], y[0], z[0]);
   glVertex3d(x[0], y[0], z[1]);
   glVertex3d(x[1], y[0], z[0]);
   glVertex3d(x[1], y[0], z[1]);
   glVertex3d(x[1], y[1], z[0]);
   glVertex3d(x[1], y[1], z[1]);
   glVertex3d(x[0], y[1], z[0]);
   glVertex3d(x[0], y[1], z[1]);
   glEnd();

   // Write the coordinates of the lower left and upper right corner.
   //   glEnable (GL_COLOR_MATERIAL);
   //   GLfloat textcol[3] = { 0, 0, 0 };
   //   glColor3fv (textcol);

   if (drawaxes == 1)
   {
      glPushAttrib (GL_LIST_BIT);
      extern GLuint fontbase;
      glListBase (fontbase);

      ostringstream buf;
      buf << "(" << x[0] << "," << y[0] << ","  << z[0] << ")" ;
      glRasterPos3d (x[0], y[0], z[0]);
      glCallLists (buf.str().size(), GL_UNSIGNED_BYTE, buf.str().c_str());

      ostringstream buf1;
      buf1 << "(" << x[1] << "," << y[1] << "," << z[1] << ")" ;
      glRasterPos3d (x[1], y[1], z[1]);
      glCallLists (buf1.str().size(), GL_UNSIGNED_BYTE, buf1.str().c_str());

      glPopAttrib();
   }

   glEndList();
}

void VisualizationSceneScalarData::DrawPolygonLevelLines(
   double * point, int n, Array<double> & level)
{
   int l,k,k1;
   double curve,t;
   double p[3];

   for (l=0; l<level.Size(); l++)
   {
      glBegin(GL_LINE_LOOP);
      // for some SUN's: glBegin(GL_LINES);

      curve = level[l];
      for(k=0; k < n; k++)
      {
         k1 = (k+1)%n;
         if ( (curve <=point[4*k+3] && curve >= point[4*k1+3]) ||
              (curve >=point[4*k+3] && curve <= point[4*k1+3]) )
         {
            if ((curve - point[4*k1+3]) == 0.)
               t = 1.;
            else if ((curve - point[4*k+3]) == 0.)
               t = 0.;
            else
               t = (curve - point[4*k+3]) / (point[4*k1+3]-point[4*k+3]);
            p[0] = (1.0-t)*point[4*k+0]+t*point[4*k1+0];
            p[1] = (1.0-t)*point[4*k+1]+t*point[4*k1+1];
            p[2] = (1.0-t)*point[4*k+2]+t*point[4*k1+2];
            glVertex3dv(p);
         }
      }
      glEnd();
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
}

const char *strings_on_off[] = { "off", "on" };

void VisualizationSceneScalarData::PrintState()
{
   cout << "\nlight " << strings_on_off[light ? 1 : 0]
        << "\nperspective " << strings_on_off[OrthogonalProjection ? 0 : 1]
        << "\nviewcenter " << ViewCenterX << ' ' << ViewCenterY
        << "\nzoom " << (OrthogonalProjection ? ViewScale :
                         tan(M_PI / 8.) / tan(ViewAngle * (M_PI / 360.0)))
        << "\nvaluerange " << minv << ' ' << maxv
        << '\n' << endl;
}

void VisualizationSceneScalarData::ShrinkPoints(DenseMatrix &pointmat,
                                                int i, int fn, int fo)
{
   int dim = mesh->Dimension();

   if (shrink != 1.0)
   {
      if (dim == 2)
      {
         int k;
         double cx=0.0, cy=0.0;

         for (k = 0; k < pointmat.Width(); k++)
         {
            cx += pointmat(0,k);
            cy += pointmat(1,k);
         }
         cx /= pointmat.Width();
         cy /= pointmat.Width();

         for (k = 0; k < pointmat.Width(); k++)
         {
            pointmat(0,k) = shrink*pointmat(0,k) + (1-shrink)*cx;
            pointmat(1,k) = shrink*pointmat(1,k) + (1-shrink)*cy;
         }
      }
      else
      {
         int attr = mesh->GetBdrAttribute(i);
         for (int k = 0; k < pointmat.Width(); k++)
            for (int d = 0; d < dim; d++)
               pointmat(d,k) = shrink*pointmat(d,k) + (1-shrink)*bdrc(d,attr-1);
      }
   }

   if (shrinkmat != 1.0)
   {
      int attr, elem1, elem2;
      if (dim == 2)
         attr = mesh->GetAttribute(i);
      else
      {
         mesh->GetFaceElements(fn, &elem1, &elem2);
         if (fo % 2 == 0)
            attr = mesh->GetAttribute(elem1);
         else
            attr = mesh->GetAttribute(elem2);
      }

      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < dim; d++)
            pointmat(d,k) = shrinkmat*pointmat(d,k) + (1-shrinkmat)*matc(d,attr-1);
   }
}

void VisualizationSceneScalarData::ComputeBdrAttrCenter()
{
   DenseMatrix pointmat;
   Vector nbdrc(mesh->bdr_attributes.Max());
   int dim = mesh->Dimension();

   bdrc.SetSize(dim,mesh->bdr_attributes.Max());
   bdrc = 0.0;
   nbdrc = 0.0;

   for (int i = 0; i < mesh -> GetNBE(); i++)
   {
      mesh->GetBdrPointMatrix(i, pointmat);
      nbdrc(mesh->GetBdrAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < dim; d++)
            bdrc(d,mesh->GetBdrAttribute(i)-1) += pointmat(d,k);
   }

   for (int i = 0; i < mesh->bdr_attributes.Max(); i++)
      if (nbdrc(i) != 0)
         for (int d = 0; d < dim; d++)
            bdrc(d,i) /= nbdrc(i);
}

void VisualizationSceneScalarData::ComputeElemAttrCenter()
{
   DenseMatrix pointmat;
   Vector nmatc(mesh->attributes.Max());
   int dim = mesh->Dimension();

   matc.SetSize(dim,mesh->attributes.Max());
   matc = 0.0;
   nmatc = 0.0;

   for (int i = 0; i < mesh -> GetNE(); i++)
   {
      mesh->GetPointMatrix(i, pointmat);
      nmatc(mesh->GetAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < dim; d++)
            matc(d,mesh->GetAttribute(i)-1) += pointmat(d,k);
   }

   for (int i = 0; i < mesh->attributes.Max(); i++)
      if (nmatc(i) != 0)
         for (int d = 0; d < dim; d++)
            matc(d,i) /= nmatc(i);
}


Plane::Plane(double A,double B,double C,double D)
{
   eqn[0] = A;
   eqn[1] = B;
   eqn[2] = C;
   eqn[3] = D;

   CartesianToSpherical();

   double x[2] = {vsdata -> x[0], vsdata -> x[1]};
   double y[2] = {vsdata -> y[0], vsdata -> y[1]};
   double z[2] = {vsdata -> z[0], vsdata -> z[1]};
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
