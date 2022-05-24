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

   // *INDENT-OFF*
   double mat[16] =
   {
      -left[0], up[0], -dir[0], 0.0,
      -left[1], up[1], -dir[1], 0.0,
      -left[2], up[2], -dir[2], 0.0,
      0.0, 0.0, 0.0, 1.0
   };
   // *INDENT-ON*
   return glm::make_mat4(mat);
}

glm::mat4 Camera::TransposeRotMatrix()
{
   GetLeft();
   // *INDENT-OFF*
   double mat_t[16] =
   {
      -left[0], -left[1], -left[2], 0.0,
         up[0],    up[1],    up[2], 0.0,
       -dir[0],  -dir[1],  -dir[2], 0.0,
      0.0, 0.0, 0.0, 1.0
   };
   // *INDENT-ON*
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

   cut_lambda = 0.0;
   cut_updated = false;

   background = BG_WHITE;
   GetAppWindow()->getRenderer().setClearColor(1.f, 1.f, 1.f, 1.f);
   _use_cust_l0_pos = false;
   light_mat_idx = 3;
   use_light = true;

   palette.Init();
}

VisualizationScene::~VisualizationScene() {}


void VisualizationScene
:: DrawCutTriangle(gl3::GlDrawable& buff,
                   const double (&p)[4][3], const double (&cv)[4],
                   const double minv, const double maxv)
{
   // element center
   double c[3];
   c[0] = c[1] = c[2] = 0.0;
   for (int j = 0; j < 3; j++)
   {
      c[0] += p[j][0]; c[1] += p[j][1]; c[2] += p[j][2];
   }
   c[0] /= 3.0; c[1] /= 3.0; c[2] /= 3.0;

   double l = cut_lambda;
   double q[3][3];

   // corners of the cut frame
   for (int j = 0; j < 3; j++)
   {
      q[j][0] = l*p[j][0] + (1.0-l)*c[0];
      q[j][1] = l*p[j][1] + (1.0-l)*c[1];
      q[j][2] = l*p[j][2] + (1.0-l)*c[2];
   }

   double d[4][3];
   double cvv[4];

   // bottom trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[0][k]; d[1][k] = p[1][k]; d[2][k] = q[1][k]; d[3][k] = q[0][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);

   // diagonal trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[1][k]; d[1][k] = p[2][k]; d[2][k] = q[2][k]; d[3][k] = q[1][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);

   // left trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[2][k]; d[1][k] = p[0][k]; d[2][k] = q[0][k]; d[3][k] = q[2][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);
}


void VisualizationScene
:: DrawCutQuad(gl3::GlDrawable& buff,
               const double (&p)[4][3], const double (&cv)[4],
               const double minv, const double maxv)
{
   // element center
   double c[3];
   c[0] = c[1] = c[2] = 0.0;
   for (int j = 0; j < 4; j++)
   {
      c[0] += p[j][0]; c[1] += p[j][1]; c[2] += p[j][2];
   }
   c[0] /= 4.0; c[1] /= 4.0; c[2] /= 4.0;

   double l = cut_lambda;
   double q[4][3];

   // corners of the cut frame
   for (int j = 0; j < 4; j++)
   {
      q[j][0] = l*p[j][0] + (1.0-l)*c[0];
      q[j][1] = l*p[j][1] + (1.0-l)*c[1];
      q[j][2] = l*p[j][2] + (1.0-l)*c[2];
   }

   double d[4][3];
   double cvv[4];

   // bottom trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[0][k]; d[1][k] = p[1][k]; d[2][k] = q[1][k]; d[3][k] = q[0][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);

   // right trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[1][k]; d[1][k] = p[2][k]; d[2][k] = q[2][k]; d[3][k] = q[1][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);

   // top trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[2][k]; d[1][k] = p[3][k]; d[2][k] = q[3][k]; d[3][k] = q[2][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);

   // left trapezoid
   for (int k = 0; k < 3; k++)
   {
      d[0][k] = p[3][k]; d[1][k] = p[0][k]; d[2][k] = q[0][k]; d[3][k] = q[3][k];
   }
   DrawQuad(buff, d, cvv, minv, maxv);
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
      Normalize(normals);
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

glTF_Builder::material_id
VisualizationScene::AddPaletteMaterial(glTF_Builder &bld)
{
#ifdef GLVIS_USE_LIBPNG
   const bool palette_smooth = palette.GetSmoothSetting();
   auto sampler =
      bld.addSampler(
         /* magFilter: */
         palette_smooth ? glTF_Builder::mag_filter::LINEAR :
         /**/             glTF_Builder::mag_filter::NEAREST,
         /* minFilter: */
         palette_smooth ? glTF_Builder::min_filter::LINEAR :
         /**/             glTF_Builder::min_filter::NEAREST,
         /* wrapS: */ glTF_Builder::wrap_type::CLAMP_TO_EDGE,
         /* wrapT: */ glTF_Builder::wrap_type::CLAMP_TO_EDGE);
   // create palette image
   const int palette_size = palette.GetNumColors();
   vector<array<float,4>> palette_data(palette_size);
#if 0
   glGetTextureImage(
      palette.GetColorTexture(), 0,
      gl3::GLDevice::useLegacyTextureFmts() ? GL_RGBA : GL_RGBA32F,
      GL_FLOAT,
      palette_size,
      palette_data.data());
#elif 0
   glGetTexImage(GL_TEXTURE_2D, 0,
                 gl3::GLDevice::useLegacyTextureFmts() ? GL_RGBA : GL_RGBA32F,
                 GL_FLOAT, palette_data.data());
#else
   const double *palette_data_raw = palette.GetData();
   for (int i = 0; i < palette_size; ++i)
   {
      for (int j = 0; j < 3; ++j)
      {
         palette_data[i][j] = (float) palette_data_raw[j + 3*i];
      }
      palette_data[i][3] = 1.0f;
   }
#endif
   auto palette_img =
      bld.addImage(
         /* imageName: */      "palette",
         /* width: */           palette_size,
         /* height: */          1,
         /* color4f *pixels: */ palette_data.data());
   auto palette_tex = bld.addTexture(sampler, palette_img);
   glTF_Builder::pbr_matallic_roughness palette_pbr_mr =
   {
      /* haveTexture: */       true,
      /* baseColorFactor: */   { 1.f, 1.f, 1.f, 1.f },
      /* baseColorTexture: */  palette_tex,
      /* metallicFactor: */    1.f,
      /* roughnessFactor: */   .3f
   };
#else // GLVIS_USE_LIBPNG
   glTF_Builder::pbr_matallic_roughness palette_pbr_mr =
   {
      /* haveTexture: */       false,
      /* baseColorFactor: */   { 0.f, 1.f, .5f, 1.f },
      /* baseColorTexture: */  {0},
      /* metallicFactor: */    1.f,
      /* roughnessFactor: */   .3f
   };
#endif // GLVIS_USE_LIBPNG
   auto palette_mat =
      bld.addMaterial(
         /* materialName: */         "Palette Material",
         /* pbrMetallicRoughness: */ palette_pbr_mr,
         /* doubleSided: */          true);

   return palette_mat;
}

glTF_Builder::material_id
VisualizationScene::AddBlackMaterial(glTF_Builder &bld)
{
   glTF_Builder::pbr_matallic_roughness black_pbr_mr =
   {
      /* haveTexture: */       false,
      /* baseColorFactor: */   GetLineColor(),
      /* baseColorTexture: */  {0},
      /* metallicFactor: */    1.f,
      /* roughnessFactor: */   1.f
   };
   auto black_mat =
      bld.addMaterial(
         /* materialName: */         "Black Material",
         /* pbrMetallicRoughness: */ black_pbr_mr,
         /* doubleSided: */          true);

   return black_mat;
}

glTF_Builder::material_id
VisualizationScene::AddPaletteLinesMaterial(
   glTF_Builder &bld, glTF_Builder::material_id palette_mat)
{
   glTF_Builder::pbr_matallic_roughness palette_pbr_mr_copy;
   bld.getMaterialPBRMR(palette_mat, palette_pbr_mr_copy);
   palette_pbr_mr_copy.metallicFactor = 1.f;
   palette_pbr_mr_copy.roughnessFactor = 1.f;
   auto palette_lines_mat =
      bld.addMaterial(
         /* materialName: */         "PaletteLines Material",
         /* pbrMetallicRoughness: */ palette_pbr_mr_copy,
         /* doubleSided: */          true);
   return palette_lines_mat;
}

glTF_Builder::node_id
VisualizationScene::AddModelNode(glTF_Builder &bld, const string &nodeName)
{
   auto new_node = bld.addNode(nodeName);
   // Coordinate system switch: (x,y,z) -> (x,z,-y).
   // In glTF, the translation (T), rotation (R), and scale (S) properties are
   // applied in the order: T * R * S.
   // In GLVis, we apply them in the order: R * S * T.
   bld.addNodeScale(
      new_node, { (float)xscale, (float)zscale, (float)yscale });
   bld.addNodeTranslation(
      new_node, { float(-xscale*(bb.x[0]+bb.x[1])/2),
                  float(-zscale*bb.z[0]),
                  float(yscale*(bb.y[0]+bb.y[1])/2)
                });
   return new_node;
}

// Used in VisualizationScene::AddTriangles() below.
void minmax(const float *data, size_t components, size_t stride, size_t count,
            vector<float> &mins, vector<float> &maxs)
{
   if (count == 0)
   {
      mins.assign(components, +numeric_limits<float>::infinity());
      maxs.assign(components, -numeric_limits<float>::infinity());
      return;
   }
   mins.resize(components);
   maxs.resize(components);
   for (size_t c = 0; c != components; ++c)
   {
      const auto entry = data[c];
      mins[c] = entry;
      maxs[c] = entry;
   }
   if (count == 1) { return; }
   for (size_t i = 1; i != count; ++i)
   {
      data += stride;
      for (size_t c = 0; c != components; ++c)
      {
         const auto entry = data[c];
         if (entry < mins[c]) { mins[c] = entry; }
         else if (entry > maxs[c]) { maxs[c] = entry; }
      }
   }
}

int VisualizationScene::AddTriangles(glTF_Builder &bld,
                                     glTF_Builder::mesh_id mesh,
                                     glTF_Builder::buffer_id buffer,
                                     glTF_Builder::material_id material,
                                     const gl3::GlDrawable &gl_drawable)
{
   int num_buf = 0, buf_layout = -1;
   for (int layout = 0; layout < gl3::NUM_LAYOUTS; ++layout)
   {
      const gl3::IVertexBuffer *buf = gl_drawable.buffers[layout][1].get();
      if (buf && buf->count() != 0)
      {
         num_buf++;
         buf_layout = layout;
         cout << "triangles: layout = " << layout << ", # vertices = "
              << buf->count() << '\n';
      }
   }
   int num_ibuf = 0, ibuf_layout = -1;
   for (int layout = 0; layout < gl3::NUM_LAYOUTS; ++layout)
   {
      const gl3::IIndexedBuffer *ibuf =
         gl_drawable.indexed_buffers[layout][1].get();
      if (ibuf && ibuf->getIndices().size() != 0)
      {
         num_ibuf++;
         ibuf_layout = layout;
         cout << "indexed triangles: layout = " << layout << ", # vertices = "
              << ibuf->count() << ", # indices = " << ibuf->getIndices().size()
              << '\n';
      }
   }
   if (num_buf + num_ibuf == 0) { return 0; }

   if (num_buf + num_ibuf > 1)
   {
      cout << "glTF export: skipping" << num_buf + num_ibuf - 1
           << " triangle buffer(s).\n";
   }
   const gl3::IVertexBuffer *surf_buf = nullptr;
   const gl3::IIndexedBuffer *surf_ibuf = nullptr;
   const vector<int> *surf_indices = nullptr;
   if (num_ibuf)
   {
      surf_buf = surf_ibuf = gl_drawable.indexed_buffers[ibuf_layout][1].get();
      surf_indices = &surf_ibuf->getIndices();
      buf_layout = ibuf_layout;
   }
   else
   {
      surf_buf = gl_drawable.buffers[buf_layout][1].get();
   }
   const size_t surf_vertices_count = surf_buf->count();
   const size_t surf_vertices_stride = surf_buf->getStride(); // in bytes
   const float *surf_vertices_data =
      reinterpret_cast<const float *>(surf_buf->getData());
   vector<float> vmins, vmaxs;
   int components = surf_vertices_stride/sizeof(float);
   switch (buf_layout)
   {
      case gl3::LAYOUT_VTX:
      case gl3::LAYOUT_VTX_COLOR: components = 3; break;
      case gl3::LAYOUT_VTX_TEXTURE0: components = 5; break;
      case gl3::LAYOUT_VTX_NORMAL:
      case gl3::LAYOUT_VTX_NORMAL_COLOR: components = 6; break;
      case gl3::LAYOUT_VTX_NORMAL_TEXTURE0: components = 8; break;
   }
   minmax(surf_vertices_data, components, surf_vertices_stride/sizeof(float),
          surf_vertices_count, vmins, vmaxs);

   glTF_Builder::buffer_view_id surf_indices_buf_view;
   if (surf_indices)
   {
      surf_indices_buf_view =
         bld.addBufferView(
            /* buffer: */     buffer,
            /* data: */       surf_indices->data(),
            /* byteLength: */ surf_indices->size()*sizeof(int),
            /* byteStride: */ sizeof(int),
            /* byteAlign: */  sizeof(int),
            /* target: */ glTF_Builder::target_type::ELEMENT_ARRAY_BUFFER);
   }
#if 0
   auto surf_vertices_buf_view =
      bld.addBufferView(
         /* buffer: */     buffer,
         /* data: */       surf_vertices_data,
         /* byteLength: */ surf_vertices_count*surf_vertices_stride,
         /* byteStride: */ surf_vertices_stride,
         /* byteAlign: */  sizeof(float),
         /* target: */     glTF_Builder::target_type::ARRAY_BUFFER);
#else
   auto surf_vertices_buf_view =
      bld.addBufferView(
         /* buffer: */     buffer,
         /* data: */       surf_vertices_data,
         /* byteLength: */ 0,
         /* byteStride: */ surf_vertices_stride,
         /* byteAlign: */  sizeof(float),
         /* target: */     glTF_Builder::target_type::ARRAY_BUFFER);
   // Coordinate system switch: (x,y,z) -> (x,z,-y), see:
   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#coordinate-system-and-units
   switch (buf_layout)
   {
      case gl3::Vertex::layout:
      {
         auto surf_vert_data =
            reinterpret_cast<const gl3::Vertex *>(surf_vertices_data);
         for (size_t i = 0; i != surf_vertices_count; ++i)
         {
            auto &v = surf_vert_data[i];
            gl3::Vertex vtx = { { v.coord[0], v.coord[2], -v.coord[1] } };
            bld.appendToBufferView(surf_vertices_buf_view, &vtx, sizeof(vtx));
         }
         break;
      }
      case gl3::VertexNorm::layout:
      {
         auto surf_vert_data =
            reinterpret_cast<const gl3::VertexNorm *>(surf_vertices_data);
         for (size_t i = 0; i != surf_vertices_count; ++i)
         {
            auto &v = surf_vert_data[i];
            gl3::VertexNorm vn =
            {
               { v.coord[0], v.coord[2], -v.coord[1] },
               { v.norm[0], v.norm[2], -v.norm[1] }
            };
            bld.appendToBufferView(surf_vertices_buf_view, &vn, sizeof(vn));
         }
         break;
      }
      case gl3::VertexNormTex::layout:
      {
         auto surf_vert_data =
            reinterpret_cast<const gl3::VertexNormTex *>(surf_vertices_data);
         for (size_t i = 0; i != surf_vertices_count; ++i)
         {
            auto &v = surf_vert_data[i];
            gl3::VertexNormTex tv =
            {
               { v.coord[0], v.coord[2], -v.coord[1] },
               { v.norm[0], v.norm[2], -v.norm[1] },
               v.texCoord
            };
            bld.appendToBufferView(surf_vertices_buf_view, &tv, sizeof(tv));
         }
         break;
      }
      default:
      {
         cout << "glTF export: coorditate switch for layout " << buf_layout
              << " is not implemented here:" << MFEM_LOCATION;
         bld.appendToBufferView(surf_vertices_buf_view,
                                surf_vertices_data,
                                surf_vertices_count*surf_vertices_stride);
         break;
      }
   }
#endif
   glTF_Builder::accessor_id surf_indices_acc {glTF_Builder::INVALID_ID};
   if (surf_indices)
   {
      surf_indices_acc =
         bld.addAccessor(
            /* bufferView: */    surf_indices_buf_view,
            /* byteOffset: */    0,
            /* componentType: */ glTF_Builder::component_type::UNSIGNED_INT,
            /* count: */         surf_indices->size(),
            /* tensorType: */    glTF_Builder::tensor_type::SCALAR);
   }
   auto surf_vertices_pos_acc =
      bld.addAccessorVec3f(
         /* bufferView: */ surf_vertices_buf_view,
         /* byteOffset: */ 0,
         /* count: */      surf_vertices_count,
         // Coordinate system switch: (x,y,z) -> (x,z,-y), see above
         /* vec3f min: */  { vmins[0], vmins[2], -vmaxs[1] },
         /* vec3f max: */  { vmaxs[0], vmaxs[2], -vmins[1] });
   unsigned floatOffset = 3;
   glTF_Builder::accessor_id surf_vertices_nor_acc{glTF_Builder::INVALID_ID};
   const bool have_normals =
      (buf_layout == gl3::LAYOUT_VTX_NORMAL ||
       buf_layout == gl3::LAYOUT_VTX_NORMAL_COLOR ||
       buf_layout == gl3::LAYOUT_VTX_NORMAL_TEXTURE0);
   if (have_normals)
   {
      surf_vertices_nor_acc =
         bld.addAccessor(
            /* bufferView: */    surf_vertices_buf_view,
            /* byteOffset: */    floatOffset*sizeof(float),
            /* componentType: */ glTF_Builder::component_type::FLOAT,
            /* count: */         surf_vertices_count,
            /* tensorType: */    glTF_Builder::tensor_type::VEC3);
      floatOffset += 3;
   }
   glTF_Builder::accessor_id surf_vertices_tex_acc{glTF_Builder::INVALID_ID};
   const bool have_texcoords =
      (buf_layout == gl3::LAYOUT_VTX_TEXTURE0 ||
       buf_layout == gl3::LAYOUT_VTX_NORMAL_TEXTURE0);
   if (have_texcoords)
   {
      surf_vertices_tex_acc =
         bld.addAccessorVec2f(
            /* bufferView: */ surf_vertices_buf_view,
            /* byteOffset: */ floatOffset*sizeof(float),
            /* count: */      surf_vertices_count,
            /* vec2f min: */  { vmins[floatOffset], vmins[floatOffset+1] },
            /* vec2f max: */  { vmaxs[floatOffset], vmaxs[floatOffset+1] });
      // floatOffset += 2;
   }
   bld.addMeshTriangles(
      /* mesh: */             mesh,
      /* vertexPositions: */  surf_vertices_pos_acc,
      /* vertexNormals: */    surf_vertices_nor_acc,
      /* vertexTexCoords0: */ surf_vertices_tex_acc,
      /* vertexIndices: */    surf_indices_acc,
      /* material: */         material);

   return surf_indices ? surf_indices->size()/3 : surf_vertices_count/3;
}

int VisualizationScene::AddLines(glTF_Builder &bld,
                                 glTF_Builder::mesh_id mesh,
                                 glTF_Builder::buffer_id buffer,
                                 glTF_Builder::material_id material,
                                 const gl3::GlDrawable &gl_drawable)
{
   int num_buf = 0, buf_layout = -1;
   for (int layout = 0; layout < gl3::NUM_LAYOUTS; ++layout)
   {
      const gl3::IVertexBuffer *buf = gl_drawable.buffers[layout][0].get();
      if (buf && buf->count() != 0)
      {
         num_buf++;
         buf_layout = layout;
         cout << "lines: layout = " << layout << ", # vertices = "
              << buf->count() << '\n';
      }
   }
   int num_ibuf = 0 /* , ibuf_layout = -1 */;
   for (int layout = 0; layout < gl3::NUM_LAYOUTS; ++layout)
   {
      const gl3::IIndexedBuffer *ibuf =
         gl_drawable.indexed_buffers[layout][0].get();
      if (ibuf && ibuf->getIndices().size() != 0)
      {
         num_ibuf++;
         /* ibuf_layout = layout; */
         cout << "indexed lines: layout = " << layout << ", # vertices = "
              << ibuf->count() << ", # indices = " << ibuf->getIndices().size()
              << '\n';
      }
   }
   if (num_buf + num_ibuf == 0) { return 0; }
   if (num_buf == 0)
   {
      cout << "glTF export: indexed lines are not implemented.\n";
      return 0;
   }
   if (num_buf + num_ibuf > 1)
   {
      cout << "glTF export: skipping" << num_buf + num_ibuf - 1
           << " line buffer(s).\n";
   }
   const gl3::IVertexBuffer *lines_buf =
      gl_drawable.buffers[buf_layout][0].get();
   const size_t lines_vertices_count = lines_buf->count();
   const size_t lines_vertices_stride = lines_buf->getStride(); // in bytes
   const float *lines_vertices_data =
      reinterpret_cast<const float *>(lines_buf->getData());
   vector<float> vmins, vmaxs;
   int components = lines_vertices_stride/sizeof(float);
   switch (buf_layout)
   {
      case gl3::LAYOUT_VTX:
      case gl3::LAYOUT_VTX_COLOR: components = 3; break;
      case gl3::LAYOUT_VTX_TEXTURE0: components = 5; break;
      case gl3::LAYOUT_VTX_NORMAL:
      case gl3::LAYOUT_VTX_NORMAL_COLOR: components = 6; break;
      case gl3::LAYOUT_VTX_NORMAL_TEXTURE0: components = 8; break;
   }
   minmax(lines_vertices_data, components, lines_vertices_stride/sizeof(float),
          lines_vertices_count, vmins, vmaxs);

#if 0
   auto lines_vertices_buf_view =
      bld.addBufferView(
         /* buffer: */     buffer,
         /* data: */       lines_vertices_data,
         /* byteLength: */ lines_vertices_count*lines_vertices_stride,
         /* byteStride: */ lines_vertices_stride,
         /* byteAlign: */  sizeof(float),
         /* target: */     glTF_Builder::target_type::ARRAY_BUFFER);
#else
   auto lines_vertices_buf_view =
      bld.addBufferView(
         /* buffer: */     buffer,
         /* data: */       lines_vertices_data,
         /* byteLength: */ 0,
         /* byteStride: */ lines_vertices_stride,
         /* byteAlign: */  sizeof(float),
         /* target: */     glTF_Builder::target_type::ARRAY_BUFFER);
   // Coordinate system switch: (x,y,z) -> (x,z,-y), see:
   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#coordinate-system-and-units
   switch (buf_layout)
   {
      case gl3::Vertex::layout:
      {
         auto lines_vert_data =
            reinterpret_cast<const gl3::Vertex *>(lines_vertices_data);
         for (size_t i = 0; i != lines_vertices_count; ++i)
         {
            auto &v = lines_vert_data[i];
            gl3::Vertex vtx = { { v.coord[0], v.coord[2], -v.coord[1] } };
            bld.appendToBufferView(lines_vertices_buf_view, &vtx, sizeof(vtx));
         }
         break;
      }
      case gl3::VertexColor::layout:
      {
         auto lines_vert_data =
            reinterpret_cast<const gl3::VertexColor *>(lines_vertices_data);
         for (size_t i = 0; i != lines_vertices_count; ++i)
         {
            auto &v = lines_vert_data[i];
            gl3::VertexColor vc =
            {
               { v.coord[0], v.coord[2], -v.coord[1] }, v.color
            };
            bld.appendToBufferView(lines_vertices_buf_view, &vc, sizeof(vc));
         }
         break;
      }
      case gl3::VertexTex::layout:
      {
         auto lines_vert_data =
            reinterpret_cast<const gl3::VertexTex *>(lines_vertices_data);
         for (size_t i = 0; i != lines_vertices_count; ++i)
         {
            auto &v = lines_vert_data[i];
            gl3::VertexTex vt =
            {
               { v.coord[0], v.coord[2], -v.coord[1] }, v.texCoord
            };
            bld.appendToBufferView(lines_vertices_buf_view, &vt, sizeof(vt));
         }
         break;
      }
      default:
      {
         cout << "glTF export: coorditate switch for layout " << buf_layout
              << " is not implemented here:" << MFEM_LOCATION;
         bld.appendToBufferView(lines_vertices_buf_view,
                                lines_vertices_data,
                                lines_vertices_count*lines_vertices_stride);
         break;
      }
   }
#endif
   auto lines_vertices_pos_acc =
      bld.addAccessorVec3f(
         /* bufferView: */ lines_vertices_buf_view,
         /* byteOffset: */ 0,
         /* count: */      lines_vertices_count,
         // Coordinate system switch: (x,y,z) -> (x,z,-y), see above
         /* vec3f min: */  { vmins[0], vmins[2], -vmaxs[1] },
         /* vec3f max: */  { vmaxs[0], vmaxs[2], -vmins[1] });
   unsigned floatOffset = 3;
   glTF_Builder::accessor_id lines_vertices_col_acc{glTF_Builder::INVALID_ID};
   const bool have_colors = (buf_layout == gl3::LAYOUT_VTX_COLOR);
   if (have_colors)
   {
      lines_vertices_col_acc =
         bld.addAccessor(
            /* bufferView: */    lines_vertices_buf_view,
            /* byteOffset: */    floatOffset*sizeof(float),
            /* componentType: */ glTF_Builder::component_type::UNSIGNED_BYTE,
            /* count: */         lines_vertices_count,
            /* tensorType: */    glTF_Builder::tensor_type::VEC4);
      // floatOffset += 2;
   }
   glTF_Builder::accessor_id lines_vertices_tex_acc{glTF_Builder::INVALID_ID};
   const bool have_texcoords = (buf_layout == gl3::LAYOUT_VTX_TEXTURE0);
   if (have_texcoords)
   {
      lines_vertices_tex_acc =
         bld.addAccessorVec2f(
            /* bufferView: */ lines_vertices_buf_view,
            /* byteOffset: */ floatOffset*sizeof(float),
            /* count: */      lines_vertices_count,
            /* vec2f min: */  { vmins[floatOffset], vmins[floatOffset+1] },
            /* vec2f max: */  { vmaxs[floatOffset], vmaxs[floatOffset+1] });
      // floatOffset += 2;
   }
   bld.addMeshLines(
      /* mesh: */             mesh,
      /* vertexPositions: */  lines_vertices_pos_acc,
      /* vertexTexcoords0: */ lines_vertices_tex_acc,
      /* vertexColors0: */    lines_vertices_col_acc,
      /* material: */         material);

   return lines_vertices_count/2;
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
      GetAppWindow()->getRenderer().setClearColor(1.f, 1.f, 1.f, 1.f);
   }
   else
   {
      background = BG_BLK;
      GetAppWindow()->getRenderer().setClearColor(0.f, 0.f, 0.f, 1.f);
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
   modelView.translate(-(bb.x[0]+bb.x[1])/2, -(bb.y[0]+bb.y[1])/2,
                       -(bb.z[0]+bb.z[1])/2);
   return modelView.mtx;
}
