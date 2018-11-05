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

#include "gl3print.hpp"
#include "aux_vis.hpp"
#include <cmath>

namespace gl3
{

inline GL2PSvertex VertFBtoGL2PS(FeedbackVertex v, int (&vp)[4], float half_w, float half_h) {
    //clip coords -> ndc
    float x = v.pos[0] / v.pos[3];
    float y = v.pos[1] / v.pos[3];
    float z = v.pos[2] / v.pos[3];
    return {
        //GL2PSvertex.xyz
        //ndc -> device coords
        {
            half_w * x + vp[0] + half_w,
            half_h * y + vp[1] + half_h,
            z
        },
        //GL2PSvertex.rgba
        { v.color[0], v.color[1], v.color[2], v.color[3] }
    };
}

void processTriangleTransformFeedback(FeedbackVertex * buf, size_t numVerts) {
    int vp[4];
    GetGlState()->getViewport(vp);
    float half_w = (vp[2] - vp[0]) * 0.5f,
          half_h = (vp[3] - vp[1]) * 0.5f;

    for (size_t i = 0; i < numVerts; i += 3) {
        GL2PSvertex tri_vtx[3];
        if (!GetGlState()->isClipPlaneEnabled()
            || (buf[i].clipCoord >= 0.f
                && buf[i+1].clipCoord >= 0.f
                && buf[i+2].clipCoord >= 0.f)) {
            for (int j = 0; j < 3; j++) {
                tri_vtx[j] = VertFBtoGL2PS(buf[i+j], vp, half_w, half_h);
            }
            gl2psAddPolyPrimitive(GL2PS_TRIANGLE, 3, tri_vtx, 0, 0.f, 0.f, 0xffff, 1, 1, 0, 0, 0);
        } else if (buf[i].clipCoord < 0.f
                    && buf[i+1].clipCoord < 0.f
                    && buf[i+2].clipCoord < 0.f) {
            continue;
        } else {
            //clip through middle of triangle
            for (int j = 0; j < 3; j++) {
                int i_a = i+j;
                int i_b = i+((j+1) % 3);
                int i_c = i+((j+2) % 3);
                if ((buf[i_a].clipCoord < 0.f) == (buf[i_b].clipCoord < 0.f)) {
                    //pts a, b are on same side of clip plane
                    //compute clip pts
                    FeedbackVertex n[2];
                    //perspective-correct interpolation factors, needed for colors
                    float c_w_a = buf[i_a].clipCoord / buf[i_a].pos[3];
                    float c_w_b = buf[i_b].clipCoord / buf[i_b].pos[3];
                    float c_w_c = buf[i_c].clipCoord / buf[i_c].pos[3];
                    for (int ci = 0; ci < 4; ci++) {
                        //a --- n_0 --- c
                        n[0].pos[ci] = (buf[i_a].pos[ci] * buf[i_c].clipCoord
                                        - buf[i_c].pos[ci] * buf[i_a].clipCoord);
                        n[0].color[ci] = (buf[i_a].color[ci] * c_w_c
                                        - buf[i_c].color[ci] * c_w_a)
                                        / (c_w_c - c_w_a);
                        //b --- n_1 --- c
                        n[1].pos[ci] = (buf[i_b].pos[ci] * buf[i_c].clipCoord
                                        - buf[i_c].pos[ci] * buf[i_b].clipCoord);
                        n[1].color[ci] = (buf[i_b].color[ci] * c_w_c
                                        - buf[i_c].color[ci] * c_w_b)
                                        / (c_w_c - c_w_b);
                    }
                    if (buf[i_c].clipCoord < 0.f) {
                        //pts a, b are in clip plane
                        //create quadrilateral a -- n_0 -- n_1 -- b
                        GL2PSvertex quad[4] = {
                            VertFBtoGL2PS(buf[i_a], vp, half_w, half_h),
                            VertFBtoGL2PS(n[0], vp, half_w, half_h),
                            VertFBtoGL2PS(n[1], vp, half_w, half_h),
                            VertFBtoGL2PS(buf[i_b], vp, half_w, half_h)
                        };
                        // split along vtx 0-2
                        gl2psAddPolyPrimitive(GL2PS_QUADRANGLE, 4, quad, 0, 0.f, 0.f, 0xffff, 1, 1, 0, 0, 0);
                    } else {
                        //pt c is in clip plane
                        //add triangle c -- n_0 -- n_1
                        tri_vtx[0] = VertFBtoGL2PS(buf[i_c], vp, half_w, half_h);
                        tri_vtx[1] = VertFBtoGL2PS(n[0], vp, half_w, half_h);
                        tri_vtx[2] = VertFBtoGL2PS(n[1], vp, half_w, half_h);
                        gl2psAddPolyPrimitive(GL2PS_TRIANGLE, 3, tri_vtx, 0, 0.f, 0.f, 0xffff, 1, 1, 0, 0, 0);
                    }
                    break;
                }
            }
        }
    }
}

void processLineTransformFeedback(FeedbackVertex * buf, size_t numVerts) {
    int vp[4];
    GetGlState()->getViewport(vp);
    float half_w = (vp[2] - vp[0]) * 0.5f,
          half_h = (vp[3] - vp[1]) * 0.5f;

    for (size_t i = 0; i < numVerts; i += 2) {
        GL2PSvertex line_vtx[2];
        if (!GetGlState()->isClipPlaneEnabled() ||
            (buf[i].clipCoord >= 0.f && buf[i+1].clipCoord >= 0.f)) {
            line_vtx[0] = VertFBtoGL2PS(buf[i], vp, half_w, half_h);
            line_vtx[1] = VertFBtoGL2PS(buf[i+1], vp, half_w, half_h);
        } else if (buf[i].clipCoord < 0.f && buf[i+1].clipCoord < 0.f) {
            //outside of clip plane;
            continue;
        } else {
            int i_a, i_b;
            if (buf[i].clipCoord < 0.f) {
                //vertex i lies in clipped region
                i_a = i+1;
                i_b = i;
            } else { //buf[i+1].clipCoord < 0.f
                //vertex i+1 lies in clipped region
                i_a = i;
                i_b = i+1;
            }
            //compute new vertex (CbVa - CaVb), where Vb lies in the clipped region
            FeedbackVertex clip_vert;
            //perspective-correct interpolation factors for color
            float c_w_a = buf[i_a].clipCoord / buf[i_a].pos[3];
            float c_w_b = buf[i_b].clipCoord / buf[i_b].pos[3];
            for (int j = 0; j < 4; j++) {
                clip_vert.pos[j] = buf[i_a].pos[j] * buf[i_b].clipCoord
                                    - buf[i_b].pos[j] * buf[i_a].clipCoord;
                clip_vert.color[j] = (buf[i_a].color[j] * c_w_b
                                    - buf[i_b].color[j] * c_w_a)
                                    / (c_w_b - c_w_a);
            }
            line_vtx[0] = VertFBtoGL2PS(clip_vert, vp, half_w, half_h);
            line_vtx[1] = VertFBtoGL2PS(buf[i_a], vp, half_w, half_h);
        }
        gl2psAddPolyPrimitive(GL2PS_LINE, 2, line_vtx, 0, 0.f, 0.f,
                              0xFFFF, 1, 0.2, 0, 0, 0);  
    }
}

}
