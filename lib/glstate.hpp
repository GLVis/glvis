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
#ifndef GLSTATE_HPP
#define GLSTATE_HPP
#include <SDL2/SDL.h>
#include "platform_gl.hpp"
#include <SDL2/SDL_opengl.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_transform.hpp>

struct GlMatrix {
    glm::mat4 mtx;

    GlMatrix()
        : mtx(4) {}

    void rotate(double angle, double x, double y, double z) {
        mtx = glm::rotate(mtx, glm::radians(angle), glm::vec3(x,y,z));
    }

    void mult(glm::mat4 rhs) {
        mtx = mtx * rhs;
    }

    void translate(double x, double y, double z) {
        mtx = glm::translate(mtx, glm::vec3(x, y, z));
    }

    void scale(double x, double y, double z) {
        mtx = glm::scale(mtx, glm::vec3(x, y, z));
    }

    void identity() {
        mtx = glm::mat4();
    }
};

class GlState
{
    GLuint program;
public:
    GlMatrix modelView;
    GlMatrix projection;

    GlState()
        : program(0) {
    }

    ~GlState() {
    }

    /**
     * Compiles the rendering pipeline shaders.
     */
    bool compileShaders();

    void loadMatrixUniforms() {
        GLuint locModelView = glGetUniformLocation(program, "modelViewMatrix");
        GLuint locProject = glGetUniformLocation(program, "projectionMatrix");

    }
};

#endif
