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
#include <string>

#include "platform_gl.hpp"

#include "material.hpp"

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

struct GlMatrix {
    glm::mat4 mtx;

    /**
     * Applies a rotation transform to the matrix.
     */
    void rotate(float angle, double x, double y, double z) {
        mtx = glm::rotate(mtx, glm::radians(angle), glm::vec3(x,y,z));
    }

    void mult(glm::mat4 rhs) {
        mtx = mtx * rhs;
    }

    /**
     * Applies a translation transform to the matrix.
     */
    void translate(double x, double y, double z) {
        mtx = glm::translate(mtx, glm::vec3(x, y, z));
    }

    /**
     * Applies a scale transform to the matrix.
     */
    void scale(double x, double y, double z) {
        mtx = glm::scale(mtx, glm::vec3(x, y, z));
    }

    /**
     * Sets the matrix to an orthographic projection.
     */
    void ortho(double left,
               double right,
               double bottom,
               double top,
               double z_near,
               double z_far) {
        mtx = glm::ortho(left, right, bottom, top, z_near, z_far);
    }

    /**
     * Sets the matrix to a perspective projection.
     */
    void perspective(double fov, double aspect, double z_near, double z_far) {
        mtx = glm::perspective(glm::radians(fov), aspect, z_near, z_far);
    }

    /**
     * Sets the matrix to the identity matrix.
     */
    void identity() {
        mtx = glm::mat4(1.0);
    }
};

class GlState
{
public:
    enum render_type {
        RENDER_TEXT,
        RENDER_COLOR,
        RENDER_COLOR_TEX
    };

    enum shader_attrib {
        ATTR_VERTEX = 0,
        ATTR_TEXT_VERTEX,
        ATTR_NORMAL,
        ATTR_COLOR,
        ATTR_TEXCOORD0,
        ATTR_TEXCOORD1,
        NUM_ATTRS
    };

protected:
    render_type _shaderMode;
    bool _render_feedback = false;

    GLuint default_program;
    GLuint feedback_program;
    GLuint global_vao;

    int _w;
    int _h;

    //client state variables
    bool gl_lighting = false,
         gl_depth_test = false,
         gl_blend = false,
         gl_clip_plane = false;
    float _static_color[4];

    static const int MAX_LIGHTS = 3;
    //cached uniforms
    int _num_lights;
    float _ambient[4];
    Light _pt_lights[3];
    Material _mat;

    glm::vec4 _clip_plane;
    GlMatrix _projection_cp;

    //shader attribs
    bool _attr_enabled[NUM_ATTRS];

    //shader uniform locations
    GLuint locUseColorTex, locContainsText;
    GLuint locUseClipPlane, locClipPlane;
    GLuint locSpec, locShin;
    GLuint locModelView, locProject, locProjectText, locNormal;
    GLuint locNumLights, locGlobalAmb;
    GLuint locPosition[MAX_LIGHTS], locDiffuse[MAX_LIGHTS], locSpecular[MAX_LIGHTS];

    void initShaderState(GLuint program);
public:
    GlMatrix modelView;
    GlMatrix projection;

    GlState()
        : _shaderMode(RENDER_COLOR)
        , default_program(0)
        , feedback_program(0)
        , global_vao(0)
        , _ambient{0.2, 0.2, 0.2, 1.0}
        , _clip_plane(0.0, 0.0, 0.0, 0.0)
        , _attr_enabled{false} {
        modelView.identity();
        projection.identity();
    }

    ~GlState() {
        if (global_vao != 0) {
            glBindVertexArray(0);
            glDeleteVertexArrays(1, &global_vao);
        }
        if (default_program)
            glDeleteProgram(default_program);
        if (feedback_program)
            glDeleteProgram(feedback_program);
    }

    /**
     * Compiles the rendering pipeline shaders.
     */
    bool compileShaders();

    /**
     * Switches to the transform feedback rendering pipeline.
     */
    bool renderToFeedback() {
        if (feedback_program == 0) {
            return false;
        }
        if (_render_feedback != true) {
            glUseProgram(feedback_program);
            initShaderState(feedback_program);
            _render_feedback = true;
        }
        return true;
    }

    /**
     * Switches to the default rendering pipeline.
     */
    void renderToDefault() {
        if (_render_feedback != false) {
            glUseProgram(default_program);
            initShaderState(default_program);
            _render_feedback = false;
        }
    }

    void enableDepthTest() {
        if (!gl_depth_test) {
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            gl_depth_test = true;
        }
    }

    void disableDepthTest() {
        if (gl_depth_test) {
            glDisable(GL_DEPTH_TEST);
            gl_depth_test = false;
        }
    }

    void enableBlend() {
        if (!gl_blend) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            gl_blend = true;
        }
    }

    void disableBlend() {
        if (gl_blend) {
            glDisable(GL_BLEND);
            gl_blend = false;
        }
    }

    void enableClipPlane() {
        if (!gl_clip_plane) {
            glUniform1i(locUseClipPlane, GL_TRUE); 
            gl_clip_plane = true;
        }
    }

    void disableClipPlane() {
        if (gl_clip_plane) {
            glUniform1i(locUseClipPlane, GL_FALSE);
            gl_clip_plane = false;
        }
    }

    void setClipPlane(const double * eqn) {
        _clip_plane = glm::vec4(eqn[0], eqn[1], eqn[2], eqn[3]);
        _clip_plane = glm::inverseTranspose(modelView.mtx) * _clip_plane;
        glUniform4fv(locClipPlane, 1, glm::value_ptr(_clip_plane));
    }

    void setStaticColor(float r, float g, float b, float a = 1.0) {
        _static_color[0] = r;
        _static_color[1] = g;
        _static_color[2] = b;
        _static_color[3] = a;
        if (!_attr_enabled[ATTR_COLOR]) {
            glVertexAttrib4fv(ATTR_COLOR, _static_color);
        }
    }

    void enableAttribArray(GlState::shader_attrib attr) {
        if (!_attr_enabled[attr]) {
            glEnableVertexAttribArray(attr);
            _attr_enabled[attr] = true;
        }
    }

    void disableAttribArray(GlState::shader_attrib attr) {
        if (_attr_enabled[attr]) {
            glDisableVertexAttribArray(attr);
            _attr_enabled[attr] = false;
        }
        if (attr == ATTR_COLOR) {
            glVertexAttrib4fv(ATTR_COLOR, _static_color);
        } else if (attr == ATTR_NORMAL) {
            glVertexAttrib3f(ATTR_NORMAL, 0.f, 0.f, 1.f);
        }
    }

    /**
     * Loads the current model-view and projection matrix into the shader.
     */
    void loadMatrixUniforms(bool viewportChange = false) {
        if (viewportChange) {
            glm::mat4 proj2d = glm::ortho<float>(0, _w, 0, _h, -5.0, 5.0);
            glUniformMatrix4fv(locProjectText, 1, GL_FALSE, glm::value_ptr(proj2d));
        }
        glUniformMatrix4fv(locModelView, 1, GL_FALSE, glm::value_ptr(modelView.mtx));
        glUniformMatrix4fv(locProject, 1, GL_FALSE, glm::value_ptr(projection.mtx));
        glm::mat3 normal(modelView.mtx);
        normal = glm::inverseTranspose(normal);
        glUniformMatrix3fv(locNormal, 1, GL_FALSE, glm::value_ptr(normal));
    }

    /**
     * Sets the material parameters to use in lighting calculations.
     */
    void setMaterial(Material mat) {
        _mat = mat;
        glUniform4fv(locSpec, 1, mat.specular);
        glUniform1f(locShin, mat.shininess);
    }

    /**
     * Sets an individual point light's parameters.
     */
    void setLight(int i, Light lt) {
        if (i >= MAX_LIGHTS)
            return;
        _pt_lights[i] = lt;
        glUniform3fv(locPosition[i], 1, lt.position);
        glUniform4fv(locDiffuse[i], 1, lt.diffuse);
        glUniform4fv(locSpecular[i], 1, lt.specular);

    }

    void setLightPosition(int i, float * pos) {
        if (i >= MAX_LIGHTS)
            return;
        _pt_lights[i].position[0] = pos[0];
        _pt_lights[i].position[1] = pos[1];
        _pt_lights[i].position[2] = pos[2];
        glUniform3fv(locPosition[i], 1, pos); 
    }

    /**
     * Sets the number of point lights to include in the calculation.
     */
    void setNumLights(int n) {
        if (n >= MAX_LIGHTS)
            return;
        _num_lights = n;
        if (gl_lighting) {
            glUniform1i(locNumLights, _num_lights);
        }
    }

    /**
     * Sets the global ambient light intensity to use.
     */
    void setGlobalAmbLight(float * amb) {
        for (int i = 0; i < 4; i++) {
            _ambient[i] = amb[i];
        }
        if (gl_lighting) {
            glUniform4fv(locGlobalAmb, 1, amb);
        }
    }

    /**
     * Disables lighting in the shader, passing through colors directly.
     */
    bool disableLight() {
        if (gl_lighting) {
            glUniform1i(locNumLights, 0);
            gl_lighting = false;
            return true;
        }
        return false;
    }

    /**
     * Enables lighting in the shader, using the pre-set parameters.
     */
    void enableLight() {
        if (!gl_lighting) {
            glUniform1i(locNumLights, _num_lights);
            gl_lighting = true;
        }
    }

    /**
     * Prepares the shader pipeline for text rendering.
     */
    void setModeRenderText() {
        if (_shaderMode != RENDER_TEXT) {
            glDepthMask(GL_FALSE);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glUniform1i(locContainsText, GL_TRUE);
            _shaderMode = RENDER_TEXT;
        }
    }

    /**
     * Prepares the shader pipeline for rendering with color textures.
     */
    void setModeColorTexture(bool setUniforms = true) {
        if (_shaderMode == RENDER_TEXT) {
            if (gl_depth_test)
                glDepthMask(GL_TRUE);
            if (!gl_blend)
                glDisable(GL_BLEND);
            glUniform1i(locContainsText, GL_FALSE);
        }
        if (_shaderMode != RENDER_COLOR_TEX) {
            glUniform1i(locUseColorTex, GL_TRUE);
            _shaderMode = RENDER_COLOR_TEX;
        }
    }

    /**
     * Prepares the shader pipeline for rendering with color values.
     */
    void setModeColor(bool setUniforms = true) {
        if (_shaderMode == RENDER_TEXT) {
            if (gl_depth_test)
                glDepthMask(GL_TRUE);
            if (!gl_blend)
                glDisable(GL_BLEND);
            glUniform1i(locContainsText, GL_FALSE);
        }
        if (_shaderMode != RENDER_COLOR) {
            glUniform1i(locUseColorTex, GL_FALSE);
            _shaderMode = RENDER_COLOR;
        }
    }

    /**
     * Sets the viewport of the OpenGL context.
     */
    void setViewport(GLsizei w, GLsizei h) {
        this->_w = w;
        this->_h = h;
        glViewport(0, 0, w, h);
        loadMatrixUniforms(true);
    }

    void getViewport(GLint (&vp)[4]) {
        vp[0] = vp[1] = 0;
        vp[2] = _w;
        vp[3] = _h;
    }

    render_type getRenderMode() { return _shaderMode; }

    bool isClipPlaneEnabled() { return gl_clip_plane; }
};

#endif
