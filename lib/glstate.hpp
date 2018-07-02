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
        ATTR_VERTEX,
        ATTR_TEXT_VERTEX,
        ATTR_COLOR,
        ATTR_NORMAL,
        ATTR_TEXCOORD0,
        ATTR_TEXCOORD1,
        NUM_ATTRS
    };

protected:
    render_type _shaderMode;

    GLuint program;

    int _w;
    int _h;

    //client state variables
    bool _lighting = false,
         _depthTest = false,
         _blend = false,
         _colorMat = false,
         _clipPlane = false;
    float _staticColor[4];

    //lighting uniforms
    int _numLights;
    float _ambient[4];

    glm::vec4 _clip_plane;
    GlMatrix _projection_cp;

    //shader locations
    int _attr_locs[NUM_ATTRS];
    bool _attr_enabled[NUM_ATTRS];

    GLuint locUseClipPlane, locClipPlane;
    GLuint locAmb, locDif, locSpec, locShin;
    GLuint locModelView, locProject, locProjectText, locNormal;

    glm::vec3 getRasterPoint(double x, double y, double z) {
        return glm::project(glm::vec3(x, y, z), modelView.mtx, projection.mtx, glm::vec4(0, 0, _w, _h));
    }

public:
    GlMatrix modelView;
    GlMatrix projection;

    void enableDepthTest() {
        if (!_depthTest) {
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            _depthTest = true;
        }
    }

    void disableDepthTest() {
        if (_depthTest) {
            glDisable(GL_DEPTH_TEST);
            _depthTest = false;
        }
    }

    void enableBlend() {
        if (!_blend) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            _blend = true;
        }
    }

    void disableBlend() {
        if (_blend) {
            glDisable(GL_BLEND);
            _blend = false;
        }
    }

    /**
     * Sets the ambient and diffuse lighting to full intensity, tracking the object color.
     */
    void enableColorMaterial() { _colorMat = true; }

    void disableColorMaterial() { _colorMat = false; }

    void enableClipPlane() {
        if (!_clipPlane) {
            glUniform1i(locUseClipPlane, GL_TRUE); 
            _clipPlane = true;
        }
    }

    void disableClipPlane() {
        if (_clipPlane) {
            glUniform1i(locUseClipPlane, GL_FALSE);
            _clipPlane = false;
        }
    }

    void setClipPlane(const double * eqn) {
        _clip_plane = glm::vec4(eqn[0], eqn[1], eqn[2], eqn[3]);
        _clip_plane = glm::inverseTranspose(modelView.mtx) * _clip_plane;
        glUniform4fv(locClipPlane, 1, glm::value_ptr(_clip_plane));
    }

    void setStaticColor(float r, float g, float b, float a = 1.0) {
        _staticColor[0] = r;
        _staticColor[1] = g;
        _staticColor[2] = b;
        _staticColor[3] = a;
        if (!_attr_enabled[ATTR_COLOR]) {
            glVertexAttrib4fv(_attr_locs[ATTR_COLOR], _staticColor);
        }
    }

    GlState()
        : program(0),
         _ambient{0.2, 0.2, 0.2, 1.0} {
    }

    ~GlState() {
        if (program)
            glDeleteProgram(program);
    }

    /**
     * Compiles the rendering pipeline shaders.
     */
    bool compileShaders();

    /**
     * Initializes samplers and uniforms in the shader pipeline.
     */
    void initShaderState();

    int getAttribLoc(GlState::shader_attrib attr) {
        return _attr_locs[attr];
    }

    void enableAttribArray(GlState::shader_attrib attr) {
        glEnableVertexAttribArray(_attr_locs[attr]);
        _attr_enabled[attr] = true;
    }

    //TODO: if color, use glVertexAttrib
    void disableAttribArray(GlState::shader_attrib attr) {
        glDisableVertexAttribArray(_attr_locs[attr]);
        _attr_enabled[attr] = false;
        if (attr == ATTR_COLOR) {
            glVertexAttrib4fv(_attr_locs[ATTR_COLOR], _staticColor);
        } else if (attr == ATTR_NORMAL) {
            glVertexAttrib3f(_attr_locs[ATTR_NORMAL], 0.f, 0.f, 1.f);
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
     * Set temporary model-view matrix directly in the shader.
     */
    void setImmModelView(GlMatrix modelview) {
        glUniformMatrix4fv(locModelView, 1, GL_FALSE, glm::value_ptr(modelview.mtx));
        glm::mat3 normal(modelview.mtx);
        normal = glm::inverseTranspose(normal);
        glUniformMatrix3fv(locNormal, 1, GL_FALSE, glm::value_ptr(normal));
    }

    /**
     * Set temporary projection matrix directly in the shader.
     */
    void setImmProjection(GlMatrix proj) {
        glUniformMatrix4fv(locProject, 1, GL_FALSE, glm::value_ptr(proj.mtx));
    }

    /**
     * Sets the material parameters to use in lighting calculations.
     */
    void setMaterial(Material mat) {
        if (!_colorMat) {
            glUniform4fv(locAmb, 1, mat.ambient);
            glUniform4fv(locDif, 1, mat.diffuse);
        } else {
            float fullIntens[4] = {1.0, 1.0, 1.0, 1.0};
            glUniform4fv(locAmb, 1, fullIntens);
            glUniform4fv(locDif, 1, fullIntens);
        }
        glUniform4fv(locSpec, 1, mat.specular);
        glUniform1f(locShin, mat.shininess);
    }

    /**
     * Sets an individual point light's parameters.
     */
    void setLight(int i, Light lt) {
        std::string location = "lights[" + std::to_string(i) + "].";
        GLuint locPosition = glGetUniformLocation(program, (location + "position").c_str());
        GLuint locDiffuse = glGetUniformLocation(program, (location + "diffuse").c_str());
        GLuint locSpecular = glGetUniformLocation(program, (location + "specular").c_str());
        glUniform3fv(locPosition, 1, lt.position);
        glUniform4fv(locDiffuse, 1, lt.diffuse);
        glUniform4fv(locSpecular, 1, lt.specular);
    }

    void setLightPosition(int i, float * pos) {
        std::string location = "lights[" + std::to_string(i) + "].";
        GLuint locPosition = glGetUniformLocation(program, (location + "position").c_str());
        glUniform3fv(locPosition, 1, pos); 
    }

    /**
     * Sets the number of point lights to include in the calculation.
     */
    void setNumLights(int n) {
        _numLights = n;
        if (_lighting) {
            GLuint locNumLights = glGetUniformLocation(program, "numLights");
            glUniform1i(locNumLights, _numLights);
        }
    }

    /**
     * Sets the global ambient light intensity to use.
     */
    void setGlobalAmbLight(float * amb) {
        for (int i = 0; i < 4; i++) {
            _ambient[i] = amb[i];
        }
        if (_lighting) {
            GLuint locGLight = glGetUniformLocation(program, "g_ambient");
            glUniform4fv(locGLight, 1, amb);
        }
    }

    /**
     * Disables lighting in the shader, passing through colors directly.
     */
    bool disableLight() {
        if (_lighting) {
            GLuint locNumLights = glGetUniformLocation(program, "numLights");
            glUniform1i(locNumLights, 0);
            _lighting = false;
            return true;
        }
        return false;
    }

    /**
     * Enables lighting in the shader, using the pre-set parameters.
     */
    void enableLight() {
        if (!_lighting) {
            GLuint locGLight = glGetUniformLocation(program, "g_ambient");
            GLuint locNumLights = glGetUniformLocation(program, "numLights");
            glUniform4fv(locGLight, 1, _ambient);
            glUniform1i(locNumLights, _numLights);
            _lighting = true;
        }
    }

    /**
     * Prepares the shader pipeline for text rendering.
     */
    void setModeRenderText(double x, double y, double z) {
        setModeRenderText();
        disableAttribArray(ATTR_VERTEX);
        glVertexAttrib3f(_attr_locs[ATTR_VERTEX], x, y, z);
    }

    /**
     * Prepares the shader pipeline for text rendering.
     */
    void setModeRenderText() {
        if (_shaderMode != RENDER_TEXT) {
            //glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            GLuint locContainsText = glGetUniformLocation(program, "containsText");
            glUniform1i(locContainsText, GL_TRUE);
            _shaderMode = RENDER_TEXT;
        }
    }

    /**
     * Prepares the shader pipeline for rendering with color textures.
     */
    void setModeColorTexture(bool setUniforms = true) {
        if (_shaderMode == RENDER_TEXT) {
            if (_depthTest)
                glEnable(GL_DEPTH_TEST);
            if (!_blend)
                glDisable(GL_BLEND);
            GLuint locContainsText = glGetUniformLocation(program, "containsText");
            glUniform1i(locContainsText, GL_FALSE);
        }
        if (_shaderMode != RENDER_COLOR_TEX) {
            GLuint locUseColorTex = glGetUniformLocation(program, "useColorTex");
            glUniform1i(locUseColorTex, GL_TRUE);
            _shaderMode = RENDER_COLOR_TEX;
            if (setUniforms) {
                loadMatrixUniforms();
            }
        }
    }

    /**
     * Prepares the shader pipeline for rendering with color values.
     */
    void setModeColor(bool setUniforms = true) {
        if (_shaderMode == RENDER_TEXT) {
            if (_depthTest)
                glEnable(GL_DEPTH_TEST);
            if (!_blend)
                glDisable(GL_BLEND);
            GLuint locContainsText = glGetUniformLocation(program, "containsText");
            glUniform1i(locContainsText, GL_FALSE);
        }
        if (_shaderMode != RENDER_COLOR) {
            GLuint locUseColorTex = glGetUniformLocation(program, "useColorTex");
            glUniform1i(locUseColorTex, GL_FALSE);
            _shaderMode = RENDER_COLOR;
            if (setUniforms) {
                loadMatrixUniforms();
            }
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

    render_type getRenderMode() {
        return _shaderMode;
    }
};

#endif
