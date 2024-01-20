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

#ifndef GLVIS_SHADER_HPP
#define GLVIS_SHADER_HPP

#include "renderer.hpp"

#include <unordered_map>

namespace gl3
{
class ShaderProgram
{
public:
   ShaderProgram()
   {
      program_id = glCreateProgram();
      if (program_id == 0)
      {
         std::cerr << "Failed to create an OpenGL program object." << std::endl;
      }
   }

   bool create(std::string vertexShader,
               std::string fragmentShader,
               std::unordered_map<int, std::string> inAttributes,
               int numOutputs);

   bool isCompiled() const { return is_compiled; }

   int uniform(std::string uniformName) const
   {
      auto unifId = uniform_idx.find(uniformName);
      if (unifId != uniform_idx.end())
      {
         return unifId->second;
      }
      return -1;
   }

   std::unordered_map<std::string, GLuint> getUniformMap() const
   {
      return uniform_idx;
   }

   GLuint getProgramId() const { return program_id; }

   void bind() const { glUseProgram(program_id); }

private:
   static void GetGLSLVersion();

   std::string formatShader(const std::string& inShader, GLenum shaderType);
   GLuint compileShader(const std::string& inShader, GLenum shaderType);
   bool linkShaders(const std::vector<GLuint>& shaders);
   void mapShaderUniforms();

   static int glsl_version;
   const static bool glsl_is_es;
   std::unordered_map<int, std::string> attrib_idx;
   int num_outputs = 0;

   ShaderPrgmHandle program_id = 0;
   ShaderHandle vertex_shader = 0;
   ShaderHandle fragment_shader = 0;
   bool is_compiled = false;
   std::unordered_map<std::string, GLuint> uniform_idx;
};

}

#endif // GLVIS_SHADER_HPP
