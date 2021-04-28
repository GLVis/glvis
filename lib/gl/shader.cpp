// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "shader.hpp"
#include <regex>

#ifdef __EMSCRIPTEN__
// TODO: for webgl glsl the version ends in "es": #version GLSL_VER es
const std::string GLSL_HEADER = "precision mediump float;\n";
#else
const std::string GLSL_HEADER = "#version GLSL_VER\n";
#endif

namespace gl3
{

bool ShaderProgram::create(std::string vertexShader,
                           std::string fragmentShader,
                           std::unordered_map<int, std::string> inAttributes,
                           int numOutputs)
{
   attrib_idx = inAttributes;
   num_outputs = numOutputs;
   is_compiled = false;
   GetGLSLVersion();

   std::string fmtVS = formatShader(vertexShader, GL_VERTEX_SHADER);
   vertex_shader = compileShader(fmtVS, GL_VERTEX_SHADER);
   if (vertex_shader == 0)
   {
      return false;
   }

   std::string fmtFS = formatShader(fragmentShader, GL_FRAGMENT_SHADER);
   fragment_shader = compileShader(fmtFS, GL_FRAGMENT_SHADER);
   if (fragment_shader == 0)
   {
      return false;
   }

   if (!linkShaders({vertex_shader, fragment_shader}))
   {
      std::cerr << "Failed to link shaders for program." << std::endl;
      return false;
   }

   mapShaderUniforms();

   is_compiled = true;
   return is_compiled;
}

void ShaderProgram::setOutputFramebuffer(const FBOHandle& fbo)
{
   glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
   if (num_outputs == 0)
   {
      std::cerr << "Warning: shader has no outputs." << std::endl;
      return;
   }
   vector<GLenum> output_bufs(num_outputs);
   for (int iout = 0; iout < num_outputs; iout++)
   {
      output_bufs[iout] = GL_COLOR_ATTACHMENT0 + iout;
   }
   glDrawBuffers(num_outputs, output_bufs.data());
}

void ShaderProgram::setDefaultDrawFramebuffer()
{
   if (num_outputs == 0)
   {
      std::cerr << "Warning: shader has no outputs." << std::endl;
   }
   if (num_outputs > 1)
   {
      std::cerr << "Warning: attempting to set a shader with more than one "
                << "output on the default framebuffer.";
   }
   glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
   GLenum back_buf = GL_BACK;
   glDrawBuffers(1, &back_buf);
}

int ShaderProgram::glsl_version = -1;
bool ShaderProgram::glsl_es = false;

void ShaderProgram::GetGLSLVersion()
{
   if (glsl_version == -1)
   {
      std::string verStr = (char*)glGetString(GL_VERSION);
      int ver_major, ver_minor;
      int vs_idx = verStr.find_first_of(".");
      ver_major = std::stoi(verStr.substr(vs_idx - 1, vs_idx));
      ver_minor = std::stoi(verStr.substr(vs_idx + 1, 1));
      int opengl_ver = ver_major * 100 + ver_minor * 10;

#ifndef __EMSCRIPTEN__
      // The GLSL version is the same as the OpenGL version when the OpenGL
      // version is >= 3.30, otherwise it is:
      //
      // GL Version | GLSL Version
      // -------------------------
      //    2.0     |   1.10
      //    2.1     |   1.20
      //    3.0     |   1.30
      //    3.1     |   1.40
      //    3.2     |   1.50

      if (opengl_ver < 330)
      {
         if (ver_major == 2)
         {
            glsl_version = opengl_ver - 90;
         }
         else if (ver_major == 3)
         {
            glsl_version = opengl_ver - 170;
         }
         else
         {
            std::cerr << "fatal: unsupported OpenGL version " << opengl_ver << std::endl;
            glsl_version = 100;
         }
      }
      else
      {
         glsl_version = opengl_ver;
      }
#else
      // GL Version | WebGL Version | GLSL Version
      //    2.0     |      1.0      |   1.00 ES
      //    3.0     |      2.0      |   3.00 ES
      //    3.1     |               |   3.10 ES
      if (opengl_ver < 300)
      {
         glsl_version = 100;
      }
      else
      {
         glsl_version = 300;
      }
      glsl_es = true;
#endif
      std::cerr << "Using GLSL " << glsl_version << std::endl;
   }
}

std::string ShaderProgram::formatShader(const std::string& inShader,
                                        GLenum shaderType)
{
   std::string formatted = inShader;

   // replace some identifiers depending on the version of glsl we're using
   if (glsl_version >= 130)
   {
      if (shaderType == GL_VERTEX_SHADER)
      {
         formatted = std::regex_replace(formatted, std::regex("attribute"), "in");
         formatted = std::regex_replace(formatted, std::regex("varying"), "out");
      }

      else if (shaderType == GL_FRAGMENT_SHADER)
      {
         formatted = std::regex_replace(formatted, std::regex("varying"), "in");

         // requires GL_ARB_explicit_attrib_location extension or GLSL 3.30
         // although gl_FragColor was deprecated in GLSL 1.3
         for (int i = 0; i < num_outputs; i++)
         {
            std::string indexString = "gl_FragData\\[";
            indexString += std::to_string(i) + "\\]";
            std::string outputString = "out vec4 fragColor_";
            outputString += std::to_string(i) + ";\n";
            if (glsl_version > 130)
            {
               formatted = std::regex_replace(formatted, std::regex(indexString),
                                              "fragColor_" + std::to_string(i));
            }
            if (glsl_version > 130 && glsl_version < 330)
            {
               // append output variable specifications without layout specifiers
               formatted = outputString + formatted;
            }
            else if (glsl_version >= 330 || (glsl_version == 300 && glsl_es))
            {
               std::string layoutString = "layout(location = ";
               layoutString += std::to_string(i) + ") ";
               formatted = layoutString + outputString + formatted;
            }
            if (i == 0)
            {
               formatted = std::regex_replace(formatted, std::regex("gl_FragColor"),
                                              "fragColor_0");
            }
         }
      }

      else
      {
         std::cerr << "buildShaderString: unknown shader type" << std::endl;
         return {};
      }

      formatted = std::regex_replace(formatted, std::regex("texture2D"), "texture");
   }

   if (!GLDevice::isOpenGL3())
   {
      formatted = "#define USE_ALPHA\n" + formatted;
   }
   // add the header
   formatted = std::regex_replace(GLSL_HEADER, std::regex("GLSL_VER"),
                                  std::to_string(glsl_version)) + formatted;
#ifdef __EMSCRIPTEN__
   // special prepend for WebGL 2 shaders
   if (glsl_version == 300)
   {
      formatted = "#version 300 es\n" + formatted;
   }
#endif

   return formatted;
}

GLuint ShaderProgram::compileShader(const std::string& inShader,
                                    GLenum shaderType)
{
   int shader_len = inShader.length();
   const char *shader_cstr = inShader.c_str();

   GLuint shader = glCreateShader(shaderType);
   glShaderSource(shader, 1, &shader_cstr, &shader_len);
   glCompileShader(shader);

   GLint stat;
   glGetShaderiv(shader, GL_COMPILE_STATUS, &stat);
   // glGetObjectParameteriv
   if (stat == GL_FALSE)
   {
      std::cerr << "Failed to compile shader" << std::endl;
      int err_len;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &err_len);
      char *error_text = new char[err_len];
      glGetShaderInfoLog(shader, err_len, &err_len, error_text);
      std::cerr << error_text << std::endl;
      delete[] error_text;
      return 0;
   }
   return shader;
}

bool ShaderProgram::linkShaders(const std::vector<GLuint>& shaders)
{
   program_id = glCreateProgram();
   if (program_id == 0)
   {
      std::cerr << "Failed to create an OpenGL program object." << std::endl;
   }
   // Bind all incoming attributes to their VAO indices.
   for (auto attrib_pair : attrib_idx)
   {
      glBindAttribLocation(program_id, attrib_pair.first,
                           attrib_pair.second.c_str());
   }

#ifndef __EMSCRIPTEN__
   // Bind fragment output variables to MRT indices.
   for (int i = 0; i < num_outputs; i++)
   {
      std::string fragOutVar = "fragColor_" + std::to_string(i);
      glBindFragDataLocation(program_id, i, fragOutVar.c_str());
   }
#endif

   for (GLuint i : shaders)
   {
      glAttachShader(program_id, i);
   }
   glLinkProgram(program_id);

   GLint stat;
   glGetProgramiv(program_id, GL_LINK_STATUS, &stat);
   if (stat == GL_FALSE)
   {
      cerr << "fatal: Shader linking failed" << endl;
   }
   return (stat == GL_TRUE);
}

void ShaderProgram::mapShaderUniforms()
{
   int num_unifs;
   glGetProgramiv(program_id, GL_ACTIVE_UNIFORMS, &num_unifs);
   int max_unif_len;
   glGetProgramiv(program_id, GL_ACTIVE_UNIFORM_MAX_LENGTH, &max_unif_len);

   for (int i = 0; i < num_unifs; i++)
   {
      vector<char> unif_buf(max_unif_len+1);
      GLsizei name_length;
      GLint gl_size;
      GLenum gl_type;
      glGetActiveUniform(program_id, i, max_unif_len, &name_length, &gl_size,
                         &gl_type,
                         unif_buf.data());
      std::string unif_name(unif_buf.data(), name_length);
      GLuint location = glGetUniformLocation(program_id, unif_name.c_str());
      uniform_idx[unif_name] = location;
   }
}

}
