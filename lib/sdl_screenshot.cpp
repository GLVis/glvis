// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "sdl.hpp"
#include "mfem.hpp"

#if defined(GLVIS_USE_LIBTIFF)
#include "tiffio.h"
#elif defined(GLVIS_USE_LIBPNG)
#include <png.h>
#endif

#if defined(GLVIS_USE_LIBTIFF)
const char *glvis_screenshot_ext = ".tif";
#elif defined(GLVIS_USE_LIBPNG)
const char *glvis_screenshot_ext = ".png";
#else
const char *glvis_screenshot_ext = ".bmp";
#endif

// https://wiki.libsdl.org/SDL_CreateRGBSurfaceFrom
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
Uint32 rmask = 0xff000000;
Uint32 gmask = 0x00ff0000;
Uint32 bmask = 0x0000ff00:
               Uint32 amask = 0x000000ff;
#else // little endian, like x86
Uint32 rmask = 0x000000ff;
Uint32 gmask = 0x0000ff00;
Uint32 bmask = 0x00ff0000;
Uint32 amask = 0xff000000;
#endif

// https://halfgeek.org/wiki/Vertically_invert_a_surface_in_SDL
#define SDL_LOCKIFMUST(s) (SDL_MUSTLOCK(s) ? SDL_LockSurface(s) : 0)
#define SDL_UNLOCKIFMUST(s) { if(SDL_MUSTLOCK(s)) SDL_UnlockSurface(s); }

using namespace std;

int InvertSurfaceVertical(SDL_Surface *surface)
{
   Uint8 *t, *a, *b, *last;
   Uint16 pitch;

   if ( SDL_LOCKIFMUST(surface) < 0 )
   {
      return -2;
   }

   /* do nothing unless at least two lines */
   if (surface->h < 2)
   {
      SDL_UNLOCKIFMUST(surface);
      return 0;
   }

   /* get a place to store a line */
   pitch = surface->pitch;
   t = (Uint8*)malloc(pitch);

   if (t == NULL)
   {
      SDL_UNLOCKIFMUST(surface);
      return -2;
   }

   /* get first line; it's about to be trampled */
   memcpy(t,surface->pixels,pitch);

   /* now, shuffle the rest so it's almost correct */
   a = (Uint8*)surface->pixels;
   last = a + pitch * (surface->h - 1);
   b = last;

   while (a < b)
   {
      memcpy(a,b,pitch);
      a += pitch;
      memcpy(b,a,pitch);
      b -= pitch;
   }

   /* in this shuffled state, the bottom slice is too far down */
   memmove( b, b+pitch, last-b );

   /* now we can put back that first row--in the last place */
   memcpy(last,t,pitch);

   /* everything is in the right place; close up. */
   free(t);
   SDL_UNLOCKIFMUST(surface);

   return 0;
}

int SdlWindow::screenshotHelper(bool convert)
{
#ifdef GLVIS_DEBUG
   cout << "Screenshot: glFinish() ... " << flush;
#endif
   glFinish();
#ifdef GLVIS_DEBUG
   cout << "done." << endl;
#endif
#ifndef __EMSCRIPTEN__
   if (isExposePending())
   {
      MFEM_WARNING("Expose pending, some events may not have been handled." << endl);
   }
   string filename = screenshot_file.c_str();
   string convert_name = screenshot_file.c_str();
   bool call_convert = false;
   if (convert)
   {
      // check the extension of 'fname' to see if convert is needed
      size_t ext_size = strlen(glvis_screenshot_ext);
      if (filename.size() < ext_size ||
          filename.compare(filename.size() - ext_size,
                           ext_size, glvis_screenshot_ext) != 0)
      {
         call_convert = true;
         filename += glvis_screenshot_ext;
      }
   }
   else // do not call convert
   {
      filename += glvis_screenshot_ext;
   }

   int w, h;
   getGLDrawSize(w, h);
   if (isSwapPending())
   {
#if GLVIS_DEBUG
      cerr << "Screenshot: reading image data from back buffer..." << endl;
#endif
      glReadBuffer(GL_BACK);
   }
   else
   {
#if GLVIS_DEBUG
      cerr << "Screenshot: reading image data from front buffer..." << endl;
#endif
      glReadBuffer(GL_FRONT);
   }
#if defined(GLVIS_USE_LIBTIFF)
   // Save a TIFF image. This requires the libtiff library, see www.libtiff.org
   TIFF* image;

   // MyExpose(w,h);

   unsigned char *pixels = new unsigned char[3*w];
   if (!pixels)
   {
      return 1;
   }

   image = TIFFOpen(filename.c_str(), "w");
   if (!image)
   {
      delete [] pixels;
      return 2;
   }

   TIFFSetField(image, TIFFTAG_IMAGEWIDTH, w);
   TIFFSetField(image, TIFFTAG_IMAGELENGTH, h);
   TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
   TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
   TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
   TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 3);
   TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
   TIFFSetField(image, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
   TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   for (int i = 0; i < h; i++)
   {
      glReadPixels(0, h-i-1, w, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      if (TIFFWriteScanline(image, pixels, i, 0) < 0)
      {
         TIFFClose(image);
         delete [] pixels;
         return 3;
      }
   }

   TIFFFlushData(image);
   TIFFClose(image);
   delete [] pixels;

#elif defined(GLVIS_USE_LIBPNG)
   // Save as png image. Requires libpng.

   png_byte *pixels = new png_byte[3*w];
   if (!pixels)
   {
      return 1;
   }

   png_structp png_ptr =
      png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (!png_ptr)
   {
      delete [] pixels;
      return 1;
   }
   png_infop info_ptr = png_create_info_struct(png_ptr);
   if (!info_ptr)
   {
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      delete [] pixels;
      return 1;
   }

   FILE *fp = fopen(filename.c_str(), "wb");
   if (!fp)
   {
      png_destroy_write_struct(&png_ptr, &info_ptr);
      delete [] pixels;
      return 2;
   }

   if (setjmp(png_jmpbuf(png_ptr)))
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      delete [] pixels;
      return 3;
   }

   png_uint_32 ppi = isHighDpi() ? 144 : 72; // pixels/inch
   png_uint_32 ppm = ppi/0.0254 + 0.5;            // pixels/meter
   png_set_pHYs(png_ptr, info_ptr, ppm, ppm, PNG_RESOLUTION_METER);

   png_init_io(png_ptr, fp);
   png_set_IHDR(png_ptr, info_ptr, w, h, 8, PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT);

   png_write_info(png_ptr, info_ptr);
   for (int i = 0; i < h; i++)
   {
      glReadPixels(0, h-1-i, w, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      png_write_row(png_ptr, pixels);
   }
   png_write_end(png_ptr, info_ptr);

   fclose(fp);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   delete [] pixels;

#else
   // use SDL for screenshots

   // https://stackoverflow.com/questions/20233469/how-do-i-take-and-save-a-bmp-screenshot-in-sdl-2
   unsigned char * pixels = new unsigned char[w*h*4]; // 4 bytes for RGBA
   glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

   SDL_Surface * surf = SDL_CreateRGBSurfaceFrom(pixels, w, h, 8*4, w*4, rmask,
                                                 gmask, bmask, amask);
   if (surf == nullptr)
   {
      std::cerr << "unable to take screenshot: " << SDL_GetError() << std::endl;
   }
   else
   {
      if (InvertSurfaceVertical(surf))
      {
         std::cerr << "failed to invert surface, your screenshot may be upside down" <<
                   std::endl;
      }
      SDL_SaveBMP(surf, filename.c_str());
      SDL_FreeSurface(surf);
      // automatically convert to png if not being used
      if (!call_convert)
      {
         call_convert = true;
         convert_name += ".png";
      }
   }
   delete [] pixels;
#endif

   if (call_convert)
   {
      ostringstream cmd;
      cmd << "convert " << filename << ' ' << convert_name;
      if (system(cmd.str().c_str()))
      {
         return 1;
      }
      remove(filename.c_str());
   }
   return 0;
#else
   cout << "Screenshots not yet implemented for JS" << endl;
   return 1;
#endif
}
