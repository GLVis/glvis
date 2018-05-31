#include "font.hpp"
#include "aux_vis.hpp"

int GLVisFont::SetFontFile(const char *font_file, int height)
{
  int err;

  if (init == -1) // library init failed
  {
     return -1;
  }

  if (init <= 0)
  {
     err = FT_Init_FreeType(&library);
     if (err)
     {
        cout << "GLVis: Can not initialize FreeType library!" << endl;
        return (init = -1);
     }
  }
  else
  {
     FreeSeq();
     delete [] image;
     FT_Done_Face(face);
  }

  err = FT_New_Face(library, font_file, 0, &face);
  if (err)
  {
     cout << "GLVis: Can not open font file: " << font_file << endl;
     FT_Done_FreeType(library);
     return (init = -2);
  }

  if (1)
  {
     // set font height in points
     int ppi_w, ppi_h;
     GetAppWindow()->getDpi(ppi_w, ppi_h);
     err = FT_Set_Char_Size(face, 0, height*64, ppi_w, ppi_h);
     if (err)
     {
        cout << "GLVis: Can not set font height: " << height << " pts"
             << endl;
        FT_Done_Face(face);
        FT_Done_FreeType(library);
        return (init = -3);
     }
  }
  else
  {
     // set font height in pixels
     err = FT_Set_Pixel_Sizes(face, 0, height);
     if (err)
     {
        cout << "GLVis: Can not set font height: " << height << " pixels"
             << endl;
        FT_Done_Face(face);
        FT_Done_FreeType(library);
        return (init = -3);
     }
  }

  glyph = NULL;
  use_kerning = FT_HAS_KERNING(face);

  image = NULL;

#ifdef GLVIS_OGL3
#endif

  init = 1;

  return 0;
}


