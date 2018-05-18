/*
 * (c) Copyright 1993, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <GL/glx.h>
#include <GL/glew.h>
#include "tk.h"
#include "aux_gl.hpp"

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif


static struct {
    int keyField;
    void (*KeyFunc)(void);
} keyTable[200];

static struct {
    int keyField;
    void (*KeyFunc)(GLenum);
} modkeyTable[200];

static struct {
    unsigned int mouseField;
    void (*MouseFunc)(AUX_EVENTREC *);
} mouseDownTable[20], mouseUpTable[20], mouseLocTable[20];

static int keyTableCount = 0;
static int modkeyTableCount = 0;
static int mouseDownTableCount = 0;
static int mouseUpTableCount = 0;
static int mouseLocTableCount = 0;
static GLenum displayModeType = 0;


static void DefaultHandleReshape(GLsizei w, GLsizei h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, (GLdouble)w, 0.0, (GLdouble)h, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

static void DefaultHandleExpose(GLsizei w, GLsizei h)
{
}

static GLenum MouseLoc(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseLocTableCount; i++) {
	if ((button & Button1Mask) &&
            AUX_LEFTBUTTON == mouseLocTable[i].mouseField) {

            info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2] = button;
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
	if ((button & Button3Mask) &&
            AUX_RIGHTBUTTON == mouseLocTable[i].mouseField) {

            info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2] = button;
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
	if ((button & Button2Mask) &&
            AUX_MIDDLEBUTTON == mouseLocTable[i].mouseField) {

            info.event = AUX_MOUSELOC;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2] = button;
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseLocTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum MouseUp(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseUpTableCount; i++) {
        if ((button & AUX_LEFTBUTTON) == mouseUpTable[i].mouseField) {
            info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((button & AUX_RIGHTBUTTON) == mouseUpTable[i].mouseField) {
            info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((button & AUX_MIDDLEBUTTON) == mouseUpTable[i].mouseField) {
            info.event = AUX_MOUSEUP;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseUpTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum MouseDown(int x, int y, GLenum button)
{
    AUX_EVENTREC info;
    GLenum flag;
    int i;

    flag = GL_FALSE;
    for (i = 0; i < mouseDownTableCount; i++) {
        if ((button & AUX_LEFTBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_LEFTBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((button & AUX_RIGHTBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_RIGHTBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
        if ((button & AUX_MIDDLEBUTTON) == mouseDownTable[i].mouseField) {
	    info.event = AUX_MOUSEDOWN;
	    info.data[AUX_MOUSEX] = x;
	    info.data[AUX_MOUSEY] = y;
	    info.data[2]          = (button>>3);
	    info.data[AUX_MOUSESTATUS] = AUX_MIDDLEBUTTON;
	    (*mouseDownTable[i].MouseFunc)(&info);
	    flag |= GL_TRUE;
	}
    }
    return flag;
}

static GLenum KeyDown(int key, GLenum status)
{
    GLenum flag;
    int i;

    flag = GL_FALSE;
    if (keyTableCount) {
	for (i = 0; i < keyTableCount; i++) {
	    if (key == keyTable[i].keyField) {
		(*keyTable[i].KeyFunc)();
		flag |= GL_TRUE;
	    }
	}
    }
    if (modkeyTableCount) {
	for (i = 0; i < modkeyTableCount; i++) {
	    if (key == modkeyTable[i].keyField) {
		(*modkeyTable[i].KeyFunc)(status);
		flag |= GL_TRUE;
	    }
	}
    }
    return flag;
}

void auxExposeFunc(void (*Func)(GLsizei, GLsizei))
{
    tkExposeFunc(Func);
}

void auxReshapeFunc(void (*Func)(GLsizei, GLsizei))
{
    tkExposeFunc(Func);
    tkReshapeFunc(Func);
}

void auxIdleFunc(void (*Func)(void))
{
    tkIdleFunc(Func);
}

void auxKeyFunc(int key, void (*Func)(void))
{
    keyTable[keyTableCount].keyField = key;
    keyTable[keyTableCount++].KeyFunc = Func;
}

void auxModKeyFunc(int key, void (*Func)(GLenum))
{
   modkeyTable[modkeyTableCount].keyField = key;
   modkeyTable[modkeyTableCount++].KeyFunc = Func;
}

void auxKeyFuncReplace(int key, void (*Func)(void))
{
   int n = 0;
   for (int i = 0; i < keyTableCount; i++)
      if (keyTable[i].keyField == key)
         if (++n == 1)
            keyTable[i].KeyFunc = Func;
   if (n > 1)
      printf("auxKeyFuncReplace : A key is assigned multiple functions!\n");
}

void auxCallKeyFunc(int key, GLenum status)
{
   KeyDown(key, status);
}

void auxMouseFunc(int mouse, int mode, void (*Func)(AUX_EVENTREC *))
{
    if (mode == AUX_MOUSEDOWN) {
	mouseDownTable[mouseDownTableCount].mouseField = mouse;
	mouseDownTable[mouseDownTableCount++].MouseFunc = Func;
    } else if (mode == AUX_MOUSEUP) {
	mouseUpTable[mouseUpTableCount].mouseField = mouse;
	mouseUpTable[mouseUpTableCount++].MouseFunc = Func;
    } else if (mode == AUX_MOUSELOC) {
	mouseLocTable[mouseLocTableCount].mouseField = mouse;
	mouseLocTable[mouseLocTableCount++].MouseFunc = Func;
    }
}

void auxMainLoop(void (*Func)(void))
{
    tkDisplayFunc(Func);
    tkExec();
}

void auxInitPosition(int x, int y, int width, int height)
{
    tkInitPosition(x, y, width, height);
}

void auxInitDisplayMode(GLenum type)
{
    displayModeType = type;
    tkInitDisplayMode(type);
}

GLenum auxInitWindow(const char *title)
{
    int useDoubleAsSingle = 0;

    if (tkInitWindow(title) == GL_FALSE)
    {
	if (AUX_WIND_IS_SINGLE(displayModeType))
        {
	    tkInitDisplayMode(displayModeType | AUX_DOUBLE);
	    if (tkInitWindow(title) == GL_FALSE)
            {
		return GL_FALSE;    /*  curses, foiled again	*/
	    }
	    fprintf(stderr, "Can't initialize a single buffer visual.\n");
	    fprintf(stderr, "Will use a double buffer visual instead,");
	    fprintf(stderr, "only drawing into the front buffer.\n");
	    displayModeType = displayModeType | AUX_DOUBLE;
	    useDoubleAsSingle = 1;
	}
        else
        {
           return GL_FALSE;
        }
    }
    tkReshapeFunc(DefaultHandleReshape);
    tkExposeFunc(DefaultHandleExpose);
    tkMouseUpFunc(MouseUp);
    tkMouseDownFunc(MouseDown);
    tkMouseMoveFunc(MouseLoc);
    tkKeyDownFunc(KeyDown);
    // auxKeyFunc(AUX_ESCAPE, auxQuit);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClearIndex(0);
    glLoadIdentity();
    if (useDoubleAsSingle)
	glDrawBuffer(GL_FRONT);
    return GL_TRUE;
}

void auxCloseWindow(void)
{
    tkCloseWindow();
    keyTableCount = 0;
    modkeyTableCount = 0;
    mouseDownTableCount = 0;
    mouseUpTableCount = 0;
    mouseLocTableCount = 0;
}

void auxQuit(void)
{
    tkQuit();
}

void auxSwapBuffers(void)
{
    tkSwapBuffers();
}

Display *auxXDisplay(void)
{

    return (tkGetXDisplay());
}

Window auxXWindow(void)
{

    return (tkGetXWindow());
}

void auxSetOneColor(int index, float r, float g, float b)
{
    tkSetOneColor(index, r, g, b);
}

void auxSetFogRamp(int density, int startIndex)
{
    tkSetFogRamp(density, startIndex);
}

void auxSetGreyRamp(void)
{
    tkSetGreyRamp();
}

void auxSetRGBMap(int size, float *rgb)
{
    tkSetRGBMap(size, rgb);
}

GLint auxGetColorMapSize(void)
{

    return tkGetColorMapSize();
}

void auxGetMouseLoc(int *x, int *y)
{
    tkGetMouseLoc(x, y);
}
