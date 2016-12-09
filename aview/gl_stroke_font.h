/*
 * freeglut_internal.h
 *
 * The freeglut library private include file.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 2 1999
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef GL_STROKE_FONT_H
#define GL_STROKE_FONT_H

/* The stroke font structures */

typedef struct taganti_StrokeVertex anti_StrokeVertex;
struct taganti_StrokeVertex {
  float X, Y;
};

typedef struct taganti_StrokeStrip anti_StrokeStrip;
struct taganti_StrokeStrip {
  int Number;
  const anti_StrokeVertex *Vertices;
};

typedef struct taganti_StrokeChar anti_StrokeChar;
struct taganti_StrokeChar {
  float Right;
  int Number;
  const anti_StrokeStrip *Strips;
};

typedef struct taganti_StrokeFont anti_StrokeFont;
struct taganti_StrokeFont {
  const char *Name;                   /* The source font name      */
  int Quantity;                       /* Number of chars in font   */
  float Height;                       /* Height of the characters  */
  const anti_StrokeChar **Characters; /* The characters mapping    */
};

extern const anti_StrokeFont antiStrokeRoman;
extern const anti_StrokeFont antiStrokeMonoRoman;

/*
 * Those pointers will be used by following definitions:
 */
#define ANTI_STROKE_ROMAN (&antiStrokeRoman)
#define ANTI_STROKE_MONO_ROMAN (&antiStrokeMonoRoman)

void antiStrokeCharacter(const anti_StrokeFont *font, int character);
int antiStrokeWidth(const anti_StrokeFont *font, int character);
int antiStrokeLength(const anti_StrokeFont *font, const unsigned char *string);
float antiStrokeHeight(const anti_StrokeFont *fontID);

#endif // GL_STROKE_FONT_H
