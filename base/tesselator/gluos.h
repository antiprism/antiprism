/*
** gluos.h - operating system dependencies for GLU
**
** $Header: //depot/main/gfx/lib/localglu/include/gluos.h#4 $
*/

#ifdef _WIN32

#if !defined(WIN32_LEAN_AND_MEAN)
#define WIN32_LEAN_AND_MEAN
#endif

#define NOGDI
#define NOIME
#include <windows.h>

#define localGLAPI

#else

/* Disable Microsoft-specific keywords */
#define WINGDIAPI
#define localGLAPI

#endif
