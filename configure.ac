#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(antiprism, 0.25, adrian@antiprism.com)
AC_SUBST(RELEASE_DATE, ["23 November 2017"])
AC_CANONICAL_TARGET  # For m4/fltk.m4
AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([base/antiprism.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

AC_LANG([C++])
AX_CXX_COMPILE_STDCXX(11)

gl_LD_VERSION_SCRIPT
LT_INIT([win32-dll])

AC_CHECK_FUNC([gettimeofday], [utimer=1], [utimer=0])
if test $utimer = 0 ; then
   AC_MSG_CHECKING([for timeGetTime in winmm])
   my_ac_save_LIBS=$LIBS
   LIBS="-lwinmm $LIBS"
   AC_TRY_LINK([#include <windef.h>
             #include <mmsystem.h>
             ],
             [timeGetTime();],
             [ dnl worked!
             utimer=2; AC_MSG_RESULT([yes])
             ],
             [ dnl failed, restore LIBS
             LIBS=$my_ac_save_LIBS; utimer=0; AC_MSG_RESULT(no);
             AC_MSG_ERROR([no suitable subsecond time function found])
             ],            )
fi

AC_CHECK_LIB([m], [acos])

NO_GLUT=0
GLUT=1
OPENGLUT=2
FREEGLUT=3
FLTKGLUT=4

AC_PATH_X
if test x"$no_x" = x"yes" ; then
   AC_CANONICAL_HOST
   case $host_os in
      *cygwin* ) CPPFLAGS="-I/usr/include/opengl $CPPFLAGS"
                 LDFLAGS="-L/usr/lib/w32api $LDFLAGS"
                 ;;
   esac
fi


glut_type=$NO_GLUT
AC_ARG_WITH(
   [openglut],
   [AS_HELP_STRING([--with-openglut], [build using OpenGLUT])],
   [glut_type=$OPENGLUT],
   [AC_ARG_WITH([freeglut],
         [AS_HELP_STRING([--with-freeglut], [build using Freeglut])],
         [glut_type=$FREEGLUT],
         [AC_ARG_WITH([fltkglut],
            [AS_HELP_STRING([--with-fltkglut], [build using FLTK GLUT])],
            [glut_type=$FLTKGLUT])]
         )]
)

FOUND_GLU=0
AX_CHECK_GL
if test ! x"$no_gl" = x"yes" ; then
   AX_CHECK_GLU
   if test ! x"$no_glu" = x"yes" ; then
      FOUND_GLU=1
      GLPROG_LIBS="$GLU_LIBS $GL_LIBS $LIBS"
      GLPROG_CXXFLAGS="$GLU_CFLAGS $LIBS"
      if test $glut_type = $NO_GLUT ; then
         AX_CHECK_GLUT
         if test ! x"$no_glut" = x"yes" ; then
            glut_type=$GLUT
         fi
      else
         if test $glut_type = $OPENGLUT ; then
            glut_type=$NO_GLUT
            AC_CHECK_HEADER([GL/openglut.h],     # must be in GL directory
               [AC_CHECK_LIB([openglut], [main],
                  [GLUT_LIBS="-lopenglut"; glut_type=$OPENGLUT]) ])
         else
            if test $glut_type = $FREEGLUT ; then
               glut_type=$NO_GLUT
               AC_CHECK_HEADER([GL/freeglut.h],  # must be in GL directory
                  [AC_CHECK_LIB([freeglut], [main],
                     [GLUT_LIBS="-lfreeglut"; glut_type=$FREEGLUT]) ])
            else
               if test $glut_type = $FLTKGLUT ; then
                  AM_PATH_FLTK(1.3)
                  glut_type=$NO_GLUT
                  if test ! CXXFLAGS = "" ; then
                     GLUT_CXXFLAGS="$FLTK_CXXFLAGS"
                     GLUT_LDFLAGS= "$FLTK_LDFLAGS"
                     GLUT_LIBS="$FLTK_LIBS"
                     glut_type=$FLTKGLUT
                  fi
               fi
            fi
         fi
      fi
   fi
fi

if test $FOUND_GLU = 0 ; then
   echo
   echo "*********************************************"
   echo "   missing some OpenGL or GLU library or header"
   echo "   check OpenGL and GLU development packages are installed"
   echo "*********************************************"
   echo
fi

AH_TEMPLATE([FOUND_GLU], [The GL and GLU libraries and headers were found])
AC_DEFINE_UNQUOTED([FOUND_GLU], [$FOUND_GLU])
AC_SUBST(FOUND_GLU)
AC_SUBST(GLPROG_LIBS)
AC_SUBST(GLPROG_CXXFLAGS)
AC_SUBST(HAVE_OPENGL_GL_H)
AC_SUBST(HAVE_GL_GL_H)
AC_SUBST(HAVE_OPENGL_GLU_H)
AC_SUBST(HAVE_GL_GLU_H)
AC_SUBST(HAVE_GLUT_GLUT_H)
AC_SUBST(HAVE_GL_GLUT_H)


if test $glut_type = $NO_GLUT ; then
   echo
   echo "*********************************************"
   echo "   missing GLUT library or header"
   echo "   check GLUT development package is installed"
   echo "*********************************************"
   echo
fi

if test ! $FOUND_GLU = 0 && test ! $glut_type = 0 ; then
   AM_CONDITIONAL(BUILD_ANTIVIEW, true)
   echo "OpenGL/GLUT found, antiview will be built"
else
   AM_CONDITIONAL(BUILD_ANTIVIEW, false)
   echo
   echo "*********************************************"
   echo "   antiview will not be built"
   echo "   missing OpenGL or GLUT components (details above)"
   echo "*********************************************"
   echo
fi

AC_SUBST(GLUT_LIBS)
AC_SUBST(GLUT_CFLAGS)
AC_SUBST(GLUT_CXXFLAGS)
AC_SUBST(GLUT_LDLAGS)

AH_TEMPLATE([GLUT_TYPE], [The GLUT library that will used 1=GLUT, 2=OpenGLUT])
AC_DEFINE_UNQUOTED([GLUT_TYPE], [$glut_type])


# Timer

AH_TEMPLATE([UTIMER], [usec timer to use 1=gettimeofday, 2=timeGetTime])
AC_DEFINE_UNQUOTED([UTIMER], [$utimer])
UTIMER=$utimer
AC_SUBST(UTIMER)

AC_CHECK_FUNC([usleep], [usleep=1],
   [AC_CHECK_HEADER([winbase.h], [usleep=2],
      [AC_MSG_ERROR([no suitable subsecond sleep function found])], [#include <windows.h>])
   ])
AH_TEMPLATE([USLEEP], [usec sleep function to use 1=usleep, 2=Sleep])
AC_DEFINE_UNQUOTED([USLEEP], [$usleep])
USLEEP=$usleep
AC_SUBST(USLEEP)

# Install type, neded for Debian package
AC_ARG_ENABLE(debian,
   [  --enable-debian    Install compatible with debian packaging],
   [case "${enableval}" in
       yes) debian=true ;;
       no)  debian=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
    esac],
    [debian=false])
AM_CONDITIONAL(DEBIAN_INSTALL, test x$debian = xtrue)


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset modf pow sqrt strcasecmp strchr strcspn strncasecmp strpbrk strrchr strspn strstr strtol])

AC_CONFIG_FILES([Makefile
                 Doxyfile
                 aview/Makefile
                 base/muparser/Makefile
                 base/tesselator/Makefile
                 base/qhull/Makefile
                 base/Makefile
                 src/Makefile
                 src_extra/Makefile
                 doc_src/common_defs.inc
                 ])
AC_OUTPUT
