#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(antiprism, 0.20.pre14, adrian@antiprism.com)
AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([base/antiprism.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

AC_LANG([C++])

LT_INIT([win32_dll])

# Check whether Python module should be built
#AC_ARG_ENABLE( python,
#   [  --disable-python        exclude mod_python and support for metric modules written in python],
#   [ if test x"$enableval" != xyes; then enable_python="no"; fi ], [ enable_python="yes" ] )

#AM_PATH_PYTHON(2.3, [], [enable_python="no"])
      
AM_CONDITIONAL(BUILD_PY_MOD, false)
AM_CONDITIONAL(USE_SWIG, false)

#enable_python="disabled" #disable python for this release
if test x"$enable_python" = x"yes" ;
then
   AM_CONDITIONAL(BUILD_PY_MOD, true)
   echo "will build Python Module";
# Check whether Python module wrapper will be regenerated with swig
   AC_PROG_SWIG(1.3.21)
   SWIG_ENABLE_CXX
   SWIG_PYTHON
   if test "$SWIG_LIB" ; then
      AM_CONDITIONAL(USE_SWIG, true)
      echo "SWIG is available to rebuild the Python module interface"
   else
      echo
      echo "*********************************************"
      echo "   SWIG is not available. It is not needed to build and "
      echo "   install this package. However, if the library source is "
      echo "   modified then the changes will not be included in the "
      echo "   Python module."
      echo "*********************************************"
      echo
   fi
#else
#   echo
#   echo "*********************************************"
#   echo "   the Python module will not be built"
#   echo "   either configure was run with --disable-python"
#   echo "   or a Python development package is not installed (or found)"
#   echo "*********************************************"
#   echo
fi

# Add -fno-strict-aliasing to CFLAGS if accepted
ALIASING_FLAG=
if $CC -fno-strict-aliasing -S -o /dev/null -xc /dev/null >/dev/null 2>&1;
then
   ALIASING_FLAG="-fno-strict-aliasing"
   CPPFLAGS="$ALIASING_FLAG $CPPFLAGS"
fi
AC_SUBST(ALIASING_FLAG)


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
         [glut_type=$FREEGLUT])]
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
            AC_CHECK_HEADER([GL/openglut.h],
               [AC_CHECK_LIB([openglut], [main],
                  [GLUT_LIBS="-lopenglut"; glut_type=$OPENGLUT]) ])
         else
            if test $glut_type = $FREEGLUT ; then
               glut_type=$NO_GLUT
               AC_CHECK_HEADER([GL/freeglut.h],
                  [AC_CHECK_LIB([freeglut], [main],
                     [GLUT_LIBS="-lfreeglut"; glut_type=$FREEGLUT]) ])
            fi
         fi
      fi
   fi
   if FOUND_GL=1 ; then
         # Do not use X if not required; if no_x is yes then leave it as yes
         break
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

# Python binding


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
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset pow sqrt strchr strcspn strpbrk strrchr strspn strstr strtol])

AC_CONFIG_FILES([Makefile
                 aview/Makefile
                 base/tesselator/Makefile
                 base/Makefile
                 src/Makefile
                 src_extra/Makefile
                 ])
                 #python/Makefile
                 #python/distutils/setup.py
AC_OUTPUT
