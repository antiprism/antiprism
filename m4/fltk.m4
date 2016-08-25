# Configure paths for FLTK
# Alvin Beach 26/01/2007
# Used sdl.m4 as example. Thanks SDL Team!
# stole most of the commands from:
# http://www.kdevelop.org/mediawiki/index.php/FLTK_hello_world_application_template

dnl AM_PATH_FLTK([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for FLTK, and define FLTK_CXXFLAGS, FLTK_LIBS, and FLTK_LDFLAGS
dnl
AC_DEFUN([AM_PATH_FLTK],
[dnl
dnl Get the cflags and libraries from the fltk-config script
dnl

dnl Prefix where FLTK is installed
AC_ARG_WITH(fltk-prefix,[  --with-fltk-prefix=PREFIX    Prefix where FLTK is installed (default: /usr)],
        fltk_prefix="$withval", fltk_prefix="/usr")

dnl Switch to pass to fltk-config
AC_ARG_WITH(fltk-ldopts,[  --with-fltk-ldopts=OPTS    Linker options to pass to fltk-config (e.g "--use-images --use-gl")],
        fltk_ldopts="$withval", fltk_ldopts="")

AC_REQUIRE([AC_CANONICAL_TARGET])

PATH="$fltk_prefix/bin:$prefix/bin:$prefix/usr/bin:$PATH"

AC_PATH_PROG(FLTK_CONFIG, fltk-config, no, $PATH)

dnl Assume fltk is not installed
FLTK_CXXFLAGS=""
FLTK_LIBS=""
FLTK_LDFLAGS=""

if test "$FLTK_CONFIG" != "no" ; then
        min_fltk_version=ifelse([$1], ,0.0.0,$1)

        AC_MSG_CHECKING(for FLTK - version >= $min_fltk_version)

        fltk_version="`$FLTK_CONFIG --version`"

        dnl Convert versions to an integer e.g. 1.1.7 becomes 117
        min_vers=`echo $min_fltk_version | sed 's/\.//g'`
        fltk_vers=`echo $fltk_version | sed 's/\.//g'`

        if test -n "$fltk_vers" -a "$fltk_vers" -ge $min_vers; then
                dnl All is good. Set the FLTK flags and such
                AC_MSG_RESULT(yes)
                ifelse([$2], , :, [$2])

                if test -n "$fltk_ldopts"; then
                        AC_MSG_NOTICE([FLTK linker options: $fltk_ldopts])
                fi

                FLTK_CXXFLAGS="`$FLTK_CONFIG --cxxflags`"
                dnl Antiprism: Always need --use-glut
                FLTK_LDFLAGS="`$FLTK_CONFIG $fltk_ldopts --use-glut --ldflags`"
                FLTK_LIBS="`$FLTK_CONFIG $fltk_ldopts --use-glut --ldflags`"
        else
                dnl A version of FLTK less-than $min_fltk_version is installed...fltk should be upgraded
                AC_MSG_RESULT(no)
                AC_MSG_ERROR([Installed version of FLTK is too old. Please upgrade to atleast $min_fltk_version])
                ifelse([$3], , :, [$3])
        fi
else
        dnl fltk-config could not be found. fltk must not be installed
        AC_MSG_ERROR([fltk-config could not be found. Either FLTK is not installed or FLTK's location needs to be specified with --with-fltk-prefix ])
        ifelse([$3], , :, [$3])
fi
])
