# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_BLOSC([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of BLOSC library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   former only looks for serial BLOSC installations via h5cc. The latter
#   only looks for parallel BLOSC installations via h5pcc. If the optional
#   argument is omitted, serial installations will be preferred over
#   parallel ones.
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the BLOSC library.
#     yes  - do check for BLOSC library in standard locations.
#     path - complete path to the BLOSC helper script h5cc or h5pcc.
#
#   If BLOSC is successfully found, this macro calls
#
#     AC_SUBST(BLOSC_VERSION)
#     AC_SUBST(BLOSC_CC)
#     AC_SUBST(BLOSC_CFLAGS)
#     AC_SUBST(BLOSC_CPPFLAGS)
#     AC_SUBST(BLOSC_LDFLAGS)
#     AC_SUBST(BLOSC_LIBS)
#     AC_SUBST(BLOSC_FC)
#     AC_SUBST(BLOSC_FFLAGS)
#     AC_SUBST(BLOSC_FLIBS)
#     AC_DEFINE(HAVE_BLOSC)
#
#   and sets with_hdf5="yes".  Additionally, the macro sets
#   with_hdf5_fortran="yes" if a matching Fortran wrapper script is found.
#   Note that Autconf's Fortran support is not used to perform this check.
#   H5CC and H5FC will contain the appropriate serial or parallel BLOSC
#   wrapper script locations.
#
#   If BLOSC is disabled or not found, this macros sets with_hdf5="no" and
#   with_hdf5_fortran="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. BLOSC_{C,CPP,LD}FLAGS may be used when building with C or C++.
#   BLOSC_F{FLAGS,LIBS} should be used when building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for BLOSC support
#        AX_LIB_BLOSC()
#
#     2) dnl Check for serial BLOSC support
#        AX_LIB_BLOSC([serial])
#
#     3) dnl Check for parallel BLOSC support
#        AX_LIB_BLOSC([parallel])
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "BLOSC support:  $with_hdf5"
#
#   You could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that BLOSC uses:
#
#     AX_LIB_BLOSC([parallel])
#     if test "$with_hdf5" = "yes"; then
#             CC="$BLOSC_CC"
#     else
#             AC_MSG_ERROR([Unable to find BLOSC, we need parallel BLOSC.])
#     fi
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 12

# This code was adapted by Holger Hoefling on May 2nd 2016 to also 
# substitute variables BLOSC_MAJOR_VERSION, BLOSC_MINOR_VERSION 
# as well as BLOSC_REVISION_VERSION so that the 3 parts of the
# version number are available separately for testing 
# Also changed to check for hdf5_hl header library existence


AC_DEFUN([AX_LIB_BLOSC], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    : # Recognized value
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_BLOSC within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(
    [--with-hdf5=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [location of h5cc for serial BLOSC configuration],
            [parallel], [location of h5pcc for parallel BLOSC configuration],
            [location of h5cc or h5pcc for BLOSC configuration])
  ),
  [if test "$withval" = "no"; then
     with_hdf5="no"
   elif test "$withval" = "yes"; then
     with_hdf5="yes"
   else
     with_hdf5="yes"
     H5CC="$withval"
   fi],
   [with_hdf5="yes"]
)

dnl Set defaults to blank
BLOSC_CC=""
BLOSC_VERSION=""
BLOSC_CFLAGS=""
BLOSC_CPPFLAGS=""
BLOSC_LDFLAGS=""
BLOSC_LIBS=""
BLOSC_FC=""
BLOSC_FFLAGS=""
BLOSC_FLIBS=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_hdf5" = "yes"; then
    if test -z "$H5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [H5CC],
            m4_case(m4_normalize([$1]),
                [serial],   [h5cc],
                [parallel], [h5pcc],
                [h5cc h5pcc]),
            [])
    else
        AC_MSG_CHECKING([Using provided BLOSC C wrapper])
        AC_MSG_RESULT([$H5CC])
    fi
    AC_MSG_CHECKING([for BLOSC libraries])
    if test ! -f "$H5CC" || test ! -x "$H5CC"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN(m4_case(m4_normalize([$1]),
            [serial],  [
Unable to locate serial BLOSC compilation helper script 'h5cc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc.
BLOSC support is being disabled (equivalent to --with-hdf5=no).
],            [parallel],[
Unable to locate parallel BLOSC compilation helper script 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5pcc.
BLOSC support is being disabled (equivalent to --with-hdf5=no).
],            [
Unable to locate BLOSC compilation helper scripts 'h5cc' or 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc or h5pcc.
BLOSC support is being disabled (equivalent to --with-hdf5=no).
]))
        with_hdf5="no"
        with_hdf5_fortran="no"
    else
        dnl Get the h5cc output
        BLOSC_SHOW=$(eval $H5CC -show)

        dnl Get the actual compiler used
        BLOSC_CC=$(eval $H5CC -show | $AWK '{print $[]1}')
        if test "$BLOSC_CC" = "ccache"; then
            BLOSC_CC=$(eval $H5CC -show | $AWK '{print $[]2}')
        fi

        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "BLOSC Version: X.Y.Z"
        BLOSC_VERSION=$(eval $H5CC -showconfig | $GREP 'BLOSC Version:' \
            | $AWK '{print $[]3}')
        BLOSC_MAJOR_VERSION=$(echo $BLOSC_VERSION | $AWK -F \. '{print $[]1}')
        BLOSC_MINOR_VERSION=$(echo $BLOSC_VERSION | $AWK -F \. {'print $[]2'})
        BLOSC_REVISION_VERSION=$(echo $BLOSC_VERSION | $AWK -F \. {'print $[]3'})

        dnl A ideal situation would be where everything we needed was
        dnl in the AM_* variables. However most systems are not like this
        dnl and seem to have the values in the non-AM variables.
        dnl
        dnl We try the following to find the flags:
        dnl (1) Look for "NAME:" tags
        dnl (2) Look for "H5_NAME:" tags
        dnl (3) Look for "AM_NAME:" tags
        dnl
        BLOSC_tmp_flags=$(eval $H5CC -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

        dnl Find the installation directory and append include/
        BLOSC_tmp_inst=$(eval $H5CC -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

        dnl Add this to the CPPFLAGS
        BLOSC_CPPFLAGS="-I${BLOSC_tmp_inst}/include"

        dnl Now sort the flags out based upon their prefixes
        for arg in $BLOSC_SHOW $BLOSC_tmp_flags ; do
          case "$arg" in
            -I*) echo $BLOSC_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || BLOSC_CPPFLAGS="$arg $BLOSC_CPPFLAGS"
              ;;
            -L*) echo $BLOSC_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || BLOSC_LDFLAGS="$arg $BLOSC_LDFLAGS"
              ;;
            -l*) echo $BLOSC_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || BLOSC_LIBS="$arg $BLOSC_LIBS"
              ;;
          esac
        done

        BLOSC_LIBS="$BLOSC_LIBS -lhdf5"
        AC_MSG_RESULT([yes (version $[BLOSC_VERSION])])

        dnl See if we can compile
        AC_LANG_PUSH([C])
        ax_lib_hdf5_save_CC=$CC
        ax_lib_hdf5_save_CPPFLAGS=$CPPFLAGS
        ax_lib_hdf5_save_LIBS=$LIBS
        ax_lib_hdf5_save_LDFLAGS=$LDFLAGS
        CC=$BLOSC_CC
        CPPFLAGS=$BLOSC_CPPFLAGS
        LIBS=$BLOSC_LIBS
        LDFLAGS=$BLOSC_LDFLAGS
        AC_CHECK_HEADER([hdf5.h], [ac_cv_hadf5_h=yes], [ac_cv_hadf5_h=no])
        AC_CHECK_LIB([hdf5], [H5Fcreate], [ac_cv_libhdf5=yes],
                     [ac_cv_libhdf5=no])
        if test "$ac_cv_hadf5_h" = "no" && test "$ac_cv_libhdf5" = "no" ; then
          AC_MSG_WARN([Unable to compile BLOSC test program])
        fi
        AC_CHECK_HEADER([hdf5_hl.h], [ac_cv_hadf5_hl_h=yes], [ac_cv_hadf5_hl_h=no], [#include <hdf5.h>])
        AC_CHECK_LIB([hdf5_hl], [H5LTpath_valid], [ac_cv_libhdf5_hl=yes],
                     [ac_cv_libhdf5_hl=no])
        if test "$ac_cv_hadf5_hl_h" = "no" && test "$ac_cv_libhdf5_hl" = "no"; then
          AC_MSG_WARN([Unable to compile BLOSC_HL test program])
        fi
        dnl Look for BLOSC's high level library
        AC_HAVE_LIBRARY([hdf5_hl], [BLOSC_LIBS="$BLOSC_LIBS -lhdf5_hl"], [], [])

        CC=$ax_lib_hdf5_save_CC
        CPPFLAGS=$ax_lib_hdf5_save_CPPFLAGS
        LIBS=$ax_lib_hdf5_save_LIBS
        LDFLAGS=$ax_lib_hdf5_save_LDFLAGS
        AC_LANG_POP([C])

        AC_MSG_CHECKING([for matching BLOSC Fortran wrapper])
        dnl Presume BLOSC Fortran wrapper is just a name variant from H5CC
        H5FC=$(eval echo -n $H5CC | $SED -n 's/cc$/fc/p')
        if test -x "$H5FC"; then
            AC_MSG_RESULT([$H5FC])
            with_hdf5_fortran="yes"
            AC_SUBST([H5FC])

            dnl Again, pry any remaining -Idir/-Ldir from compiler wrapper
            for arg in `$H5FC -show`
            do
              case "$arg" in #(
                -I*) echo $BLOSC_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || BLOSC_FFLAGS="$arg $BLOSC_FFLAGS"
                  ;;#(
                -L*) echo $BLOSC_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || BLOSC_FFLAGS="$arg $BLOSC_FFLAGS"
                     dnl BLOSC installs .mod files in with libraries,
                     dnl but some compilers need to find them with -I
                     echo $BLOSC_FFLAGS | $GREP -e "-I${arg#-L}" >/dev/null \
                      || BLOSC_FFLAGS="-I${arg#-L} $BLOSC_FFLAGS"
                  ;;
              esac
            done

            dnl Make Fortran link line by inserting Fortran libraries
            for arg in $BLOSC_LIBS
            do
              case "$arg" in #(
                -lhdf5_hl) BLOSC_FLIBS="$BLOSC_FLIBS -lhdf5hl_fortran $arg"
                  ;; #(
                -lhdf5)    BLOSC_FLIBS="$BLOSC_FLIBS -lhdf5_fortran $arg"
                  ;; #(
                *) BLOSC_FLIBS="$BLOSC_FLIBS $arg"
                  ;;
              esac
            done
        else
            AC_MSG_RESULT([no])
            with_hdf5_fortran="no"
        fi

	AC_SUBST([BLOSC_VERSION])
	AC_SUBST([BLOSC_MAJOR_VERSION])
	AC_SUBST([BLOSC_MINOR_VERSION])
	AC_SUBST([BLOSC_REVISION_VERSION])
	AC_SUBST([BLOSC_CC])
	AC_SUBST([BLOSC_CFLAGS])
	AC_SUBST([BLOSC_CPPFLAGS])
	AC_SUBST([BLOSC_LDFLAGS])
	AC_SUBST([BLOSC_LIBS])
	AC_SUBST([BLOSC_FC])
	AC_SUBST([BLOSC_FFLAGS])
	AC_SUBST([BLOSC_FLIBS])
	AC_DEFINE([HAVE_BLOSC], [1], [Defined if you have BLOSC support])
    fi
fi
])
