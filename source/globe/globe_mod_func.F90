#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_mod_func.F90
!> \brief Interface functions, used for short macros

!> \file globe_mod_func.F90
!> This file includes some interface function declarations to make
!! functions variable type independent. This allows the use of short
!! and compatible macros, without a lot of huge function calls in the
!! program code.\n
!! The first interface block pre-defines functions to allow the use of
!! optional parameters.\n
!! The following named interface blocks define the interface between the
!! variable type independent function call and the different variable
!! type dependent functions. All variable type dependent functions
!! are defined in the file globe_mod_func_inc.f90 and are included in
!! the \c contains part of this file. The include allows the use of
!! macro variables, whereby the functions only need a single definition
!! by replacing the different variable types by macro variables.

#ifndef __DOXYGEN

!     ******************************************************************
!     GLOBAL FUNCTION DEFINITION
!     ******************************************************************

      module globe_functions
      implicit none

!     * declaration of procedures with optional parameters

      interface
        subroutine globe_read_srvnc (kfile, f90, line, head, data, nbyte, dtype, iostat)
        implicit none
          integer,             INTENT(in)  :: kfile, line
          character(len=*),    INTENT(in)  :: f90
          integer,             INTENT(out) :: head(8)
          integer,             INTENT(in)  :: nbyte
          integer,   OPTIONAL, INTENT(in)  :: dtype
          integer*1, OPTIONAL, INTENT(out) :: data(nbyte)
          integer,   OPTIONAL, INTENT(out) :: iostat
        end subroutine globe_read_srvnc

        subroutine error_read (ioerr, sfile, f90, line)
        implicit none
          integer,                    INTENT(in) :: ioerr
          character(len=*),           INTENT(in) :: sfile
          character(len=*), OPTIONAL, INTENT(in) :: f90
          integer,          OPTIONAL, INTENT(in) :: line
        end subroutine error_read
      end interface

!     * function to convert any number to string (left-aligned)

      interface globe_string
        module procedure globe_i1_str, &
     &                   globe_i2_str, &
     &                   globe_i4_str, &
     &                   globe_i8_str, &
     &                   globe_r4_str, &
     &                   globe_r8_str, &
     &                   globe_c8_str, &
     &                   globe_c16_str
      end interface

!     * interface of procedure globe_read_srvnc to allow a one line call (macro)
!     * (This long definition is needed to prevent the declaration of
!        array size and variable type at the procedure call (macro).)
!     * Because of the multiple procedure definitions, this construct is used
!       only as interface (layer) between a short (one line) macro and the
!       real long procedure call. So, changes at the real procedure do not
!       need to be performed multiple times.

      interface globe_read_srvnc_
        module procedure globe_read_srvnc_si1, &
     &                   globe_read_srvnc_si2, &
     &                   globe_read_srvnc_si4, &
     &                   globe_read_srvnc_si8, &
     &                   globe_read_srvnc_sr4, &
     &                   globe_read_srvnc_sr8, &
     &                   globe_read_srvnc_sc8, &
     &                   globe_read_srvnc_sc16,&
     &                   globe_read_srvnc_ai1, &
     &                   globe_read_srvnc_ai2, &
     &                   globe_read_srvnc_ai4, &
     &                   globe_read_srvnc_ai8, &
     &                   globe_read_srvnc_ar4, &
     &                   globe_read_srvnc_ar8, &
     &                   globe_read_srvnc_ac8, &
     &                   globe_read_srvnc_ac16
      end interface

!     * interface of procedure globe_mpsc to allow a one line call (macro)
!     * (This long definition is needed to prevent the declaration of
!        array size and variable type at the procedure call (macro).)
!     * Because of the multiple procedure definitions, this construct is used
!       only as interface (layer) between a short (one line) macro and the
!       real long procedure call. So, changes at the real procedure do not
!       need to be performed multiple times.

      interface globe_mpsc_
        module procedure globe_mpsc_1i1, &
     &                   globe_mpsc_1i2, &
     &                   globe_mpsc_1i4, &
     &                   globe_mpsc_1i8, &
     &                   globe_mpsc_1r4, &
     &                   globe_mpsc_1r8, &
     &                   globe_mpsc_1c8, &
     &                   globe_mpsc_1c16,&
     &                   globe_mpsc_2i1, &
     &                   globe_mpsc_2i2, &
     &                   globe_mpsc_2i4, &
     &                   globe_mpsc_2i8, &
     &                   globe_mpsc_2r4, &
     &                   globe_mpsc_2r8, &
     &                   globe_mpsc_2c8, &
     &                   globe_mpsc_2c16
      end interface

!     * interface of procedure globe_mpga to allow a one line call (macro)
!     * (This long definition is needed to prevent the declaration of
!        array size and variable type at the procedure call (macro).)
!     * Because of the multiple procedure definitions, this construct is used
!       only as interface (layer) between a short (one line) macro and the
!       real long procedure call. So, changes at the real procedure do not
!       need to be performed multiple times.

      interface globe_mpga_
        module procedure globe_mpga_1i1, &
     &                   globe_mpga_1i2, &
     &                   globe_mpga_1i4, &
     &                   globe_mpga_1i8, &
     &                   globe_mpga_1r4, &
     &                   globe_mpga_1r8, &
     &                   globe_mpga_1c8, &
     &                   globe_mpga_1c16,&
     &                   globe_mpga_2i1, &
     &                   globe_mpga_2i2, &
     &                   globe_mpga_2i4, &
     &                   globe_mpga_2i8, &
     &                   globe_mpga_2r4, &
     &                   globe_mpga_2r8, &
     &                   globe_mpga_2c8, &
     &                   globe_mpga_2c16
      end interface

      contains

#define __TYPE__ integer*1
#define __FORMAT__ '(I0)'
#define __NAME1__ globe_i1_str
#define __NAME2__ globe_read_srvnc_si1
#define __NAME3__ globe_read_srvnc_ai1
#define __NAME4__ globe_mpsc_1i1
#define __NAME5__ globe_mpsc_2i1
#define __NAME6__ globe_mpga_1i1
#define __NAME7__ globe_mpga_2i1
#define __INCLUDE_IOSTAT__
#include "./globe_mod_func_inc.f90"
#undef __INCLUDE_IOSTAT__

#define __TYPE__ integer*2
#define __FORMAT__ '(I0)'
#define __NAME1__ globe_i2_str
#define __NAME2__ globe_read_srvnc_si2
#define __NAME3__ globe_read_srvnc_ai2
#define __NAME4__ globe_mpsc_1i2
#define __NAME5__ globe_mpsc_2i2
#define __NAME6__ globe_mpga_1i2
#define __NAME7__ globe_mpga_2i2
#include "./globe_mod_func_inc.f90"

#define __TYPE__ integer*4
#define __FORMAT__ '(I0)'
#define __NAME1__ globe_i4_str
#define __NAME2__ globe_read_srvnc_si4
#define __NAME3__ globe_read_srvnc_ai4
#define __NAME4__ globe_mpsc_1i4
#define __NAME5__ globe_mpsc_2i4
#define __NAME6__ globe_mpga_1i4
#define __NAME7__ globe_mpga_2i4
#include "./globe_mod_func_inc.f90"

#define __TYPE__ integer*8
#define __FORMAT__ '(I0)'
#define __NAME1__ globe_i8_str
#define __NAME2__ globe_read_srvnc_si8
#define __NAME3__ globe_read_srvnc_ai8
#define __NAME4__ globe_mpsc_1i8
#define __NAME5__ globe_mpsc_2i8
#define __NAME6__ globe_mpga_1i8
#define __NAME7__ globe_mpga_2i8
#include "./globe_mod_func_inc.f90"

#define __TYPE__ real*4
#define __FORMAT__ '(g0)'
#define __NAME1__ globe_r4_str
#define __NAME2__ globe_read_srvnc_sr4
#define __NAME3__ globe_read_srvnc_ar4
#define __NAME4__ globe_mpsc_1r4
#define __NAME5__ globe_mpsc_2r4
#define __NAME6__ globe_mpga_1r4
#define __NAME7__ globe_mpga_2r4
#include "./globe_mod_func_inc.f90"

#define __TYPE__ real*8
#define __FORMAT__ '(g0)'
#define __NAME1__ globe_r8_str
#define __NAME2__ globe_read_srvnc_sr8
#define __NAME3__ globe_read_srvnc_ar8
#define __NAME4__ globe_mpsc_1r8
#define __NAME5__ globe_mpsc_2r8
#define __NAME6__ globe_mpga_1r8
#define __NAME7__ globe_mpga_2r8
#include "./globe_mod_func_inc.f90"

#define __TYPE__ complex*8
#define __FORMAT__ '(g0,g0)'
#define __NAME1__ globe_c8_str
#define __NAME2__ globe_read_srvnc_sc8
#define __NAME3__ globe_read_srvnc_ac8
#define __NAME4__ globe_mpsc_1c8
#define __NAME5__ globe_mpsc_2c8
#define __NAME6__ globe_mpga_1c8
#define __NAME7__ globe_mpga_2c8
#include "./globe_mod_func_inc.f90"

#define __TYPE__ complex*16
#define __FORMAT__ '(g0,g0)'
#define __NAME1__ globe_c16_str
#define __NAME2__ globe_read_srvnc_sc16
#define __NAME3__ globe_read_srvnc_ac16
#define __NAME4__ globe_mpsc_1c16
#define __NAME5__ globe_mpsc_2c16
#define __NAME6__ globe_mpga_1c16
#define __NAME7__ globe_mpga_2c16
#include "./globe_mod_func_inc.f90"

      end module globe_functions
#endif /* __DOXYGEN */
