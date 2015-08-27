#define __ACTIVATE
#include "../globe/globe_macros.f90"

!     * GLOBAL INTERFACE SUB-FUNCTION DEFINITIONS
!       Its used in globe_mod_func.F90 as command line:
!       #include "./globe_mod_func_inc.f90"

!> \file globe_mod_func_inc.f90
!! \brief Interface subroutine definitions

!> \file globe_mod_func_inc.f90
!! This file includes the definitions of the subroutines which are
!! included to the file globe_mod_func.F90 in the section \c contains.
!! To avoid multiple definitions of similar subroutines, which only
!! differ by the variable type, macro variables are used.
!! The macro variables used in the file globe_mod_func_inc.f90 are
!! defined in globe_mod_func.F90 in the section \c contains.

!     ******************************************************************
!     * sub-function to convert any number to string (left-aligned)
!     ******************************************************************

!> \brief Convert any number to left-aligned string

!> \detail This function converts any number to a string (left-aligned).\n

!> \param num any numerical scalar variable

!> \param __NAME1__ is replaced by the following function names:\n
!! globe_i1_str, globe_i2_str, globe_i4_str, globe_i8_str,
!! globe_r4_str, globe_r8_str, globe_c8_str, and globe_c16_str

!> \param __TYPE__ is replaced by the appropriate variable type

!> \param __FORMAT__ is replaced by the appropriate format string of the
!! variable type

!> \return left aligned string

!     ******************************************************************

      function __NAME1__(num) result(str3)
      implicit none
        __TYPE__ , INTENT(in) :: num

        character*15 :: str1
        character*80 :: str2, str3

        write(str1, __FORMAT__ ) num
        write(str2,'(A)') adjustl(str1)
        str3 = str2
      return
      end function __NAME1__

!     ******************************************************************
!     * sub-procedures of globe_read_srvnc to allow a one line call (macro)
!       (one for scalar and one for vector data, special case for iostat)
!     ******************************************************************

!> \brief Interface procedure to globe_read_srvnc

!> \detail See description at subroutine \c __NAME3__

!     ******************************************************************

#ifdef __INCLUDE_IOSTAT__
      subroutine __NAME2__ (kfile, f90, line, head, data, iostat)
#else /* __INCLUDE_IOSTAT__ */
      subroutine __NAME2__ (kfile, f90, line, head, data)
#endif /* __INCLUDE_IOSTAT__ */
      implicit none
        integer,             intent(in)  :: kfile, line
        character(len=*),    intent(in)  :: f90
        integer,             intent(out) :: head(8)
#ifdef __INCLUDE_IOSTAT__
        __TYPE__ , optional, intent(out) :: data
        integer,   optional, intent(out) :: iostat
#else /* __INCLUDE_IOSTAT__ */
        __TYPE__ ,           intent(out) :: data
#endif /* __INCLUDE_IOSTAT__ */
        integer   :: nbyte, dtype
        integer*1 :: bdata(KIND(data))
#ifdef __INCLUDE_IOSTAT__
        if ( present(iostat) ) then
          call globe_read_srvnc(kfile, f90, line, head, nbyte=0, iostat=iostat)
          return
        endif
#endif /* __INCLUDE_IOSTAT__ */
        nbyte = KIND(data)
        dtype = RANGE(data) * KIND(data)
        call globe_read_srvnc(kfile, f90, line, head, bdata, nbyte, dtype)
        data = transfer(bdata, data)
      return
      end subroutine __NAME2__

!     ******************************************************************

!> \brief Interface procedure to globe_read_srvnc

!> \detail Interface procedures of globe_read_srvnc()
!! (globe_functions.F90), which make the call variable type independent,
!! set the field sizes, and convert the \c data to the right type.\n
!> The parameter \c iostat is only used with the variable type
!! \c integer*1 to allow the return of the status message. For all
!! other variable types the second form of \c __NAME2__ is used.
!! Both routines are similar, but \c __NAME2__ is used for scalar
!! \c data, and \c __NAME3__ for \c data arrays.

!> \param kfile  File ID of the input file.
!> \param f90    Name of the Fortran file which includes the call to
!! this routine, used in case an error message is needed.
!> \param line   Line number in the Fortran file with the call to
!! this routine, used in case an error message is needed.
!> \param head   The header of the data field in SRV-format.
!> \param data   The input data converted to the requested type.
!> \param iostat The status of the read operation of the data.
!! header. No data are returned in case \c iostat is present.

!> \param __NAME2__ is replaced by the following function names:\n
!! globe_read_srvnc_si1, globe_read_srvnc_si2, globe_read_srvnc_si4,
!! globe_read_srvnc_si8, globe_read_srvnc_sr4, globe_read_srvnc_sr8,
!! globe_read_srvnc_sc8, globe_read_srvnc_sc16

!> \param __NAME3__ is replaced by the following function names:\n
!! globe_read_srvnc_ai1, globe_read_srvnc_ai2, globe_read_srvnc_ai4,
!! globe_read_srvnc_ai8, globe_read_srvnc_ar4, globe_read_srvnc_ar8,
!! globe_read_srvnc_ac8, globe_read_srvnc_ac16

!> \param __INCLUDE_IOSTAT__ Flag to trigger the use of the first form
!! of the call to subroutine \c __NAME2__ which include the variable
!! \c iostat. This is only used with the \c __TYPE__ \c integer*1.

!> \param __TYPE__ is replaced by the appropriate variable type.

!     ******************************************************************

      subroutine __NAME3__ (kfile, f90, line, head, data)
      implicit none
        integer,          intent(in)  :: kfile, line
        character(len=*), intent(in)  :: f90
        integer,          intent(out) :: head(8)
        __TYPE__ ,        intent(out) :: data(:)

        integer :: nbyte, dtype
        integer*1, allocatable :: bdata(:)

        nbyte = SIZE(data) * KIND(data)
        dtype = RANGE(data) * KIND(data)
        allocate(bdata(nbyte))
        call globe_read_srvnc(kfile, f90, line, head, bdata, nbyte, dtype)
        data = transfer(bdata, data)
        deallocate(bdata)
      return
      end subroutine __NAME3__

!     ******************************************************************
!     * sub-procedures of globe_mpsc to allow a one line call (macro)
!       (one for vector and one for array data)
!     ******************************************************************

!> \brief Interface procedure of globe_mpsc

!> \detail See description at subroutine \c __NAME5__

!     ******************************************************************

      subroutine __NAME4__ (f90, line, pf, pp)
      implicit none
        character(len=*), INTENT(in)  :: f90
        integer,          INTENT(in)  :: line
        __TYPE__ ,        INTENT(in)  :: pf(:)
        __TYPE__ ,        INTENT(out) :: pp(:)

        integer :: nb1pf, nb1pp

        nb1pf = SIZE(pf) * KIND(pf)
        nb1pp = SIZE(pp) * KIND(pp)
        call globe_mpsc(f90, line, pf, nb1pf, 1, pp, nb1pp, 1)
      return
      end subroutine __NAME4__

!     ******************************************************************

!> \brief Interface procedure of globe_mpsc

!> \detail Interface procedures of globe_mpsc() (globe_mpimod.F90),
!! which make the call variable type independent, and set the field
!! sizes.\n
!! Both routines are similar, but \c __NAME4__ is used for one
!! dimensional arrays, and \c __NAME5__ for two dimensional arrays.

!> \param f90  Name of the Fortran file which includes the call to
!! this routine, used in case an error message is needed.
!> \param line Line number in the Fortran file with the call to
!! this routine, used in case an error message is needed.
!> \param pf   Global field to scatter to all parallel processes.
!> \param pp   Local process field part of \c pf.

!> \param __NAME4__ is replaced by the following function names:\n
!! globe_mpsc_1i1, globe_mpsc_1i2, globe_mpsc_1i4, globe_mpsc_1i8,
!! globe_mpsc_1r4, globe_mpsc_1r8, globe_mpsc_1c8, globe_mpsc_1c16

!> \param __NAME5__ is replaced by the following function names:\n
!! globe_mpsc_2i1, globe_mpsc_2i2, globe_mpsc_2i4, globe_mpsc_2i8,
!! globe_mpsc_2r4, globe_mpsc_2r8, globe_mpsc_2c8, globe_mpsc_2c16

!> \param __TYPE__ is replaced by the appropriate variable type.

!     ******************************************************************

      subroutine __NAME5__ (f90, line, pf, pp)
      implicit none
        character(len=*), INTENT(in)  :: f90
        integer,          INTENT(in)  :: line
        __TYPE__ ,        INTENT(in)  :: pf(:,:)
        __TYPE__ ,        INTENT(out) :: pp(:,:)

        integer :: nb1pf, d2pf, nb1pp, d2pp

        nb1pf = SIZE(pf,1) * KIND(pf)
        d2pf  = SIZE(pf,2)
        nb1pp = SIZE(pp,1) * KIND(pp)
        d2pp  = SIZE(pp,2)
        call globe_mpsc(f90, line, pf, nb1pf, d2pf, pp, nb1pp, d2pp)
      return
      end subroutine __NAME5__

!     ******************************************************************
!     * sub-procedures of globe_mpga to allow a one line call (macro)
!       (one for vector and one for array data)
!     ******************************************************************

!> \brief Interface procedures of globe_mpga

!> \detail See description at subroutine \c __NAME7__

!     ******************************************************************

      subroutine __NAME6__ (f90, line, pf, pp)
      implicit none
        character(len=*), INTENT(in)  :: f90
        integer,          INTENT(in)  :: line
        __TYPE__ ,        INTENT(out) :: pf(:)
        __TYPE__ ,        INTENT(in)  :: pp(:)

        integer :: nb1pf, nb1pp

        nb1pf = SIZE(pf) * KIND(pf)
        nb1pp = SIZE(pp) * KIND(pp)
        call globe_mpga(f90, line, pf, nb1pf, 1, pp, nb1pp, 1)
      return
      end subroutine __NAME6__

!     ******************************************************************

!> \brief Interface procedures of globe_mpga

!> \detail Interface procedures of globe_mpga() (globe_mpimod.F90),
!! which make the call variable type independent, and set the field
!! sizes.\n
!! Both routines are similar, but \c __NAME6__ is used for one
!! dimensional arrays, and \c __NAME7__ for two dimensional arrays.

!> \param f90  Name of the Fortran file which includes the call to
!! this routine, used in case an error message is needed.
!> \param line Line number in the Fortran file with the call to
!! this routine, used in case an error message is needed.
!> \param pf   Global field gathered from all parallel processes.
!> \param pp   Local process field to become part of \c pf.

!> \param __NAME6__ is replaced by the following function names:\n
!! globe_mpga_1i1, globe_mpga_1i2, globe_mpga_1i4, globe_mpga_1i8,
!! globe_mpga_1r4, globe_mpga_1r8, globe_mpga_1c8, globe_mpga_1c16

!> \param __NAME7__ is replaced by the following function names:\n
!! globe_mpga_2i1, globe_mpga_2i2, globe_mpga_2i4, globe_mpga_2i8,
!! globe_mpga_2r4, globe_mpga_2r8, globe_mpga_2c8, globe_mpga_2c16

!> \param __TYPE__ is replaced by the appropriate variable type.

!     ******************************************************************

      subroutine __NAME7__ (f90, line, pf, pp)
      implicit none
        character(len=*), INTENT(in)  :: f90
        integer,          INTENT(in)  :: line
        __TYPE__ ,        INTENT(out) :: pf(:,:)
        __TYPE__ ,        INTENT(in)  :: pp(:,:)

        integer :: nb1pf, d2pf, nb1pp, d2pp

        nb1pf = SIZE(pf,1) * KIND(pf)
        d2pf  = SIZE(pf,2)
        nb1pp = SIZE(pp,1) * KIND(pp)
        d2pp  = SIZE(pp,2)
        call globe_mpga(f90, line, pf, nb1pf, d2pf, pp, nb1pp, d2pp)
      return
      end subroutine __NAME7__

#undef __TYPE__
#undef __FORMAT__
#undef __NAME1__
#undef __NAME2__
#undef __NAME3__
#undef __NAME4__
#undef __NAME5__
#undef __NAME6__
#undef __NAME7__
