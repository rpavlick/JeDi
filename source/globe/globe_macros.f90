#ifndef __IGNORE
! To use this macros, add the commands:
! #define __ACTIVATE
! #include "../globe/globe_macros.f90"
! to the begin of the fortran file and rename its extension to .F90

!> \file globe_macros.f90
!> \brief Global macro definitions

!> \file globe_macros.f90
!> This file includes the global macro definitions used in all models.
!! This allows the use of short subroutine calls in the code, where the
!! macros add the calculation of the required field sizes and line numbers.

#endif /* __IGNORE */
#ifdef __ACTIVATE
#ifndef __IGNORE

! **********************************************************************

#if defined (__DEBUG) && ! defined (__NONAN) && ! defined (__IFORT) && ! defined (__XLF90)
#define __allocate(VAR,SIZE) if (.not. allocated(VAR)) allocate(VAR SIZE) ; VAR = VAR * 0.0 / 0.0
#else
#define __allocate(VAR,SIZE) if (.not. allocated(VAR)) allocate(VAR SIZE)
#endif

#define __deallocate(VAR) if (allocated(VAR)) deallocate(VAR)

#define __globe_mpbc(VAR) call globe_mpbc(VAR,product(shape(VAR))*kind(VAR))

#define __globe_mpsc_from_to(FROM,TO) call globe_mpsc_(__FILE__,__LINE__,FROM,TO)

#define __globe_mpga_to_from(TO,FROM) call globe_mpga_(__FILE__,__LINE__,TO,FROM)

#define __globe_writeoutput(FILE,VAR,CODE) call globe_mpwritearray(FILE,VAR,size(VAR,1),size(VAR),CODE)

#define __globe_read_srvnc(FILE,HEAD,DATA) call globe_read_srvnc_(FILE,__FILE__,__LINE__,HEAD,DATA)

#define __error_read(IOERR,FILE) call error_read(IOERR,FILE,__FILE__,__LINE__)

#define __date_string() (((kglobe_year + kglobe_firstyear - 1) * 100) + kglobe_month) * 100 + kglobe_day

#define __diag(FILE,MSG) call globe_write_diag(FILE,MSG," ")
#define __diag_num(FILE,MSG,NUM) call globe_write_diag(FILE,MSG,TRIM(globe_string(NUM)))
! For formated output, use:
!       character(len=135) :: msg
!       write(msg,'(A,I3,F)') 'output: ', var1, var2
!       __diag(kfile,trim(msg))

! **********************************************************************
!                           Documentation macros
! **********************************************************************

#ifdef __DOXYGEN
#define __vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv() !> \}
#define __SEC() !> \name
#define __oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo() !> \{
#else
#define __vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv() !     ******************************************************************
#define __SEC() !     *
#define __oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo() !     ******************************************************************
#endif

#endif /* __IGNORE */
#endif /* __ACTIVATE */
