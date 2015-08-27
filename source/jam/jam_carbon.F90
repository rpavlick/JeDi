#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam_carbon.F90

!     ==================================================================
!
!     J A M  ---   C A R B O N
!
!     ==================================================================

!     ******************************************************************
!     JAM_STEP_CARBON_Y
!     ******************************************************************

!> \brief Read of CO<sub>2</sub> data from text file

!> \detail This routine reads some CO<sub>2</sub> data from a text file.

!> \param koutyear Actual year to read the CO<sub>2</sub> data

!     ******************************************************************

      subroutine jam_step_carbon_y (koutyear)
      use jam_mod
      implicit none

      integer, INTENT(in) :: koutyear

      if (kmypid == kroot) then
        __diag_num(kfile_diag, 'Reading CO2 concentration for year ', koutyear)
        call jam_read_co2_list(xCO2g_a_init, koutyear)
      endif
      __globe_mpbc(xCO2g_a_init)

      return
      end subroutine jam_step_carbon_y
