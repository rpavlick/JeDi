#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_mod_stat.F90
!> \brief global status variable definitions

!> \file globe_mod_stat.F90
!> This file contains the \c module \c globe_stat_mod which defines all
!! global status variables used by the GLOBE routines.

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      module globe_stat_mod

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                       water balance
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable :: xglobe_H2Og_a(:)        !< H<sub>2</sub>O in boundary layer (rel. hum.)
      real,allocatable :: xglobe_H2Ol_v(:)        !< H<sub>2</sub>O in vegetation (water saturation)

      real,allocatable :: epglobe_H2Ol_s(:)       !< soil water redistribution

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                      carbon balance
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable :: xglobe_CO2g_a(:)        !< CO<sub>2</sub> in atmosphere (ppm)

      real,allocatable :: rglobe_Co_v(:)          !< biomass

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                       surface state
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable :: xglobe_T_a(:)           !< surface temperature
      real,allocatable :: globe_FVEG(:)           !< vegetation cover

      real             :: xglobe_pairg_as = 101325.0  !< surface pressure

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      end module globe_stat_mod
