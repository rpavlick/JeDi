#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_mod_flux.F90
!> \brief global flux variable definitions

!> \file globe_mod_flux.F90
!> This file contains the \c module \c globe_flux_mod which defines all
!! global flux variables used by the GLOBE routines.

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      module globe_flux_mod

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                      Energy balance
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable    :: fglobe_RADs_as(:)    !< short wave radiation
      real,allocatable    :: fglobe_RADl_as(:)    !< long wave radiation

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                       Water balance
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable    :: fglobe_H2Ol_as(:)    !< precipitation
      real,allocatable    :: fglobe_H2Os_as(:)    !< snowfall
      real,allocatable    :: fglobe_H2Ol_sv(:)    !< root water uptake
      real,allocatable    :: fglobe_H2Ol_sr(:)    !< saturation excess flow
      real,allocatable    :: fglobe_H2Ol_sb(:)    !< soil drainage
      real,allocatable    :: fglobe_H2Ol_br(:)    !< river discharge
      real,allocatable    :: fglobe_H2Og_sa(:)    !< bare soil evaporation
      real,allocatable    :: fglobe_H2Og_sa_pot(:)!< potential evaporation
      real,allocatable    :: fglobe_H2Og_va(:)    !< transpiration

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                      Carbon balance
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable    :: fglobe_CO2d_as(:)    !< CO<sup>2</sup> dissolved in rain
      real,allocatable    :: fglobe_CO2g_sa(:)    !< soil-atmosphere gas exchange

      real,allocatable    :: fglobe_CO2g_av(:)    !< gross carbon uptake
      real,allocatable    :: fglobe_CO2g_va(:)    !< aboveground respiration
      real,allocatable    :: fglobe_CO2g_vs(:)    !< belowground respiration
      real,allocatable    :: fglobe_CORG_vs(:)    !< litter fall
      real,allocatable    :: fglobe_Cco_vv(:)     !< NPP

      real,allocatable    :: fglobe_CO2d_sr(:)    !< dissolved carbon export
      real,allocatable    :: fglobe_CO2d_sb(:)    !< dissolved carbon flow to base
      real,allocatable    :: fglobe_CO2d_br(:)    !< dissolved carbon groundwater to river
      real,allocatable    :: fglobe_CO2d_sv(:)    !< dissolved carbon in root uptake

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      end module globe_flux_mod
