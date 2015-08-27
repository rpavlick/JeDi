#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_sub.F90
!> \brief Other subroutines of the GLOBE module

!> \file globe_sub.F90
!> This file includes initialization subroutines of the GLOBE module
!! which do not fit in the other Fortran files.

!     ******************************************************************
!     GLOBE_FIELDS_INIT
!     ******************************************************************

!> \brief Initialization of variable fields with constant values

!> \detail This routine initializes some variable fields with constant
!! values.\n

!> \b Variables \b set: xglobe_T_a, globe_FVEG, fglobe_RADs_as,
!! fglobe_RADl_as
!>
!! - \b water \b variables: xglobe_H2Og_a, xglobe_H2Ol_v, fglobe_H2Ol_as,
!! fglobe_H2Os_as, fglobe_H2Ol_sv, fglobe_H2Ol_sr, fglobe_H2Ol_sb,
!! fglobe_H2Ol_br, fglobe_H2Og_sa, fglobe_H2Og_sa_pot, fglobe_H2Og_va
!>
!! - \b carbon \b variables: xglobe_CO2g_a, fglobe_CO2d_as, fglobe_CO2g_av,
!! fglobe_CO2g_va, fglobe_CO2g_vs, fglobe_CO2g_sa, fglobe_CO2d_sr,
!! fglobe_CO2d_br, fglobe_CO2d_sv, fglobe_CORG_vs, fglobe_Cco_vv

!     ******************************************************************

      subroutine globe_fields_init ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_fields_init')

!     * set some fields

      xglobe_T_a(:)      = 0.0  ! surface temperature
      globe_FVEG(:)      = 0.0

      fglobe_RADs_as(:)  = 0.0
      fglobe_RADl_as(:)  = 0.0

!     * set water variables

      xglobe_H2Og_a(:)   = 0.7
      xglobe_H2Ol_v(:)   = 0.0

      fglobe_H2Ol_as(:)  = 0.0
      fglobe_H2Os_as(:)  = 0.0
      fglobe_H2Ol_sv(:)  = 0.0
      fglobe_H2Ol_sr(:)  = 0.0
      fglobe_H2Ol_sb(:)  = 0.0
      fglobe_H2Ol_br(:)  = 0.0
      fglobe_H2Og_sa(:)  = 0.0
      fglobe_H2Og_sa_pot(:) = 0.0
      fglobe_H2Og_va(:)  = 0.0

!     * set carbon variables

      xglobe_CO2g_a(:)   = 360.0

      fglobe_CO2d_as(:)  = 0.0
      fglobe_CO2g_av(:)  = 0.0
      fglobe_CO2g_va(:)  = 0.0
      fglobe_CO2g_vs(:)  = 0.0
      fglobe_CO2g_sa(:)  = 0.0
      fglobe_CO2d_sr(:)  = 0.0
      fglobe_CO2d_br(:)  = 0.0
      fglobe_CO2d_sv(:)  = 0.0
      fglobe_CORG_vs(:)  = 0.0
      fglobe_Cco_vv(:)   = 0.0

!     * end

      __diag(kglobe_fdiag,'globe_fields_init: done')

      return
      end subroutine globe_fields_init

!     ******************************************************************
!     GLOBE_NPRO_FIELDS
!     ******************************************************************

!> \brief Initialization of variables related to parallel computing

!> \detail This routine determines field sizes related to the parallel
!! execution of the program and initializes the surface state
!! variables.\n
!! To split a variable field for the processing of its parts
!! with different CPUs used by the program at execution time, all parts
!! need to have the same size. This means that the field size has to be
!! a multiple of the number of CPUs used. To realize that, often
!! additional elements have to be added to the fields.\n
!! The first part calculates the field sizes according to the number of
!! CPUs and the size of the data fields.\n
!! The second part allocates the size corrected surface state variables,
!! related to the number of CPUs used.\n
!! The last part fills the size corrected surface state variables, where
!! additional fields are filled with a copy of the last field value.\n

!> \b Variables \b set: nglobe_cgpt, nglobe_gpts2, kglobe_gLPIDs2,
!! pglobe_gLP_lon2, pglobe_gLP_lat2, pglobe_gLP_hght2,
   !> pglobe_gLP_paw2,
!> pglobe_gLP_len2, pglobe_gLP_wid2, pglobe_gLP_area2, pglobe_gLP_togr2

!     ******************************************************************

      subroutine globe_npro_fields ()
      use globe_mod
      implicit none

      integer :: j, kextrapoints

      __diag(kglobe_fdiag,'globe_npro_fields')

!     * compute grid field sizes given the number of CPUs

      kextrapoints = MOD(nglobe_landpts, nglobe_npro)
      if (kextrapoints .ne. 0) then
        nglobe_cgpt = (nglobe_landpts + nglobe_npro - kextrapoints) / nglobe_npro
      else
        nglobe_cgpt = nglobe_landpts / nglobe_npro
      endif
      nglobe_gpts2 = nglobe_cgpt * nglobe_npro

!     * output grid partitioning

      if (kglobe_mypid == kglobe_kroot) then
        __diag_num(kglobe_fdiag,'globe_npro_fields: kextrapoints = ',kextrapoints)
        __diag_num(kglobe_fdiag,'globe_npro_fields: nglobe_cgpt  = ',nglobe_cgpt)
        __diag_num(kglobe_fdiag,'globe_npro_fields: nglobe_gpts2 = ',nglobe_gpts2)
      endif

!     * allocate CPU-corrected lists

      allocate(kglobe_gLPIDs2(nglobe_gpts2))
      allocate(pglobe_gLP_lon2(nglobe_gpts2))
      allocate(pglobe_gLP_lat2(nglobe_gpts2))
      allocate(pglobe_gLP_hght2(nglobe_gpts2))
      allocate(pglobe_gLP_paw2(nglobe_gpts2))
      allocate(pglobe_gLP_len2(nglobe_gpts2))
      allocate(pglobe_gLP_wid2(nglobe_gpts2))
      allocate(pglobe_gLP_area2(nglobe_gpts2))
      allocate(pglobe_gLP_togr2(nglobe_gpts2))

!     * copy land grid point list to CPU-corrected lists

      if (kglobe_mypid == kglobe_kroot) then
        kglobe_gLPIDs2(1:nglobe_landpts)   = kglobe_gLPIDs(:)
        pglobe_gLP_lon2(1:nglobe_landpts)  = pglobe_gLP_lon(:)
        pglobe_gLP_lat2(1:nglobe_landpts)  = pglobe_gLP_lat(:)
        pglobe_gLP_hght2(1:nglobe_landpts) = pglobe_gLP_elevation(:)
        pglobe_gLP_paw2(1:nglobe_landpts)  = pglobe_gLP_paw(:)
        pglobe_gLP_togr2(1:nglobe_landpts) = pglobe_gLP_topograd(:)
        pglobe_gLP_len2(1:nglobe_landpts)  = pglobe_gLP_len(:)
        pglobe_gLP_wid2(1:nglobe_landpts)  = pglobe_gLP_wid(:)
        pglobe_gLP_area2(1:nglobe_landpts) = pglobe_gLP_area(:)
        if (nglobe_gpts2 .gt. nglobe_landpts) then
          do j = nglobe_landpts + 1, nglobe_gpts2
            kglobe_gLPIDs2(j)   = kglobe_gLPIDs(nglobe_landpts)
            pglobe_gLP_lon2(j)  = pglobe_gLP_lon(nglobe_landpts)
            pglobe_gLP_lat2(j)  = pglobe_gLP_lat(nglobe_landpts)
            pglobe_gLP_hght2(j) = pglobe_gLP_elevation(nglobe_landpts)
            pglobe_gLP_paw2(j)  = pglobe_gLP_paw(nglobe_landpts)
            pglobe_gLP_togr2(j) = pglobe_gLP_topograd(nglobe_landpts)
            pglobe_gLP_len2(j)  = pglobe_gLP_len(nglobe_landpts)
            pglobe_gLP_wid2(j)  = pglobe_gLP_wid(nglobe_landpts)
            pglobe_gLP_area2(j) = pglobe_gLP_area(nglobe_landpts)
          enddo
        endif
      endif

!     * end

      __diag(kglobe_fdiag,'globe_npro_fields: done')

      return
      end subroutine globe_npro_fields

!     ******************************************************************
!     GLOBE_LANDSEA_ALLOC
!     ******************************************************************

!> \brief Allocation of fields read from the surface state file

!> \detail This routine allocates the fields read from the surface state
!! file.

!     ******************************************************************

      subroutine globe_landsea_alloc ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_landsea_alloc')

      if (kglobe_mypid == kglobe_kroot) then
        __allocate(kglobe_gLPIDs,(nglobe_landpts))

        __allocate(pglobe_longitude,(nglobe_nlon*nglobe_nlat))
        __allocate(pglobe_latitude,(nglobe_nlon*nglobe_nlat))
        __allocate(pglobe_gLP_lon,(nglobe_landpts))
        __allocate(pglobe_gLP_lat,(nglobe_landpts))
        __allocate(pglobe_gLP_len,(nglobe_landpts))
        __allocate(pglobe_gLP_wid,(nglobe_landpts))
        __allocate(pglobe_gLP_area,(nglobe_landpts))

        __allocate(pglobe_gLP_elevation,(nglobe_landpts))
        __allocate(pglobe_gLP_paw,(nglobe_landpts))
        __allocate(pglobe_gLP_topograd,(nglobe_landpts))

        __allocate(pglobe_gLP_numnn,(nglobe_landpts))
        __allocate(pglobe_gLP_idnn,(nglobe_landpts,8))
        __allocate(pglobe_gLP_distnn,(nglobe_landpts,8))
        __allocate(kglobe_lpt_nn,(nglobe_landpts,8))
      endif

!     * allocate the temporary sub-procedure variables

      __allocate(globe_tmp_lonlat_4,(nglobe_nlon*nglobe_nlat))
      __allocate(globe_tmp_lonlat,(nglobe_nlon*nglobe_nlat))
      __allocate(globe_tmp_nLPs_4,(nglobe_landpts))

      __diag(kglobe_fdiag,'globe_landsea_alloc: done')

      return
      end subroutine globe_landsea_alloc

!     ******************************************************************
!     GLOBE_ALLOC
!     ******************************************************************

!> \brief Allocation of fields of the GLOBE module

!> \detail This routine allocates most of the fields of the GLOBE module.

!     ******************************************************************

      subroutine globe_alloc ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_alloc')

!     * allocate surface state fields

      __allocate(xglobe_T_a,(nglobe_cgpt))
      __allocate(globe_FVEG,(nglobe_cgpt))

!     * allocate energy balance fields

      __allocate(fglobe_RADs_as,(nglobe_cgpt))
      __allocate(fglobe_RADl_as,(nglobe_cgpt))

!     * allocate water fields

      __allocate(xglobe_H2Og_a,(nglobe_cgpt))
      __allocate(xglobe_H2Ol_v,(nglobe_cgpt))

      __allocate(fglobe_H2Ol_as,(nglobe_cgpt))
      __allocate(fglobe_H2Os_as,(nglobe_cgpt))
      __allocate(fglobe_H2Ol_sv,(nglobe_cgpt))
      __allocate(fglobe_H2Ol_sr,(nglobe_cgpt))
      __allocate(fglobe_H2Ol_sb,(nglobe_cgpt))
      __allocate(fglobe_H2Ol_br,(nglobe_cgpt))
      __allocate(fglobe_H2Og_sa,(nglobe_cgpt))
      __allocate(fglobe_H2Og_sa_pot,(nglobe_cgpt))
      __allocate(fglobe_H2Og_va,(nglobe_cgpt))

!     * allocate carbon fields

      __allocate(xglobe_CO2g_a,(nglobe_cgpt))
      __allocate(rglobe_Co_v,(nglobe_cgpt))

      __allocate(fglobe_CO2d_as,(nglobe_cgpt))
      __allocate(fglobe_CO2d_sr,(nglobe_cgpt))
      __allocate(fglobe_CO2d_br,(nglobe_cgpt))
      __allocate(fglobe_CO2d_sv,(nglobe_cgpt))
      __allocate(fglobe_CO2d_sb,(nglobe_cgpt))

      __allocate(fglobe_CO2g_sa,(nglobe_cgpt))
      __allocate(fglobe_CO2g_av,(nglobe_cgpt))
      __allocate(fglobe_CO2g_va,(nglobe_cgpt))
      __allocate(fglobe_CO2g_vs,(nglobe_cgpt))
      __allocate(fglobe_CORG_vs,(nglobe_cgpt))
      __allocate(fglobe_Cco_vv,(nglobe_cgpt))

!     * allocate the temporary subprocedure variables
      __allocate(globe_tmp_nLPs,(nglobe_landpts))
      __allocate(globe_tmp_nGPts2,(nglobe_gpts2))

!     * end

      __diag(kglobe_fdiag,'globe_alloc: done')

      return
      end subroutine globe_alloc

!     ******************************************************************
!     GLOBE_DEALLOC
!     ******************************************************************

!> \brief Deallocation of all GLOBE module fields

!> \detail This routine deallocates all fields of the GLOBE module.

!     ******************************************************************

      subroutine globe_dealloc ()
      use globe_mod
      implicit none

!     ------------------------------------------------------------------
!     * grid description
!     ------------------------------------------------------------------

      __deallocate(kglobe_gLPIDs)
      __deallocate(pglobe_latitude)
      __deallocate(pglobe_longitude)
      __deallocate(pglobe_gLP_lat)
      __deallocate(pglobe_gLP_lon)
      __deallocate(pglobe_gLP_elevation)
      __deallocate(pglobe_gLP_paw)
      __deallocate(pglobe_gLP_topograd)
      __deallocate(pglobe_gLP_len)
      __deallocate(pglobe_gLP_wid)
      __deallocate(pglobe_gLP_area)

      __deallocate(kglobe_gLPIDs2)
      __deallocate(pglobe_gLP_lat2)
      __deallocate(pglobe_gLP_lon2)
      __deallocate(pglobe_gLP_hght2)
      __deallocate(pglobe_gLP_paw2)
      __deallocate(pglobe_gLP_togr2)
      __deallocate(pglobe_gLP_len2)
      __deallocate(pglobe_gLP_wid2)
      __deallocate(pglobe_gLP_area2)

#ifndef __JEDI
!     ------------------------------------------------------------------
!     * array for saving nearest neighbors and their indices
!     ------------------------------------------------------------------

      __deallocate(pglobe_gLP_numnn)
      __deallocate(pglobe_gLP_idnn)
      __deallocate(pglobe_gLP_distnn)
      __deallocate(kglobe_lpt_nn)
#endif /* __JEDI */

!     ------------------------------------------------------------------
!     * exchange fluxes
!     ------------------------------------------------------------------

!     * surface state

      __deallocate(xglobe_T_a)
      __deallocate(globe_FVEG)

!     * energy balance

      __deallocate(fglobe_RADs_as)
      __deallocate(fglobe_RADl_as)

!     * water balance

      __deallocate(xglobe_H2Og_a)
      __deallocate(xglobe_H2Ol_v)

      __deallocate(fglobe_H2Ol_as)
      __deallocate(fglobe_H2Os_as)
      __deallocate(fglobe_H2Ol_sv)
      __deallocate(fglobe_H2Ol_sr)
      __deallocate(fglobe_H2Ol_sb)
      __deallocate(fglobe_H2Ol_br)
      __deallocate(fglobe_H2Og_sa)
      __deallocate(fglobe_H2Og_sa_pot)
      __deallocate(fglobe_H2Og_va)

!     * carbon balance

      __deallocate(xglobe_CO2g_a)
      __deallocate(rglobe_Co_v)

      __deallocate(fglobe_CO2d_as)
      __deallocate(fglobe_CO2d_sb)
      __deallocate(fglobe_CO2g_sa)

      __deallocate(fglobe_CO2g_av)
      __deallocate(fglobe_CO2g_va)
      __deallocate(fglobe_CO2g_vs)
      __deallocate(fglobe_CORG_vs)
      __deallocate(fglobe_Cco_vv)

      __deallocate(fglobe_CO2d_sr)
      __deallocate(fglobe_CO2d_br)
      __deallocate(fglobe_CO2d_sv)

!     * temporary sub-procedure variables

      __deallocate(globe_tmp_lonlat_4)
      __deallocate(globe_tmp_lonlat)
      __deallocate(globe_tmp_nLPs_4)
      __deallocate(globe_tmp_nLPs)
      __deallocate(globe_tmp_nGPts2)

#ifndef __PLASIM
!     * variables of plasim

      __deallocate(globe_plasim_dalb)
      __deallocate(globe_plasim_dcsoil)
      __deallocate(globe_plasim_dcveg)
      __deallocate(globe_plasim_devap)
      __deallocate(globe_plasim_dpet)
      __deallocate(globe_plasim_dflux)
      __deallocate(globe_plasim_dforest)
      __deallocate(globe_plasim_dglac)
      __deallocate(globe_plasim_dgpp)
      __deallocate(globe_plasim_dgppl)
      __deallocate(globe_plasim_dgppw)
      __deallocate(globe_plasim_dlai)
      __deallocate(globe_plasim_dlhfl)
      __deallocate(globe_plasim_dlitter)
      __deallocate(globe_plasim_dls)
      __deallocate(globe_plasim_dnogrow)
      __deallocate(globe_plasim_dnpp)
      __deallocate(globe_plasim_dp)
      __deallocate(globe_plasim_dprc)
      __deallocate(globe_plasim_dprl)
      __deallocate(globe_plasim_dprs)
      __deallocate(globe_plasim_dq)
      __deallocate(globe_plasim_dres)
      __deallocate(globe_plasim_drhs)
      __deallocate(globe_plasim_drunoff)
      __deallocate(globe_plasim_dshfl)
      __deallocate(globe_plasim_dsmelt)
      __deallocate(globe_plasim_dsndch)
      __deallocate(globe_plasim_dsnow)
      __deallocate(globe_plasim_dswfl)
      __deallocate(globe_plasim_dt)
      __deallocate(globe_plasim_dtd2)
      __deallocate(globe_plasim_dtd3)
      __deallocate(globe_plasim_dtd4)
      __deallocate(globe_plasim_dtd5)
      __deallocate(globe_plasim_dtsoil)
      __deallocate(globe_plasim_dveg)
      __deallocate(globe_plasim_dwatc)
      __deallocate(globe_plasim_dwmax)
      __deallocate(globe_plasim_dz0)

!     * variables of landmod

      __deallocate(globe_plasim_doro)
      __deallocate(globe_plasim_dts)
      __deallocate(globe_plasim_dtsm)
      __deallocate(globe_plasim_dqs)
      __deallocate(globe_plasim_driver)
      __deallocate(globe_plasim_duroff)
      __deallocate(globe_plasim_dvroff)
      __deallocate(globe_plasim_darea)
      __deallocate(globe_plasim_dsoilz)
      __deallocate(globe_plasim_dsoilt)
      __deallocate(globe_plasim_dsnowt)
      __deallocate(globe_plasim_dtclsoil)
      __deallocate(globe_plasim_dsnowz)
      __deallocate(globe_plasim_dwater)
      __deallocate(globe_plasim_pgrow)
      __deallocate(globe_plasim_plai)
      __deallocate(globe_plasim_pgs)
      __deallocate(globe_plasim_pz0_max)
      __deallocate(globe_plasim_dtcl)
      __deallocate(globe_plasim_dwcl)
      __deallocate(globe_plasim_dtclim)
      __deallocate(globe_plasim_dwclim)
      __deallocate(globe_plasim_dz0clim)
      __deallocate(globe_plasim_dz0climo)
      __deallocate(globe_plasim_dalbclim)
      __deallocate(globe_plasim_dh2ol_added)
      __deallocate(globe_plasim_geomask)
#endif /* __PLASIM */

      return
      end subroutine globe_dealloc

#ifndef __PLASIM
!     ******************************************************************
!     GLOBE_PLASIM_ALLOC
!     ******************************************************************

!> \brief Allocation of fields related to the PLASIM module

!> \detail This routine allocates all fields which are used to transfer
!! data between PLASIM, SIMBA, and the land module of GLOBE, extracted
!! from PLASIM.

!     ******************************************************************

      subroutine globe_plasim_alloc ()
      use globe_mod
      implicit none

!     * variables of plasim

      allocate(globe_plasim_dalb(nglobe_cgpt))
      allocate(globe_plasim_dcsoil(nglobe_cgpt))
      allocate(globe_plasim_dcveg(nglobe_cgpt))
      allocate(globe_plasim_devap(nglobe_cgpt))
      allocate(globe_plasim_dpet(nglobe_cgpt))
      allocate(globe_plasim_dflux(nglobe_cgpt,globe_plasim_NLEP))
      allocate(globe_plasim_dforest(nglobe_cgpt))
      allocate(globe_plasim_dglac(nglobe_cgpt))
      allocate(globe_plasim_dgpp(nglobe_cgpt))
      allocate(globe_plasim_dgppl(nglobe_cgpt))
      allocate(globe_plasim_dgppw(nglobe_cgpt))
      allocate(globe_plasim_dlai(nglobe_cgpt))
      allocate(globe_plasim_dlhfl(nglobe_cgpt))
      allocate(globe_plasim_dlitter(nglobe_cgpt))
      allocate(globe_plasim_dls(nglobe_cgpt))
      allocate(globe_plasim_dnogrow(nglobe_cgpt))
      allocate(globe_plasim_dnpp(nglobe_cgpt))
      allocate(globe_plasim_dp(nglobe_cgpt))
      allocate(globe_plasim_dprc(nglobe_cgpt))
      allocate(globe_plasim_dprl(nglobe_cgpt))
      allocate(globe_plasim_dprs(nglobe_cgpt))
      allocate(globe_plasim_dq(nglobe_cgpt,globe_plasim_NLEP))
      allocate(globe_plasim_dres(nglobe_cgpt))
      allocate(globe_plasim_drhs(nglobe_cgpt))
      allocate(globe_plasim_drunoff(nglobe_cgpt))
      allocate(globe_plasim_dshfl(nglobe_cgpt))
      allocate(globe_plasim_dsmelt(nglobe_cgpt))
      allocate(globe_plasim_dsndch(nglobe_cgpt))
      allocate(globe_plasim_dsnow(nglobe_cgpt))
      allocate(globe_plasim_dswfl(nglobe_cgpt,globe_plasim_NLEP))
      allocate(globe_plasim_dt(nglobe_cgpt,globe_plasim_NLEP))
      allocate(globe_plasim_dtd2(nglobe_cgpt))
      allocate(globe_plasim_dtd3(nglobe_cgpt))
      allocate(globe_plasim_dtd4(nglobe_cgpt))
      allocate(globe_plasim_dtd5(nglobe_cgpt))
      allocate(globe_plasim_dtsoil(nglobe_cgpt))
      allocate(globe_plasim_dveg(nglobe_cgpt))
      allocate(globe_plasim_dwatc(nglobe_cgpt))
      allocate(globe_plasim_dwmax(nglobe_cgpt))
      allocate(globe_plasim_dz0(nglobe_cgpt))

      return
      end subroutine globe_plasim_alloc

!     ******************************************************************
!     GLOBE_LANDMOD_ALLOC
!     ******************************************************************

!> \brief Allocation of fields related to the land module of PLASIM

!> \detail This routine allocates all fields related to the data
!! transfer to the GLOBE land module, which was extracted from PLASIM.

!     ******************************************************************

      subroutine globe_landmod_alloc ()
      use globe_mod
      implicit none

!     * variables of landmod

      allocate(globe_plasim_doro(nglobe_cgpt))
      allocate(globe_plasim_dts(nglobe_cgpt))
      allocate(globe_plasim_dtsm(nglobe_cgpt))
      allocate(globe_plasim_dqs(nglobe_cgpt))
      allocate(globe_plasim_driver(nglobe_cgpt))
      allocate(globe_plasim_duroff(nglobe_cgpt))
      allocate(globe_plasim_dvroff(nglobe_cgpt))
      allocate(globe_plasim_darea(nglobe_cgpt))
      allocate(globe_plasim_dsoilz(globe_plasim_NLSOIL))
      allocate(globe_plasim_dsoilt(nglobe_cgpt,globe_plasim_NLSOIL))
      allocate(globe_plasim_dsnowt(nglobe_cgpt))
      allocate(globe_plasim_dtclsoil(nglobe_cgpt))
      allocate(globe_plasim_dsnowz(nglobe_cgpt))
      allocate(globe_plasim_dwater(nglobe_cgpt))
      allocate(globe_plasim_pgrow(nglobe_cgpt))
      allocate(globe_plasim_plai(nglobe_cgpt))
      allocate(globe_plasim_pgs(nglobe_cgpt))
      allocate(globe_plasim_pz0_max(nglobe_cgpt))
      allocate(globe_plasim_dtcl(nglobe_cgpt,0:13))
      allocate(globe_plasim_dwcl(nglobe_cgpt,0:13))
      allocate(globe_plasim_dtclim(nglobe_cgpt))
      allocate(globe_plasim_dwclim(nglobe_cgpt))
      allocate(globe_plasim_dz0clim(nglobe_cgpt))
      allocate(globe_plasim_dz0climo(nglobe_cgpt))
      allocate(globe_plasim_dalbclim(nglobe_cgpt))
      allocate(globe_plasim_dh2ol_added(nglobe_cgpt))
      allocate(globe_plasim_geomask(nglobe_cgpt))

      return
      end subroutine globe_landmod_alloc
#endif /* __PLASIM */
