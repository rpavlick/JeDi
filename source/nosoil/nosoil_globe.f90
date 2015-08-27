!     ------------------------------------------------------------------
!     GLOBE INTERFACE FOR A SETUP WITHOUT AN EXPLICIT LAND MODEL
!     ------------------------------------------------------------------

      subroutine soil_init ()
      use globe_mod
      implicit none

      if (kglobe_mypid == kglobe_kroot)                                &
     &  write(*,*) 'nosoil:  no soil simulated'

      return
      end subroutine soil_init

!     ------------------------------------------------------------------

      subroutine soil_step (globe_step_id)
      implicit none

      integer :: globe_step_id

      return
      end subroutine soil_step

!     ------------------------------------------------------------------

      subroutine soil_stop ()
      implicit none

      return
      end subroutine soil_stop

!     ------------------------------------------------------------------

      subroutine soil_output ()
      implicit none

      return
      end subroutine soil_output

!     ------------------------------------------------------------------

      subroutine soil_diag ()
      implicit none

      return
      end subroutine soil_diag

!     ------------------------------------------------------------------

      subroutine soil_findfield (kFldCode)
      implicit none

      integer       :: kFldCode

      return
      end subroutine soil_findfield

!     ------------------------------------------------------------------

      subroutine soil_skip_reset ()
      implicit none

      return
      end subroutine soil_skip_reset

!     ------------------------------------------------------------------

      subroutine soil_skip_max (pDTskip)
      implicit none

      real :: pDTskip

      return
      end subroutine soil_skip_max

!     ------------------------------------------------------------------

      subroutine soil_skip (nskip_yrs)
      implicit none

      integer :: nskip_yrs

      return
      end subroutine soil_skip
