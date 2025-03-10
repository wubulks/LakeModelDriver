#include <define.h>

MODULE MOD_Lake_1DAccVars
! ---------------------------------- code history -------------------------------------
! Description:
!     Define the lake accumulation variables and subroutines for the lake module.
!
! Subroutines:
!     AllocLakeAccVars  : Allocate memory for lake variables.
!     ReleaseLakeAccVars: Release memory for lake variables.
!     FlushLakeAccVars  : Flush lake variables.
!     LakeVarsAcc       : Accumulate lake variables.
!     LakeVarsSaveHist  : Save lake variables to history file.
!s
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
!================================================================================
   USE MOD_Precision
   USE MOD_Vars_Global
!================================================================================

   real(r8), allocatable :: a_dplak          (:)  !lake depth [m]
   real(r8), allocatable :: a_zlake        (:,:)  !Lake layer node depth [m]
   real(r8), allocatable :: a_zilak        (:,:)  !Lake layer interface depth [m]
   real(r8), allocatable :: a_dzlake       (:,:)  !Lake layer thickness [m]
   real(r8), allocatable :: a_ziarea       (:,:)  !Lake layerinterface area [m2], only for Simstrat
   real(r8), allocatable :: a_lkz0m          (:)  !Roughness length for momentum [m]
   real(r8), allocatable :: a_lkz0h          (:)  !Roughness length for sensible heat  [m]
   real(r8), allocatable :: a_lkz0q          (:)  !Roughness length for latent heat [m]
   real(r8), allocatable :: a_felak          (:)  !Lake fetch length [m]
   real(r8), allocatable :: a_gamma          (:)  !Mixing enhancement factor, the meaning is different for each model [-]
   real(r8), allocatable :: a_etal           (:)  !Lake extinction coefficient [1/m]
   real(r8), allocatable :: a_btpri          (:)  !Beta prime in Monin-Obukhov theory [-]
   real(r8), allocatable :: a_frlak          (:)  !Lake fraction [-] 
   real(r8), allocatable :: a_tmsno          (:)  !Mean snow temperature, only for FLake [K], only for FLake
   real(r8), allocatable :: a_tmice          (:)  !Mean ice temperature, only for FLake [K], only for FLake
   real(r8), allocatable :: a_tmmnw          (:)  !Mean temperature of the water column [K]
   real(r8), allocatable :: a_tmwml          (:)  !Mixed-layer temperature [K], only for FLake
   real(r8), allocatable :: a_tmbot          (:)  !Temperature at the water-bottom sediment interface [K], only for FLake
   real(r8), allocatable :: a_tmups          (:)  !Temperature at the bottom of the upper layer of the sediments [K], only for FLake
   real(r8), allocatable :: a_mldp           (:)  !Mixed layer depth [m], only for FLake
   real(r8), allocatable :: a_upsdp          (:)  !Bottom of the upper layer of the sediments [m], only for FLake
   real(r8), allocatable :: a_icedp          (:)  !Mean temperature of the lake [K], for FLake and Simstrat
   real(r8), allocatable :: a_bicedp         (:)  !black ice depth [m], only for Simstrat
   real(r8), allocatable :: a_wicedp         (:)  !white ice depth [m], only for Simstrat
   real(r8), allocatable :: a_CTfrac         (:)  !Shape factor (thermocline)
   real(r8), allocatable :: a_rhosnw         (:)  !snow density [kg/m3], only for Simstrat
   real(r8), allocatable :: a_uwatv        (:,:)  !Water velocity in x-direction [m/s], only for Simstrat
   real(r8), allocatable :: a_vwatv        (:,:)  !Water velocity in y-direction [m/s], only for Simstrat
   real(r8), allocatable :: a_lksal        (:,:)  !Salinity [‰], only for Simstrat
   real(r8), allocatable :: a_tke          (:,:)  !Turbulent kinetic energy (TKE) [J/kg], only for Simstrat
   real(r8), allocatable :: a_etke           (:)  !Seiche energy [J], only for Simstrat
   real(r8), allocatable :: a_eps          (:,:)  !TKE dissipation rate [W/kg], only for Simstrat
   real(r8), allocatable :: a_num          (:,:)  !Turbulent viscosity (momentum) [m2/s], only for Simstrat
   real(r8), allocatable :: a_nuh          (:,:)  !Turbulent diffusivity (heat) [m2/s], only for Simstrat
   real(r8), allocatable :: a_lkrho        (:,:)  !Density of water [kg/m3]

   PUBLIC :: allocate_LakeAccVars
   PUBLIC :: deallocate_LakeAccVars
   PUBLIC :: Flush_LakeAccVars
   PUBLIC :: accumulate_LakeTimeVars

   PRIVATE :: acc1d_lake, acc2d_lake

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_LakeAccVars ()
! ---------------------------------- code history -------------------------------------
! Description:
!     Allocate memory for lake variables. Will be called in the subroutine allocate_acc_fluxes.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE
!==============================================================================
!==============================================================================
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (a_dplak               (numpatch)) 
            allocate (a_zlake       (nl_lake,numpatch)) 
            allocate (a_zilak     (nl_lake+1,numpatch)) 
            allocate (a_dzlake      (nl_lake,numpatch))
            allocate (a_ziarea    (nl_lake+1,numpatch)) 
            allocate (a_lkz0m               (numpatch)) 
            allocate (a_lkz0h               (numpatch)) 
            allocate (a_lkz0q               (numpatch)) 
            allocate (a_felak               (numpatch)) 
            allocate (a_gamma               (numpatch)) 
            allocate (a_etal                (numpatch)) 
            allocate (a_btpri               (numpatch)) 
            allocate (a_frlak               (numpatch)) 
            allocate (a_tmsno               (numpatch)) 
            allocate (a_tmice               (numpatch)) 
            allocate (a_tmmnw               (numpatch)) 
            allocate (a_tmwml               (numpatch)) 
            allocate (a_tmbot               (numpatch)) 
            allocate (a_tmups               (numpatch)) 
            allocate (a_mldp                (numpatch)) 
            allocate (a_upsdp               (numpatch)) 
            allocate (a_icedp               (numpatch)) 
            allocate (a_bicedp              (numpatch)) 
            allocate (a_wicedp              (numpatch))
            allocate (a_CTfrac              (numpatch)) 
            allocate (a_rhosnw              (numpatch)) 
            allocate (a_uwatv       (nl_lake,numpatch)) 
            allocate (a_vwatv       (nl_lake,numpatch)) 
            allocate (a_lksal       (nl_lake,numpatch)) 
            allocate (a_tke       (nl_lake+1,numpatch)) 
            allocate (a_etke                (numpatch)) 
            allocate (a_eps       (nl_lake+1,numpatch)) 
            allocate (a_num       (nl_lake+1,numpatch)) 
            allocate (a_nuh       (nl_lake+1,numpatch)) 
            allocate (a_lkrho       (nl_lake,numpatch)) 
         ENDIF
      ENDIF

   END SUBROUTINE allocate_LakeAccVars



   SUBROUTINE deallocate_LakeAccVars ()
! ---------------------------------- code history -------------------------------------
! Description:
!     Release memory for lake variables. Will be called in the subroutine deallocate_acc_fluxes.
!
! Original author:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE
!==============================================================================
!==============================================================================
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate ( a_dplak  )  
            deallocate ( a_zlake  )  
            deallocate ( a_zilak  )  
            deallocate ( a_dzlake ) 
            deallocate ( a_ziarea ) 
            deallocate ( a_lkz0m  )  
            deallocate ( a_lkz0h  )  
            deallocate ( a_lkz0q  )  
            deallocate ( a_felak  )  
            deallocate ( a_gamma  )  
            deallocate ( a_etal   )  
            deallocate ( a_btpri  )  
            deallocate ( a_frlak  )  
            deallocate ( a_tmsno  )  
            deallocate ( a_tmice  )  
            deallocate ( a_tmmnw  )  
            deallocate ( a_tmwml  )  
            deallocate ( a_tmbot  )  
            deallocate ( a_tmups  )  
            deallocate ( a_mldp   )  
            deallocate ( a_upsdp  )  
            deallocate ( a_icedp  )  
            deallocate ( a_bicedp )  
            deallocate ( a_wicedp )
            deallocate ( a_CTfrac )  
            deallocate ( a_rhosnw )  
            deallocate ( a_uwatv  )  
            deallocate ( a_vwatv  )  
            deallocate ( a_lksal  )  
            deallocate ( a_tke    )  
            deallocate ( a_etke   )  
            deallocate ( a_eps    )  
            deallocate ( a_num    )  
            deallocate ( a_nuh    )  
            deallocate ( a_lkrho  ) 
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_LakeAccVars



   SUBROUTINE Flush_LakeAccVars ()
! ---------------------------------- code history -------------------------------------
! Description:
!     flush lake variables. Will be called in the subroutine FLUSH_acc_fluxes.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      USE MOD_Vars_Global, only: spval
      IMPLICIT NONE
!==============================================================================
!==============================================================================
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            a_dplak         (:) = spval 
            a_zlake       (:,:) = spval 
            a_zilak       (:,:) = spval 
            a_dzlake      (:,:) = spval
            a_ziarea      (:,:) = spval 
            a_lkz0m         (:) = spval 
            a_lkz0h         (:) = spval 
            a_lkz0q         (:) = spval 
            a_felak         (:) = spval 
            a_gamma         (:) = spval 
            a_etal          (:) = spval 
            a_btpri         (:) = spval 
            a_frlak         (:) = spval 
            a_tmsno         (:) = spval 
            a_tmice         (:) = spval 
            a_tmmnw         (:) = spval 
            a_tmwml         (:) = spval 
            a_tmbot         (:) = spval 
            a_tmups         (:) = spval 
            a_mldp          (:) = spval 
            a_upsdp         (:) = spval 
            a_icedp         (:) = spval 
            a_bicedp        (:) = spval 
            a_wicedp        (:) = spval
            a_CTfrac        (:) = spval 
            a_rhosnw        (:) = spval 
            a_uwatv       (:,:) = spval 
            a_vwatv       (:,:) = spval 
            a_lksal       (:,:) = spval 
            a_tke         (:,:) = spval 
            a_etke          (:) = spval 
            a_eps         (:,:) = spval 
            a_num         (:,:) = spval 
            a_nuh         (:,:) = spval 
            a_lkrho       (:,:) = spval 
         ENDIF
      ENDIF

   END SUBROUTINE Flush_LakeAccVars



   SUBROUTINE accumulate_LakeTimeVars ()
! ---------------------------------- code history -------------------------------------
! Description:
!     Accumulate lake variables. Will be called in the subroutine accumulate_fluxes.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_Lake_TimeVars
   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE
   integer :: i
!==============================================================================
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            CALL acc1d_lake (dplak       , a_dplak        )   
            CALL acc2d_lake (zlake       , a_zlake        )   
            CALL acc2d_lake (zilak       , a_zilak        )   
            CALL acc2d_lake (dzlake      , a_dzlake       )  
            CALL acc2d_lake (ziarea      , a_ziarea       )  
            CALL acc1d_lake (lkz0m       , a_lkz0m        )   
            CALL acc1d_lake (lkz0h       , a_lkz0h        )   
            CALL acc1d_lake (lkz0q       , a_lkz0q        )   
            CALL acc1d_lake (felak       , a_felak        )   
            CALL acc1d_lake (gamma       , a_gamma        )   
            CALL acc1d_lake (etal        , a_etal         )   
            CALL acc1d_lake (btpri       , a_btpri        )   
            CALL acc1d_lake (frlak       , a_frlak        )   
            CALL acc1d_lake (tmsno       , a_tmsno        )   
            CALL acc1d_lake (tmice       , a_tmice        )   
            CALL acc1d_lake (tmmnw       , a_tmmnw        )   
            CALL acc1d_lake (tmwml       , a_tmwml        )   
            CALL acc1d_lake (tmbot       , a_tmbot        )   
            CALL acc1d_lake (tmups       , a_tmups        )   
            CALL acc1d_lake (mldp        , a_mldp         )   
            CALL acc1d_lake (upsdp       , a_upsdp        )   
            CALL acc1d_lake (icedp       , a_icedp        )   
            CALL acc1d_lake (bicedp      , a_bicedp       )   
            CALL acc1d_lake (wicedp      , a_wicedp       )
            CALL acc1d_lake (CTfrac      , a_CTfrac       )
            CALL acc1d_lake (rhosnw      , a_rhosnw       )   
            CALL acc2d_lake (uwatv       , a_uwatv        )   
            CALL acc2d_lake (vwatv       , a_vwatv        )   
            CALL acc2d_lake (lksal       , a_lksal        )   
            CALL acc2d_lake (tke         , a_tke          )   
            CALL acc1d_lake (etke        , a_etke         )   
            CALL acc2d_lake (eps         , a_eps          )   
            CALL acc2d_lake (num         , a_num          )   
            CALL acc2d_lake (nuh         , a_nuh          )   
            CALL acc2d_lake (lkrho       , a_lkrho        )  
            
         ENDIF
      ENDIF


   END SUBROUTINE accumulate_LakeTimeVars


   !------
   SUBROUTINE acc1d_lake (var, s)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:)
   real(r8), intent(inout) :: s  (:)
   ! Local variables
   integer :: i

      DO i = lbound(var,1), ubound(var,1)
         IF (var(i) /= spval) THEN
            IF (s(i) /= spval) THEN
               s(i) = s(i) + var(i)
            ELSE
               s(i) = var(i)
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE acc1d_lake

   !------
   SUBROUTINE acc2d_lake (var, s)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:,:)
   real(r8), intent(inout) :: s  (:,:)
   ! Local variables
   integer :: i1, i2

      DO i2 = lbound(var,2), ubound(var,2)
         DO i1 = lbound(var,1), ubound(var,1)
            IF (var(i1,i2) /= spval) THEN
               IF (s(i1,i2) /= spval) THEN
                  s(i1,i2) = s(i1,i2) + var(i1,i2)
               ELSE
                  s(i1,i2) = var(i1,i2)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE acc2d_lake


END MODULE MOD_Lake_1DAccVars

