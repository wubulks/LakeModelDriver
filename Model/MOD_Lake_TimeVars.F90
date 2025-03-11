#include <define.h>

MODULE MOD_Lake_TimeVars
! ---------------------------------- code history -------------------------------------
! Description:
!     Define the lake time variables and subroutines for the lake module.
!
! Subroutines:
!     allocate_LakeTimeVars  : Allocate memory for lake variables.
!     deallocate_LakeTimeVars: Release memory for lake variables.
!     WRITE_LakeTimeVars : Write lake variables to restart file.
!     READ_LakeTimeVars  : Read lake variables from restart file.
!     CHECK_LakeTimeVars : Check lake variables.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
!================================================================================
   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
!================================================================================

   real(r8), allocatable :: dplak            (:)  !lake depth [m]
   real(r8), allocatable :: zlake          (:,:)  !Lake layer node depth [m]
   real(r8), allocatable :: zilak          (:,:)  !Lake layer interface depth [m]
   real(r8), allocatable :: dzlake          (:,:)  !Lake layer thickness [m]
   real(r8), allocatable :: ziarea         (:,:)  !Lake layerinterface area [m2], only for Simstrat
   real(r8), allocatable :: lkz0m            (:)  !Roughness length for momentum [m]
   real(r8), allocatable :: lkz0h            (:)  !Roughness length for sensible heat  [m]
   real(r8), allocatable :: lkz0q            (:)  !Roughness length for latent heat [m]
   real(r8), allocatable :: felak            (:)  !Lake fetch length [m]
   real(r8), allocatable :: gamma            (:)  !Mixing enhancement factor, the meaning is different for each model [-]
   real(r8), allocatable :: etal             (:)  !Lake extinction coefficient [1/m]
   real(r8), allocatable :: btpri            (:)  !Beta prime in Monin-Obukhov theory [-]
   real(r8), allocatable :: frlak            (:)  !Lake fraction [-] 
   real(r8), allocatable :: tmsno            (:)  !Mean snow temperature, only for FLake [K], only for FLake
   real(r8), allocatable :: tmice            (:)  !Mean ice temperature, only for FLake [K], only for FLake
   real(r8), allocatable :: tmmnw            (:)  !Mean temperature of the water column [K]
   real(r8), allocatable :: tmwml            (:)  !Mixed-layer temperature [K], only for FLake
   real(r8), allocatable :: tmbot            (:)  !Temperature at the water-bottom sediment interface [K], only for FLake
   real(r8), allocatable :: tmups            (:)  !Temperature at the bottom of the upper layer of the sediments [K], only for FLake
   real(r8), allocatable :: mldp             (:)  !Mixed layer depth [m], only for FLake
   real(r8), allocatable :: upsdp            (:)  !Bottom of the upper layer of the sediments [m], only for FLake
   real(r8), allocatable :: icedp            (:)  !Mean temperature of the lake [K], for FLake and Simstrat
   real(r8), allocatable :: bicedp           (:)  !black ice depth [m], only for Simstrat
   real(r8), allocatable :: wicedp           (:)  !white ice depth [m], only for Simstrat
   real(r8), allocatable :: CTfrac           (:)  !Shape factor (thermocline)
   real(r8), allocatable :: rhosnw           (:)  !snow density [kg/m3], only for Simstrat
   real(r8), allocatable :: uwatv          (:,:)  !Water velocity in x-direction [m/s], only for Simstrat
   real(r8), allocatable :: vwatv          (:,:)  !Water velocity in y-direction [m/s], only for Simstrat
   real(r8), allocatable :: lksal          (:,:)  !Salinity [‰], only for Simstrat
   real(r8), allocatable :: tke            (:,:)  !Turbulent kinetic energy (TKE) [J/kg], only for Simstrat
   real(r8), allocatable :: etke             (:)  !Seiche energy [J], only for Simstrat
   real(r8), allocatable :: eps            (:,:)  !TKE dissipation rate [W/kg], only for Simstrat
   real(r8), allocatable :: num            (:,:)  !Turbulent viscosity (momentum) [m2/s], only for Simstrat
   real(r8), allocatable :: nuh            (:,:)  !Turbulent diffusivity (heat) [m2/s], only for Simstrat
   real(r8), allocatable :: lkrho          (:,:)  !Density of water [kg/m3]


!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_LakeTimeVars ()
! ---------------------------------- code history -------------------------------------
! Description:
!     Allocate memory for lake variables. Will be called in the subroutine allocate_acc_fluxes.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE
!==============================================================================
!================================================================================
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (dplak            (numpatch)); dplak     (:) = spval 
            allocate (zlake    (nl_lake,numpatch)); zlake   (:,:) = spval 
            allocate (zilak  (nl_lake+1,numpatch)); zilak   (:,:) = spval 
            allocate (dzlake    (nl_lake,numpatch)); dzlake   (:,:) = spval
            allocate (ziarea (nl_lake+1,numpatch)); ziarea  (:,:) = spval 
            allocate (lkz0m            (numpatch)); lkz0m     (:) = spval 
            allocate (lkz0h            (numpatch)); lkz0h     (:) = spval 
            allocate (lkz0q            (numpatch)); lkz0q     (:) = spval 
            allocate (felak            (numpatch)); felak     (:) = spval 
            allocate (gamma            (numpatch)); gamma     (:) = spval 
            allocate (etal             (numpatch)); etal      (:) = spval 
            allocate (btpri            (numpatch)); btpri     (:) = spval 
            allocate (frlak            (numpatch)); frlak     (:) = spval 
            allocate (tmsno            (numpatch)); tmsno     (:) = spval 
            allocate (tmice            (numpatch)); tmice     (:) = spval 
            allocate (tmmnw            (numpatch)); tmmnw     (:) = spval 
            allocate (tmwml            (numpatch)); tmwml     (:) = spval 
            allocate (tmbot            (numpatch)); tmbot     (:) = spval 
            allocate (tmups            (numpatch)); tmups     (:) = spval 
            allocate (mldp             (numpatch)); mldp      (:) = spval 
            allocate (upsdp            (numpatch)); upsdp     (:) = spval 
            allocate (icedp            (numpatch)); icedp     (:) = spval 
            allocate (bicedp           (numpatch)); bicedp    (:) = spval 
            allocate (wicedp           (numpatch)); wicedp    (:) = spval 
            allocate (CTfrac           (numpatch)); CTfrac    (:) = spval
            allocate (rhosnw           (numpatch)); rhosnw    (:) = spval 
            allocate (uwatv    (nl_lake,numpatch)); uwatv   (:,:) = spval 
            allocate (vwatv    (nl_lake,numpatch)); vwatv   (:,:) = spval 
            allocate (lksal    (nl_lake,numpatch)); lksal   (:,:) = spval 
            allocate (tke    (nl_lake+1,numpatch)); tke     (:,:) = spval 
            allocate (etke             (numpatch)); etke      (:) = spval 
            allocate (eps    (nl_lake+1,numpatch)); eps     (:,:) = spval 
            allocate (num    (nl_lake+1,numpatch)); num     (:,:) = spval 
            allocate (nuh    (nl_lake+1,numpatch)); nuh     (:,:) = spval 
            allocate (lkrho    (nl_lake,numpatch)); lkrho   (:,:) = spval

         ENDIF
      ENDIF

   END SUBROUTINE allocate_LakeTimeVars



   SUBROUTINE deallocate_LakeTimeVars()
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
!================================================================================
     
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate ( dplak  )  
            deallocate ( zlake  )  
            deallocate ( zilak  )  
            deallocate ( dzlake ) 
            deallocate ( ziarea )  
            deallocate ( lkz0m  )  
            deallocate ( lkz0h  )  
            deallocate ( lkz0q  )  
            deallocate ( felak  )  
            deallocate ( gamma  )  
            deallocate ( etal   )  
            deallocate ( btpri  )  
            deallocate ( frlak  )  
            deallocate ( tmsno  )  
            deallocate ( tmice  )  
            deallocate ( tmmnw  )  
            deallocate ( tmwml  )  
            deallocate ( tmbot  )  
            deallocate ( tmups  )  
            deallocate ( mldp   )  
            deallocate ( upsdp  )  
            deallocate ( icedp  )  
            deallocate ( bicedp )  
            deallocate ( wicedp )
            deallocate ( CTfrac )  
            deallocate ( rhosnw )  
            deallocate ( uwatv  )  
            deallocate ( vwatv  )  
            deallocate ( lksal  )  
            deallocate ( tke    )  
            deallocate ( etke   )  
            deallocate ( eps    )  
            deallocate ( num    )  
            deallocate ( nuh    )  
            deallocate ( lkrho  )  
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_LakeTimeVars



   SUBROUTINE WRITE_LakeTimeVars(idate, lc_year, site, dir_restart)
! ---------------------------------- code history -------------------------------------
! Description:
!     Write lake variables to restart file. Will be called in the subroutine WRITE_TimeVariables.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFVector
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   IMPLICIT NONE
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: idate(3)
   integer, intent(in) :: lc_year      !year of land cover type data
   character(len=*), intent(in) :: site
   character(len=*), intent(in) :: dir_restart

!  ------------------------- Local variables ---------------------------
   character(len=256) :: file_restart
   character(len=14)  :: cdate
   character(len=256) :: cyear         !character for lc_year
   integer :: compress
!==============================================================================
      compress = DEF_REST_CompressLevel

      ! land cover type year
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      ! Define dimensions
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lakeI',    nl_lake+1     )

      ! Write variables
      CALL ncio_write_vector (file_restart, 'dplak      ', 'patch', landpatch, dplak, compress) 
      CALL ncio_write_vector (file_restart, 'zlake      ', 'lake' , nl_lake  , 'patch'  , landpatch, zlake     , compress) 
      CALL ncio_write_vector (file_restart, 'zilak      ', 'lakeI', nl_lake+1, 'patch'  , landpatch, zilak     , compress) 
      CALL ncio_write_vector (file_restart, 'dzlake      ', 'lake' , nl_lake  , 'patch'  , landpatch, dzlake     , compress)
      CALL ncio_write_vector (file_restart, 'ziarea     ', 'lakeI', nl_lake+1, 'patch'  , landpatch, ziarea    , compress) 
      CALL ncio_write_vector (file_restart, 'lkz0m      ', 'patch', landpatch, lkz0m    , compress) 
      CALL ncio_write_vector (file_restart, 'lkz0h      ', 'patch', landpatch, lkz0h    , compress) 
      CALL ncio_write_vector (file_restart, 'lkz0q      ', 'patch', landpatch, lkz0q    , compress) 
      CALL ncio_write_vector (file_restart, 'felak      ', 'patch', landpatch, felak    , compress) 
      CALL ncio_write_vector (file_restart, 'gamma      ', 'patch', landpatch, gamma    , compress) 
      CALL ncio_write_vector (file_restart, 'etal       ', 'patch', landpatch, etal     , compress) 
      CALL ncio_write_vector (file_restart, 'btpri      ', 'patch', landpatch, btpri    , compress) 
      CALL ncio_write_vector (file_restart, 'frlak      ', 'patch', landpatch, frlak    , compress) 
      CALL ncio_write_vector (file_restart, 'tmsno      ', 'patch', landpatch, tmsno    , compress) 
      CALL ncio_write_vector (file_restart, 'tmice      ', 'patch', landpatch, tmice    , compress) 
      CALL ncio_write_vector (file_restart, 'tmmnw      ', 'patch', landpatch, tmmnw    , compress) 
      CALL ncio_write_vector (file_restart, 'tmwml      ', 'patch', landpatch, tmwml    , compress) 
      CALL ncio_write_vector (file_restart, 'tmbot      ', 'patch', landpatch, tmbot    , compress) 
      CALL ncio_write_vector (file_restart, 'tmups      ', 'patch', landpatch, tmups    , compress) 
      CALL ncio_write_vector (file_restart, 'mldp       ', 'patch', landpatch, mldp     , compress) 
      CALL ncio_write_vector (file_restart, 'upsdp      ', 'patch', landpatch, upsdp    , compress) 
      CALL ncio_write_vector (file_restart, 'icedp      ', 'patch', landpatch, icedp    , compress) 
      CALL ncio_write_vector (file_restart, 'bicedp     ', 'patch', landpatch, bicedp   , compress) 
      CALL ncio_write_vector (file_restart, 'wicedp     ', 'patch', landpatch, wicedp   , compress)
      CALL ncio_write_vector (file_restart, 'CTfrac     ', 'patch', landpatch, CTfrac   , compress) 
      CALL ncio_write_vector (file_restart, 'rhosnw     ', 'patch', landpatch, rhosnw   , compress) 
      CALL ncio_write_vector (file_restart, 'uwatv      ', 'lake' , nl_lake  , 'patch'  , landpatch, uwatv     , compress) 
      CALL ncio_write_vector (file_restart, 'vwatv      ', 'lake' , nl_lake  , 'patch'  , landpatch, vwatv     , compress) 
      CALL ncio_write_vector (file_restart, 'lksal      ', 'lake' , nl_lake  , 'patch'  , landpatch, lksal     , compress) 
      CALL ncio_write_vector (file_restart, 'tke        ', 'lakeI', nl_lake+1, 'patch'  , landpatch, tke       , compress) 
      CALL ncio_write_vector (file_restart, 'etke       ', 'patch', landpatch, etke     , compress) 
      CALL ncio_write_vector (file_restart, 'eps        ', 'lakeI', nl_lake+1, 'patch'  , landpatch, eps       , compress) 
      CALL ncio_write_vector (file_restart, 'num        ', 'lakeI', nl_lake+1, 'patch'  , landpatch, num       , compress) 
      CALL ncio_write_vector (file_restart, 'nuh        ', 'lakeI', nl_lake+1, 'patch'  , landpatch, nuh       , compress) 
      CALL ncio_write_vector (file_restart, 'lkrho      ', 'lake' , nl_lake  , 'patch'  , landpatch, lkrho     , compress) 

   END SUBROUTINE WRITE_LakeTimeVars



   SUBROUTINE READ_LakeTimeVars(idate, lc_year, site, dir_restart)
! ---------------------------------- code history -------------------------------------
! Description:
!     Read lake variables. Will be called in the subroutine READ_TimeVariables.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_NetCDFVector
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   IMPLICIT NONE
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: idate(3)
   integer, intent(in) :: lc_year      !year of land cover type data
   character(len=*), intent(in) :: site
   character(len=*), intent(in) :: dir_restart
!  ------------------------- Local variables ---------------------------
   character(len=256) :: file_restart
   character(len=14)  :: cdate, cyear
!==============================================================================

      ! land cover type year
      write(cyear,'(i4.4)') lc_year

      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
      
      CALL ncio_read_vector (file_restart, 'dplak    '  ,            landpatch, dplak     )   !lake depth
      CALL ncio_read_vector (file_restart, 'zlake    '  , nl_lake  , landpatch, zlake     )   !Lake layer node depth
      CALL ncio_read_vector (file_restart, 'zilak    '  , nl_lake+1, landpatch, zilak     )   !Lake layer interface depth
      CALL ncio_read_vector (file_restart, 'dzlake    '  , nl_lake  , landpatch, dzlake     )   !Lake layer thickness
      CALL ncio_read_vector (file_restart, 'ziarea   '  , nl_lake+1, landpatch, ziarea    )   !Lake layerinterface area [m2], only for Simstrat
      CALL ncio_read_vector (file_restart, 'lkz0m    '  ,            landpatch, lkz0m     )   !Roughness length for momentum [m]
      CALL ncio_read_vector (file_restart, 'lkz0h    '  ,            landpatch, lkz0h     )   !Roughness length for sensible heat  [m]
      CALL ncio_read_vector (file_restart, 'lkz0q    '  ,            landpatch, lkz0q     )   !Roughness length for latent heat [m]
      CALL ncio_read_vector (file_restart, 'felak    '  ,            landpatch, felak     )   !Lake fetch length
      CALL ncio_read_vector (file_restart, 'gamma    '  ,            landpatch, gamma     )   !Mixing enhancement factor, the meaning is different for each model [-]
      CALL ncio_read_vector (file_restart, 'etal     '  ,            landpatch, etal      )   !Lake extinction coefficient [1/m]
      CALL ncio_read_vector (file_restart, 'btpri    '  ,            landpatch, btpri     )   !Beta prime in Monin-Obukhov theory
      CALL ncio_read_vector (file_restart, 'frlak    '  ,            landpatch, frlak     )   !Lake fraction [-] 
      CALL ncio_read_vector (file_restart, 'tmsno    '  ,            landpatch, tmsno     )   !Mean snow temperature, only for FLake [K], only for FLake
      CALL ncio_read_vector (file_restart, 'tmice    '  ,            landpatch, tmice     )   !Mean ice temperature, only for FLake [K], only for FLake
      CALL ncio_read_vector (file_restart, 'tmmnw    '  ,            landpatch, tmmnw     )   !Mean temperature of the water column [K]
      CALL ncio_read_vector (file_restart, 'tmwml    '  ,            landpatch, tmwml     )   !Mixed-layer temperature [K], only for FLake
      CALL ncio_read_vector (file_restart, 'tmbot    '  ,            landpatch, tmbot     )   !Temperature at the water-bottom sediment interface [K], only for FLake
      CALL ncio_read_vector (file_restart, 'tmups    '  ,            landpatch, tmups     )   !Temperature at the bottom of the upper layer of the sediments [K], only for FLake
      CALL ncio_read_vector (file_restart, 'mldp     '  ,            landpatch, mldp      )   !Mixed layer depth [m], only for FLake
      CALL ncio_read_vector (file_restart, 'upsdp    '  ,            landpatch, upsdp     )   !Bottom of the upper layer of the sediments [m], only for FLake
      CALL ncio_read_vector (file_restart, 'icedp    '  ,            landpatch, icedp     )   !Mean temperature of the lake [K], for FLake and Simstrat
      CALL ncio_read_vector (file_restart, 'bicedp   '  ,            landpatch, bicedp    )   !black ice depth (m), only for Simstrat
      CALL ncio_read_vector (file_restart, 'wicedp   '  ,            landpatch, wicedp    )   !white ice depth (m), only for Simstrat
      CALL ncio_read_vector (file_restart, 'CTfrac   '  ,            landpatch, CTfrac    )   !Shape factor (thermocline)
      CALL ncio_read_vector (file_restart, 'rhosnw   '  ,            landpatch, rhosnw    )   !snow density (kg/m3), only for Simstrat
      CALL ncio_read_vector (file_restart, 'uwatv    '  , nl_lake  , landpatch, uwatv     )   !Water velocity in x-direction [m/s], only for Simstrat
      CALL ncio_read_vector (file_restart, 'vwatv    '  , nl_lake  , landpatch, vwatv     )   !Water velocity in y-direction [m/s], only for Simstrat
      CALL ncio_read_vector (file_restart, 'lksal    '  , nl_lake  , landpatch, lksal     )   !Salinity [‰], only for Simstrat
      CALL ncio_read_vector (file_restart, 'tke      '  , nl_lake+1, landpatch, tke       )   !Turbulent kinetic energy (TKE) [J/kg], only for Simstrat
      CALL ncio_read_vector (file_restart, 'etke     '  ,            landpatch, etke      )   !Seiche energy [J], only for Simstrat
      CALL ncio_read_vector (file_restart, 'eps      '  , nl_lake+1, landpatch, eps       )   !TKE dissipation rate [W/kg], only for Simstrat
      CALL ncio_read_vector (file_restart, 'num      '  , nl_lake+1, landpatch, num       )   !Turbulent viscosity (momentum), only for Simstrat
      CALL ncio_read_vector (file_restart, 'nuh      '  , nl_lake+1, landpatch, nuh       )   !Turbulent diffusivity (heat), only for Simstrat
      CALL ncio_read_vector (file_restart, 'lkrho   '  , nl_lake  , landpatch, lkrho    )   !Density of water [kg/m3]

   END SUBROUTINE READ_LakeTimeVars


#ifdef RangeCheck
   SUBROUTINE CHECK_LakeTimeVars()
! ---------------------------------- code history -------------------------------------
! Description:
!     Check lake variables. Will be called in the subroutine check_TimeVariables.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   IMPLICIT NONE
!==============================================================================

      CALL check_vector_data ('dplak       [m]    ', dplak      ) ! lake depth
      CALL check_vector_data ('zlake       [m]    ', zlake      ) ! Lake layer node depth
      CALL check_vector_data ('zilak       [m]    ', zilak      ) ! Lake layer interface depth
      CALL check_vector_data ('dzlake      [m]    ', dzlake     ) ! Lake layer thickness
      CALL check_vector_data ('ziarea      [m2]   ', ziarea     ) ! Lake layerinterface area [m2], only for Simstrat
      CALL check_vector_data ('lkz0m       [m]    ', lkz0m      ) ! Roughness length for momentum [m]
      CALL check_vector_data ('lkz0h       [m]    ', lkz0h      ) ! Roughness length for sensible heat  [m]
      CALL check_vector_data ('lkz0q       [m]    ', lkz0q      ) ! Roughness length for latent heat [m]
      CALL check_vector_data ('felak       [m]    ', felak      ) ! Lake fetch length
      CALL check_vector_data ('gamma       [-]    ', gamma      ) ! Mixing enhancement factor, the meaning is different for each model [-]
      CALL check_vector_data ('etal        [1/m]  ', etal       ) ! Lake extinction coefficient [1/m]
      CALL check_vector_data ('btpri       [-]    ', btpri      ) ! Beta prime in Monin-Obukhov theory
      CALL check_vector_data ('frlak       [-]    ', frlak      ) ! Lake fraction [-] 
      CALL check_vector_data ('tmsno       [K]    ', tmsno      ) ! Mean snow temperature, only for FLake [K], only for FLake
      CALL check_vector_data ('tmice       [K]    ', tmice      ) ! Mean ice temperature, only for FLake [K], only for FLake
      CALL check_vector_data ('tmmnw       [K]    ', tmmnw      ) ! Mean temperature of the water column [K]
      CALL check_vector_data ('tmwml       [K]    ', tmwml      ) ! Mixed-layer temperature [K], only for FLake
      CALL check_vector_data ('tmbot       [K]    ', tmbot      ) ! Temperature at the water-bottom sediment interface [K], only for FLake
      CALL check_vector_data ('tmups       [K]    ', tmups      ) ! Temperature at the bottom of the upper layer of the sediments [K], only for FLake
      CALL check_vector_data ('mldp        [m]    ', mldp       ) ! Mixed layer depth [m], only for FLake
      CALL check_vector_data ('upsdp       [m]    ', upsdp      ) ! Bottom of the upper layer of the sediments [m], only for FLake
      CALL check_vector_data ('icedp       [m]    ', icedp      ) ! Mean temperature of the lake [K], for FLake and Simstrat
      CALL check_vector_data ('bicedp      [m]    ', bicedp     ) ! black ice depth (m), only for Simstrat
      CALL check_vector_data ('wicedp      [m]    ', wicedp     ) ! white ice depth (m), only for Simstrat
      CALL check_vector_data ('CTfrac      [-]    ', CTfrac     ) ! Shape factor (thermocline)
      CALL check_vector_data ('rhosnw      [kg/m3]', rhosnw     ) ! snow density (kg/m3), only for Simstrat
      CALL check_vector_data ('uwatv       [m/s]  ', uwatv      ) ! Water velocity in x-direction [m/s], only for Simstrat
      CALL check_vector_data ('vwatv       [m/s]  ', vwatv      ) ! Water velocity in y-direction [m/s], only for Simstrat
      CALL check_vector_data ('lksal       [‰]    ', lksal      ) ! Salinity [‰], only for Simstrat
      CALL check_vector_data ('tke         [J/kg] ', tke        ) ! Turbulent kinetic energy (TKE) [J/kg], only for Simstrat
      CALL check_vector_data ('etke        [J]    ', etke       ) ! Seiche energy [J], only for Simstrat
      CALL check_vector_data ('eps         [W/kg] ', eps        ) ! TKE dissipation rate [W/kg], only for Simstrat
      CALL check_vector_data ('num         [m2/s] ', num        ) ! Turbulent viscosity (momentum), only for Simstrat
      CALL check_vector_data ('nuh         [m2/s] ', nuh        ) ! Turbulent diffusivity (heat), only for Simstrat
      CALL check_vector_data ('lkrho       [kg/m3]', lkrho       ) ! Density of water [kg/m3]

   END SUBROUTINE CHECK_LakeTimeVars
#endif


   SUBROUTINE InitLakeTimeVars(ipatch, lakedepth, t_lake, lake_icefrac, savedtke1)
! ---------------------------------- code history -------------------------------------
! Description:
!     Initialize lake variables. Will be called in the subroutine initialize.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_LandPatch, only: numpatch
   USE MOD_Lake_Utils, only: LakeIni
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_Lake_Namelist
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_Lake_Const, only: nlice
   IMPLICIT NONE
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer , intent(in) :: ipatch
   real(r8), intent(in) :: lakedepth
   real(r8), intent(inout) :: t_lake(nl_lake)
   real(r8), intent(inout) :: lake_icefrac(nl_lake)
   real(r8), intent(inout) :: savedtke1
   integer :: i
!==============================================================================
      i = ipatch
      !  Initialize lake variables
      lkz0m(i) = DEF_External_Lake%DEF_LAKE_Z0M 
      lkz0h(i) = DEF_External_Lake%DEF_LAKE_Z0H
      lkz0q(i) = DEF_External_Lake%DEF_LAKE_Z0Q
      felak(i) = DEF_External_Lake%DEF_LAKE_FETCH
      btpri(i) = DEF_External_Lake%DEF_LAKE_BETAPRIME
      gamma(i) = DEF_External_Lake%DEF_LAKE_GAMMA
      etal(i) = DEF_External_Lake%DEF_LAKE_ETA
      dplak(i) = lakedepth
      if (dplak(i) < 1.0) dplak(i) = 1.0
      CALL LakeIni( &
                  ! "in" arguments
                  ! -------------------
                  nlake = nl_lake     ,  nsnow = -maxsnl      ,  nsoil = nl_soil           ,  nlice = nlice        ,& 
                  dplak = dplak(i)    ,  tskin = 285.         ,&
                  ! "out" arguments
                  ! -------------------
                  zlake = zlake(:,i)  ,  zilak = zilak(:,i)   ,  dzlake = dzlake(:,i)      ,  lktmp = t_lake       ,& 
                  rhosnw = rhosnw(i)  ,  lkrho = lkrho(:,i)   ,  icefr = lake_icefrac      ,  stke1 = savedtke1    ,&
                  tmsno = tmsno(i)    ,  tmice = tmice(i)     ,  tmmnw = tmmnw(i)          ,  tmwml = tmwml(i)     ,&
                  tmbot = tmbot(i)    ,  tmups = tmups(i)     ,  mldp = mldp(i)            ,  upsdp = upsdp(i)     ,&
                  CTfrac = CTfrac(i)  ,  icedp = icedp(i)     ,  bicedp = bicedp(i)        ,  wicedp = wicedp(i)   ,&
                  uwatv = uwatv(:,i)  ,  vwatv = vwatv(:,i)   ,  lksal = lksal(:,i)        ,  tke = tke(:,i)       ,&
                  eps = eps(:,i)      ,  etke = etke(i)       ,  num = num(:,i)            ,  nuh = nuh(:,i)       ,&
                  ziarea = ziarea(:,i))

   END SUBROUTINE InitLakeTimeVars


END MODULE MOD_Lake_TimeVars


