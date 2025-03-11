#include <define.h>

MODULE MOD_Lake_Namelist
! ---------------------------------- code history -------------------------------------
! Description:
!     Define the namelist and history variables for the lake module.
!
! Type:
!     lake_namelist_type : Define the namelist variables for the lake module
!     lake_histvars_type : Define the history variables for the lake module
!
! Subroutines:
!     read_lake_namelist : read lake namelist variables
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
!================================================================================
   USE MOD_Precision
   IMPLICIT NONE
!================================================================================
   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: read_lake_namelist
   ! PRIVATE :: sync_hist_vars_one

   TYPE lake_namelist_type

      !- Lake Model Options
      ! [1] : CoLM-Lake(CoLML)
      ! [2] : FLake
      ! [3] : Simstrat
      ! [4] : XOML (not supported yet)
      integer :: DEF_LAKE_MODEL = 1          

      !- Calculate the surface flux using a unified surface flux layer
      ![.true.]  : Calculate surface fluxes using a unified surface layer
      ![.false.] : Calculate surface fluxes in each lake model
      logical :: DEF_USE_SURFACE_LAYER = .false.

      !- Lake Surface Layer Flux Scheme Options
      ! [1] : CoLM-Lake(CoLML)
      ! [2] : FLake
      ! [3] : Simstrat
      integer :: DEF_LAKE_SURFACE_FLUX_SCHEME = 1

      !- Lake Surface Flux Scheme Options
      ! [.true.]  : Use CoLM-Lake Surface flux scheme
      ! [.false.] : Use Original Surface flux scheme in each lake model
      logical :: DEF_USE_COLML_FLUX_SCHEME = .true.  
      
      !- Lake Layer Scheme Options
      ! [1] : CoLML layer scheme
      ! [2] : CSSPL layer scheme
      ! [3] : equal layer scheme
      ! [4] : MIN XU, log layer scheme
      ! [5] : WMEJ, log layer scheme
      integer :: DEF_LAKE_LAYER_SCHEME = 1           
      
      !- Lake Surface Roughness Scheme Options
      ! [0] : Default scheme
      ! [1] : Externally constant
      integer :: DEF_LAKE_ROUGHNESSS_SCHEME = 0      
      real(r8):: DEF_LAKE_Z0M = 0.001               ! roughness length over ground, momentum, only used when DEF_LAKE_ROUGHNESSS_SCHEME = 1
      real(r8):: DEF_LAKE_Z0H = 0.0001              ! roughness length over ground,sensible heat, only used when DEF_LAKE_ROUGHNESSS_SCHEME = 1
      real(r8):: DEF_LAKE_Z0Q = 0.0001              ! roughness length over ground, latent heat, only used when DEF_LAKE_ROUGHNESSS_SCHEME = 1
      
      !- Lake Betaprime Scheme Options
      ! [0] : Default scheme
      ! [1] : Externally constant
      integer :: DEF_LAKE_BETAPRIME_SCHEME = 0      
      real(r8):: DEF_LAKE_BETAPRIME = 0.3           ! betaprime in MOST, only used when DEF_LAKE_BETAPRIME_SCHEME = 1
      
      !- Lake light extinction coefficient Scheme Options
      ! [0] : Default scheme
      ! [1] : Externally constant
      integer :: DEF_LAKE_ETA_SCHEME = 0           
      real(r8):: DEF_LAKE_ETA = 0.5                 ! light extinction coefficient, only used when DEF_LAKE_ETA_SCHEME = 1
      
      !- Lake fetch length Scheme Options
      ! [0] : Default scheme
      ! [1] : Externally constant
      integer :: DEF_LAKE_FETCH_SCHEME = 0         
      real(r8):: DEF_LAKE_FETCH = 1000.0            ! lake fetch length, only used when DEF_LAKE_FETCH_SCHEME = 1

      !- Lake mxing enhencement factor
      ! [1.0] : Default scheme, no enhancement
      ! [>1.0]: Enhancement lake mixing
      ! [<1.0]: Suppression lake mixing
      ! [<0.0]: Error
      ! Simstrat: Control the α_Seiche (Fraction of wind energy which goes into seiche energy [-])
      ! CoLM-Lake: Control the mixfact (Mixing enhancement factor [-])
      ! FLake: Control the c_relax_C (Constant in the relaxation equation for the shape factor)
      real(r8):: DEF_LAKE_GAMMA = 1.0  ! lake mixing enhancement factor

   END TYPE lake_namelist_type


   TYPE lake_histvars_type

      logical :: dplak        = .true.
      logical :: zlake        = .true.
      logical :: zilak        = .true.
      logical :: dzlake        = .true.
      logical :: ziarea       = .true.
      logical :: lkz0m        = .true.
      logical :: lkz0h        = .true.
      logical :: lkz0q        = .true.
      logical :: felak        = .true.
      logical :: gamma        = .true.
      logical :: etal         = .true.
      logical :: btpri        = .true.
      logical :: frlak        = .true.
      logical :: tmsno        = .true.
      logical :: tmice        = .true.
      logical :: tmmnw        = .true.
      logical :: tmwml        = .true.
      logical :: tmbot        = .true.
      logical :: tmups        = .true.
      logical :: mldp         = .true.
      logical :: upsdp        = .true.
      logical :: icedp        = .true.
      logical :: bicedp       = .true.
      logical :: wicedp       = .true.
      logical :: CTfrac       = .true.
      logical :: rhosnw       = .true.
      logical :: uwatv        = .true.
      logical :: vwatv        = .true.
      logical :: lksal        = .true.
      logical :: tke          = .true.
      logical :: etke         = .true.
      logical :: eps          = .true.
      logical :: num          = .true.
      logical :: nuh          = .true.
      logical :: lkrho        = .true.

   END TYPE lake_histvars_type

   type (lake_namelist_type) :: DEF_External_Lake
   type (lake_histvars_type) :: DEF_External_Lake_hist

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE read_lake_namelist (nlfile)

   USE MOD_SPMD_Task
   USE MOD_Lake_Const, only: COLML, FLAKE, SIMSTRAT
   IMPLICIT NONE

   character(len=*), intent(in) :: nlfile

   ! Local variables
   logical :: fexists, set_defaults
   integer :: ivar
   integer :: ierr
   character(len=256) :: lakedir

   namelist /External_Lake/ DEF_External_Lake, DEF_External_Lake_hist

      ! ----- open the namelist file -----
      IF (p_is_master) THEN
         open(10, status='OLD', file=nlfile, form="FORMATTED")
         read(10, nml=External_Lake, iostat=ierr)
         IF (ierr /= 0) THEN
            CALL CoLM_Stop (' ***** ERROR: Problem reading namelist: '// trim(nlfile))
         ENDIF
         close(10)
      ENDIF

   IF (DEF_External_Lake%DEF_USE_SURFACE_LAYER) THEN
      write(*,*) '                  *****                  '
      SELECT CASE (DEF_External_Lake%DEF_LAKE_SURFACE_FLUX_SCHEME)
      CASE (COLML)
         write(*,*) '[External Lake] Using CoLM-Lake Surface Flux Scheme'
      CASE (FLAKE)
         write(*,*) '[External Lake] Using FLake Surface Flux Scheme'
      CASE (SIMSTRAT)
         write(*,*) '[External Lake] Using Simstrat Surface Flux Scheme'
      CASE DEFAULT
         write(*,*) '[External Lake] Wrong lake surface flux scheme option, Stop!'
         CALL CoLM_stop ()
      END SELECT
      write(*,*) '                  *****                  '
   ENDIF


#ifdef USEMPI
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_MODEL,           1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_USE_SURFACE_LAYER,           1, mpi_logical , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_SURFACE_FLUX_SCHEME,    1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_USE_COLML_FLUX_SCHEME,       1, mpi_logical , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_LAYER_SCHEME,           1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_ROUGHNESSS_SCHEME,      1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_Z0M,                    1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_Z0H,                    1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_Z0Q,                    1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_BETAPRIME_SCHEME,       1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_BETAPRIME,              1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_ETA_SCHEME,             1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_ETA,                    1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_FETCH_SCHEME,           1, mpi_integer , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_FETCH,                  1, mpi_real8   , p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_External_Lake%DEF_LAKE_GAMMA,                  1, mpi_real8   , p_address_master, p_comm_glb, p_err)
#endif

      CALL sync_hist_lakevars (set_defaults = .true.) 

   END SUBROUTINE read_lake_namelist



   SUBROUTINE sync_hist_lakevars (set_defaults)
   USE MOD_Namelist, only: sync_hist_vars_one
   IMPLICIT NONE

   logical, intent(in) :: set_defaults

   CALL sync_hist_vars_one (DEF_External_Lake_hist%dplak       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%zlake       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%zilak       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%dzlake       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%ziarea      ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%lkz0m       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%lkz0h       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%lkz0q       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%felak       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%gamma       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%etal        ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%btpri       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%frlak       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmsno       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmice       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmmnw       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmwml       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmbot       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tmups       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%mldp        ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%upsdp       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%icedp       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%bicedp      ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%wicedp      ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%rhosnw      ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%uwatv       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%vwatv       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%lksal       ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%tke         ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%etke        ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%eps         ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%num         ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%nuh         ,  set_defaults)
   CALL sync_hist_vars_one (DEF_External_Lake_hist%lkrho       ,  set_defaults)

   END SUBROUTINE sync_hist_lakevars


END MODULE MOD_Lake_Namelist