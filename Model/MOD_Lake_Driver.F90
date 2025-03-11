#include <define.h>

MODULE  MOD_Lake_Driver

! ================================== code history =====================================
! Description:
!     This is the new version of the CWRF lake simulation driver, which includes a total of 7 lake models.
!     The CALL of the lake model is controlled by <lkopt>:
!         (1) CoLML, Yongjiu Dai et al.       
!         (2) FLake, Dmitrii Mironov et al.   
!         (3) Simstrat, Goudsmit et al.       
!         (4) XOML, Tiejun Ling et al.        
!
!     -*- list of key output variables from Lake Driver -------------------------------------------------
!     ---------------------------------------------------------------------------------------------------      
!              Lake model:    COLML     FLake     Simstrat      XOML      note
!                   tmsno:        x        io           io         x
!                   tmice:        d        io           io         x
!                   lktmp:       io         o           io        io
!                   tmmnw:        x        io            x         x
!                   tmwml:        x        io            x         x                  
!                   tmbot:        x        io            x         x
!                   tmups:        x        io            x         x
!                   t_ssb:       io        io          sio         x      
!                   xwliq:       io         x            x         x
!                   xwice:       io         x            x         x                                                 
!                   icefr:       io         o            o         x                
!                   stke1:       io         i            i         i              
!                   icedp:       io        io           io         x                           
!                  bicedp:        x         x           io         x
!                  wicedp:        x         x           io         x
!                   snwdp:       io        io           io         x
!                   snwml:        o         o            o         x
!                   uwatv:        x         x           io        io     
!                   vwatv:        x         x           io        io
!                   lksal:        x         x           io        io
!                     tke:        x         x           io        io
!                    etke:       io         x           io         x        
!                     eps:        x         x           io        io
!                     num:        x         x           io         x
!                     nuh:        x         x           io         x
!                      Km:        x         x            x        io
!                      Kh:        x         x            x        io
!                      Ke:        x         x            x        io
!                   lkrho:        x         x           io        io
!                  rhosnw:        x         x           io         x
!                    mldp:        x        io            x         x
!                   upsdp:        x        io            x         x
!                  CTfrac:        x        io            x         x
!                   -----:   o=output    d=diagnostic output
!                            i=input     x=not set nor used (=spv_water)
!                            s=only snow layers otherwise [snow + sediment] layers
!                            r=required when restarting
!     ---------------------------------------------------------------------------------------------------      
!     -*- list of key output variables from Lake Driver -------------------------------------------------
!
!-lxz  TODO: make the whole lakectl to consider ifdef FrcICE
!-lxz  TODO: combine with seaice using xice to link with icefr(1)
!-WMEJ TODO: check all variable comments.
!-WMEJ TODO: organize variables in the order of input, output, and inout.
!-WMEJ TODO: add the necessary variables.
!-WMEJ TODO: review the variable names and comments.
! * WMEJICE: Options for calculating the ice fraction in the lake model. [OFF by default]
!            IF defined, the mass fraction will be used to calculate the ice fraction.
!
! Original author: 
!     Yongjiu Dai, 2000
!
! Revisions:
!     Min Xu, 02/2012
!     Xin-Zhong Liang,
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 
!
! =====================================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Lake_Driver
   PUBLIC :: Lake_Interface

   ! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: LakeCtl

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------

   SUBROUTINE Lake_Interface( &
         ! "in" arguments    
         ! -------------------
         dtlak     , rlat      , rlon      , bifall    ,&
         hwref     , htref     , hqref     , usurf     ,&
         vsurf     , tmref     , qmref     , arhos     ,&
         psurf     , lwdns     , sabg      , hpbl      ,& 
         sols      , soll      , solsd     , solld     ,&
         crain     , lrain     , csnow     , lsnow     ,&
         tprec     , ipatch    ,&
         ! "inout" arguments
         ! -------------------
         tskin     , lktmp     , t_ssb     , snlay     ,&
         zssb      , zissb     , dzssb     , snwcv     ,&
         stke1     , snwag     , snwdp     , icefr     ,&
         xwliq     , xwice     , snout     ,&
! SNICAR model variables
         forc_aer  , sabg_lyr  , snofrz    ,&
         mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
         mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
         ! "out" arguments
         ! -------------------
         fsena     , fevpa     , lfevpa    , fseng     ,&
         fevpg     , lwups     , fgrnd     , trad      ,&
         qseva     , qsubl     , qsdew     , qfros     ,&
         taux      , tauy      , ustar     , qstar     ,&
         tstar     , lkems     , snwml     , zol       ,&
         t2m       , q2m       , fm        , fq        ,&
         rib       , fh        , z0m       , urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!    Serves as the primary interface between COLMMAIN driver program and external routines,
!    explicitly handling dynamic parameters through argument passing while
!    implicitly accessing other variables via module association.
!    The associated module exists purely for aesthetic coherence in source hierarchy
!    rather than functional necessity - a concession to coding standard formalities
!
! Called: (* means optional)
!    -> Lake_Driver         : Lake Model Driver
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Vars_Global
   USE MOD_Lake_TimeVars
   USE MOD_Lake_Namelist
   USE MOD_Vars_TimeInvariants
   USE MOD_Lake_Const, only: nlice
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                    ! patch index

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      rlat                    ,&! latitude [radians]
      rlon                    ,&! longitude [radians]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sabg                    ,&! solar absorbed by ground, the reflection has been removed [W/m2]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! atmospheric boundary layer height, only USE in CoLML_SurFlux [m]
      tprec                   ,&! snowfall/rainfall temperature [kelvin]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      lktmp(nl_lake)          ,&! lake temperature (K)
      icefr(nl_lake)            ! lake ice fraction [-]

   real(r8), intent(inout)  :: &
      dzssb(maxsnl+1:nl_soil) ,&! snow + soil + bedrock layer thickness [m]
      zssb (maxsnl+1:nl_soil) ,&! snow + soil + bedrock node layer depth [m]
      zissb(maxsnl:nl_soil)   ,&! snow + soil + bedrock layer interface level [m]
      t_ssb(maxsnl+1:nl_soil) ,&! snow + soil + bedrock layer temperature [K]
      xwliq(maxsnl+1:nl_soil) ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice(maxsnl+1:nl_soil)   ! snow + soil + bedrock layer ice lens (kg/m2)
        
   real(r8), intent(inout)  :: &
      snwcv                   ,&! snow mass (kg/m2)
      snwdp                   ,&! snow depth [m]
      snwag                   ,&! non dimensional snow age [-]
      snout                   ,&! rate of water out of snow bottom (mm/s)
      tskin                   ,&! ground surface temperature [k]
      stke1                     ! top level eddy conductivity [W/m/K]

   ! ----- SNICAR variables -----
   real(r8), intent(in)     :: &
      forc_aer(14)              ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

   real(r8), intent(in)     :: &
      sabg_lyr(maxsnl+1:1)      ! solar radiation absorbed by ground [W/m2]
   
   real(r8), intent(inout)  :: &
      mss_bcpho(maxsnl+1:0)   ,&! mass of hydrophobic BC in snow  (col,lyr) [kg]
      mss_bcphi(maxsnl+1:0)   ,&! mass of hydrophillic BC in snow (col,lyr) [kg]
      mss_ocpho(maxsnl+1:0)   ,&! mass of hydrophobic OC in snow  (col,lyr) [kg]
      mss_ocphi(maxsnl+1:0)   ,&! mass of hydrophillic OC in snow (col,lyr) [kg]
      mss_dst1 (maxsnl+1:0)   ,&! mass of dust species 1 in snow  (col,lyr) [kg]
      mss_dst2 (maxsnl+1:0)   ,&! mass of dust species 2 in snow  (col,lyr) [kg]
      mss_dst3 (maxsnl+1:0)   ,&! mass of dust species 3 in snow  (col,lyr) [kg]
      mss_dst4 (maxsnl+1:0)     ! mass of dust species 4 in snow  (col,lyr) [kg]
   
   real(r8), intent(out)    :: &
      snofrz(maxsnl+1:0)        ! snow freezing flag (col,lyr) [0/1]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      lkems                   ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      z0m                     ,&! roughness length over ground, momentum [m]
      trad                      ! radiative temperature [K]
        
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

   integer :: i
!==============================================================================
   i = ipatch

   CALL Lake_Driver( &
         ! "in" arguments    
         ! -------------------
         nl_lake      , -maxsnl         , nl_soil         , nlice      ,& 
         dtlak        , rlat            , rlon            ,&
         hwref        , htref           , hqref           , usurf      ,&
         vsurf        , tmref           , qmref           , arhos      ,&
         psurf        , lwdns           , sabg            , 1000.      ,& 
         crain        , lrain           , csnow           , lsnow      ,&
         hpbl         , vf_quartz(1:,i) , vf_gravels(1:,i), vf_om(1:,i),&
         vf_sand(1:,i), wf_gravels(1:,i), wf_sand(1:,i)   , porsl(1:,i),&
         csol(1:,i)   , k_solids(1:,i)  , dksatf(1:,i)    , dkdry(1:,i),&
         dksatu(1:,i) , BA_alpha(1:,i)  , BA_beta(1:,i)   , tprec      ,&
         sols         , soll            , solsd     , solld     ,&
         ipatch       , bifall          ,&
         ! "inout" arguments
         ! -------------------
         tskin        , lktmp           , t_ssb           ,&
         zlake(:,i)   , zilak(:,i)      , dzlake(:,i)     , dplak(i)   ,&
         zssb         , zissb           , dzssb           ,&
         stke1        , snwcv           , snwag           ,&
         snwdp        , xwliq           , xwice           , lkz0m(i)   ,&
         lkz0h(i)     , lkz0q(i)        , felak(i)        , gamma(i)   ,&
         etal(i)      , btpri(i)        , icefr           , snlay      ,&
         tmsno(i)     , tmice(i)        , tmmnw(i)        , tmwml(i)   ,&
         tmbot(i)     , tmups(i)        , icedp(i)        , mldp(i)    ,&
         upsdp(i)     , CTfrac(i)       , frlak(i)        , ziarea(:,i),&
         uwatv(:,i)   , vwatv(:,i)      , lksal(:,i)      , tke(:,i)   ,&
         eps(:,i)     , etke(i)         , num(:,i)        , nuh(:,i)   ,&
         bicedp(i)    , wicedp(i)       , lkrho(:,i)      , rhosnw(i)  ,&
         snout        ,&
! SNICAR model variables
         forc_aer     , sabg_lyr        , snofrz          ,&
         mss_bcpho    , mss_bcphi       , mss_ocpho       , mss_ocphi  ,&
         mss_dst1     , mss_dst2        , mss_dst3        , mss_dst4   ,&
! END SNICAR model variables
         ! "out" arguments
         ! -------------------
         fsena        , fevpa           , lfevpa          , fseng      ,&
         fevpg        , lwups           , fgrnd           , trad       ,&
         qseva        , qsubl           , qsdew           , qfros      ,&
         taux         , tauy            , ustar           , qstar      ,&
         tstar        , lkems           , snwml           , zol        ,&
         t2m          , q2m             , fm              , fq         ,&
         rib          , fh              )

      z0m = lkz0m(i)

   END SUBROUTINE Lake_Interface



   SUBROUTINE Lake_Driver( &
            ! "in" arguments    
            ! -------------------
            nlake     , nsnow     , nsoil     , nlice     ,& 
            dtlak     , rlat      , rlon      ,&
            hwref     , htref     , hqref     , usurf     ,&
            vsurf     , tmref     , qmref     , arhos     ,&
            psurf     , lwdns     , sabg      , zcbcv     ,& 
            crain     , lrain     , csnow     , lsnow     ,&
            hpbl      , vf_quartz , vf_gravels, vf_om     ,&
            vf_sand   , wf_gravels, wf_sand   , porsl     ,&
            csol      , k_solids  , dksatf    , dkdry     ,&
            dksatu    , BA_alpha  , BA_beta   , tprec     ,&
            sols      , soll      , solsd     , solld     ,&
            ipatch    , bifall    ,&
            ! "inout" arguments
            ! -------------------
            tskin     , lktmp     , t_ssb     ,&
            zlake     , zilak     , dzlake     , dplak     ,&
            zssb      , zissb     , dzssb     ,&
            stke1     , snwcv     , snwag     ,&
            snwdp     , xwliq     , xwice     , z0m       ,&
            z0h       , z0q       , felak     , gamma     ,&
            etal      , btpri     , icefr     , snlay     ,&
            tmsno     , tmice     , tmmnw     , tmwml     ,&
            tmbot     , tmups     , icedp     , mldp      ,&
            upsdp     , CTfrac    , frlak     , ziarea    ,&
            uwatv     , vwatv     , lksal     , tke       ,&
            eps       , etke      , num       , nuh       ,&
            bicedp    , wicedp    , lkrho     , rhosnw    ,&
            snout     ,&
! SNICAR model variables
            forc_aer  , sabg_lyr  , snofrz    ,&
            mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
            mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
            ! "out" arguments
            ! -------------------
            fsena     , fevpa     , lfevpa    , fseng     ,&
            fevpg     , lwups     , fgrnd     , trad      ,&
            qseva     , qsubl     , qsdew     , qfros     ,&
            taux      , tauy      , ustar     , qstar     ,&
            tstar     , lkems     , snwml     , zol       ,&
            t2m       , q2m       , fm        , fq        ,&
            rib       , fh        , urban_call)
                  
! ---------------------------------- code history -------------------------------------
! Description:
!     This SUBROUTINE is responsible for managing the energy exchange between
!     lakes and the atmosphere, and for preparing certain variables for the lake models.
!     All variable processing related to the lakes should be completed here.
!
! Called: (* means optional)
!    -> LakeCtl             : Lake Model Controller
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: rflk_depth_bs_ref, d2r, ShallowDeepPartition, SHALLOW, DEEP, PI
   USE MOD_Lake_Namelist, only: DEF_External_Lake
   USE MOD_Vars_Global, only: spval
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer (positive value)
      nsoil                   ,&! Maximum number of soil layer
      nlice                     ! Maximum number of ice layer

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      rlat                    ,&! latitude [radians]
      rlon                    ,&! longitude [radians]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sabg                    ,&! solar absorbed by ground, the reflection has been removed [W/m2]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! atmospheric boundary layer height, only USE in CoLML_SurFlux [m]
      zcbcv                   ,&! convective boundary height [m]
      tprec                   ,&! snowfall/rainfall temperature [kelvin]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

   real(r8), intent(in)     :: &
      porsl     (nsoil)       ,&! soil porosity [-]
      dksatu    (nsoil)       ,&! thermal conductivity of saturated unfrozen soil [W/m-K]
      k_solids  (nsoil)       ,&! thermal conductivity of mineralssoil [W/m-K]
      dkdry     (nsoil)       ,&! thermal conductivity of dry soil [W/m-K]
      vf_quartz (nsoil)       ,&! volumetric fraction of quartz within mineral soil [-]
      vf_gravels(nsoil)       ,&! volumetric fraction of gravels [-]
      vf_om     (nsoil)       ,&! volumetric fraction of organic matter [-]
      vf_sand   (nsoil)       ,&! volumetric fraction of sand [-]
      wf_gravels(nsoil)       ,&! gravimetric fraction of gravels [-]
      wf_sand   (nsoil)       ,&! gravimetric fraction of sand [-]
      csol      (nsoil)       ,&! heat capacity of soil solids [J/(m3 K)]
      dksatf    (nsoil)       ,&! thermal conductivity of saturated frozen soil [W/m-K]
      BA_alpha  (nsoil)       ,&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nsoil)         ! beta in Balland and Arp(2005) thermal conductivity scheme

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlake(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      ziarea(nlake+1)         ,&! lake layer interface area [m2]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)


   real(r8), intent(inout)  :: &
      dzssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer thickness [m]
      zssb  (-nsnow+1:nsoil)  ,&! snow + soil + bedrock node layer depth [m]
      zissb (-nsnow:nsoil)    ,&! snow + soil + bedrock layer interface level [m]
      t_ssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer temperature [K]
      xwliq (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice (-nsnow+1:nsoil)    ! snow + soil + bedrock layer ice lens (kg/m2)
        
   real(r8), intent(inout)  :: &
      dplak                   ,&! lake depth [m]
      snwcv                   ,&! snow mass (kg/m2)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth [m]
      snout                   ,&! rate of water out of snow bottom (mm/s)
      tmsno                   ,&! snow temperature [K]
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      mldp                    ,&! mixed layer depth [m]
      upsdp                   ,&! bottom of the upper layer of the sediments [m]
      icedp                   ,&! ice depth [m]
      bicedp                  ,&! black snow depth (m)
      wicedp                  ,&! white snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      CTfrac                  ,&! Shape factor (thermocline)
      tskin                   ,&! ground surface temperature [k]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient, 1 means no enhancement [-]
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory [-]
      frlak                     ! lake fraction [-]

   ! ----- SNICAR variables -----
   real(r8), intent(in)     :: &
      forc_aer(14)              ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

   real(r8), intent(in)     :: &
      sabg_lyr(-nsnow+1:1)      ! solar radiation absorbed by ground [W/m2]
   
   real(r8), intent(inout)  :: &
      mss_bcpho(-nsnow+1:0)   ,&! mass of hydrophobic BC in snow  (col,lyr) [kg]
      mss_bcphi(-nsnow+1:0)   ,&! mass of hydrophillic BC in snow (col,lyr) [kg]
      mss_ocpho(-nsnow+1:0)   ,&! mass of hydrophobic OC in snow  (col,lyr) [kg]
      mss_ocphi(-nsnow+1:0)   ,&! mass of hydrophillic OC in snow (col,lyr) [kg]
      mss_dst1 (-nsnow+1:0)   ,&! mass of dust species 1 in snow  (col,lyr) [kg]
      mss_dst2 (-nsnow+1:0)   ,&! mass of dust species 2 in snow  (col,lyr) [kg]
      mss_dst3 (-nsnow+1:0)   ,&! mass of dust species 3 in snow  (col,lyr) [kg]
      mss_dst4 (-nsnow+1:0)     ! mass of dust species 4 in snow  (col,lyr) [kg]
   
   real(r8), intent(out)    :: &
      snofrz(-nsnow+1:0)        ! snow freezing flag (col,lyr) [0/1]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      lkems                   ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      trad                      ! radiative temperature [K]
        
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   integer                  :: &
      sflopt                  ,&! snow fraction option: 1) external constant, 2) snow depth-dependent(Subin 2012)
      lkopt                   ,&! lake models: 1) CoLML, 2) FLake, 3) Simstrat, 4) XOML
      zopt                    ,&! lake surface roughness : 1) external constant, 2) wind-dependent(Subin 2012)
      betaopt                 ,&! fraction of solar radiation in the NIR : 1) external constant, 2) equation(CoLM)
      fetchopt                ,&! lake fetch : 1) external constant, 2) lake depth-dependent(Subin)
      etaopt                  ,&! extinction coefficient : 1) external constant, 2) lake depth-dependent(Subin)
      scwat                   ,&! surface category of water characteristics: (1) shallow lake, (2) deep lake
      j                       ,&! loop counter
      i                         ! loop counter

   real(r8)                 :: &
      htvp                    ,&
      fgrnd1                  ,&! ground heat flux [W/m2]
      hfice                   ,&! heat flux from ice to water [W/m2]
      hfsnow                  ,&! heat flux from snow to water [W/m2]
      hfsnowice               ,&! heat flux from snow to ice [W/m2]
      stateC10                ,&! Drag coefficient
      u_taus                  ,&! friction velocity [m/s]
      xlat                    ,&! latitude [degrees]
      xlon                    ,&! longitude [degrees]
      dpsed                   ,&! sediment depth [m]
      T_bs                    ,&! bottom sediment temperature condition [K]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      u2m                     ,&! 2 m height wind speed [m/s]
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      rnird                   ,&! reflected direct beam nir solar radiation (W/m**2)
      rniri                   ,&! reflected diffuse nir solar radiation (W/m**2]
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
!================================================================================
      !- Solar radiation related (not used)
      rnird = 0. ! reflected direct beam nir solar radiation (W/m**2) 
      rniri = 0. ! reflected diffuse nir solar radiation (W/m**2)

      T_bs  = t_ssb(nsoil)         ! +WMEJ Assume the temperature is consistent with the soil layers of other models. 
      ! dpsed = zissb(nsoil+1)       ! +WMEJ Assume the thickness is consistent with the soil layers of other models.
      dpsed = rflk_depth_bs_ref  ! +WMEJ The default recommended thickness in the model is 10m., only used in CoLML
      IF (dplak <= ShallowDeepPartition) THEN
         scwat = SHALLOW
      ELSE
         scwat = DEEP
      ENDIF

      !- Convert latitude and longitude to degrees
      xlat = rlat/d2r
      xlon = rlon/d2r

      !- Lake model controller
      !   Avoid using the USE MODULE  method to import, 
      !   as it may need to be adjusted for different lakes in the future,
      !   and is not conducive to subsequent testing.
      lkopt = DEF_External_Lake%DEF_LAKE_MODEL
      zopt = DEF_External_Lake%DEF_LAKE_ROUGHNESSS_SCHEME
      betaopt = DEF_External_Lake%DEF_LAKE_BETAPRIME_SCHEME
      fetchopt = DEF_External_Lake%DEF_LAKE_FETCH_SCHEME
      etaopt = DEF_External_Lake%DEF_LAKE_ETA_SCHEME
      sflopt = DEF_External_Lake%DEF_LAKE_SURFACE_FLUX_SCHEME

      IF (DEF_External_Lake%DEF_USE_SURFACE_LAYER) THEN 
         CALL LakeFluxes ( &
               ! "in" arguments    
               ! -------------------
               lkopt   , sflopt  , scwat     , nsnow    ,& 
               nlake   , zopt    , betaopt   , fetchopt ,& 
               etaopt  , xlat    , dtlak     , dplak    ,& 
               hwref   , htref   , hqref     , tmref    ,&
               qmref   , usurf   , vsurf     , arhos    ,&
               psurf   , sols    , soll      , solsd    ,&
               solld   , sabg    , rnird     , rniri    ,&
               crain   , lrain   , csnow     , lsnow    ,&
               lwdns   , icefr   , hpbl      , zcbcv    ,&
               dzlake  , zlake   , zssb      , dzssb    ,&
               t_ssb   , xwliq   , xwice     , zissb    ,&
               icedp   , bicedp  , wicedp    , snwdp    ,&
               stke1   , ipatch  , nsoil     ,&
               ! "inout" arguments
               ! -------------------
               tskin   , z0m     , z0h       , z0q      ,&
               felak   , btpri   , rhosnw    , lktmp    ,&
               ! "out" arguments
               ! -------------------
               fseng   , fsena   , fevpg     , fevpa    ,&
               lfevpa  , lwups   , fgrnd     , fgrnd1   ,&
               zol     , rib     , trad      , htvp     ,&
               lkems   , wdm     , ram       , rah      ,&
               raw     , shfdt   , taux      , tauy     ,&
               t2m     , q2m     , u10m      , v10m     ,&
               fh2m    , fq2m    , fm10m     , fq10m    ,&
               fm      , fh      , fq        , ustar    ,&
               qstar   , tstar   , u2m       , u_taus   ,&
               hfice   , hfsnow  , hfsnowice , stateC10 )
      ENDIF

      CALL LakeCtl ( &
            ! "in" arguments    
            ! -------------------
            nlake     , nsnow     , nsoil     , nlice     ,& 
            scwat     , dtlak     , xlat      , xlon      , dplak     ,& 
            lkopt     , zopt      , betaopt   , fetchopt  , etaopt    ,& 
            hwref     , htref     , hqref     , usurf     , vsurf     ,& 
            tmref     , qmref     , arhos     , psurf     , lwdns     ,&
            sols      , soll      , solsd     , solld     , sabg      ,&
            rnird     , rniri     , zcbcv     , hpbl      , tprec     ,& 
            crain     , lrain     , csnow     , lsnow     , dpsed     ,&           
            vf_quartz , vf_gravels, vf_om     , vf_sand   , wf_gravels,&
            wf_sand   , porsl     , csol      , k_solids  , dksatf    ,&
            dkdry     , dksatu    , BA_alpha  , BA_beta   , T_bs      ,&
            ipatch    , bifall    ,&
            ! "inout" arguments
            ! -------------------
            zlake     , zilak     , dzlake     , lktmp     , tskin     ,& 
            dzssb     , zssb      , zissb     , t_ssb     , stke1     ,& 
            snwcv     , snwag     , snwdp     , xwliq     , xwice     ,& 
            icefr     , snlay     , icedp     , etal      , btpri     ,&
            z0m       , z0h       , z0q       , felak     , gamma     ,&
            tmsno     , tmice     , tmmnw     , tmwml     , tmbot     ,&
            tmups     , mldp      , upsdp     , CTfrac    , frlak     ,&
            ziarea    , uwatv     , vwatv     , lksal     , tke       ,&
            eps       , etke      , num       , nuh       , bicedp    ,&
            wicedp    , lkrho     , rhosnw    , snout     ,&
! SNICAR model variables
            forc_aer  , sabg_lyr  , snofrz    ,&
            mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
            mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
            ! "out" arguments
            ! -------------------
            fsena     , fevpa     , lfevpa    , fseng     , fevpg     ,&
            qseva     , qsubl     , qsdew     , qfros     ,&
            lwups     , fgrnd     , trad      , taux      , tauy      ,&
            ustar     , qstar     , tstar     , lkems     , snwml     ,&
            zol       , rib       , ram       , rah       , raw       ,&
            wdm       , t2m       , q2m       , u10m      , v10m      ,&
            fm10m     , fq10m     , fh2m      , fq2m      , fm        ,&
            fq        , fh        , shfdt     , hfice     , hfsnow    ,&
            hfsnowice , stateC10  , u_taus    , htvp      , u2m       ,&
            fgrnd1    )  


   END SUBROUTINE Lake_Driver



   SUBROUTINE LakeCtl ( &
               ! "in" arguments    
               ! -------------------
               nlake     , nsnow     , nsoil     , nlice     ,& 
               scwat     , dtlak     , xlat      , xlon      , dplak     ,& 
               lkopt     , zopt      , betaopt   , fetchopt  , etaopt    ,& 
               hwref     , htref     , hqref     , usurf     , vsurf     ,& 
               tmref     , qmref     , arhos     , psurf     , lwdns     ,&
               sols      , soll      , solsd     , solld     , sabg      ,&
               rnird     , rniri     , zcbcv     , hpbl      , tprec     ,& 
               crain     , lrain     , csnow     , lsnow     , dpsed     ,&
               vf_quartz , vf_gravels, vf_om     , vf_sand   , wf_gravels,&
               wf_sand   , porsl     , csol      , k_solids  , dksatf    ,&
               dkdry     , dksatu    , BA_alpha  , BA_beta   , T_bs      ,&
               ipatch    , bifall    ,&
               ! "inout" arguments
               ! -------------------
               zlake     , zilak     , dzlake     , lktmp     , tskin     ,& 
               dzssb     , zssb      , zissb     , t_ssb     , stke1     ,& 
               snwcv     , snwag     , snwdp     , xwliq     , xwice     ,& 
               icefr     , snlay     , icedp     , etal      , btpri     ,&
               z0m       , z0h       , z0q       , felak     , gamma     ,&
               tmsno     , tmice     , tmmnw     , tmwml     , tmbot     ,&
               tmups     , mldp      , upsdp     , CTfrac    , frlak     ,&
               ziarea    , uwatv     , vwatv     , lksal     , tke       ,&
               eps       , etke      , num       , nuh       , bicedp    ,&
               wicedp    , lkrho     , rhosnw    , snout     ,&
! SNICAR model variables
               forc_aer  , sabg_lyr  , snofrz    ,&
               mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
               mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
               ! "out" arguments
               ! -------------------
               fsena     , fevpa     , lfevpa    , fseng     , fevpg     ,&
               qseva     , qsubl     , qsdew     , qfros     ,&
               lwups     , fgrnd     , trad      , taux      , tauy      ,&
               ustar     , qstar     , tstar     , lkems     , snwml     ,&
               zol       , rib       , ram       , rah       , raw       ,&
               wdm       , t2m       , q2m       , u10m      , v10m      ,&
               fm10m     , fq10m     , fh2m      , fq2m      , fm        ,&
               fq        , fh        , shfdt     , hfice     , hfsnow    ,&
               hfsnowice , stateC10  , u_taus    , htvp      , u2m       ,&
               fgrnd1    , urban_call)  

! ---------------------------------- code history -------------------------------------
! Description:
!     Lake Model Controller. This should be used as an isolator between the lake models
!
! *************************************************************************************  
! TO ENSURE CODE READABILITY, DO NOT PERFORM ANY CALCULATIONS OR MACRO DEFINITIONS HERE.
! *************************************************************************************  
!
! Called: (* means optional)
!    *-> Lake_CoLML        : CoLM-Lake 
!    *-> Lake_FLake        : FLake 
!    *-> Lake_SIMSTRAT     : Simstrat
!    *-> Lake_XOML         : XOML
!
! Original author:
!     Min Xu, 02/2012
!
! Revisions:
!     Xin-Zhong Liang, 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reconstructed the lake model driver
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: COLML, FLAKE, SIMSTRAT, XOML
   USE MOD_Lake_CoLML, only: Lake_CoLML
   USE MOD_Lake_FLake, only: Lake_FLake
   USE MOD_Lake_Simstrat, only: Lake_Simstrat
   USE MOD_SPMD_Task
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      lkopt                   ,&! lake models: 1) CoLML, 2) FLake, 3) Simstrat, 4) XOML
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer
      nsoil                   ,&! Maximum number of soil layer
      nlice                   ,&! Maximum number of ice layer, not used in this version
      scwat                   ,&! water body type 
      zopt                    ,&! lake surface roughness : 1) external constant, 2) wind-dependent(Subin 2012)
      betaopt                 ,&! fraction of solar radiation in the NIR : 1) external constant, 2) equation(CoLM)
      fetchopt                ,&! lake fetch : 1) external constant, 2) lake depth-dependent(Subin)
      etaopt                    ! lake option : 1) external constant, 2) lake depth-dependent(Subin)
        

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      xlat                    ,&! latitude [degree]
      xlon                    ,&! longitude [degree]
      dplak                   ,&! lake depth [m]
      dpsed                   ,&! sediment depth [m]
      T_bs                    ,&! bottom sediment temperature condition [K]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      rnird                   ,&! reflected direct beam nir solar radiation (W/m**2)
      rniri                   ,&! reflected diffuse nir solar radiation (W/m**2)
      sabg                    ,&! solar absorbed by ground  [W/m2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! downward longwave radiation at surface [wm-2]
      zcbcv                   ,&! convective boundary height [m]
      tprec                   ,&! snowfall/rainfall temperature [kelvin]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

   real(r8), intent(in)     :: &
      porsl     (nsoil)       ,&! soil porosity [-]
      dksatu    (nsoil)       ,&! thermal conductivity of saturated unfrozen soil [W/m-K]
      k_solids  (nsoil)       ,&! thermal conductivity of mineralssoil [W/m-K]
      dkdry     (nsoil)       ,&! thermal conductivity of dry soil [W/m-K]
      vf_quartz (nsoil)       ,&! volumetric fraction of quartz within mineral soil [-]
      vf_gravels(nsoil)       ,&! volumetric fraction of gravels [-]
      vf_om     (nsoil)       ,&! volumetric fraction of organic matter [-]
      vf_sand   (nsoil)       ,&! volumetric fraction of sand [-]
      wf_gravels(nsoil)       ,&! gravimetric fraction of gravels [-]
      wf_sand   (nsoil)       ,&! gravimetric fraction of sand [-]
      csol      (nsoil)       ,&! heat capacity of soil solids [J/(m3 K)]
      dksatf    (nsoil)       ,&! thermal conductivity of saturated frozen soil [W/m-K]
      BA_alpha  (nsoil)       ,&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nsoil)         ! beta in Balland and Arp(2005) thermal conductivity scheme

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlake(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      ziarea(nlake+1)         ,&! lake layer interface area [m2]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)

   real(r8), intent(inout)  :: &
      dzssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer thickness [m]
      zssb  (-nsnow+1:nsoil)  ,&! snow + soil + bedrock node layer depth [m]
      zissb (-nsnow:nsoil)    ,&! snow + soil + bedrock layer interface level [m]
      t_ssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer temperature [K]
      xwliq (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice (-nsnow+1:nsoil)    ! snow + soil + bedrock layer ice lens (kg/m2)

   real(r8), intent(inout)  :: &
      snwcv                   ,&! snow mass (kg/m2)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth [m]
      snout                   ,&! rate of water out of snow bottom (mm/s)
      tmsno                   ,&! snow temperature [K]
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      mldp                    ,&! mixed layer depth [m]
      upsdp                   ,&! bottom of the upper layer of the sediments [m]
      icedp                   ,&! ice depth [m]
      bicedp                  ,&! black snow depth (m)
      wicedp                  ,&! white snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      CTfrac                  ,&! Shape factor (thermocline)
      tskin                   ,&! ground surface temperature [k]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      frlak                     ! lake fraction [-]

   ! ----- SNICAR variables -----
   real(r8), intent(in)     :: &
      forc_aer(14)              ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

   real(r8), intent(in)     :: &
      sabg_lyr(-nsnow+1:1)      ! solar radiation absorbed by ground [W/m2]
   
   real(r8), intent(inout)  :: &
      mss_bcpho(-nsnow+1:0)   ,&! mass of hydrophobic BC in snow  (col,lyr) [kg]
      mss_bcphi(-nsnow+1:0)   ,&! mass of hydrophillic BC in snow (col,lyr) [kg]
      mss_ocpho(-nsnow+1:0)   ,&! mass of hydrophobic OC in snow  (col,lyr) [kg]
      mss_ocphi(-nsnow+1:0)   ,&! mass of hydrophillic OC in snow (col,lyr) [kg]
      mss_dst1 (-nsnow+1:0)   ,&! mass of dust species 1 in snow  (col,lyr) [kg]
      mss_dst2 (-nsnow+1:0)   ,&! mass of dust species 2 in snow  (col,lyr) [kg]
      mss_dst3 (-nsnow+1:0)   ,&! mass of dust species 3 in snow  (col,lyr) [kg]
      mss_dst4 (-nsnow+1:0)     ! mass of dust species 4 in snow  (col,lyr) [kg]
   
   real(r8), intent(out)    :: &
      snofrz(-nsnow+1:0)        ! snow freezing flag (col,lyr) [0/1]

!  ------------------------- output variables ---------------------------
   real(r8), intent(inout)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      htvp                    ,&! latent heat of vapor of water (or sublimation) [j/kg]
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      fgrnd1                  ,&! ground heat flux [W/m2]
      trad                    ,&! radiative temperature [K]
      hfice                   ,&! heat flux from ice to water [W/m2]
      hfsnow                  ,&! heat flux from snow to water [W/m2]
      hfsnowice               ,&! heat flux from snow to ice [W/m2]
      stateC10                ,&! Drag coefficient
      u_taus                  ,&! friction velocity [m/s]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      lkems                   ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      u2m                     ,&! 2 m height wind speed in eastward direction [m/s]
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]

   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL
!================================================================================

      ! *************************************************************************************
      SELECT CASE (lkopt)   !- Select the lake model  !
      ! *************************************************************************************
         !-----------------------
         CASE (COLML) !- CoLM-Lake
         !-----------------------
            ! print *, 'Lake Model: CoLM-Lake'
            CALL Lake_CoLML ( &
                  ! "in" arguments
                  ! ---------------------------
                  nlake     , nsnow     , nsoil     , scwat     ,&
                  zopt      , betaopt   , fetchopt  , etaopt    ,&
                  dtlak     , xlat      , xlon      , dplak     ,&
                  hwref     , htref     , hqref     , tmref     ,&
                  usurf     , vsurf     , qmref     , arhos     ,&
                  sols      , soll      , solsd     , solld     ,&
                  sabg      , lwdns     , psurf     , hpbl      ,&
                  crain     , csnow     , lrain     , lsnow     ,&
                  porsl     , dksatu    , dkdry     , k_solids  ,&
                  vf_quartz , vf_gravels, vf_om     , vf_sand   ,&
                  wf_gravels, wf_sand   , csol      , dksatf    ,&
                  BA_alpha  , BA_beta   , zcbcv     , tprec     ,&
                  ipatch    , bifall    ,&
                  ! "inout" arguments
                  ! ---------------------------
                  dzlake     , zlake     , zilak     , lktmp     ,&
                  dzssb     , zssb      , zissb     , t_ssb     ,&
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  icefr     , tmice     , snout     , lkrho     ,& 
! SNICAR model variables
                  forc_aer  , sabg_lyr  , snofrz    ,&
                  mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
                  mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&  
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     )

         !-----------------------
         CASE (FLAKE) !- FLake
         !-----------------------
            ! print *, 'Lake Model: FLake'
            CALL Lake_FLake ( &
                  ! "in" arguments
                  ! ---------------------------
                  nlake     , nsnow     , nsoil     , scwat     ,&
                  zopt      , betaopt   , fetchopt  , etaopt    ,&
                  dtlak     , xlat      , xlon      , dplak     ,&
                  hwref     , htref     , hqref     , tmref     ,&
                  usurf     , vsurf     , qmref     , arhos     ,&
                  sols      , soll      , solsd     , solld     ,&
                  sabg      , lwdns     , psurf     , hpbl      ,&
                  crain     , csnow     , lrain     , lsnow     ,&
                  dpsed     , T_bs      , zcbcv     , ipatch    ,&
                  ! "inout" arguments
                  ! ---------------------------
                  dzlake     , zlake     , zilak     , lktmp     ,&
                  dzssb     , zssb      , zissb     , t_ssb     ,&
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  tmsno     , tmice     , tmmnw     , tmwml     ,&  
                  tmbot     , tmups     , icedp     , mldp      ,&  
                  upsdp     , CTfrac    , icefr     , &              
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     )

         !-----------------------
         CASE (SIMSTRAT) !- Simstrat
         !-----------------------
            ! print *, 'Lake Model: Simstrat'
            CALL Lake_Simstrat( &
                  ! "in" arguments    
                  ! -------------------
                  nlake     , nsnow     , scwat     , nsoil     ,&
                  zopt      , betaopt   , fetchopt  , etaopt    ,&
                  dtlak     , xlat      , xlon      , dplak     ,&
                  hwref     , htref     , hqref     , tmref     ,&
                  usurf     , vsurf     , qmref     , arhos     ,&
                  sols      , soll      , solsd     , solld     ,&
                  sabg      , lwdns     , psurf     , hpbl      ,&
                  crain     , csnow     , lrain     , lsnow     ,&
                  zcbcv     , ipatch    ,&
                  ! "inout" arguments
                  ! ---------------------------
                  dzlake     , zlake     , zilak     , lktmp     ,&
                  dzssb     , zssb      , zissb     , t_ssb     ,&
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  icefr     , icedp     , tmice     , ziarea    ,&
                  uwatv     , vwatv     , lksal     , tke       ,&
                  eps       , etke      , num       , nuh       ,&
                  bicedp    , wicedp    , lkrho     , rhosnw    ,&
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     , hfice     , hfsnow    , hfsnowice ,&
                  stateC10  , u_taus    , fgrnd1    )   

         !-----------------------
         CASE (XOML) !- XOML
         !-----------------------
            print *, '[Error] XOML is not implemented yet'
            CALL CoLM_stop ()
        END SELECT 

   END SUBROUTINE LakeCtl



   subroutine LakeFluxes ( &
                        ! "in" arguments    
                        ! -------------------
                        lkopt   , sflopt  , scwat     , nsnow    ,& 
                        nlake   , zopt    , betaopt   , fetchopt ,& 
                        etaopt  , xlat    , dtlak     , dplak    ,& 
                        hwref   , htref   , hqref     , tmref    ,&
                        qmref   , usurf   , vsurf     , arhos    ,&
                        psurf   , sols    , soll      , solsd    ,&
                        solld   , sabg    , rnird     , rniri    ,&
                        crain   , lrain   , csnow     , lsnow    ,&
                        lwdns   , icefr   , hpbl      , zcbcv    ,&
                        dzlake  , zlake   , zssb      , dzssb    ,&
                        t_ssb   , xwliq   , xwice     , zissb    ,&
                        icedp   , bicedp  , wicedp    , snwdp    ,&
                        stke1   , ipatch  , nsoil     ,&
                        ! "inout" arguments
                        ! -------------------
                        tskin   , z0m     , z0h       , z0q      ,&
                        felak   , btpri   , rhosnw    , lktmp    ,&
                        ! "out" arguments
                        ! -------------------
                        fseng   , fsena   , fevpg     , fevpa    ,&
                        lfevpa  , lwups   , fgrnd     , fgrnd1   ,&
                        zol     , rib     , trad      , htvp     ,&
                        lkems   , wdm     , ram       , rah      ,&
                        raw     , shfdt   , taux      , tauy     ,&
                        t2m     , q2m     , u10m      , v10m     ,&
                        fh2m    , fq2m    , fm10m     , fq10m    ,&
                        fm      , fh      , fq        , ustar    ,&
                        qstar   , tstar   , u2m       , u_taus   ,&
                        hfice   , hfsnow  , hfsnowice , stateC10 ,&
                        urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!     Lake surface flux driver. The solution found so far is incorrect when there is ice and snow.
!     The amount of solar radiation penetrating each dielectric layer should be determined
!     before calling surface flux, otherwise the total energy entering the lake will be unbalanced.
!
! Called: (* means optional)
!    *-> CoLML_SurfFlux      : CoLM-Lake surface flux scheme (Monin-Obukhov)
!    *-> FLake_SurfFlux      : FLake surface flux scheme (Monin-Obukhov)
!    *-> Simstrat_SurfFlux   : Simstrat surface flux scheme (bulk aerodynamic)
!    *-> LakStaParms         : Calculate the lake state parameters missing in Simstrat and FLake solutions
!    *-> Simstrat_FluxAdd    : Supplementary calculations for special variables of Simstrat models
!
! Original author:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 10/2024
!
! -------------------------------------------------------------------------------------
   use MOD_Lake_CoLML, only : CoLML_SurfFlux
   use MOD_Lake_FLake, only : FLake_SurfFlux
   use strat_forcing, only : Simstrat_SurfFlux, Simstrat_FluxAdd
   USE MOD_Lake_Const, only : COLML, FLAKE, SIMSTRAT, tpsf_L_evap,&
                                    beta_sol, tpl_rho_w_r, d2r, c_lwrad_emis
   use MOD_Lake_Subs, only : LakStaParms
   USE MOD_SPMD_Task
   !==============================================================================
   ! ----- input variables -----
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      lkopt                   ,&! lake models:  1) CoLML, 2) Simstrat, 3) FLake, 4) XOML, 5) CSSPL, 6) LISSS, 7) GRELR
      sflopt                  ,&! surface layer scheme option: 1) CoLML, 2) Simstrat, 3) FLake
      scwat                   ,&! surface category of water characteristics:
      nsnow                   ,&! number of snow layers
      nsoil                   ,&! number of soil layers
      nlake                   ,&! number of lake layers
      zopt                    ,&! option for the lake model
      betaopt                 ,&! option for the beta function
      fetchopt                ,&! option for the fetch function
      etaopt                    ! option for the eta function
      
   real(r8), intent(in)     :: &
      xlat                    ,&! latitude in degrees
      dtlak                   ,&! time step [s]
      dplak                   ,&! lake depth [m]
      stke1                   ,&! top level eddy conductivity [W/m/K]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      sabg                    ,&! solar absorbed by ground  [W/m2]
      rnird                   ,&! reflected direct beam nir solar radiation (W/m**2)
      rniri                   ,&! reflected diffuse nir solar radiation (W/m**2)
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! downward longwave radiation at surface [wm-2]
      zcbcv                   ,&! lake bottom depth [m]
      wicedp                  ,&! white ice depth [m]
      bicedp                  ,&! black ice depth [m]
      icedp                   ,&! totle ice depth [m]
      snwdp                   ,&! snow depth [m]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

   real(r8), intent(in)     :: &
      dzlake(nlake)           ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer depth [m]
      icefr(nlake)            ,&! ice fraction [-]
      zssb(-nsnow+1:nsoil)    ,&! soil + snow layer depth [m]
      zissb(-nsnow:nsoil )    ,&! snow + soil + bedrock layer interface level [m]
      dzssb(-nsnow+1:nsoil)   ,&! layer thickness [m]
      t_ssb(-nsnow+1:nsoil)   ,&! soil + snow layer temperature [K]
      xwliq(-nsnow+1:nsoil)   ,&! liquid water (kg/m2)
      xwice(-nsnow+1:nsoil)     ! ice lens (kg/m2)

   ! ----- input/output variables -----
   real(r8), intent(inout)  :: &
      lktmp(nlake)            ,&! lake temperature (K)
      tskin                   ,&! ground surface temperature [k]
      rhosnw                  ,&! snow density [kg/m3]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise use emperical eq by dplak [m]
      btpri                     ! beta prime in Monin-Obukhov theory

   ! ----- output variables -----
   real(r8), intent(out)    :: &
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      fgrnd1                  ,&! ground heat flux [W/m2]
      hfice                   ,&! heat flux from ice to water [W/m2]
      hfsnow                  ,&! heat flux from snow to water [W/m2]
      hfsnowice               ,&! heat flux from snow to ice [W/m2]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      rib                     ,&! bulk Richardson number in surface layer
      trad                    ,&! radiative temperature [K]
      htvp                    ,&! heat transfer velocity [m/s]
      lkems                   ,&! emissivity of the surface
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      stateC10                ,&! Drag coefficient
      u_taus                  ,&! friction velocity [m/s]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      t2m                     ,&! reference temperature [K]
      q2m                     ,&! reference specific humidity [kg/kg]
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      u2m                     ,&! 2 m height wind speed in eastward direction [m/s]
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      fm                      ,&! integral of profile function for momentum
      fh                      ,&! integral of profile function for heat
      fq                      ,&! integral of
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]

   !+WMEJ urban_call is optional argument for subroutine laketem, it is not used in this version
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

   ! ----- local variables -----
   integer                  :: &
      j                       ,&! loop index
      snlay                   ,&! snow layer on ice lake
      lb                        ! snow layer on ice lake
   
   real(r8)                 :: &
      U_a                     ,&! wind speed at reference height [m/s]
      rlat                    ,&! latitude in radians
      rad0                    ,&! remaining radiation [W/m2]
      ustar_                  ,&! friction velocity [m/s]
      mmtum                   ,&! momentum flux [kg/m/s**2]
      extsn                   ,&! extinction coefficient for snow
      precip                    ! total precipitation [mm/s]
   !==============================================================================
      ! define snow layer on ice lake
      snlay = 0
      DO j=-nsnow+1,0
         if(xwliq(j)+xwice(j)>0.) snlay=snlay-1
      ENDDO 
      lb = snlay+1
      U_a = max(0.1,sqrt(usurf*usurf+vsurf*vsurf))     ! limit set to 0.1
      precip = max(crain + csnow + lrain + lsnow, 0.0_r8) ! total precipitation [mm/s]
      rlat = xlat * d2r

      ! *************************************************************************************
      SELECT CASE (sflopt)   !- Select the lake model
      ! *************************************************************************************
         !-----------------------
         CASE (COLML) !- COLML
         !-----------------------
               ! print *, '    *-*    Surface Layer Scheme: CoLM-Lake    *-*    '
               !- Lake surface flux calculations
               CALL CoLML_SurfFlux ( &
                  ! "in" arguments
                  ! -------------------
                  scwat     , snlay     , zopt         , betaopt   ,&
                  fetchopt  , rlat      , dplak        , dtlak     ,& 
                  stke1     , hwref     , htref        , hqref     ,&
                  usurf     , vsurf     , tmref        , qmref     ,&
                  arhos     , psurf     , sols         , soll      ,&
                  solsd     , solld     , lwdns        , sabg      ,&
                  hpbl      , dzlake(1) , zlake(nlake) , dzssb(lb) ,&
                  lktmp(1)  , t_ssb(lb) , xwliq(lb)    , xwice(lb) ,&
                  icefr(1)  , zcbcv     , ipatch       ,&
                  ! "inout" arguments
                  ! -------------------
                  tskin     , z0m       , z0h          , z0q       ,&
                  felak     , btpri     ,&
                  ! "out" arguments
                  ! -------------------
                  fseng     , fevpg     , fsena        , fevpa     ,&
                  lfevpa    , lwups     , fgrnd        , fgrnd1    ,&
                  zol       , rib       , trad         , htvp      ,&
                  lkems     , wdm       , ram          , rah       ,&
                  raw       , shfdt     , taux         , tauy      ,&
                  t2m       , q2m       , u10m         , v10m      ,&
                  fh2m      , fq2m      , fm10m        , fq10m     ,&
                  fm        , fh        , fq           , ustar     ,&
                  qstar     , tstar     , rhosnw       , u2m       )

            ! print *, ''
            ! print *, '+WMEJ LakeFluxes: ipatch    = ', ipatch
            ! print *, '+WMEJ LakeFluxes:  scwat    = ', scwat
            ! print *, '+WMEJ LakeFluxes:  snlay    = ', snlay
            ! print *, '+WMEJ LakeFluxes:  zopt     = ', zopt
            ! print *, '+WMEJ LakeFluxes:  betaopt  = ', betaopt
            ! print *, '+WMEJ LakeFluxes:  fetchopt = ', fetchopt
            ! print *, '+WMEJ LakeFluxes:  rlat     = ', rlat
            ! print *, '+WMEJ LakeFluxes:  dplak    = ', dplak
            ! print *, '+WMEJ LakeFluxes:  dtlak    = ', dtlak
            ! print *, '+WMEJ LakeFluxes:  stke1    = ', stke1
            ! print *, '+WMEJ LakeFluxes:  hwref    = ', hwref
            ! print *, '+WMEJ LakeFluxes:  htref    = ', htref
            ! print *, '+WMEJ LakeFluxes:  hqref    = ', hqref
            ! print *, '+WMEJ LakeFluxes:  tmref    = ', tmref
            ! print *, '+WMEJ LakeFluxes:  qmref    = ', qmref
            ! print *, '+WMEJ LakeFluxes:  arhos    = ', arhos
            ! print *, '+WMEJ LakeFluxes:  psurf    = ', psurf
            ! print *, '+WMEJ LakeFluxes:  sols     = ', sols
            ! print *, '+WMEJ LakeFluxes:  soll     = ', soll
            ! print *, '+WMEJ LakeFluxes:  solsd    = ', solsd
            ! print *, '+WMEJ LakeFluxes:  solld    = ', solld
            ! print *, '+WMEJ LakeFluxes:  lwdns    = ', lwdns
            ! print *, '+WMEJ LakeFluxes:  sabg     = ', sabg
            ! print *, '+WMEJ LakeFluxes:  hpbl     = ', hpbl
            ! print *, '+WMEJ LakeFluxes:  dzlake   = ', dzlake(1)
            ! print *, '+WMEJ LakeFluxes:  zlake    = ', zlake(nlake)
            ! print *, '+WMEJ LakeFluxes:  dzssb    = ', dzssb(lb)
            ! print *, '+WMEJ LakeFluxes:  lktmp    = ', lktmp(1)
            ! print *, '+WMEJ LakeFluxes:  t_ssb    = ', t_ssb(lb)
            ! print *, '+WMEJ LakeFluxes:  xwliq    = ', xwliq(lb)
            ! print *, '+WMEJ LakeFluxes:  xwice    = ', xwice(lb)
            ! print *, '+WMEJ LakeFluxes:  icefr    = ', icefr(1)
            ! print *, '+WMEJ LakeFluxes:  zcbcv    = ', zcbcv
            ! print *, '+WMEJ LakeFluxes:  tskin    = ', tskin
            ! print *, '+WMEJ LakeFluxes:  z0m      = ', z0m
            ! print *, '+WMEJ LakeFluxes:  z0h      = ', z0h
            ! print *, '+WMEJ LakeFluxes:  z0q      = ', z0q
            ! print *, '+WMEJ LakeFluxes:  felak    = ', felak
            ! print *, '+WMEJ LakeFluxes:  btpri    = ', btpri
            ! print *, ''
            ! print *, ''

         !-----------------------
         CASE (FLAKE) !- FLake
         !-----------------------
               ! print *, '    *-*    Surface Layer Scheme: FLake    *-*    '
               CALL FLake_SurfFlux ( &
                  ! "in" arguments
                  ! ---------------------------
                  scwat     , zopt      , fetchopt     , dplak     ,&
                  tmref     , qmref     , arhos        , U_a       ,&
                  icedp     , hwref     , hqref        , tskin     ,&
                  lwdns     , sabg      , psurf        ,&
                  ! "inout" arguments
                  ! -------------------
                  z0m       , z0h       , z0q          , felak     ,&
                  ! "out" arguments
                  ! ---------------------------
                  mmtum     , ustar_    , fsena        , lfevpa    ,&
                  fevpa     , lwups     , fgrnd        )
               if (betaopt /= 1) then
                  btpri = beta_sol
               endif

         !-----------------------
         CASE (SIMSTRAT) !- Simstrat
         !-----------------------
               ! print *, '    *-*    Surface Layer Scheme: Simstrat    *-*    '
               CALL Simstrat_SurfFlux( &
                  ! "in" arguments
                  ! ---------------------------
                  snwdp   , rhosnw       , wicedp      , bicedp    ,&
                  usurf   , vsurf        , tmref       , qmref     ,&
                  psurf   , sabg         , lwdns       ,&
                  csnow   , lsnow        , crain       , lrain     ,&
                  ! "inout" arguments
                  ! -------------------
                  tskin   , lktmp(nlake) ,&
                  ! "out" arguments
                  ! ---------------------------
                  hfsnow   , hfsnowice   , hfice       , stateC10  ,&
                  u_taus   , taux        , tauy        , htvp      ,&
                  lwups    , fsena       , fseng       , lfevpa    ,&
                  fevpa    , fevpg       , sabg        , fgrnd     ,&
                  fgrnd1   , rad0        , trad        , lkems      )
               btpri = beta_sol

         !-----------------------
         CASE DEFAULT !- error
         !-----------------------
               print *, '    *-*    Error: Lake model not defined    *-*    '
               CALL CoLM_stop ()

      ! *************************************************************************************
      END SELECT
      ! *************************************************************************************

      !+WMEJ additonal output variables for coupling
      IF (sflopt == FLAKE .or. sflopt == SIMSTRAT) THEN 
         tskin = lktmp(1)
         CALL LakStaParms ( &
               ! "in" arguments
               ! -------------------
               hwref     , htref     , hqref        , usurf     ,&
               vsurf     , tmref     , qmref        , tskin     ,&
               arhos     , psurf     , zcbcv        ,& 
               ! "out" arguments
               ! -------------------
               zol       , rib       , wdm          , shfdt     ,&
               taux      , tauy      , t2m          , q2m       ,&
               u10m      , v10m      , fh2m         , fq2m      ,&
               fm10m     , fq10m     , fm           , fh        ,&
               fq        , ustar     , qstar        , tstar     ,&
               ram       , rah       , raw          , lkems     ,&
               u2m       )
      ENDIF

      !+WMEJ additonal output variables 
      IF (sflopt == FLAKE .or. sflopt == COLML) THEN
         IF (lkopt == SIMSTRAT) THEN
            CALL Simstrat_FluxAdd( &
                  ! "in" arguments
                  ! ---------------------------
                  btpri     , usurf     , vsurf        , tmref     ,&
                  snwdp     , wicedp    , bicedp       , lwdns     ,&
                  csnow     , lsnow     , crain        , lrain     ,&
                  fsena     , lfevpa    , lwups        , sabg      ,&
                  ! "out" arguments
                  ! ---------------------------
                  hfsnow    , hfsnowice , hfice        , stateC10  ,&
                  rad0      , u_taus    , taux         , tauy       )
                  !-WMEJ Lake surface shear stresses (taux, tauy) should probably not be recalculated
         ENDIF
         IF (sflopt == FLAKE) THEN
                  lkems = c_lwrad_emis
                  htvp = tpsf_L_evap
                  fseng = fsena
                  fevpg = fevpa
                  ustar = ustar_
                  fgrnd1 = btpri * sabg + lwdns - lwups - fsena - lfevpa
         ENDIF
      ENDIF

   end subroutine LakeFluxes


END MODULE  MOD_Lake_Driver

