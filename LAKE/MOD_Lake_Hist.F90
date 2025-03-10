#include <define.h>

MODULE MOD_Lake_Hist
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
   USE MOD_Lake_1DAccVars
   IMPLICIT NONE
!================================================================================

   character(len=10), PRIVATE :: HistForm ! 'Gridded', 'Vector', 'Single'

   PUBLIC :: LakeVarsSaveHist
   PRIVATE :: write_history_variable_2d
   PRIVATE :: write_history_variable_3d

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE LakeVarsSaveHist (nl_lake, file_hist, HistForm_in, itime_infile, sumarea, filter)
! ---------------------------------- code history -------------------------------------
! Description:
!     Save lake variable outputs to history file. Will be called in the subroutine hist_out.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Namelist, only: DEF_External_Lake_hist
   USE MOD_DataType
!==============================================================================
!  ------------------------- input variables ---------------------------
   type(block_data_real8_2d), intent(in) :: sumarea
   character(len=256)       , intent(in) :: file_hist
   character(len=10)        , intent(in) :: HistForm_in ! 'Gridded', 'Vector', 'Single'
   integer                  , intent(in) :: itime_infile
   integer                  , intent(in) :: nl_lake
   logical                  , intent(in) :: filter(:)

!================================================================================
      HistForm = HistForm_in
     
      CALL write_history_variable_2d ( DEF_External_Lake_hist%dplak, &
         a_dplak, file_hist, 'f_dplak', itime_infile, sumarea, filter, &
         'lake depth ','m')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%zlake, &
         a_zlake, file_hist, 'f_zlake', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Lake layer node depth','m')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%zilak, &
         a_zilak, file_hist, 'f_zilak', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'Lake layer interface depth','m')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%dzlake, &
         a_dzlake, file_hist, 'f_dzlake', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Lake layer thickness','m')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%ziarea, &
         a_ziarea, file_hist, 'f_ziarea', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'Lake layer interface area','m2')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%lkz0m, &
         a_lkz0m, file_hist, 'f_lkz0m', itime_infile, sumarea, filter, &
         'Roughness length for sensible heat','m')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%lkz0h, &
         a_lkz0h, file_hist, 'f_lkz0h', itime_infile, sumarea, filter, &
         'Roughness length for sensible heat','m')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%lkz0q, &
         a_lkz0q, file_hist, 'f_lkz0q', itime_infile, sumarea, filter, &
         'Roughness length for latent heat','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%felak, &
         a_felak, file_hist, 'f_felak', itime_infile, sumarea, filter, &
         'Lake fetch length','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%gamma, &
         a_gamma, file_hist, 'f_gamma', itime_infile, sumarea, filter, &
         'Mixing enhancement factor','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%etal, &
         a_etal, file_hist, 'f_etal', itime_infile, sumarea, filter, &
         'Lake extinction coefficient','1/m')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%btpri, &
         a_btpri, file_hist, 'f_btpri', itime_infile, sumarea, filter, &
         'Beta prime in Monin-Obukhov theory','1')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%frlak, &
         a_frlak, file_hist, 'f_frlak', itime_infile, sumarea, filter, &
         'Lake fraction','1')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmsno, &
         a_tmsno, file_hist, 'f_tmsno', itime_infile, sumarea, filter, &
         'Mean temperature of the snow column','K')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmice, &
         a_tmice, file_hist, 'f_tmice', itime_infile, sumarea, filter, &
         'Mean temperature of the ice column','K')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmmnw, &
         a_tmmnw, file_hist, 'f_tmmnw', itime_infile, sumarea, filter, &
         '!Mean temperature of the water column','K')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmwml, &
         a_tmwml, file_hist, 'f_tmwml', itime_infile, sumarea, filter, &
         'Mixed-layer temperature','K')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmbot, &
         a_tmbot, file_hist, 'f_tmbot', itime_infile, sumarea, filter, &
         'Temperature at the water-bottom sediment interface','K')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%tmups, &
         a_tmups, file_hist, 'f_tmups', itime_infile, sumarea, filter, &
         'Temperature at the bottom of the upper layer of the sediments','K')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%mldp, &
         a_mldp, file_hist, 'f_mldp', itime_infile, sumarea, filter, &
         'Mixed layer depth','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%upsdp, &
         a_upsdp, file_hist, 'f_upsdp', itime_infile, sumarea, filter, &
         'Bottom of the upper layer of the sediments','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%icedp, &
         a_icedp, file_hist, 'f_icedp', itime_infile, sumarea, filter, &
         'Mean temperature of the lake','K')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%bicedp, & 
         a_bicedp, file_hist, 'f_bicedp', itime_infile, sumarea, filter, &
         'black ice depth','m')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%wicedp, &
         a_wicedp, file_hist, 'f_wicedp', itime_infile, sumarea, filter, &
         'white ice depth','m')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%CTfrac, &
         a_CTfrac, file_hist, 'f_CTfrac', itime_infile, sumarea, filter, &
         'Shape factor (thermocline)','-')

      CALL write_history_variable_2d ( DEF_External_Lake_hist%rhosnw, &
         a_rhosnw, file_hist, 'f_rhosnw', itime_infile, sumarea, filter, &
         'snow density','kg/m3')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%uwatv, &
         a_uwatv, file_hist, 'f_uwatv', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Water velocity in x-direction','m/s')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%vwatv, &
         a_vwatv, file_hist, 'f_vwatv', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Water velocity in y-direction','m/s')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%lksal, &
         a_lksal, file_hist, 'f_lksal', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Salinity','‰')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%tke, &
         a_tke, file_hist, 'f_tke', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'Turbulent kinetic energy','J/kg')
      
      CALL write_history_variable_2d ( DEF_External_Lake_hist%etke, &
         a_etke, file_hist, 'f_etke', itime_infile, sumarea, filter, &
         'Seiche energy','J')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%eps, &
         a_eps, file_hist, 'f_eps', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'TKE dissipation rate','W/kg')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%num, &
         a_num, file_hist, 'f_num', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'Turbulent viscosity (momentum)','m2/s')

      CALL write_history_variable_3d ( DEF_External_Lake_hist%nuh, &
         a_nuh, file_hist, 'f_nuh', itime_infile, 'lakeI', 1, nl_lake+1, sumarea, filter, &
         'Turbulent diffusivity (heat)','m2/s')
      
      CALL write_history_variable_3d ( DEF_External_Lake_hist%lkrho, &
         a_lkrho, file_hist, 'f_lkrho', itime_infile, 'lake', 1, nl_lake, sumarea, filter, &
         'Density of water','kg/m3')

   END SUBROUTINE LakeVarsSaveHist



   SUBROUTINE write_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_infile, sumarea, filter, &
         longname, units)
   USE MOD_DataType
   USE MOD_HistGridded
   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_infile
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_2d ( &
            acc_vec, file_hist, varname, itime_infile, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_2d ( &
            acc_vec, file_hist, varname, itime_infile, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_infile, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_2d



   SUBROUTINE write_history_variable_3d ( is_hist, &
         acc_vec, file_hist, varname, itime_infile, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units)
   USE MOD_DataType
   USE MOD_HistGridded
   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_infile
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units

   ! Local variables
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_3d ( &
            acc_vec, file_hist, varname, itime_infile, dim1name, lb1, ndim1, &
            sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_3d ( &
            acc_vec, file_hist, varname, itime_infile, dim1name, lb1, ndim1, &
            filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_3d (acc_vec, file_hist, varname, itime_infile, &
            dim1name, ndim1, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_3d


END MODULE MOD_Lake_Hist