! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! Simulated series for calibration targets
!
! A calibration target scores a simulated series against observations.  Routed streamflow comes from
! mizuRoute, but every other target - basin storage change against GRACE, stream temperature, the
! water table under a coupled MODFLOW 6 model - reads an ordinary SUMMA variable that until now only
! reached the output file.
!
! This module collects any such variable, by name, as a basin time series: the name is resolved once
! against SUMMA's own metadata, and the value is then accumulated each model time step.  Resolving
! against the metadata rather than a list kept here means a target can name any SUMMA variable,
! and names stay honest - an unknown one is refused when the calibration starts rather than
! silently scored against something else.
! **************************************************************************************************
module simulated_series

  USE nr_type,    only: i4b,rkind,lgt
  USE summa_type, only: summa1_type_dec

  USE globalData, only: realMissing
  USE globalData, only: integerMissing
  USE globalData, only: gru_struc              ! gru-hru mapping, with the domain information

  USE var_lookup, only: iLookATTR              ! named variables for the local attributes
  USE var_lookup, only: iLookBVAR              ! named variables for the basin-average variables

  implicit none
  private

  ! which SUMMA structure a simulated series is read from
  integer(i4b), parameter, public :: ix_series_bvar  = 1   ! basin-average variables (GRU)
  integer(i4b), parameter, public :: ix_series_diag  = 2   ! diagnostic variables   (HRU/domain)
  integer(i4b), parameter, public :: ix_series_prog  = 3   ! prognostic variables   (HRU/domain)
  integer(i4b), parameter, public :: ix_series_flux  = 4   ! model fluxes           (HRU/domain)
  integer(i4b), parameter, public :: ix_series_param = 5   ! model parameters       (HRU/domain)

  ! One simulated series: a named SUMMA variable, resolved to its structure and index, and the
  ! basin value of that variable at every model time step.
  type, public :: sim_series_type
    character(len=64)        :: name = ''                 ! variable name, as the target asked for it
    character(len=64)        :: units = ''                ! its units, from SUMMA's own metadata
    integer(i4b)             :: ix_struct = integerMissing ! structure the variable lives in
    integer(i4b)             :: ix_var    = integerMissing ! index of the variable within it
    real(rkind), allocatable :: values(:)                  ! basin value at each model time step
  end type sim_series_type

  public :: init_simulated_series
  public :: collect_simulated_series
  public :: find_simulated_series
  public :: is_routed_streamflow

contains

  ! **************************************************************************************************
  ! Report whether a variable name means the routed streamflow mizuRoute produces.
  !
  ! Streamflow is the one simulated series that is not a SUMMA variable, so it is named rather than
  ! resolved, and every caller asks the same way.
  ! **************************************************************************************************
  pure function is_routed_streamflow(name) result(isFlow)
    character(*), intent(in) :: name
    logical(lgt)             :: isFlow

    select case(trim(name))
      case ('streamflow','discharge'); isFlow = .true.
      case default;                    isFlow = .false.
    end select

  end function is_routed_streamflow

  ! **************************************************************************************************
  ! Resolve a list of variable names and allocate a series for each.
  !
  ! Each name is looked up in SUMMA's metadata, in turn, as a basin-average variable, a diagnostic, a
  ! prognostic, a flux, and a parameter.  A name that matches none of those is refused here, at
  ! start-up, naming the variable that could not be found.
  !
  ! Routed streamflow is not a SUMMA variable and is not collected here; callers ask for it with
  ! is_routed_streamflow and take it from mizuRoute.
  ! **************************************************************************************************
  subroutine init_simulated_series(names,numtim,series,err,message)
    USE get_ixname_module, only: get_ixBvar,get_ixDiag,get_ixProg,get_ixFlux,get_ixParam
    USE globalData, only: bvar_meta,diag_meta,prog_meta,flux_meta,mpar_meta
    implicit none
    character(*),                           intent(in)  :: names(:)   ! variable names to collect
    integer(i4b),                           intent(in)  :: numtim     ! number of model time steps
    type(sim_series_type), allocatable,     intent(out) :: series(:)  ! the resolved series
    integer(i4b),                           intent(out) :: err        ! error code
    character(*),                           intent(out) :: message    ! error message
    integer(i4b) :: iSeries
    integer(i4b) :: ix

    err=0
    message='init_simulated_series/'

    allocate(series(size(names)),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate the simulated series'
      return
    endif

    do iSeries=1,size(names)
      series(iSeries)%name=trim(names(iSeries))

      ! resolve the name against SUMMA's metadata, most likely structure first
      ix=get_ixBvar(trim(names(iSeries)))
      if(ix/=integerMissing)then
        series(iSeries)%ix_struct=ix_series_bvar
      else
        ix=get_ixDiag(trim(names(iSeries)))
        if(ix/=integerMissing)then
          series(iSeries)%ix_struct=ix_series_diag
        else
          ix=get_ixProg(trim(names(iSeries)))
          if(ix/=integerMissing)then
            series(iSeries)%ix_struct=ix_series_prog
          else
            ix=get_ixFlux(trim(names(iSeries)))
            if(ix/=integerMissing)then
              series(iSeries)%ix_struct=ix_series_flux
            else
              ix=get_ixParam(trim(names(iSeries)))
              if(ix/=integerMissing) series(iSeries)%ix_struct=ix_series_param
            endif
          endif
        endif
      endif

      if(ix==integerMissing)then
        message=trim(message)//'"'//trim(names(iSeries))//'" is not a SUMMA variable: a calibration '// &
                'target must name "streamflow" or a variable SUMMA knows'
        err=20; return
      endif
      series(iSeries)%ix_var=ix

      ! units come from the same metadata as the variable, so the check against the observation
      ! file's units means something rather than being waved through
      select case(series(iSeries)%ix_struct)
        case (ix_series_bvar);  series(iSeries)%units=trim(bvar_meta(ix)%varunit)
        case (ix_series_diag);  series(iSeries)%units=trim(diag_meta(ix)%varunit)
        case (ix_series_prog);  series(iSeries)%units=trim(prog_meta(ix)%varunit)
        case (ix_series_flux);  series(iSeries)%units=trim(flux_meta(ix)%varunit)
        case (ix_series_param); series(iSeries)%units=trim(mpar_meta(ix)%varunit)
      end select

      allocate(series(iSeries)%values(numtim),source=realMissing,stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate values for simulated series "'//trim(names(iSeries))//'"'
        return
      endif
    enddo

  end subroutine init_simulated_series

  ! **************************************************************************************************
  ! Return the index of a named series, or integerMissing when it was not collected.
  ! **************************************************************************************************
  pure function find_simulated_series(series,name) result(iSeries)
    type(sim_series_type), intent(in) :: series(:)
    character(*),          intent(in) :: name
    integer(i4b)                      :: iSeries
    integer(i4b) :: i

    iSeries=integerMissing
    do i=1,size(series)
      if(trim(series(i)%name)==trim(name))then
        iSeries=i
        return
      endif
    enddo

  end function find_simulated_series

  ! **************************************************************************************************
  ! Record the basin value of every collected series for one model time step.
  !
  ! Basin-average variables are already one value per GRU; the rest are per HRU and domain.  Both are
  ! reduced to a single basin number the same way, as an area-weighted mean, so a target compares a
  ! basin series against a basin observation whichever structure its variable came from.  Storage and
  ! flux densities (per unit area) are what this suits, which is what the targets driving it use.
  ! **************************************************************************************************
  subroutine collect_simulated_series(modelTimeStep,summa_struct,series,err,message)
    implicit none
    integer(i4b),          intent(in)    :: modelTimeStep  ! index of the model time step
    type(summa1_type_dec), intent(in)    :: summa_struct   ! top-level SUMMA data structure
    type(sim_series_type), intent(inout) :: series(:)      ! series to record into
    integer(i4b),          intent(out)   :: err            ! error code
    character(*),          intent(out)   :: message        ! error message
    integer(i4b) :: iSeries
    integer(i4b) :: iGRU,jHRU,iDOM
    real(rkind)  :: total                                  ! area-weighted sum over the basin
    real(rkind)  :: area                                   ! total weighting area
    real(rkind)  :: wt                                     ! weight of one contribution

    err=0
    message='collect_simulated_series/'
    if(modelTimeStep < 1) return

    do iSeries=1,size(series)
      if(modelTimeStep > size(series(iSeries)%values)) cycle
      total=0._rkind
      area =0._rkind

      select case(series(iSeries)%ix_struct)

        ! basin-average variables: one value per GRU, weighted by basin area
        case (ix_series_bvar)
          do iGRU=1,summa_struct%nGRU_local
            wt=summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
            if(wt <= 0._rkind) wt=1._rkind
            total=total+wt*summa_struct%bvarStruct%gru(iGRU)%var(series(iSeries)%ix_var)%dat(1)
            area =area +wt
          enddo

        ! everything else is per HRU and domain, weighted by HRU area
        case default
          do iGRU=1,summa_struct%nGRU_local
            do jHRU=1,gru_struc(iGRU)%hruCount
              wt=summa_struct%attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
              if(wt <= 0._rkind) wt=1._rkind
              do iDOM=1,gru_struc(iGRU)%hruInfo(jHRU)%domCount
                total=total+wt*series_value(summa_struct,series(iSeries),iGRU,jHRU,iDOM)
                area =area +wt
              enddo
            enddo
          enddo

      end select

      if(area > 0._rkind)then
        series(iSeries)%values(modelTimeStep)=total/area
      else
        series(iSeries)%values(modelTimeStep)=realMissing
      endif
    enddo

  end subroutine collect_simulated_series

  ! **************************************************************************************************
  ! Read one HRU/domain value of a series from the structure it lives in.
  ! **************************************************************************************************
  function series_value(summa_struct,series,iGRU,jHRU,iDOM) result(value)
    implicit none
    type(summa1_type_dec), intent(in) :: summa_struct
    type(sim_series_type), intent(in) :: series
    integer(i4b),          intent(in) :: iGRU,jHRU,iDOM
    real(rkind)                       :: value

    select case(series%ix_struct)
      case (ix_series_diag)
        value=summa_struct%diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_prog)
        value=summa_struct%progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_flux)
        value=summa_struct%fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case (ix_series_param)
        value=summa_struct%mparStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(series%ix_var)%dat(1)
      case default
        value=realMissing
    end select

  end function series_value

end module simulated_series
