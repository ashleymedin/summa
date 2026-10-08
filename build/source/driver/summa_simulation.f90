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
module summa_simulation

USE nr_type, only: i4b, i8b, rkind, lgt
USE summa_type, only: config_info
USE summa_type, only: summa1_type_dec
USE summa_type, only: parallel_context_type

USE summa_init, only: summa_initialize
USE summa_init, only: summa_initStreamNetwork
USE summa_setup, only: summa_paramSetup
USE summa_restart, only: summa_readRestart
USE summa_forcing, only: summa_readForcing
USE summa_modelRun, only: summa_runPhysics
USE summa_writeOutput, only: summa_writeOutputFiles
USE globalData, only: integerMissing
USE globalData, only: realMissing
USE globalData, only: iulog

USE build_options, only: mizuroute_active
USE build_options, only: openwq_active

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling,        only: get_mizuroute_streamflow
USE mizuroute_coupling,        only: get_mizuroute_reach_index
USE mizuroute_coupling,        only: get_mizuroute_stream_temp
USE mizuroute_coupling,        only: mizuroute_has_stream_temp
USE finalize_mizuroute_module, only: finalize_mizuroute
#endif

#ifdef MODFLOW_ACTIVE
! Coupled MODFLOW 6: the calibration driver runs SUMMA directly rather than through its BMI,
! so it drives the coupling here instead of in summa_modflow6.f90.  mf6_coupling is the same
! MODFLOW side that the standalone couplers use, and summa_mf6_exchange the same SUMMA side.
USE mf6_coupling,       only: mf6_coupler_type
USE mf6_coupling,       only: mf6_prepare_run_dir
USE mf6_parameters,     only: write_mf6_parameters
USE summa_mf6_exchange, only: mf6x_hru_count
USE summa_mf6_exchange, only: mf6x_hru_longitude, mf6x_hru_latitude, mf6x_hru_elevation
USE summa_mf6_exchange, only: mf6x_soil_thickness
USE summa_mf6_exchange, only: mf6x_hru_area
USE summa_mf6_exchange, only: mf6x_root_reach, mf6x_root_zone_depth
USE summa_mf6_exchange, only: mf6x_get_drainage
USE summa_mf6_exchange, only: mf6x_put_lower_bound_head
USE summa_mf6_exchange, only: mf6x_put_aquifer_storage
USE summa_mf6_exchange, only: mf6x_put_aquifer_baseflow
USE summa_mf6_exchange, only: mf6x_put_surface_discharge
USE summa_mf6_exchange, only: mf6x_get_aquifer_transpire, mf6x_put_aquifer_transpire
USE summa_mf6_exchange, only: mf6x_put_transpire_lim_aqfr
USE summa_mf6_exchange, only: mf6x_put_infil_lim_aqfr
USE summa_mf6_exchange, only: mf6x_put_aquifer_reject
USE summa_mf6_exchange, only: mf6x_get_drainage_temp
USE summa_mf6_exchange, only: mf6x_get_base_nrg_flux
USE summa_mf6_exchange, only: mf6x_get_base_conductance
USE summa_mf6_exchange, only: mf6x_put_aquifer_temp
#endif

#ifdef OPENWQ_ACTIVE
USE summa_openwq, only: openwq_init
USE summa_openwq, only: openwq_run_time_start
USE summa_openwq, only: openwq_run_space_step
USE summa_openwq, only: openwq_run_time_end
USE summa_openwq, only: openwq_finalize
#endif

! module-level data structure to share configurations
implicit none
private

public :: run_simulation
public :: evaluate_objective
public :: evaluate_objectives
public :: n_calibration_targets
public :: scalarize_objectives
public :: oriented_objectives
public :: mf6_spinup_phase

! .true. only while the shared one-year cold-start spin-up is running.  start_modflow reads it to
! decide whether this coupled run WRITES the spun-up aquifer head field or READS it: without that,
! every parameter sample restarted MODFLOW from IC/STRT while SUMMA restarted from the spun-up
! state, so the soil column and the aquifer were equilibrated to different things (section 8.5).
logical(lgt), save :: mf6_spinup_phase = .false.

contains

  ! **************************************************************************************************
  ! Run a complete SUMMA simulation and return the simulated streamflow time series.
  ! The interface is model agnostic: model-specific initialization, parameter updates,
  ! simulation, and finalization are handled internally.
  ! **************************************************************************************************
  subroutine run_simulation(config,                 & ! SUMMA configuration structure
                            domain_parallel,        & ! MPI context for domain parallelism
                            instance_parallel,      & ! MPI context for model-instance parallelism
                            timeSim,flowSim,        & ! simulated time and streamflow
                            timeUnits,flowUnits,    & ! time and streamflow units
                            param_name,param_value, & ! parameter names and values
                            err, message)             ! error code and message
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    real(rkind), allocatable, intent(out) :: timeSim(:)
    real(rkind), allocatable, intent(out) :: flowSim(:)
    character(len=:), allocatable, intent(out) :: timeUnits
    character(len=:), allocatable, intent(out) :: flowUnits
    character(*), intent(in) :: param_name(:)
    real(rkind),  intent(in) :: param_value(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    ! locals
    type(summa1_type_dec), allocatable :: summa1_struc(:)
    integer(i4b), parameter            :: n=1
    character(len=512)                 :: cmessage
  
    err=0
    message='run_simulation/'
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating top-level summa structure'
      return
    endif
  
    ! populate domain and model-instance parallel contexts
    summa1_struc(n)%domain_parallel=domain_parallel
    summa1_struc(n)%instance_parallel=instance_parallel
    call initialize_summa(config, summa1_struc(n), param_name,param_value, err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call run_summa(summa1_struc(n), timeSim,flowSim, timeUnits,flowUnits, err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine run_simulation
  
  ! **************************************************************************************************
  ! Report the number of calibration targets, which is the number of objectives per parameter trial.
  !
  ! A configuration that names no target is a single-objective calibration of the one observed series
  ! in [observations], which the configuration reader has already expressed as one target.  A run with
  ! nothing to calibrate against reports zero.
  ! **************************************************************************************************
  pure function n_calibration_targets(config) result(nTarget)
    type(config_info), intent(in) :: config
    integer(i4b)                  :: nTarget

    nTarget = 0
    if(allocated(config%calib%targets)) nTarget = size(config%calib%targets)

  end function n_calibration_targets

  ! **************************************************************************************************
  ! Collapse the per-target objectives into the single value a single-objective search maximizes.
  !
  ! Targets are oriented before they are combined: efficiencies count as they stand and error metrics
  ! count negatively, so a larger scalar is a better fit whatever mix of metrics the targets use.  A
  ! single target of weight one therefore returns exactly the metric itself, which is what DDS
  ! maximized before calibration grew a target list.
  ! **************************************************************************************************
  pure function scalarize_objectives(config,objective) result(F)
    use metrics, only: metric_is_maximized
    type(config_info), intent(in) :: config
    real(rkind),       intent(in) :: objective(:)   ! one value per calibration target
    real(rkind)                   :: F              ! scalarized objective
    integer(i4b) :: iTarget

    F = 0._rkind
    if(.not.allocated(config%calib%targets)) return

    do iTarget=1,min(size(config%calib%targets),size(objective))
      if(metric_is_maximized(config%calib%targets(iTarget)%metric))then
        F = F + config%calib%targets(iTarget)%weight*objective(iTarget)
      else
        F = F - config%calib%targets(iTarget)%weight*objective(iTarget)
      endif
    enddo

  end function scalarize_objectives

  ! **************************************************************************************************
  ! Orient each target's objective so that smaller is better, as a Pareto search compares them:
  ! efficiencies are negated, error metrics kept. A value that is not finite is the worst possible.
  ! **************************************************************************************************
  pure function oriented_objectives(config,objective) result(f)
    use metrics, only: metric_is_maximized
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(config_info), intent(in) :: config
    real(rkind),       intent(in) :: objective(:)   ! one value per calibration target
    real(rkind)                   :: f(size(objective))
    integer(i4b) :: iTarget

    f = objective
    if(.not.allocated(config%calib%targets)) return

    do iTarget=1,min(size(config%calib%targets),size(objective))
      if(metric_is_maximized(config%calib%targets(iTarget)%metric)) f(iTarget) = -objective(iTarget)
      if(.not.ieee_is_finite(f(iTarget))) f(iTarget) = 0.1_rkind*huge(1._rkind)
    enddo

  end function oriented_objectives

  ! **************************************************************************************************
  ! Evaluate the objective function for a specified parameter vector.
  !
  ! Retained for callers that want one number: it evaluates every calibration target and returns the
  ! scalarized objective.  With a single target, which is what a configuration naming no target has,
  ! that number is the metric itself.
  ! **************************************************************************************************
  subroutine evaluate_objective(config,                            & ! SUMMA configuration structure
                                domain_parallel,                   & ! MPI context for domain parallelism
                                instance_parallel,                 & ! MPI context for model-instance parallelism
                                sample_id, param_name,param_value, & ! sample ID + parameter names and values
                                metric,                            & ! objective function value
                                err, message)                        ! error code and message
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    integer(i4b), intent(in)  :: sample_id
    character(*), intent(in)  :: param_name(:)
    real(rkind),  intent(in)  :: param_value(:)
    real(rkind),  intent(out) :: metric
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    ! locals
    real(rkind), allocatable :: objective(:)   ! one value per calibration target
    integer(i4b)             :: nTarget        ! number of calibration targets
    character(len=1024)      :: cmessage       ! error message of downwind routine

    err=0
    message='evaluate_objective/'
    metric=realMissing

    nTarget=max(n_calibration_targets(config),1)
    allocate(objective(nTarget),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the objective vector'
      return
    endif

    call evaluate_objectives(config,                              &
                             domain_parallel,instance_parallel,   &
                             sample_id, param_name,param_value,   &
                             objective,                           &
                             err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    metric=scalarize_objectives(config,objective)

  end subroutine evaluate_objective

  ! **************************************************************************************************
  ! Return the position of a variable over a spatial unit in a list, or integerMissing when it is not
  ! there.
  !
  ! Two targets scoring the same variable over the same unit - one storage series against two
  ! products, say - collect it once.  The same variable over two units is two series.
  ! **************************************************************************************************
  pure function find_series(names,units,unitIds,name,unitKind,unitId) result(iName)
    use simulated_series, only: ix_unit_domain
    character(*), intent(in) :: names(:)
    integer(i4b), intent(in) :: units(:)
    integer(i8b), intent(in) :: unitIds(:)
    character(*), intent(in) :: name
    integer(i4b), intent(in) :: unitKind
    integer(i8b), intent(in) :: unitId
    integer(i4b)             :: iName
    integer(i4b) :: i

    iName=integerMissing
    do i=1,size(names)
      if(trim(names(i))/=trim(name)) cycle
      if(units(i)/=unitKind) cycle
      if(unitKind/=ix_unit_domain .and. unitIds(i)/=unitId) cycle
      iName=i
      return
    enddo

  end function find_series

  ! **************************************************************************************************
  ! Evaluate every calibration target for a specified parameter vector.
  !
  ! The routine initializes SUMMA, runs the model once, and then scores that one simulation against
  ! each configured target in turn, returning one objective value per target.  Running the model once
  ! for all targets is the point: the targets are different views of the same simulation, so they must
  ! come from the same run to be comparable.
  !
  ! Objectives are returned as the metrics themselves, in the order the targets were configured.  It
  ! is the caller that decides what to do with them - maximize a weighted sum, or rank a population by
  ! Pareto dominance - so nothing here assumes which way any of them points.
  ! **************************************************************************************************
  subroutine evaluate_objectives(config,                            & ! SUMMA configuration structure
                                 domain_parallel,                   & ! MPI context for domain parallelism
                                 instance_parallel,                 & ! MPI context for model-instance parallelism
                                 sample_id, param_name,param_value, & ! sample ID + parameter names and values
                                 objective,                         & ! objective value for each target
                                 err, message,                      & ! error code and message
                                 trial_failed)                        ! .true. if the error came from the parameter set's own run
    use iso_fortran_env, only: output_unit, error_unit
    use globalData, only: ncid
    USE globalData, only: output_fileSuffix
    use var_lookup, only: iLookFREQ
    use globalData,              only: numtim
    use read_flowobs_module,     only: read_observations
    use timeseries_alignment,    only: align_timeseries
    use metrics,                 only: compute_metric
    use write_evaluation_module, only: write_evaluation
    use series_transform,        only: accumulate_series
    use series_transform,        only: remove_baseline_mean
    use series_transform,        only: balance_between_extremes
    use simulated_series,        only: sim_series_type
    use simulated_series,        only: init_simulated_series
    use simulated_series,        only: find_simulated_series
    use simulated_series,        only: is_routed_streamflow
    use simulated_series,        only: spatial_unit_index
    use simulated_series,        only: ix_unit_domain, ix_unit_gru, ix_unit_reach
    use mf6_observations,        only: is_mf6_obs
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    integer(i4b), intent(in)  :: sample_id
    character(*), intent(in)  :: param_name(:)
    real(rkind),  intent(in)  :: param_value(:)
    real(rkind),  intent(out) :: objective(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    logical(lgt), optional, intent(out) :: trial_failed
    ! locals
    logical(lgt)                       :: failed             ! the model refused or failed on this parameter set
    logical(lgt)                       :: unscored           ! no calibration target, so nothing to score
    integer(i4b)                       :: errFinal           ! error code of finalize_summa
    type(summa1_type_dec), allocatable :: summa1_struc(:)    ! top-level SUMMA data structure
    integer(i4b), parameter            :: n=1                ! number of SUMMA data structures
    integer(i4b)                       :: i                  ! looping
    integer(i4b)                       :: iTarget            ! calibration target index
    integer(i4b)                       :: nTarget            ! number of calibration targets
    character(len=4)                   :: rankString         ! include rank in the output filename
    character(len=6)                   :: sampleString       ! include sample index in the output filename
    character(len=:), allocatable      :: outputFileSuffix_orig  ! orig suffix (to restore output suffix)
    real(rkind), allocatable           :: timeSim(:)         ! simulated time
    real(rkind), allocatable           :: flowSim(:)         ! simulated streamflow
    real(rkind), allocatable           :: timeObs(:)         ! observed time
    real(rkind), allocatable           :: flowObs(:)         ! observed streamflow
    character(len=:), allocatable      :: timeSimUnits       ! simulated time units
    character(len=:), allocatable      :: flowSimUnits       ! simulated flow units
    character(len=:), allocatable      :: timeObsUnits       ! observed time units
    character(len=:), allocatable      :: flowObsUnits       ! observed flow units
    real(rkind), allocatable           :: timeAligned(:)     ! common time vector
    real(rkind), allocatable           :: flowSimAligned(:)  ! flow simulations aligned to the common time period
    real(rkind), allocatable           :: flowObsAligned(:)  ! flow observations aligned to the common time period
    type(sim_series_type), allocatable :: series(:)          ! simulated series the targets asked for
    character(len=64),     allocatable :: seriesName(:)      ! names of their variables
    integer(i4b),          allocatable :: seriesUnit(:)      ! spatial unit each is reduced over
    integer(i8b),          allocatable :: seriesUnitId(:)    ! id of that unit
    integer(i4b),          allocatable :: targetUnit(:)      ! spatial unit of each target
    integer(i4b)                       :: nSeries            ! number of them
    integer(i4b)                       :: iSeries            ! index of one of them
    real(rkind), allocatable           :: valSim(:)          ! the simulated series for the current target
    character(len=:), allocatable      :: valSimUnits        ! its units
    real(rkind), allocatable           :: valAccum(:)        ! it, integrated from a rate
    real(rkind), allocatable           :: valAnom(:)         ! a series as departures from its baseline mean
    real(rkind), allocatable           :: timeVal(:)         ! time coordinate of the simulated series for the current target
    character(len=:), allocatable      :: timeValUnits       ! its units
    integer(i4b)                       :: iMember,nMember    ! series of a target pooling several, and their number
    character(len=:), allocatable      :: vnameMember        ! observed variable of one of them
    real(rkind), allocatable           :: poolTime(:),poolSim(:),poolObs(:) ! their aligned series, end to end
    character(len=1024)                :: cmessage           ! error message of downwind routine

    err=0
    message='evaluate_objectives/'
    objective=realMissing
    failed=.false.; unscored=.false.
    if(present(trial_failed)) trial_failed=.false.

    nTarget=n_calibration_targets(config)
    if(nTarget > size(objective))then
      message=trim(message)//'the objective vector is too short for the configured calibration targets'
      err=20; return
    endif

    ! allocate top-level SUMMA structure
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating top-level summa structure'
      return
    endif

    ! populate domain and model-instance parallel contexts
    summa1_struc(n)%domain_parallel=domain_parallel
    summa1_struc(n)%instance_parallel=instance_parallel

    ! define unique output filenames for each rank and sample
    outputFileSuffix_orig=trim(output_fileSuffix)
    if(instance_parallel%size > 1)then
      write(rankString,'(I4.4)') instance_parallel%rank
      output_fileSuffix=trim(output_fileSuffix)//'_rank'//rankString
    endif
    if(sample_id > 0)then
      write(sampleString,'(I6.6)') sample_id
      output_fileSuffix=trim(output_fileSuffix)//'_sample'//sampleString
    endif

    ! the model is released after this block however it is left, so a worker can take the next trial
    trial: block

    ! initialize SUMMA; a parameter set the model refuses is a failed trial
    call initialize_summa(config, summa1_struc(n), param_name,param_value, err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); failed=.true.; exit trial; endif

    ! an objective function needs observations to compare against; without them this is an
    ! ordinary SUMMA run, so run the model and return rather than failing
    if(n_calibration_targets(summa1_struc(n)%config) == 0)then
      unscored=.true.
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, physics_failed=failed)
      if(err/=0) message=trim(message)//trim(cmessage)
      exit trial
    endif

    ! calibration run: send model chatter to stderr so stdout carries only the metric, unless a caller owns iulog
    if(iulog == output_unit) iulog = error_unit

    ! the simulated series the targets ask for, resolved before the run so an unknown variable or
    ! spatial unit is refused here rather than after a simulation has been paid for.  Routed streamflow
    ! at the network outlet is not among them: run_summa returns it for every run.
    nSeries=0
    allocate(seriesName(nTarget),seriesUnit(nTarget),seriesUnitId(nTarget),targetUnit(nTarget),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the simulated-series list'
      exit trial
    endif
    do iTarget=1,nTarget
      associate(calTarget => summa1_struc(n)%config%calib%targets(iTarget))
      targetUnit(iTarget)=spatial_unit_index(calTarget%spatial_unit)
      if(targetUnit(iTarget) == integerMissing)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for spatial unit "'// &
                trim(calTarget%spatial_unit)//'"; use gru, hru or reach'
        err=20; exit trial
      endif
      ! with mizuRoute running, streamflow comes from the network, so a GRU's outflow is named by its reach
      if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_gru .and. &
         mizuroute_active .and. summa1_struc(n)%config%use_mizuroute)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for streamflow out of a gru '// &
                'while mizuRoute routes it; name the reach instead'
        err=20; exit trial
      endif
      if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_domain) cycle
      if(is_mf6_obs(calTarget%variable)) cycle
      if(find_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                     calTarget%variable,targetUnit(iTarget),calTarget%spatial_id) == integerMissing)then
        nSeries=nSeries+1
        seriesName(nSeries)  =calTarget%variable
        seriesUnit(nSeries)  =targetUnit(iTarget)
        seriesUnitId(nSeries)=calTarget%spatial_id
      endif
      end associate
    enddo

    call init_simulated_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                               numtim,series,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! a reach belongs to the river network, so it is resolved against the mizuRoute topology
    call resolve_reach_series(summa1_struc(n),series,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! run SUMMA once; every target scores this same simulation
    if(nSeries > 0)then
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, series=series, physics_failed=failed)
    else
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, physics_failed=failed)
    endif
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! score the simulation against each calibration target
    do iTarget=1,nTarget
      associate(calTarget => summa1_struc(n)%config%calib%targets(iTarget))

      ! a modflow_obs target may pool several observations, each paired with its observed variable, and
      ! is scored on all of them together
      nMember=1
      if(is_mf6_obs(calTarget%variable)) nMember=size(calTarget%obs_name)
      do iMember=1,nMember
        vnameMember=calTarget%vname_obs
        if(allocated(calTarget%vname_pool)) vnameMember=trim(calTarget%vname_pool(iMember))

        ! the simulated series this target is compared against: routed streamflow at the network
        ! outlet, or the series collected for its variable over its spatial unit
        if(allocated(valSim)) deallocate(valSim)
        timeVal=timeSim
        timeValUnits=timeSimUnits
        if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_domain)then
          valSim=flowSim
          valSimUnits=flowSimUnits
        else if(is_mf6_obs(calTarget%variable))then
          call mf6_obs_series(summa1_struc(n), trim(calTarget%obs_name(iMember)), valSim, err, cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
            exit trial
          endif
          valSimUnits='m'
        else
          iSeries=find_simulated_series(series,calTarget%variable,targetUnit(iTarget),calTarget%spatial_id)
          if(iSeries == integerMissing)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)// &
                    '" asks for simulated variable "'//trim(calTarget%variable)//'", which was not collected'
            err=20; exit trial
          endif
          valSim=series(iSeries)%values
          valSimUnits=trim(series(iSeries)%units)
        endif

        ! read this target's observations
        if(allocated(timeObs)) deallocate(timeObs)
        if(allocated(flowObs)) deallocate(flowObs)
        call read_observations(trim(calTarget%obs_path),trim(calTarget%obs_file),trim(vnameMember), &
                               timeObs,flowObs,timeObsUnits,flowObsUnits,err,cmessage)
        if(err/=0)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          exit trial
        endif

        ! units the target states override the file's
        if(len_trim(calTarget%obs_units) > 0) flowObsUnits=trim(calTarget%obs_units)

        ! -------------------------------------------------------------------------------------------
        ! Put the two series on the same footing, as far as this target asks for
        ! -------------------------------------------------------------------------------------------

        ! integrate a simulated rate into the quantity the observations report
        if(calTarget%accumulate)then
          call accumulate_series(timeVal,valSim,timeValUnits,valAccum,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAccum,valSim)
          ! integrating kg m-2 s-1 over seconds leaves kg m-2, which is millimetres of water
          if(trim(valSimUnits)=='kg m-2 s-1') valSimUnits='mm'
        endif

        ! express both sides as departures from the same baseline, which is what an anomaly product
        ! reports and what removes the constant of integration an accumulated series carries
        if(allocated(calTarget%baseline_start) .and. allocated(calTarget%baseline_end))then
          call remove_baseline_mean(timeVal,valSim,timeValUnits,                              &
                                    calTarget%baseline_start,calTarget%baseline_end,          &
                                    valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'" (simulated): '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,valSim)

          call remove_baseline_mean(timeObs,flowObs,timeObsUnits,                             &
                                    calTarget%baseline_start,calTarget%baseline_end,          &
                                    valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'" (observed): '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,flowObs)
        endif

        ! a balance is the simulated change between seasonal extremes, one value at each observed time
        if(len_trim(calTarget%balance) > 0)then
          call balance_between_extremes(timeVal,valSim,timeValUnits,timeObs,timeObsUnits,            &
                                        calTarget%balance,calTarget%balance_window,valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,valSim)
          timeVal=timeObs
          timeValUnits=timeObsUnits
        endif

        ! align simulated and observed series
        if(allocated(timeAligned))    deallocate(timeAligned)
        if(allocated(flowSimAligned)) deallocate(flowSimAligned)
        if(allocated(flowObsAligned)) deallocate(flowObsAligned)
        call align_timeseries(timeVal,valSim,timeValUnits,valSimUnits,   &
                              timeObs,flowObs,timeObsUnits,flowObsUnits, &
                              summa1_struc(n)%config%calib%start_date,   &
                              summa1_struc(n)%config%calib%end_date,     &
                              timeAligned,flowSimAligned,flowObsAligned, &
                              err,cmessage)
        if(err/=0)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          exit trial
        endif

        if(iMember == 1)then
          call move_alloc(timeAligned,poolTime)
          call move_alloc(flowSimAligned,poolSim)
          call move_alloc(flowObsAligned,poolObs)
        else
          poolTime=[poolTime,timeAligned]
          poolSim=[poolSim,flowSimAligned]
          poolObs=[poolObs,flowObsAligned]
        endif
      enddo
      call move_alloc(poolTime,timeAligned)
      call move_alloc(poolSim,flowSimAligned)
      call move_alloc(poolObs,flowObsAligned)

      ! compute this target's metric
      call compute_metric(flowObsAligned,flowSimAligned,        &
                          calTarget%metric,calTarget%obs_transform,   &
                          objective(iTarget),err,cmessage)
      if(err/=0)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
        ! observations in the scored period that this simulation gives nothing to pair with: the parameter set failed
        failed = any(ieee_is_finite(flowObsAligned)) .and. .not.any(ieee_is_finite(flowSimAligned))
        exit trial
      endif

      ! write the aligned evaluation time series and objective value.  The SUMMA output file holds one
      ! such series, so it holds the first target's; the trials file carries every target's value.
      if(summa1_struc(n)%config%write_timeseries .and. iTarget == 1)then
        call write_evaluation(ncid(iLookFREQ%timestep),                    &
                              summa1_struc(n)%config%calib%write_aligned,  &
                              timeAligned,                                 &
                              flowObsAligned,flowSimAligned,               &
                              timeObsUnits,flowObsUnits,                   &
                              calTarget%metric,calTarget%obs_transform,          &
                              objective(iTarget),                          &
                              err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif
      endif

      end associate
    enddo

    end block trial
    if(present(trial_failed)) trial_failed=failed

    ! finalize SUMMA and release model resources, keeping the first error
    call finalize_summa(summa1_struc(n),errFinal,cmessage)
    if(err==0 .and. errFinal/=0)then; err=errFinal; message=trim(message)//trim(cmessage); endif

    ! release top-level SUMMA data structure
    ! NOTE: Deallocate here because finalize_summa operates on a single array element
    if(allocated(summa1_struc)) deallocate(summa1_struc)

    ! restore output file suffix
    output_fileSuffix=outputFileSuffix_orig
    if(err/=0 .or. unscored) return

    ! write objective function to standard output
    if(instance_parallel%size == 1)then
      write(output_unit,'(ES24.16)') objective(1)
    else
      write(output_unit,'(A,A,A,I0,A,*(1X,F12.9))') &
           'case=',trim(config%case_name),', rank=',instance_parallel%rank,', objective=', &
           (objective(iTarget), iTarget=1,nTarget)
    endif

  end subroutine evaluate_objectives

  ! ---- PRIVATE SUBROUTINES --------------------------------------------------------------------------

  ! **************************************************************************************************
  ! initialize SUMMA
  ! **************************************************************************************************
  subroutine initialize_summa(config, summa_struct, param_name, param_value, err, message)
    type(config_info),       intent(inout)    :: config
    type(summa1_type_dec)  , intent(inout)    :: summa_struct
    character(*)           , intent(in)       :: param_name(:)
    real(rkind)            , intent(in)       :: param_value(:)
    integer(i4b)           , intent(out)      :: err
    character(*)           , intent(out)      :: message
    character(len=256) :: cmessage
    character(len=64), allocatable :: summaName(:)
    real(rkind),       allocatable :: summaValue(:)

    err = 0
    message = 'initialize_summa/'

    ! declare and allocate SUMMA data structures and initialize model state
    call summa_initialize(config, summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! MODFLOW parameters go to the coupler, the rest to SUMMA
    call split_mf6_parameters(summa_struct, param_name, param_value, summaName, summaValue)

    ! initialize parameter data structures
    call summa_paramSetup(summa_struct, summaName, summaValue, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! read restart data and reset model state
    call summa_readRestart(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! stream temperature: which reach of the river network each stream HRU stands for
    ! NOTE: after the restart read, since the domain types and the attributes are known by then
    call summa_initStreamNetwork(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! initialize OpenWQ
    if(openwq_active)then
      call openwq_init(summa_struct, err)
      if(err/=0)then; message=trim(message)//'problem initializing OpenWQ'; return; endif
    endif

  end subroutine initialize_summa

  ! **************************************************************************************************
  ! Take the MODFLOW 6 parameters out of a run's parameter overrides (the caller's, or else the
  ! command line's), keeping their multipliers in summa_struct%mf6_multiplier, one if not given.
  ! **************************************************************************************************
  subroutine split_mf6_parameters(summa_struct, param_name, param_value, summaName, summaValue)
    type(summa1_type_dec),          intent(inout) :: summa_struct
    character(*),                   intent(in)    :: param_name(:)
    real(rkind),                    intent(in)    :: param_value(:)
    character(len=64), allocatable, intent(out)   :: summaName(:)
    real(rkind),       allocatable, intent(out)   :: summaValue(:)
    character(len=64), allocatable :: names(:)
    real(rkind),       allocatable :: values(:)
    logical(lgt),      allocatable :: isMf6(:)
    logical(lgt) :: fromCli
    integer(i4b) :: i, ix

    summaName  = param_name
    summaValue = param_value
    if(.not.allocated(summa_struct%config%calib%mf6_params)) return
    if(size(summa_struct%config%calib%mf6_params) == 0) return

    ! summa_paramSetup takes the command line's overrides when the caller gives none
    fromCli = size(param_name) == 0 .and. allocated(summa_struct%config%param_name)
    if(fromCli)then
      names  = summa_struct%config%param_name
      values = summa_struct%config%param_value
    else
      names  = param_name
      values = param_value
    endif

    allocate(summa_struct%mf6_multiplier(size(summa_struct%config%calib%mf6_params)))
    summa_struct%mf6_multiplier = 1._rkind
    allocate(isMf6(size(names)))
    do i=1,size(names)
      ix = findloc(summa_struct%config%calib%mf6_params(:)%name, trim(names(i)), dim=1)
      isMf6(i) = ix > 0
      if(isMf6(i)) summa_struct%mf6_multiplier(ix) = values(i)
    enddo

    summaName  = pack(names,  .not.isMf6)
    summaValue = pack(values, .not.isMf6)
    if(fromCli)then
      summa_struct%config%param_name  = summaName
      summa_struct%config%param_value = summaValue
    endif

  end subroutine split_mf6_parameters

  ! **************************************************************************************************
  ! A MODFLOW 6 observation at every model time step, from the CSV in this rank's MODFLOW directory:
  ! its last numtim rows, one per data step after any leading steady-state steps.
  ! **************************************************************************************************
  subroutine mf6_obs_series(summa_struct, obs_name, values, err, message)
    USE globalData,       only: numtim
    USE globalData,       only: data_step
    USE summaFileManager, only: OUTPUT_PATH
    USE mf6_observations, only: find_mf6_obs, read_mf6_obs
    type(summa1_type_dec),    intent(in)  :: summa_struct
    character(*),             intent(in)  :: obs_name
    real(rkind), allocatable, intent(out) :: values(:)
    integer(i4b),             intent(out) :: err
    character(*),             intent(out) :: message
    character(len=:), allocatable :: csv_file
    real(rkind),      allocatable :: times(:), rows(:)
    character(len=4)              :: rankString
    real(rkind)                   :: span
    integer(i4b)                  :: nRow
    character(len=1024)           :: cmessage

    err=0
    message='mf6_obs_series/'

    call find_mf6_obs(trim(summa_struct%config%modflow_run_dir), obs_name, csv_file, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    write(rankString,'(I4.4)') summa_struct%instance_parallel%rank
    call read_mf6_obs(trim(OUTPUT_PATH)//'modflow_rank'//rankString//'/'//csv_file, obs_name, times, rows, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    nRow=size(rows)
    if(nRow < numtim)then
      write(message,'(a,i0,a,i0,a)') trim(message)//trim(csv_file)//' holds ',nRow,' rows for ',numtim,' model time steps'
      err=20; return
    endif
    span=times(nRow)-times(nRow-numtim+1)
    if(abs(span-(numtim-1)*data_step) > 1.e-6_rkind*max(span,data_step))then
      message=trim(message)//trim(csv_file)//' does not hold one row per data step in seconds'
      err=20; return
    endif
    values=rows(nRow-numtim+1:nRow)

  end subroutine mf6_obs_series

  ! **************************************************************************************************
  ! Resolve the reach of every series scored on one, against the mizuRoute river network.
  !
  ! A reach id the network does not hold is refused here, before the run, as an unknown SUMMA variable
  ! or GRU is; a reach asked for in a run that does not route is refused the same way, and so is a
  ! reach temperature in a network no stream domain solves.
  ! **************************************************************************************************
  subroutine resolve_reach_series(summa_struct, series, err, message)
    use simulated_series, only: sim_series_type
    use simulated_series, only: ix_unit_reach
    use simulated_series, only: is_routed_stream_temp
    type(summa1_type_dec), intent(in)    :: summa_struct
    type(sim_series_type), intent(inout) :: series(:)
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message
    integer(i4b)                         :: iSeries
    character(len=256)                   :: cmessage

    err=0
    message='resolve_reach_series/'

    do iSeries=1,size(series)
      if(series(iSeries)%ix_unit /= ix_unit_reach) cycle
      if(.not.(mizuroute_active .and. summa_struct%config%use_mizuroute))then
        message=trim(message)//'a calibration target names a reach, which needs mizuRoute routing'
        err=20; return
      endif
      if(mizuroute_active)then ! build-time capability
        call get_mizuroute_reach_index(summa_struct, series(iSeries)%unit_id, series(iSeries)%ix_seg, err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
#ifdef MIZUROUTE_ACTIVE
      if(is_routed_stream_temp(series(iSeries)%name) .and. .not.mizuroute_has_stream_temp(summa_struct))then
        message=trim(message)//'a calibration target asks for "'//trim(series(iSeries)%name)//'", but no GRU has '// &
                'a stream HRU, so the network carries no temperature'
        err=20; return
      endif
#endif
    enddo

  end subroutine resolve_reach_series

  ! **************************************************************************************************
  ! run SUMMA
  ! **************************************************************************************************
  subroutine run_summa(summa_struct, timeSim,flowSim, timeUnits,flowUnits, err,message, series, physics_failed)
    USE var_lookup, only: iLookFORCE
    USE globalData, only: forc_meta
    USE globalData, only: numtim
    USE simulated_series, only: sim_series_type
    USE simulated_series, only: collect_simulated_series
    USE simulated_series, only: ix_unit_reach
    USE simulated_series, only: is_routed_stream_temp
#ifdef MODFLOW_ACTIVE
    USE globalData, only: data_step                           ! length of a SUMMA data step (s)
#endif
    ! dummy arguments
    type(summa1_type_dec), intent(inout)       :: summa_struct  ! top-level SUMMA data structure
    real(rkind), allocatable, intent(out)      :: timeSim(:)    ! simulation time
    real(rkind), allocatable, intent(out)      :: flowSim(:)    ! simulated streamflow
    character(len=:), allocatable, intent(out) :: timeUnits     ! units and reference time for simulation time
    character(len=:), allocatable, intent(out) :: flowUnits     ! units for simulated streamflow
    integer(i4b), intent(out)                  :: err           ! error code
    character(*), intent(out)                  :: message       ! error message
    ! SUMMA variables a calibration target asked for, recorded as the run proceeds
    type(sim_series_type), optional, intent(inout) :: series(:)
    ! .true. if the error came from the model physics rather than from reading or writing
    logical(lgt), optional, intent(out)            :: physics_failed
    ! locals
    integer(i4b)                               :: modelTimeStep ! index of model time step
    integer(i4b)                               :: iSeries       ! index of a simulated series
    character(len=1024)                        :: cmessage      ! error message of downwind routine
#ifdef MODFLOW_ACTIVE
    ! coupled MODFLOW 6 state, live only while summa_struct%config%use_modflow is set
    type(mf6_coupler_type)                     :: coupler       ! the MODFLOW 6 side of the coupling
    logical(lgt)                               :: coupled       ! .true. if this run is coupled to MODFLOW 6
    real, allocatable                          :: drain_hru(:)  ! per-HRU soil drainage (m s-1), SUMMA -> MODFLOW
    real, allocatable                          :: head_hru(:)   ! per-HRU prescribed head (m), MODFLOW -> SUMMA
    real, allocatable                          :: stor_hru(:)   ! per-HRU aquifer storage (m), MODFLOW -> SUMMA
    real, allocatable                          :: bflow_hru(:)  ! per-HRU aquifer baseflow (m s-1), MODFLOW -> SUMMA
    real, allocatable                          :: rej_hru(:)    ! per-HRU recharge UZF rejected (m s-1), MODFLOW -> SUMMA
    real, allocatable                          :: cond_hru(:)   ! per-HRU base conductance (s-1), SUMMA -> MODFLOW
    real, allocatable                          :: surfdis_hru(:)   ! per-HRU groundwater discharge at the surface (m s-1), MODFLOW -> SUMMA
    real, allocatable                          :: gwet_dem_hru(:)  ! per-HRU aquifer transpiration demand (m s-1), SUMMA -> MODFLOW
    real, allocatable                          :: gwet_hru(:)      ! per-HRU groundwater ET MODFLOW took (m s-1), MODFLOW -> SUMMA
    real, allocatable                          :: gwet_lim_hru(:)  ! per-HRU aquifer transpiration limiting factor (-), MODFLOW -> SUMMA
    real, allocatable                          :: infil_lim_hru(:) ! per-HRU aquifer control on the infiltrating area (-), MODFLOW -> SUMMA
    real, allocatable                          :: rtemp_hru(:)     ! per-HRU drainage temperature (K), SUMMA -> GWE
    real, allocatable                          :: atemp_hru(:)     ! per-HRU aquifer temperature (K, <= 0 unknown), GWE -> SUMMA
    real, allocatable                          :: bnrg_hru(:)      ! per-HRU conduction out the soil base (W m-2), SUMMA -> GWE
    logical(lgt)                               :: gwe_feedback  ! deepTherml = aquiferTemp: GWE sets SUMMA's aquifer temperature
    logical(lgt)                               :: atemp_known   ! atemp_hru holds a GWE temperature
    integer(i4b)                               :: errFinal      ! error code of the MODFLOW 6 shutdown
#endif

    err=0
    message='run_summa/'
    if(present(physics_failed)) physics_failed=.false.

#ifdef MODFLOW_ACTIVE
    ! start the coupled MODFLOW 6 model, if this case is configured for one
    coupled = summa_struct%config%use_modflow
    if(coupled)then
      call start_modflow(summa_struct, coupler,                        &
                         drain_hru, head_hru, stor_hru, bflow_hru,     &
                         gwet_lim_hru, infil_lim_hru, atemp_hru,       &
                         gwe_feedback, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      allocate(rej_hru(size(drain_hru)), cond_hru(size(drain_hru)), surfdis_hru(size(drain_hru)), &
               gwet_dem_hru(size(drain_hru)), gwet_hru(size(drain_hru)), rtemp_hru(size(drain_hru)), &
               bnrg_hru(size(drain_hru)), source=0.0)
      atemp_known = any(atemp_hru > 0.0)
    endif
#endif

    ! define units for time and flow
    timeUnits = trim(forc_meta(iLookFORCE%time)%varunit) ! time since reference (varies)
    flowUnits = 'm3/s' ! always in mizuroute

    ! routed streamflow time series
    allocate(timeSim(numtim), source=realMissing)
    allocate(flowSim(numtim), source=realMissing)

    ! loop through time
    do modelTimeStep=1,numtim

#ifdef MODFLOW_ACTIVE
      ! push the previous step's MODFLOW 6 state into SUMMA, before the physics reads it
      ! (explicit, one-step lag; see mf6_coupling.f90 for the full exchange)
      if(coupled)then
        if(coupler%feedback .and. (modelTimeStep > 1 .or. coupler%restarted))then
          call mf6x_put_lower_bound_head(summa_struct, head_hru)
          if(coupler%have_sy)      call mf6x_put_aquifer_storage(summa_struct, stor_hru)
          if(coupler%have_bflow)   call mf6x_put_aquifer_baseflow(summa_struct, bflow_hru)
          if(coupler%have_surfdis) call mf6x_put_surface_discharge(summa_struct, surfdis_hru)
          if(coupler%have_gwet)    call mf6x_put_aquifer_transpire(summa_struct, gwet_hru)
          if(coupler%have_evt)     call mf6x_put_transpire_lim_aqfr(summa_struct, gwet_lim_hru)
          call mf6x_put_infil_lim_aqfr(summa_struct, infil_lim_hru)
        endif
        ! every step, so a sample never starts from the last one's rejected recharge
        call mf6x_put_aquifer_reject(summa_struct, rej_hru)
        if(gwe_feedback .and. atemp_known) call mf6x_put_aquifer_temp(summa_struct, atemp_hru)
      endif
#endif

      ! read model forcing data
      call summa_readForcing(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

      ! initialize OpenWQ time step
      if(openwq_active) call openwq_run_time_start(summa_struct)

      ! run SUMMA physics and mizuRoute
      call summa_runPhysics(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(present(physics_failed)) physics_failed=.true.
        exit
      endif

      ! transfer SUMMA fluxes to OpenWQ
      if(openwq_active) call openwq_run_space_step(summa_struct)

      ! the simulation time axis, which every calibration target shares.  It comes from the forcing,
      ! not from the routing, so a target that does not need routed flow still has a time coordinate.
      timeSim(modelTimeStep) = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)

      ! save streamflow time series at the outlet, and flow or temperature at any reach a target named
      ! (unavailable when mizuRoute is not active)
#ifdef MIZUROUTE_ACTIVE
      if(mizuroute_active)then ! build-time capability
       if(summa_struct%config%use_mizuroute)then
        call get_mizuroute_streamflow(modelTimeStep, summa_struct, flowSim(modelTimeStep))
        if(present(series))then
          do iSeries=1,size(series)
            if(series(iSeries)%ix_unit /= ix_unit_reach) cycle
            if(is_routed_stream_temp(series(iSeries)%name))then
              call get_mizuroute_stream_temp(summa_struct, series(iSeries)%ix_seg, series(iSeries)%values(modelTimeStep))
            else
              call get_mizuroute_streamflow(modelTimeStep, summa_struct, series(iSeries)%values(modelTimeStep), &
                                            series(iSeries)%ix_seg)
            endif
          enddo
        endif
       endif
      endif
#endif

      ! record the SUMMA variables any calibration target asked for
      if(present(series))then
        call collect_simulated_series(modelTimeStep, summa_struct, series, err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif
      endif

      ! write the model output
      call summa_writeOutputFiles(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

#ifdef MODFLOW_ACTIVE
      ! SUMMA drainage -> MODFLOW recharge, advance MODFLOW one step, read the new water table
      ! back ready for the next iteration
      if(coupled)then
        call mf6x_get_drainage(summa_struct, drain_hru)
        if(coupler%have_evt)  call mf6x_get_aquifer_transpire(summa_struct, gwet_dem_hru)
        if(coupler%have_gwe)  call mf6x_get_drainage_temp(summa_struct, rtemp_hru)
        if(coupler%have_esl)  call mf6x_get_base_nrg_flux(summa_struct, bnrg_hru)
        if(coupler%have_base) call mf6x_get_base_conductance(summa_struct, cond_hru)
        call coupler%step(modelTimeStep, dble(data_step),           &
                          drain_hru, head_hru, stor_hru, bflow_hru, &
                          err, cmessage,                            &
                          surfdis_hru=surfdis_hru, gwet_demand_hru=gwet_dem_hru, gwet_hru=gwet_hru, &
                          gwet_lim_hru=gwet_lim_hru, rtemp_hru=rtemp_hru, atemp_hru=atemp_hru,     &
                          bnrg_hru=bnrg_hru, infil_lim_hru=infil_lim_hru, rej_hru=rej_hru, cond_hru=cond_hru)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          if(present(physics_failed)) physics_failed=.true.
          exit
        endif
        atemp_known = coupler%have_gwe .and. coupler%feedback
      endif
#endif

      ! finalize OpenWQ time step
      if(openwq_active) call openwq_run_time_end(summa_struct)
    enddo

#ifdef MODFLOW_ACTIVE
    ! shut MODFLOW down: the next parameter sample starts its own simulation, from the same
    ! aquifer initial condition this one started from
    if(coupled)then
      call coupler%finalize(errFinal, cmessage)
      if(err==0 .and. errFinal/=0)then; err=errFinal; message=trim(message)//trim(cmessage); endif
    endif
#endif

  end subroutine run_summa

#ifdef MODFLOW_ACTIVE
  ! **************************************************************************************************
  ! Start the coupled MODFLOW 6 model for this simulation.
  !
  ! Each model instance gets its own MODFLOW run directory under the case's output path, because
  ! libmf6 writes its listing and budget files into the working directory and the calibration
  ! driver runs one instance per MPI rank at the same time.  The directory is populated from the
  ! configured model directory on the first sample and reused by the rest.
  ! **************************************************************************************************
  subroutine start_modflow(summa_struct, coupler, drain_hru, head_hru, stor_hru, bflow_hru, &
                           gwet_lim_hru, infil_lim_hru, atemp_hru, gwe_feedback, err, message)
    USE globalData,       only: numtim               ! number of SUMMA data steps
    USE globalData,       only: data_step            ! length of a SUMMA data step (s)
    USE globalData,       only: model_decisions      ! SUMMA model decision structure
    USE summaFileManager, only: OUTPUT_PATH          ! path for this case's output
    USE var_lookup,       only: iLookDECISIONS       ! named indices into model_decisions
    USE mDecisions_module,only: modflowCpl           ! MODFLOW coupled groundwater parameterization
    USE mDecisions_module,only: modLatflow           ! as modflowCpl, plus lateral flow in the soil above
    USE mDecisions_module,only: prescribedHead       ! prescribed head lower boundary condition
    USE mDecisions_module,only: aquiferTempState     ! the aquifer temperature, from MODFLOW 6 GWE
    ! dummy arguments
    type(summa1_type_dec),  intent(inout) :: summa_struct
    type(mf6_coupler_type), intent(inout) :: coupler
    real, allocatable,      intent(out)   :: drain_hru(:), head_hru(:), stor_hru(:), bflow_hru(:)
    real, allocatable,      intent(out)   :: gwet_lim_hru(:), infil_lim_hru(:), atemp_hru(:)
    logical(lgt),           intent(out)   :: gwe_feedback   ! deepTherml = aquiferTemp: GWE sets the aquifer temperature
    integer(i4b),           intent(out)   :: err
    character(*),           intent(out)   :: message
    ! locals
    integer(i4b)                  :: nHRU
    double precision, allocatable :: hru_x(:), hru_y(:), hru_z(:), soil_thk(:), hru_area(:)
    double precision, allocatable :: root_reach(:), root_zone(:)
    character(len=256)            :: run_dir
    character(len=512)            :: head_file   ! per-rank spun-up MODFLOW head field (coupled restart)
    character(len=4)              :: rankString
    character(len=256)            :: cmessage

    err=0
    message='start_modflow/'

    ! the coupled-groundwater decisions must be active, exactly as for the standalone couplers
    if(model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modflowCpl .and. &
       model_decisions(iLookDECISIONS%groundwatr)%iDecision /= modLatflow)then
      message=trim(message)//'SUMMA model decision groundwatr must be "modflow" or "modLatflow" '// &
                             'when simulation.use_modflow is set'
      err=20; return
    endif
    if(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= prescribedHead)then
      message=trim(message)//'SUMMA model decision bcLowrSoiH must be "presHead" '// &
                             'when simulation.use_modflow is set'
      err=20; return
    endif

    ! per-HRU exchange buffers, in the HRU order the coupler's cell map is built in
    nHRU = mf6x_hru_count()
    allocate(drain_hru(nHRU), head_hru(nHRU), stor_hru(nHRU), bflow_hru(nHRU))
    head_hru = 0.0; stor_hru = 0.0; bflow_hru = 0.0
    allocate(gwet_lim_hru(nHRU), source=0.0)
    allocate(infil_lim_hru(nHRU), source=1.0)
    allocate(atemp_hru(nHRU), source=-1.0)
    allocate(hru_x(nHRU), hru_y(nHRU), hru_z(nHRU), soil_thk(nHRU), hru_area(nHRU), root_reach(nHRU), root_zone(nHRU))
    call mf6x_hru_longitude(summa_struct, hru_x)
    call mf6x_hru_latitude(summa_struct, hru_y)
    call mf6x_hru_elevation(summa_struct, hru_z)
    call mf6x_soil_thickness(summa_struct, soil_thk)
    call mf6x_hru_area(summa_struct, hru_area)
    call mf6x_root_reach(summa_struct, root_reach)
    call mf6x_root_zone_depth(summa_struct, root_zone)

    ! this instance's own MODFLOW directory (one per rank; sequential samples on a rank share it)
    write(rankString,'(I4.4)') summa_struct%instance_parallel%rank
    run_dir = trim(OUTPUT_PATH)//'modflow_rank'//rankString
    call mf6_prepare_run_dir(trim(summa_struct%config%modflow_run_dir), trim(run_dir), err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! this run's MODFLOW parameters, rewritten from the model directory every run so none carries over
    if(allocated(summa_struct%mf6_multiplier))then
      call write_mf6_parameters(summa_struct%config%calib%mf6_params, summa_struct%mf6_multiplier, &
                                trim(summa_struct%config%modflow_run_dir), trim(run_dir), err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    ! per-rank coupled restart: the shared spin-up writes the aquifer head field, samples read it
    head_file = trim(OUTPUT_PATH)//'modflow_spinup_heads_rank'//rankString//'.bin'

    if(mf6_spinup_phase)then
      call coupler%init(trim(summa_struct%config%modflow_config), trim(run_dir), &
                        nHRU, hru_x, hru_y, hru_z, soil_thk,                     &
                        numtim, dble(data_step), err, cmessage, hru_area=hru_area, &
                        root_reach=root_reach, root_zone=root_zone, restart_write=trim(head_file))
    else
      call coupler%init(trim(summa_struct%config%modflow_config), trim(run_dir), &
                        nHRU, hru_x, hru_y, hru_z, soil_thk,                     &
                        numtim, dble(data_step), err, cmessage, hru_area=hru_area, &
                        root_reach=root_reach, root_zone=root_zone, restart_read=trim(head_file))
    endif
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! deepTherml = aquiferTemp takes the aquifer temperature from GWE, so there must be one
    gwe_feedback = model_decisions(iLookDECISIONS%deepTherml)%iDecision == aquiferTempState
    if(gwe_feedback .and. .not. coupler%have_gwe)then
      message=trim(message)//'deepTherml = aquiferTemp needs a MODFLOW 6 GWE model; set gwe_model_name in '// &
                             trim(summa_struct%config%modflow_config)
      err=20; return
    endif

    ! a restarted sample starts from the spun-up aquifer state, not lowerBoundHead
    call coupler%restart_state(head_hru, stor_hru, gwet_lim_hru, atemp_hru, infil_lim_hru)

  end subroutine start_modflow
#endif

  ! **************************************************************************************************
  ! finalize SUMMA
  ! **************************************************************************************************
  subroutine finalize_summa(summa_struct, err, message)
    ! SUMMA global data
    use globalData, only: forcNcid                ! netcdf id for current netcdf forcing file
    use globalData, only: ncid                    ! vector of file ids of netcdf output files
    ! SUMMA buffered output structures
    use globalData, only: fullIndxSave
    use globalData, only: fullForcSave
    use globalData, only: fullProgSave
    use globalData, only: fullDiagSave
    use globalData, only: fullFluxSave
    use globalData, only: fullBvarSave
    use netcdf_util_module, only: nc_file_close   ! module to handle netcdf stuff for inputs and outputs
    type(summa1_type_dec), intent(inout) :: summa_struct
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message
    integer(i4b)                         :: iFreq
    character(len=256)                   :: cmessage

    err = 0
    message = 'finalize_summa/'

    ! deallocate SUMMA buffered output structures
    if(allocated(fullIndxSave)) deallocate(fullIndxSave)
    if(allocated(fullForcSave)) deallocate(fullForcSave)
    if(allocated(fullProgSave)) deallocate(fullProgSave)
    if(allocated(fullDiagSave)) deallocate(fullDiagSave)
    if(allocated(fullFluxSave)) deallocate(fullFluxSave)
    if(allocated(fullBvarSave)) deallocate(fullBvarSave)

    ! close NetCDF forcing file
    if(forcNcid/=integerMissing)then
      call nc_file_close(forcNcid, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      forcNcid = integerMissing
    endif

    ! close SUMMA NetCDF output files
    do iFreq=1,size(ncid)
      if(ncid(iFreq)/=integerMissing)then
        call nc_file_close(ncid(iFreq), err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ncid(iFreq) = integerMissing
      endif
    enddo

    ! report the water check of the OpenWQ coupling
    if(openwq_active) call openwq_finalize()

    ! deallocate mizuroute structures
    if(mizuroute_active)then ! build-time capability
     if(summa_struct%config%use_mizuroute) then
      call finalize_mizuroute(err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     endif
    endif

    ! more cleanup operations can be added here as required
 
    ! Allow output libraries to complete file closure
    call sleep(2)

  end subroutine finalize_summa

end module summa_simulation
