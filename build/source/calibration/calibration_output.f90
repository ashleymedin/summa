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
module calibration_output_module

  USE nr_type, only: i4b,rkind,lgt
  USE netcdf

  USE data_types,       only: target_info
  USE parameter_search, only: parameter_spec

  implicit none
  private

  public :: create_calibration_output
  public :: write_calibration_output
  public :: close_calibration_output

contains

  ! **************************************************************************************************
  ! Create calibration output file.
  !
  ! Creates a calibration NetCDF file containing one variable for each parameter represented in the
  ! parameter specification and an objective-function variable. Parameter metadata are stored as
  ! variable attributes, including native bounds, parameter transformation, sampling status, and
  ! ordered-constraint information.
  !
  ! Parameter and objective values are written by global trial index along the fixed sample dimension.
  ! **************************************************************************************************
  subroutine create_calibration_output(filename,spec,nSamples,nWorkers,case,targets,ncid,ierr,message)
    implicit none
    character(*),         intent(in)  :: filename
    type(parameter_spec), intent(in)  :: spec
    integer(i4b),         intent(in)  :: nSamples
    integer(i4b),         intent(in)  :: nWorkers
    character(*),         intent(in)  :: case
    type(target_info),    intent(in)  :: targets(:)
    integer(i4b),         intent(out) :: ncid
    integer(i4b),         intent(out) :: ierr
    character(*),         intent(out) :: message
    integer(i4b) :: dim_sample,dim_time,dim_target,dim_name
    integer(i4b) :: varid_sample,varid_objective,varid_param
    integer(i4b) :: varid_worker_rank,varid_target_name
    integer(i4b) :: varid_start_time,varid_end_time
    integer(i4b), dimension(2) :: time_dims
    integer(i4b), dimension(2) :: objective_dims
    integer(i4b), dimension(2) :: name_dims
    integer(i4b) :: iTarget
    integer(i4b) :: nTarget
    character(len=len(targets(1)%name))  :: target_name
    integer(i4b) :: iParam
    integer(i4b) :: iConstraint
    integer(i4b) :: iOrdered
    integer(i4b) :: sampled
    integer(i4b) :: ordered_set
    integer(i4b) :: ordered_index
    logical(lgt) :: file_open
    logical(lgt) :: constrained
    integer(i4b) :: ierr_close

    ierr=0
    message='create_calibration_output/'
    file_open=.false.
    nTarget=size(targets)
    if(nTarget < 1)then
      message=trim(message)//'a calibration output file needs at least one calibration target'
      ierr=20; return
    endif
    netcdf_block: block

      ! create calibration output file
      ierr=nf90_create(trim(filename),NF90_CLOBBER,ncid)
      if(ierr/=nf90_noerr) exit netcdf_block
      file_open=.true.

      ! fixed sample dimension
      ierr=nf90_def_dim(ncid,'sample',nSamples,dim_sample)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! sample coordinate
      ierr=nf90_def_var(ncid,'sample',NF90_INT,(/dim_sample/),varid_sample)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_sample,'long_name','parameter trial')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_sample,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! time-component dimension
      ierr=nf90_def_dim(ncid,'time_component',8,dim_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      time_dims=(/dim_time,dim_sample/)

      ! calibration-target dimension: one objective value per target, per trial
      ierr=nf90_def_dim(ncid,'target',nTarget,dim_target)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_def_dim(ncid,'name_length',len(target_name),dim_name)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! Parameter variables and metadata
      ! -----------------------------------------------------------------------------------------------
      do iParam=1,size(spec%params)
        ierr=nf90_def_var(ncid,trim(spec%params(iParam)%name),NF90_DOUBLE, (/dim_sample/),varid_param)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'long_name', trim(spec%params(iParam)%long_name))
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'units', trim(spec%params(iParam)%units))
        if(ierr/=nf90_noerr) exit netcdf_block

        ! native physical parameter bounds
        ierr=nf90_put_att(ncid,varid_param,'lower_bound',spec%params(iParam)%lower)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(ncid,varid_param,'upper_bound',spec%params(iParam)%upper)
        if(ierr/=nf90_noerr) exit netcdf_block

        ! parameter transformation used during search
        ierr=nf90_put_att(ncid,varid_param,'transformation', trim(spec%params(iParam)%transformation))
        if(ierr/=nf90_noerr) exit netcdf_block

        ! parameter sampling status
        if(spec%params(iParam)%sampled)then
          sampled=1
        else
          sampled=0
        endif
        ierr=nf90_put_att(ncid,varid_param,'sampled',sampled)
        if(ierr/=nf90_noerr) exit netcdf_block

        ! scalar trial value supplied by the model-specific parameter specification
        ierr=nf90_put_att(ncid,varid_param,'trial_value',spec%params(iParam)%trial_value)
        if(ierr/=nf90_noerr) exit netcdf_block

        ! ---------------------------------------------------------------------------------------------
        ! Ordered-constraint metadata
        ! ---------------------------------------------------------------------------------------------
        constrained=.false.
        ordered_set=0
        ordered_index=0
        do iConstraint=1,size(spec%ordered)
          do iOrdered=1,size(spec%ordered(iConstraint)%param_index)
            if(spec%ordered(iConstraint)%param_index(iOrdered) == iParam)then
              constrained=.true.
              ordered_set=iConstraint
              ordered_index=iOrdered
              exit
            endif
          enddo
          if(constrained) exit
        enddo
        if(constrained)then
          ierr=nf90_put_att(ncid,varid_param,'ordered_set',ordered_set)
          if(ierr/=nf90_noerr) exit netcdf_block
          ierr=nf90_put_att(ncid,varid_param,'ordered_index',ordered_index)
          if(ierr/=nf90_noerr) exit netcdf_block
          ierr=nf90_put_att(ncid,varid_param,'gap_fraction', spec%ordered(ordered_set)%gap_fraction)
          if(ierr/=nf90_noerr) exit netcdf_block
        endif
      enddo

      ! -----------------------------------------------------------------------------------------------
      ! Rank/timing information
      ! -----------------------------------------------------------------------------------------------

      ! worker rank
      ierr=nf90_def_var(ncid,'worker_rank',NF90_INT,(/dim_sample/),varid_worker_rank)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_worker_rank,'long_name', 'MPI worker rank responsible for parameter trial')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_worker_rank,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter-trial start time
      ierr=nf90_def_var(ncid,'dispatch_time',NF90_INT,time_dims,varid_start_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_start_time,'long_name','parameter trial dispatch time')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_start_time,'components', &
                        'year month day UTC_offset hour minute second millisecond')
      if(ierr/=nf90_noerr) exit netcdf_block
      
      ! parameter-trial end time
      ierr=nf90_def_var(ncid,'completion_time',NF90_INT,time_dims,varid_end_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_end_time,'long_name','parameter trial completion time')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_end_time,'components', &
                        'year month day UTC_offset hour minute second millisecond')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! Objective function: one value per calibration target, for every trial
      ! -----------------------------------------------------------------------------------------------
      objective_dims=(/dim_target,dim_sample/)
      ierr=nf90_def_var(ncid,'objective',NF90_DOUBLE,objective_dims,varid_objective)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'long_name','calibration objective function')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'units','-')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! the metric and transformation behind each target's value, in target order
      ierr=nf90_put_att(ncid,varid_objective,'metric',trim(joined_field(targets,'metric')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'obs_transform',trim(joined_field(targets,'obs_transform')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'target',trim(joined_field(targets,'name')))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_objective,'variable',trim(joined_field(targets,'variable')))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! target names, so a reader can label the objectives without parsing an attribute
      name_dims=(/dim_name,dim_target/)
      ierr=nf90_def_var(ncid,'target_name',NF90_CHAR,name_dims,varid_target_name)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,varid_target_name,'long_name','name of each calibration target')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! -----------------------------------------------------------------------------------------------
      ! Global metadata
      ! -----------------------------------------------------------------------------------------------
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'title','SUMMA calibration parameter trials')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'case_name',trim(case))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'mpi_workers',nWorkers)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(ncid,NF90_GLOBAL,'parameter_space','physical')
      if(ierr/=nf90_noerr) exit netcdf_block

      ! leave define mode
      ierr=nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! target names, written once now that the file holds them
      do iTarget=1,nTarget
        target_name=targets(iTarget)%name
        ierr=nf90_put_var(ncid,varid_target_name,target_name, &
                          start=(/1,iTarget/),count=(/len(target_name),1/))
        if(ierr/=nf90_noerr) exit netcdf_block
      enddo

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      if(file_open) ierr_close=nf90_close(ncid)
      return
    endif
    ierr=0

  end subroutine create_calibration_output

  ! **************************************************************************************************
  ! Join one field of every calibration target into a comma-separated list, in target order.
  !
  ! Written as a variable attribute so a reader can see what each objective is without opening the
  ! configuration that produced it.
  ! **************************************************************************************************
  function joined_field(targets,field) result(joined)
    implicit none
    type(target_info), intent(in) :: targets(:)
    character(*),      intent(in) :: field
    character(len=:), allocatable :: joined
    integer(i4b) :: iTarget

    joined=''
    do iTarget=1,size(targets)
      if(iTarget > 1) joined=joined//', '
      select case(trim(field))
        case ('name');          joined=joined//trim(targets(iTarget)%name)
        case ('variable');      joined=joined//trim(targets(iTarget)%variable)
        case ('metric');        joined=joined//trim(targets(iTarget)%metric)
        case ('obs_transform'); joined=joined//trim(targets(iTarget)%obs_transform)
        case default;           joined=joined//'unknown'
      end select
    enddo

  end function joined_field

  ! **************************************************************************************************
  ! Write one calibration trial.
  !
  ! Writes the complete physical-space parameter override vector and the objective value of every
  ! calibration target for one parameter trial. Parameter names must correspond to variables defined
  ! when the calibration output file was created.
  ! **************************************************************************************************
  subroutine write_calibration_output(ncid,isample,worker_rank,param_names,param_values, &
                                      objective,start_time,end_time,ierr,message)
    implicit none
    integer(i4b), intent(in) :: ncid
    integer(i4b), intent(in) :: isample
    integer(i4b), intent(in) :: worker_rank
    character(*), intent(in) :: param_names(:)
    real(rkind),  intent(in) :: param_values(:)
    real(rkind),  intent(in) :: objective(:)
    integer(i4b), intent(in) :: start_time(8)
    integer(i4b), intent(in) :: end_time(8)
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: varid_sample, varid_worker_rank
    integer(i4b) :: varid_start_time,varid_end_time
    integer(i4b) :: varid_param,varid_objective
    integer(i4b) :: iParam
    integer(i4b), dimension(1) :: start1,count1
    integer(i4b), dimension(2) :: start2,count2

    ierr=0
    message='write_calibration_output/'
    if(size(param_names) /= size(param_values))then
      message=trim(message)//'parameter name and value vectors have different sizes'
      ierr=20; return
    endif
    start1=(/isample/)
    count1=(/1/)
    start2=(/1,isample/)
    count2=(/8,1/)
    netcdf_block: block

      ! sample coordinate
      ierr=nf90_inq_varid(ncid,'sample',varid_sample)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_sample,(/isample/),start=start1,count=count1)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter values
      do iParam=1,size(param_names)
        ierr=nf90_inq_varid(ncid,trim(param_names(iParam)),varid_param)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(ncid,varid_param,(/param_values(iParam)/), start=start1,count=count1)
        if(ierr/=nf90_noerr) exit netcdf_block
      enddo

      ! worker rank
      ierr=nf90_inq_varid(ncid,'worker_rank',varid_worker_rank)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_worker_rank,(/worker_rank/),start=start1,count=count1)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! parameter-trial start times
      ierr=nf90_inq_varid(ncid,'dispatch_time',varid_start_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_start_time,start_time, start=start2,count=count2)
      if(ierr/=nf90_noerr) exit netcdf_block
     
      ! parameter-trial end times
      ierr=nf90_inq_varid(ncid,'completion_time',varid_end_time)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_end_time,end_time, start=start2,count=count2)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! objective function: every target's value for this trial
      ierr=nf90_inq_varid(ncid,'objective',varid_objective)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(ncid,varid_objective,objective, &
                        start=(/1,isample/),count=(/size(objective),1/))
      if(ierr/=nf90_noerr) exit netcdf_block

    end block netcdf_block
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine write_calibration_output

  ! **************************************************************************************************
  ! Close calibration output file.
  ! **************************************************************************************************
  subroutine close_calibration_output(ncid,ierr,message)
    implicit none
    integer(i4b), intent(in)  :: ncid
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    ierr=0
    message='close_calibration_output/'
    ierr=nf90_close(ncid)
    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine close_calibration_output

end module calibration_output_module
