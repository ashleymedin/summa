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

module glacAreaChange_module

! data types
USE nr_type

! access missing values
USE globalData,only:integerMissing     ! missing integer number
USE globalData,only:realMissing        ! missing real number

USE globalData,only:verySmall          ! a very small number used as an additive constant to check if substantial difference among real numbers
USE globalData,only:thick4area         ! an arbitrary small threshold for glacier thickness to be considered as glacier area
USE globalData,only:dJulianStart       ! julian day of start time of simulation
USE globalData,only:data_step          ! length of time steps for the outermost timeloop

! define data types
USE var_lookup,only:iLookGRID          ! named variables for the glacier grid information
USE data_types,only:&
                    glac_info,       & ! glacier information data structure
                    grid_info,       & ! glacier grid info data structure
                    grid_double        ! x%gru(:)%grid(:)%var(:)%dat2(:,:) (dp)

USE multiconst,only:&
                    secprday,        & ! seconds per day
                    gravity,         & ! gravitational acceleration    (m s-2)
                    iden_water,      & ! intrinsic density of water    (kg m-3)
                    iden_ice           ! intrinsic density of ice      (kg m-3)

! access domain types
USE globalData,only:upland             ! domain type for upland areas
USE globalData,only:glacCln1           ! first domain type for glacier clean areas
USE globalData,only:glacCln2           ! second domain type for glacier clean areas
USE globalData,only:glacDbr            ! domain type for glacier debris areas
USE globalData,only:wetland            ! domain type for wetland areas                

implicit none
! define solver settings
character(len=32),parameter :: limiter='superbee' !'minmod'
character(len=32),parameter :: method='MUSCL' !'upwind'

! privacy
private::run_flowModel,run_debrisModel,diffusion_MUSCL,diffusion_upwind,advection_upwind
private::minmod,superbee,flux,SIA,midpt,pluss,minus,make_subgrid
public::glacAreaChange
public::time_updateGlacArea
public::updateGlacDomain
contains
! ************************************************************************************************
! public subroutine glacAreaChange: get new glacier area and elevation
! ************************************************************************************************
subroutine glacAreaChange(&
                    ! model control
                    t_total,                 & ! intent(in):    total time to run (s), time since last update
                    nHRU,                    & ! intent(in):    number of HRUs that have a glacier domain
                    nDOM,                    & ! intent(in):    number of glacier domains
                    ndebris,                 & ! intent(in):    number of debris domains in each HRU
                    nclean,                  & ! intent(in):    number of clean domains in each HRU
                    hruInd,                  & ! intent(in):    hruInd of each glacier domain
                    ! glacier topography      
                    nGlac,                   & ! intent(inout): number of glaciers in GRU
                    glacInfo,                & ! intent(inout): information for each glacier
                    gridInfo,                & ! intent(in):    information for each grid
                    gridData,                & ! intent(inout): grid data for each glacier
                    ! mass balance
                    massChange,              & ! intent(in):    rchange in glacier water equivalent (kg m-2) in each glacier domain over the nYears
                    elev,                    & ! intent(inout): elevation in each glacier domain (m)
                    ! debris
                    debris_thick_dom,        & ! intent(inout): debris thickness in each glacier domain (m)
                    iden_soil_mean,          & ! intent(in):    mean intrinsic density of soil in each glacier domain (kg m-3)
                    theta_sat_mean,          & ! intent(in):    mean soil porosity in each glacier domain (-)
                    dbr_conc,              & ! intent(in):    englacial debris concentration (kg m-3)
                    dbr_crit,                & ! intent(in):    critical debris thickness start debris-free terminal wedge (m)
                    lat_moraine_wid,         & ! intent(inout): lateral moraine width (rockfall length) (m)
                    ! area
                    glacierAblArea,          & ! intent(inout): per glacier ablation area (m2)
                    glacierAccArea,          & ! intent(inout): per glacier accumulation area (m2)
                    area,                    & ! intent(inout): area of each domain (m2)
                    ablFrac,                 & ! intent(out):   per domain ablation fraction (-)
                    ! error handling
                    err, message)              ! intent(out):   error control
   ! ---------------------------------------------------------------------------------------------
  implicit none
  ! model control
  real(rkind), intent(in)            :: t_total                         ! total time to run (s), time since last update
  integer(i4b), intent(in)           :: nHRU                            ! number of HRUs that have a glacier domain
  integer(i4b), intent(in)           :: nDOM                            ! number of glacier domains
  integer(i4b), intent(in)           :: ndebris(:)                      ! number of debris domains in each HRU
  integer(i4b), intent(in)           :: nclean(:)                       ! number of clean domains in each HRU
  integer(i8b), intent(in)           :: hruInd(:)                       ! hruInd of each glacier domain
  ! glacier topograpy
  integer(i4b), intent(inout)        :: nGlac                           ! number of glaciers in GRU
  type(glac_info), intent(inout)     :: glacInfo(:)                     ! information for each glacier
  type(grid_info), intent(in)        :: gridInfo(:)                     ! information for each grid
  type(grid_double), intent(inout)   :: gridData                        ! data for each grid
  ! mass balance, realMissing value if domain is missing (i.e. glacier does not have one of ablation or accumulation)
  real(rkind), intent(in)            :: massChange(:)                   ! change in glacier water equivalent (kg m-2) in each glacier domain
  real(rkind), intent(inout)         :: elev(:)                         ! elevation of each glacier domain (m)
  ! debris
  real(rkind), intent(inout)         :: debris_thick_dom(:)             ! debris thickness in glacier domain (m)
  real(rkind), intent(in)            :: iden_soil_mean(:)               ! mean intrinsic density of soil (kg m-3)
  real(rkind), intent(in)            :: theta_sat_mean(:)               ! mean soil porosity (-)
  real(rkind), intent(in)            :: dbr_conc                        ! englacial debris concentration (kg m-3)
  real(rkind), intent(in)            :: dbr_crit                        ! critical debris thickness start debris-free terminal wedge (m)
  real(rkind), intent(in)            :: lat_moraine_wid                 ! lateral moraine width (rockfall length) (m) 
  ! area 
  real(rkind), intent(inout)         :: glacierAblArea(nGlac)           ! per glacier ablation area (m2)
  real(rkind), intent(inout)         :: glacierAccArea(nGlac)           ! per glacier accumulation area (m2)
  real(rkind), intent(inout)         :: area(:)                         ! area of each glacier domain (m2)
  real(rkind), intent(out)           :: ablFrac(nDOM)                   ! per domain ablation fraction (-)
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message 
  ! locals
  real(rkind)                        :: elev0(nDOM)                     ! initial elevation of each glacier domain (m)
  real(rkind)                        :: area0(nDOM)                     ! initial area of each glacier domain (m2)
  real(rkind), allocatable           :: surface(:,:)                    ! surface elevation of each glacier domain (m)
  real(rkind), allocatable           :: bed(:,:)                        ! bed elevation of each glacier domain (m)
  integer(i4b), allocatable          :: cell2hru(:,:)                   ! index mapping from grid cells to HRUs
  integer(i4b), allocatable          :: glacierMask(:,:)                ! 1-0 mask of glacier domain
  real(rkind), allocatable           :: debris(:,:)                     ! debris thickness of each glacier cell (m)
  integer(i4b), allocatable          :: glacClnMask(:,:)                ! mask for clean glacier area
  integer(i4b), allocatable          :: glacHiMask(:,:)                 ! mask for high elevation clean glacier area
  integer(i4b), allocatable          :: glacLoMask(:,:)                 ! mask for low elevation clean glacier area
  integer(i4b), allocatable          :: glacDbrMask(:,:)                ! mask for debris ablation area
  integer(i4b), allocatable          :: glacAblMask(:,:)                ! mask for ablation area
  integer(i4b)                       :: i,j,k,n,iGlac,iGrid,iDOM,iHRU   ! loop indices
  integer(i4b)                       :: dbr                             ! debris loop index 0 or 1
  integer(i4b)                       :: ind                             ! indice for mass balance interpolation
  real(rkind), allocatable           :: slope(:,:), intercept(:,:)      ! slope and intercept for linear extrapolation of mass balance
  real(rkind)                        :: sum_mass                        ! sum of mass balance values for averaging
  integer(i4b)                       :: nx, ny                          ! number of grid cells in x and y directions
  real(rkind)                        :: dx, dy                          ! grid cell size in x and y directions
  real(rkind)                        :: volume                          ! volume of each glacier km3
  real(rkind)                        :: totVolume                       ! total volume of all glaciers km3, might want to send out of routine
  real(rkind)                        :: ELA_elev(2)                     ! ELA elevation for debris and clean ablation (m)
  real(rkind)                        :: ELA_use                         ! ELA used for domain calculations
  real(rkind)                        :: ELA_use_glac                    ! ELA used for current glacier (make it within glacier surface range so can have a glacier debris emergence area)
  real(rkind), allocatable           :: hgt(:,:)                        ! height of each glacier domain (m)
  integer(i4b)                       :: maxCount                        ! maximum number of valid points
  integer(i4b)                       :: validCount(2)                   ! number of valid points for each debris condition
  real(rkind)                        :: iden_soil, theta_sat            ! mean values for soil properties
  real(rkind), allocatable           :: validElev(:,:)                  ! filter out points where equal to realMissing
  real(rkind), allocatable           :: validMassChange(:,:)            ! filter out points where elev equal to realMissing
  integer(i4b), allocatable          :: glacid_to_index(:)              ! mapping array from glacier id to index in gridInfo
  integer(i4b)                       :: hiInd, loInd                    ! indices for high and low elevation domains
  real(rkind), allocatable           :: sortedElev(:)                   ! sorted elevations for clean glacier area
  integer(i4b), allocatable          :: sortedIndices(:)                ! sorted indices for clean glacier area
  real(rkind)                        :: cumulativeArea                  ! cumulative area for clean glacier area
  integer(i4b)                       :: numCells                        ! number of cells in glacier grid
  real(rkind)                        :: areaAdd                         ! area to add to each clean high glacier domain
  real(rkind),parameter              :: min_thickness=0.02_rkind        ! minimum thickness of debris cover to be considered as debris cover (m)
  logical(lgt),parameter             :: printFlag=.true.               ! flag to print details for debugging

  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='glacAreaChange/'
  totVolume = 0._rkind
  ELA_elev = realMissing
  ELA_use = realMissing
  validCount = 0_i4b
  iden_soil = 0._rkind
  theta_sat = 0._rkind

  ! Make one mass balance regression for debris 0 and less than 0 (accumulation and clean ablation)
  !  and another mass balance regression for debris 0 and greater than 0 (accumulation and debris ablation)
  do dbr = 0,1
    if(dbr==0)then
      if(sum(nclean)>0)  validCount(dbr+1) = count(elev/=realMissing .and. debris_thick_dom==0._rkind)
    else 
      if(sum(ndebris)>0) validCount(dbr+1) = count(elev/=realMissing .and. debris_thick_dom >0._rkind)
    endif
  enddo
  maxCount = maxval(validCount)
  allocate(validElev(maxCount,2),validMassChange(maxCount,2))
  allocate(slope(maxCount+1,2),intercept(maxCount+1,2))
  validElev = realMissing
  validMassChange = realMissing
  slope = 0._rkind
  intercept = 0._rkind
  elev0 = elev
  area0 = area

  ! Get elevation relationships with mass balance for each debris condition
  do dbr = 0,1
    ! Filter out points where no area and elevation is missing and debris as above
    j = 1
    if(validCount(dbr+1) == 0) cycle ! no valid points
    if(sum(nclean)==0 .and. dbr==0) cycle ! no clean, skip 
    if(sum(ndebris)==0 .and. dbr==1) cycle ! no debris, skip 
    do iDOM = 1, nDOM
      if(dbr==0)then ! no debris, dbr = 0, index = 1
        if(elev(iDOM)/=realMissing .and. debris_thick_dom(iDOM)==0._rkind)then
          validElev(j,dbr+1) = elev(iDOM)
          validMassChange(j,dbr+1) = massChange(iDOM)
          j = j + 1
        endif
      else ! debris > 0, dbr = 1, index = 2
        if(elev(iDOM)/=realMissing .and. debris_thick_dom(iDOM) >0._rkind)then
          validElev(j,dbr+1) = elev(iDOM)
          validMassChange(j,dbr+1) = massChange(iDOM)
          if(debris_thick_dom(iDOM) > 0._rkind)then ! for now, just take mean of means since is usually HRU independent
            iden_soil = iden_soil + iden_soil_mean(iDOM)/validCount(dbr+1)
            theta_sat = theta_sat + theta_sat_mean(iDOM)/validCount(dbr+1)
          endif
          j = j + 1
        endif
      endif
    enddo
  
    ! Sort the valid elevations in descending order
    do i = 1, validCount(dbr+1)-1
      do j = 1, validCount(dbr+1)-i
        if(validElev(j,dbr+1) < validElev(j+1,dbr+1))then
          call swap_real(validElev(j,dbr+1), validElev(j+1,dbr+1))
          call swap_real(validMassChange(j,dbr+1), validMassChange(j+1,dbr+1))
        endif
      enddo
    enddo

    ! Combine points with the same elevation into a single point that is the average of those points
    i = 1
    do while (i <= validCount(dbr+1))
      n = 1
      sum_mass = validMassChange(i,dbr+1)
      do j = i+1, validCount(dbr+1)
        if(validElev(i,dbr+1) == validElev(j,dbr+1))then
          sum_mass = sum_mass + validMassChange(j,dbr+1)
          n = n + 1
        else
          exit
        endif
      enddo
      if(n > 1)then
        validMassChange(i,dbr+1) = sum_mass / n
        ! Shift remaining elements to the left
        do k = i+1, validCount(dbr+1)-n
          validElev(k,dbr+1) = validElev(k+n,dbr+1)
          validMassChange(k,dbr+1) = validMassChange(k+n,dbr+1)
        enddo
        validCount(dbr+1) = validCount(dbr+1) - n + 1
      endif
      i = i + 1
    enddo
    
    ! Calculate piecewise linear regression for mass balance as a function of elevation 
    !   also find ELA elevation, assuming monotonically decreasing mass balance with elevation decrease
    !   NOTE: function should probably be capped above and below, but for now just extrapolate
    if(validCount(dbr+1) == 1)then
      slope(1,dbr+1) = 0._rkind
      intercept(1,dbr+1) = validMassChange(1,dbr+1)
      if(validMassChange(1,dbr+1) > 0.0)then
        ELA_elev(dbr+1) = -1.e6_rkind ! all domains are accumulation
      else
        ELA_elev(dbr+1) = 1.e6_rkind ! all domains are ablation
      endif
    else
      ind = 0
      do i = 1, validCount(dbr+1)-1
        slope(i+1,dbr+1)= (validMassChange(i+1,dbr+1) - validMassChange(i,dbr+1)) / (validElev(i+1,dbr+1) - validElev(i,dbr+1))
        intercept(i+1,dbr+1) = validMassChange(i,dbr+1) - slope(i+1,dbr+1) * validElev(i,dbr+1)
        if(validMassChange(i,dbr+1) >= 0._rkind) ind = i
      enddo
      if(ind == validCount(dbr+1)-1) ind = validCount(dbr+1)+1 ! all domains accumulation, extrapolate below lowest point
      if(ind == 0) ind = 1 ! all domains ablation, extrapolate above highest point
      slope(1,dbr+1) = slope(2,dbr+1)
      intercept(1,dbr+1) = intercept(2,dbr+1)
      if(slope(1,dbr+1)<0._rkind)then ! don't propogate mass balance inversions
        slope(1,dbr+1) = 1.e-6_rkind
        intercept(1,dbr+1) = validMassChange(1,dbr+1)
      endif
      slope(validCount(dbr+1)+1,dbr+1) = slope(validCount(dbr+1),dbr+1)
      intercept(validCount(dbr+1)+1,dbr+1) = intercept(validCount(dbr+1),dbr+1)
      if(slope(validCount(dbr+1)+1,dbr+1)<0._rkind)then ! don't propogate mass balance inversions
        slope(validCount(dbr+1)+1,dbr+1) = 1.e-6_rkind
        intercept(validCount(dbr+1)+1,dbr+1) = validMassChange(validCount(dbr+1),dbr+1)
      endif
      ELA_elev(dbr+1) = -intercept(ind,dbr+1) / slope(ind,dbr+1)
    endif

  enddo ! end of loop for debris elevation relationships 
  ELA_use = ELA_elev(1) ! default to clean ablation ELA
  if(ELA_use<0._rkind)then ! no clean ablation
    ELA_use = ELA_elev(2)
  endif

  ! debugging print
  if(printFlag)then
    do i = 1, nDOM
      write(*,'(a,2(1x,f7.1),1x,f5.2,1x,f8.1)') "Original domain elevation (m), area (km2), debris depth (m), mass change (kg m-2) =",&
            elev(i), area(i)*1.e-6_rkind, debris_thick_dom(i), massChange(i)
    enddo
    write(*,'(a,f6.1)') "ELA used (m) = ",ELA_use
  endif

  ! Initialize new domain vars (would have to have some glacier area to get here so okay to reset, will be recalculated)
  do i = 1,nDOM
    area(i) = 0._rkind
    ablFrac(i) = 0._rkind
    elev(i) = 0._rkind
    debris_thick_dom(i) = 0._rkind
  enddo

  ! Allocate the mapping array from glacier id to index in gridInfo
  allocate(glacid_to_index(nGlac))

  ! Populate the mapping array if glaciers exist
  do iGlac = 1, nGlac
    glacid_to_index(iGlac) = -1  ! Initialize with an invalid index
    do iGrid = 1, size(gridInfo(:)%grid_id)
      if(gridInfo(iGrid)%grid_id==glacInfo(iGlac)%glac_id )then
        glacid_to_index(iGlac) = iGrid
        exit
      endif
    enddo
  enddo

  ! run flow for each glacier
  do iGlac = 1, nGlac
    iGrid = glacid_to_index(iGlac)

    ! set up glacier grid
    ny = gridInfo(iGrid)%ny
    nx = gridInfo(iGrid)%nx
    dx = gridInfo(iGrid)%dx
    dy = gridInfo(iGrid)%dy
    
    ! set height arrays and masks
    allocate(hgt(nx,ny), glacClnMask(nx,ny), glacHiMask(nx,ny), glacLoMask(nx,ny), glacDbrMask(nx,ny), glacAblMask(nx,ny))

    ! set up grid data
    allocate(surface(nx,ny), bed(nx,ny), cell2hru(nx,ny), glacierMask(nx,ny), debris(nx,ny))
    surface = gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny)
    bed = gridData%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)
    debris = gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny)
    cell2hru = int(gridData%grid(iGrid)%var(iLookGRID%cell2hru)%dat2(1:nx,1:ny))
    glacierMask = int(gridData%grid(iGrid)%var(iLookGRID%glacierMask)%dat2(1:nx,1:ny))
    if(printFlag)then 
      volume = sum(surface - bed) * dx * dy * 1.e-9_rkind ! km3
      write(*,'(a,i2,a,4(1x,f8.3))') "<GLACIER ",iGlac, " START: ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume
    endif
    ELA_use_glac = min(ELA_use, maxval(surface)+verySmall) ! don't let ELA be above max surface elevation of glacier
    ELA_use_glac = max(ELA_use_glac, minval(surface)-verySmall) ! don't let ELA be below min surface elevation of glacier

    ! run flow model if glacier has area
    if (glacierAblArea(iGlac) + glacierAccArea(iGlac)>0._rkind) then
      call run_flowModel(t_total, debris, surface, bed, glacierMask, slope, intercept, validElev, validCount, maxCount, dbr_conc, &
                         dbr_crit, lat_moraine_wid, iden_soil, theta_sat, ELA_use_glac, nx, ny, dx, dy, volume, printFlag)
    else
      if(printFlag) write(*,'(a,i2,a)') ">GLACIER ",iGlac, " SKIP: no glacier area"
      volume = 0._rkind
      deallocate(hgt, surface, bed, cell2hru, glacierMask, debris, glacClnMask, glacHiMask, glacLoMask, glacDbrMask, glacAblMask)
      cycle
    endif

    ! update basin variables
    totVolume = totVolume+volume ! add volume to total volume, includes debris, not currently used
    hgt = surface - bed
    glacierAccArea(iGlac) = sum(merge(1_i4b, 0_i4b, hgt>thick4area .and. surface>= ELA_use_glac))*dx*dy
    glacierAblArea(iGlac) = sum(merge(1_i4b, 0_i4b, hgt>thick4area .and. surface<  ELA_use_glac))*dx*dy
    ! debugging print
    if(printFlag) write(*,'(a,i2,a,4(1x,f8.3))') ">GLACIER ",iGlac, " END  : ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume

    ! Loop through HRUs and calculate new domain areas and elevations for each HRU
    ! Order of domains will go HRU 1: clean1, clean2, debris; HRU 2: clean1, clean2, debris etc.
    ! It is possible one of the domains is missing if there is not two clean or debris cover 
    n = 0 ! initialize domain counter
    do iHRU = 1, nHRU
      ! start at last index + 1
      glacClnMask = merge(1_i4b, 0_i4b, cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris< min_thickness)
      glacDbrMask = merge(1_i4b, 0_i4b, cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris>=min_thickness)
      ! if no domain to go into, then give all area to the other domain
      if(nclean(iHRU)==0 .and. sum(glacClnMask)>0) glacDbrMask = glacDbrMask + glacClnMask ! give all area to debris domain (will reduce average domain debris thickness)
      if(ndebris(iHRU)==0 .and. sum(glacDbrMask)>0) glacClnMask = glacClnMask + glacDbrMask ! give all area to clean domain (will make debris thickness zero)

      if(nclean(iHRU)>0)then
        if(nclean(iHRU)>1)then 
          ! two clean domains, split area between the two based on elevation
          glacHiMask = 0_i4b
          glacLoMask = 0_i4b
          n = n+2 ! two clean domains, now n is last index of clean domains
          ! give higher cells to higher elevation domain
          hiInd = n
          loInd = n-1
          if(elev0(n-1)>elev0(n))then
            hiInd = n-1
            loInd = n
          endif
          ! keep area ratio the same between the two clean domains (ensures will keep 2 clean domains)
          areaAdd = area0(hiInd)/(area0(n-1)+area0(n)) * sum(glacClnMask) *dx*dy

          if(areaAdd > 0._rkind)then
            ! Flatten the surface array and sort it in descending order
            numCells = nx * ny
            allocate(sortedElev(numCells), sortedIndices(numCells))
            sortedElev = reshape(surface, [numCells])
            sortedIndices = [(i, i=1, numCells)]

            ! Sort the elevations in descending order
            do i = 1, numCells - 1
              do j = i + 1, numCells
                if(sortedElev(i) < sortedElev(j))then
                  ! Swap elevations
                  call swap_real(sortedElev(i), sortedElev(j))
                  ! Swap indices
                call swap_integer(sortedIndices(i), sortedIndices(j))
                endif
              enddo
            enddo

            ! Add cells to the mask until the cumulative area matches areaAdd
            cumulativeArea = 0._rkind
            do i = 1, numCells
              ! Get the 2D indices from the 1D index
              j = mod(sortedIndices(i) - 1, nx) + 1_i4b
              k = (sortedIndices(i) - 1) / nx + 1_i4b
            
              ! Add the cell to the mask if in clean area
              if(glacClnMask(j, k) == 1_i4b)then
                glacHiMask(j, k) = 1_i4b
                cumulativeArea = cumulativeArea + dx * dy
              endif
              if(cumulativeArea >= areaAdd) exit
            enddo
            deallocate(sortedElev, sortedIndices)
          endif

          ! Set the mask for the other domain
          glacLoMask = merge(glacClnMask, 0_i4b, glacHiMask == 0_i4b)

          ! Calculate areas, elevations, and debris thickness for 2 clean domains
          area(hiInd) = area(hiInd) + sum(glacHiMask)*dx*dy
          elev(hiInd) = elev(hiInd) + sum(surface * glacHiMask) *dx*dy
          glacAblMask = merge(glacHiMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use_glac)
          ablFrac(hiInd) = ablFrac(hiInd)+ sum(glacAblMask) *dx*dy
          debris_thick_dom(hiInd) = 0._rkind

          area(loInd) = area(loInd) + sum(glacLoMask)*dx*dy
          elev(loInd) = elev(loInd) + sum(surface * glacLoMask) *dx*dy
          glacAblMask = merge(glacLoMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use_glac)
          ablFrac(loInd) = ablFrac(loInd)+ sum(glacAblMask) *dx*dy
          debris_thick_dom(loInd) = 0._rkind
        else ! only one clean domain
          n = n+1 ! now n is last index of clean domains
          area(n) = area(n) + sum(glacClnMask)*dx*dy
          elev(n) = elev(n) + sum(surface * glacClnMask) *dx*dy
          debris_thick_dom(n) = 0._rkind
          glacAblMask = merge(glacClnMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use_glac)
          ablFrac(n) = ablFrac(n)+ sum(glacAblMask) *dx*dy
        endif
      endif

      if(ndebris(iHRU)>0)then
        n = n+1 ! currently only one debris domain possible, now n is last index of debris domains
        area(n) = area(n) + sum(glacDbrMask)*dx*dy
        elev(n) = elev(n) + sum(surface * glacDbrMask) *dx*dy
        debris_thick_dom(n) = debris_thick_dom(n) + sum(debris * glacDbrMask) *dx*dy
        glacAblMask = merge(glacDbrMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use_glac)
        ablFrac(n) = ablFrac(n)+ sum(glacAblMask) *dx*dy ! should be 1.0
      endif
    enddo

    ! update gridData and deallocate
    gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) = surface
    gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny) = debris
    deallocate(hgt, surface, bed, cell2hru, glacierMask, debris, glacClnMask, glacHiMask, glacLoMask, glacDbrMask, glacAblMask)

  enddo ! end of glacier loop

  ! Set elevations to realMissing if no area in domain
  do iDOM = 1,nDOM
    if(elev(iDOM)==0._rkind) elev(iDOM)=realMissing
    if(area(iDOM)>0._rkind)then 
      elev(iDOM) = elev(iDOM) / area(iDOM)
      ablFrac(iDOM) = ablFrac(iDOM) / area(iDOM)
      debris_thick_dom(iDOM) = debris_thick_dom(iDOM) / area(iDOM)
    else
      area(iDOM) = 0._rkind
      ablFrac(iDOM) = 0._rkind
      debris_thick_dom(iDOM) = 0._rkind
    endif
  enddo
  if(printFlag)then
    do i = 1, nDOM
      write(*,'(a,2(1x,f7.1),1x,f5.2)') "After area change domain elevation (m), area (km2), debris depth (m) =", &
            elev(i), area(i)*1.e-6_rkind, debris_thick_dom(i)
    enddo
  endif

  deallocate(glacid_to_index, validElev, validMassChange, slope, intercept)

end subroutine glacAreaChange

! ************************************************************************************************
! private subroutine run_flowModel: set up Shallow Ice Approximation (SIA) diffusion flow model 
!   and run for time period
!   Follows the implementation of Jarosch et al., 2013
! ************************************************************************************************
subroutine run_flowModel(t_total, debris, S, B, glacierMask, slope, intercept, validElev, validCount, maxCount, dbr_conc, &
                         dbr_crit, lat_moraine_wid, iden_soil, theta_sat, ELA, nx, ny, dx, dy, volume, printFlag)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: t_total, dx, dy, B(nx,ny)
  real(rkind), intent(inout) :: debris(nx,ny), S(nx,ny)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
  real(rkind), intent(in) :: dbr_conc, dbr_crit, lat_moraine_wid, iden_soil, theta_sat, ELA
  integer(i4b), intent(in) :: validCount(2), maxCount
  integer(i4b), intent(in) :: ny, nx, glacierMask(nx,ny)
  logical(lgt), intent(in) :: printFlag
  real(rkind), intent(out) :: volume
  ! Local variables
  real(rkind) :: dt, t, max_dt, min_dt, deltat, div_q(nx,ny), dt_cfl, meanS
  real(rkind) :: gamma, m_dot(nx,ny), H(nx,ny), zeroMask(nx,ny), distance(nx,ny)
  integer(i4b),parameter :: n=3 ! Glen's flow law exponent
  real(rkind),parameter :: A=2.4e-24 ! Modern Glen parameter
  real(rkind),parameter :: cfl= 0.124 ! Courant-Friedrichs-Lewy condition
  integer(i4b) :: i, j, isteps

  gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
  max_dt = 31._rkind * secprday! max timestep in seconds, a month
  min_dt = 0._rkind ! min timestep in seconds 
  t = 0._rkind
  isteps = 0 ! counter for debugging print, only print once a max_dt
  volume = 0._rkind

  ! calculation of spatial mask distance to side for lateral moraine rockfall, use same each time step for efficiency
  !  use a two-pass chamfer distance transform algorithm, approximate Euclidean distance
  if(sum(debris)>0._rkind)then
    distance = 1.e6_rkind ! initialize to large value
    zeroMask = merge(1_i4b, 0_i4b, S==B)
    ! First pass: top-left to bottom-right
    do i = 1, nx
      do j = 1, ny
        if (i > 1) distance(i, j) = min(distance(i, j), distance(i-1, j) + dx)
        if (j > 1)  distance(i, j) = min(distance(i, j), distance(i, j-1) + dy)
      end do
    end do
    ! Second pass: bottom-right to top-left
    do i = nx, 1, -1
      do j = ny, 1, -1
        if (i < nx) distance(i, j) = min(distance(i, j), distance(i+1, j) + dx)
        if (j < ny) distance(i, j) = min(distance(i, j), distance(i, j+1) + dy)
      end do
    end do
    distance = merge(distance, 0._rkind, zeroMask==1_i4b) ! set distance to zero where no glacier (may grow into these areas, include them in 0 distance)
  endif

  ! time loop
  do while (t < t_total)
    dt = t_total - t
    H = merge(S-B, 0._rkind, glacierMask==1_i4b)
    ! debugging print
    if(printFlag)then 
      ! only print once a max_dt
      if(t > isteps*max_dt)then
        isteps = isteps + 1
        if(count(H>thick4area)>0_i4b)&
            write(*,'(a,f4.2,a,2(1x,f6.1))') "  time (yr) = ", t/secprday/365.25_rkind, " mean max H (m) =", &
                  sum(H)/count(H>thick4area), maxval(H)
      endif
    endif
    if(count(H>thick4area)==0_i4b)then
      if(printFlag) write(*,'(a,f4.2,a)') "  time (yr) = ", t/secprday/365.25_rkind, " no glacier present"
      debris = 0._rkind
      S = B
      return
    endif

    ! get mass balance rate over surface (m s-1)
    m_dot = massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
    
    S = S - debris ! remove debris from glacier surface for flow calculation

    ! Select diffusion method and call a step
    select case (method)
      ! MUSCL scheme, mass conserving, but more expensive
      case ("MUSCL");   call diffusion_MUSCL(S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl)
      ! Upwinding scheme, not mass conserving, but stable and quicker
      case ("upwind"); call diffusion_upwind(S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl)
      case default; stop 'Error: method not recognized in glacier flow model'
    end select

    ! update time step
    deltat = min(dt_cfl, dt)
    if(deltat > max_dt) deltat = max_dt
    if(deltat < min_dt) deltat = min_dt
    t = t + deltat

    ! Update S
    S = S + (m_dot + div_q) * deltat
    S = merge(S, B, S > B)

    ! Check that the glacier is in boundaries, fix small violations, how small is arbitrary
    if(any((S - B) > 0._rkind .and. glacierMask==0_i4b))then
      if(any((S - B) > 10._rkind .and. glacierMask==0_i4b)) stop 'Glacier exceeds boundaries in flow model'
      S = merge(B, S, (S - B) > 0._rkind .and. glacierMask==0_i4b)
    endif

    ! Check that glacier surface is not infinite (unstable), bring down to mean glacier height
    ! This is a temporary fix, should be replaced with a more sophisticated method
    if(any((S - B) > 1.e6_rkind .and. glacierMask==1_i4b))then
      meanS = sum(merge(S, 0._rkind, glacierMask==1_i4b .and. S-B<1.e6_rkind)) / count(S-B<=1.e6_rkind)
      S = merge(S, meanS, (S - B) <= 1.e6_rkind)
    endif

    ! Update debris thickness if there is debris, using englacial debris advection transport model on sub-grid sub-step scale
    if(sum(debris)>0._rkind)then
      call run_debrisModel(S, B, debris, dbr_crit, lat_moraine_wid, gamma, n, deltat, m_dot, distance,&
                           dbr_conc, theta_sat, iden_soil, ELA, nx, ny, dx, dy) 
    endif
    S = S + debris ! add debris back to surface  
  enddo ! end of time loop

  H = merge(S-B, 0._rkind, glacierMask==1_i4b)
  ! debugging print
  if(printFlag)then 
    if(count(H>thick4area)>0_i4b)&
        write(*,'(a,f4.2,a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = ", t/secprday/365.25_rkind, " mean max H (m) =", &
              sum(H)/count(H>thick4area), maxval(H), ", mean max debris (m) =", sum(debris)/count(H>thick4area), maxval(debris)
  endif
  if(count(H>thick4area)==0_i4b)then
    if(printFlag) write(*,'(a,f4.2,a)') "  time (yr) = ", t/secprday/365.25_rkind, " no glacier present"
    debris = 0._rkind
    S = B
    return
  endif

  ! Calculate volume of glacier (includes debris)
  volume = sum(S-B) * dx * dy * 1.e-9_rkind ! km3

end subroutine run_flowModel


! ************************************************************************************************
! private subroutine diffusion_upwind: upwinding diffusion scheme for SIA flow
!   not mass conserving, but stable and quicker
! ************************************************************************************************
subroutine diffusion_upwind(S, B, mask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, n
  real(rkind), intent(out) :: div_q(nx,ny), dt_cfl
  ! Local variables
  integer(i4b) :: l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny)
  real(rkind) :: Skl(nx,ny), Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkml(nx,ny), Hkl(nx,ny), Hklp(nx,ny), Hklm(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_dn(nx,ny)
  real(rkind) :: H_l_upwind_up(nx,ny), H_l_upwind_dn(nx,ny)
  real(rkind) :: f_l_p(nx,ny), f_l_m(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_up(nx,ny), H_k_dn(nx,ny)
  real(rkind) :: H_k_upwind_up(nx,ny), H_k_upwind_dn(nx,ny)
  real(rkind) :: f_k_p(nx,ny), f_k_m(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: divisor, div_k(nx,ny), div_l(nx,ny)
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography  

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lm = [1, (i, i=1, ny-1)]
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  km = [1, (i, i=1, nx-1)]
  
  ! indices
  Sklp  = S(k ,lp)
  Sklm  = S(k ,lm)
  Skplp = S(kp,lp)
  Skplm = S(kp,lm)
  Skpl  = S(kp,l )
  Skl   = S(k ,l )
  Skmlp = S(km,lp)
  Skmlm = S(km,lm)
  Skml  = S(km,l )
  !
  Hkpl  = H(kp,l )
  Hkml  = H(km,l )
  Hkl   = H(k ,l )
  Hklp  = H(k ,lp)
  Hklm  = H(k ,lm)

  ! calculate l upstream
  H_l_up = midpt(Hklp, Hkl)
  H_l_dn = midpt(Hkl, Hklm)
  H_l_upwind_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Sklp(j,i) > Skl(j,i))then
        H_l_upwind_up(j,i) = Hklp(j,i)
      else
        H_l_upwind_up(j,i) = Hkl(j,i)
      endif
    enddo
  enddo
  H_l_upwind_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) > Sklm(j,i))then
        H_l_upwind_dn(j,i) = Hkl(j,i)
      else
        H_l_upwind_dn(j,i) = Hklm(j,i)
      endif
    enddo
  enddo

  ! calculate l flux
  f_l_p = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_m = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up = gamma * H_l_up**(n+1_i4b) * H_l_upwind_up * f_l_p
  D_l_dn = gamma * H_l_dn**(n+1_i4b) * H_l_upwind_dn * f_l_m
  ! Enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k upstream 
  H_k_up = midpt(Hkpl, Hkl)
  H_k_dn = midpt(Hkl, Hkml)
  H_k_upwind_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skpl(j,i) > Skl(j,i))then
        H_k_upwind_up(j,i) = Hkpl(j,i)
      else
        H_k_upwind_up(j,i) = Hkl(j,i)
      endif
    enddo
  enddo
  H_k_upwind_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) > Skml(j,i))then
        H_k_upwind_dn(j,i) = Hkl(j,i)
      else
        H_k_upwind_dn(j,i) = Hkml(j,i)
      endif
    enddo
  enddo

  ! calculate k flux
  f_k_p = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_m = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up = gamma * H_k_up**(n+1_i4b) * H_k_upwind_up * f_k_p
  D_k_dn = gamma * H_k_dn**(n+1_i4b) * H_k_upwind_dn * f_k_m
  ! Enforce zero diffusion outside the mask
  D_k_up = merge(0._rkind, D_k_up, mask==0_i4b)
  D_k_dn = merge(0._rkind, D_k_dn, mask==0_i4b)

  ! calculate delta t and t
  divisor = max(maxval(abs(D_k_up)), maxval(abs(D_k_dn)), maxval(abs(D_l_up)), maxval(abs(D_l_dn)))

  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dx**2_i4b, dy**2_i4b) / divisor
  endif

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy)
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx)
  div_q = div_k + div_l

end subroutine diffusion_upwind


! ************************************************************************************************
! subroutine diffusion_MUSCL: MUSCL scheme for SIA flow
!   mass conserving, more accurate, but less stable and slower
! ************************************************************************************************
subroutine diffusion_MUSCL(S, B, mask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, n
  real(rkind), intent(out) :: div_q(nx,ny), dt_cfl
  ! Local variables
  integer(i4b) :: l(ny), lp(ny), lm(ny), lpp(ny), lmm(ny), k(nx), kp(nx), km(nx), kpp(nx), kmm(nx)
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny)
  real(rkind) :: Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkppl(nx,ny), Hkml(nx,ny), Hkmml(nx,ny), Hkl(nx,ny), Hklp(nx,ny)
  real(rkind) :: Hklpp(nx,ny), Hklm(nx,ny), Hklmm(nx,ny)
  real(rkind) :: H_l_up_m(nx,ny), H_l_up_p(nx,ny)
  real(rkind) :: H_l_dn_m(nx,ny), H_l_dn_p(nx,ny)
  real(rkind) :: f_l_p(nx,ny), f_l_m(nx,ny)
  real(rkind) :: D_l_up_m(nx,ny), D_l_up_p(nx,ny), D_l_up_min(nx,ny), D_l_up_max(nx,ny)
  real(rkind) :: D_l_dn_m(nx,ny), D_l_dn_p(nx,ny), D_l_dn_min(nx,ny), D_l_dn_max(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_up_m(nx,ny), H_k_up_p(nx,ny), H_k_dn_m(nx,ny), H_k_dn_p(nx,ny)
  real(rkind) :: f_k_p(nx,ny), f_k_m(nx,ny)
  real(rkind) :: D_k_up_m(nx,ny), D_k_up_p(nx,ny), D_k_up_min(nx,ny), D_k_up_max(nx,ny)
  real(rkind) :: D_k_dn_m(nx,ny), D_k_dn_p(nx,ny), D_k_dn_min(nx,ny), D_k_dn_max(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: divisor, div_k(nx,ny), div_l(nx,ny)
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lpp = [(i, i=3, ny), ny, ny]
  lm = [1, (i, i=1, ny-1)]
  lmm = [1, 1, (i, i=1, ny-2)]
  if(ny == 1)then
    lpp = [1]
    lmm = [1]
  endif
  
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  kpp = [(i, i=3, nx), nx, nx]
  km = [1, (i, i=1, nx-1)]
  kmm = [1, 1, (i, i=1, nx-2)]
  if(nx == 1)then
    kpp = [1]
    kmm = [1]
  endif

  ! indices
  Sklp  = S(k  ,lp )
  Sklm  = S(k  ,lm )
  Skplp = S(kp ,lp )
  Skplm = S(kp ,lm )
  Skpl  = S(kp ,l  )
  Skl   = S(k  ,l  )
  Skmlp = S(km ,lp )
  Skmlm = S(km ,lm )
  Skml  = S(km ,l  )
  !
  Hkpl  = H(kp ,l  )
  Hkppl = H(kpp,l  )
  Hkml  = H(km ,l  )
  Hkmml = H(kmm,l  )
  Hkl   = H(k  ,l  )
  Hklp  = H(k  ,lp )
  Hklpp = H(k  ,lpp)
  Hklm  = H(k  ,lm )
  Hklmm = H(k  ,lmm)

  ! calculate l+1/2 index
  H_l_up_m = minus(Hklm, Hkl, Hklp)
  H_l_up_p = pluss(Hkl, Hklp, Hklpp)

  ! calculate l-1/2 index
  H_l_dn_m = minus(Hklmm, Hklm, Hkl)
  H_l_dn_p = pluss(Hklm, Hkl, Hklp)

  ! calculate l flux
  f_l_p = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_m  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up_m = gamma * H_l_up_m**(n+2_i4b) * f_l_p  ! equation 30 Jarosh 2013
  D_l_up_p = gamma * H_l_up_p**(n+2_i4b) * f_l_p  ! equation 30 Jarosh 2013
  D_l_up_min = min(D_l_up_m, D_l_up_p)            ! equation 31 Jarosh 2013
  D_l_up_max = max(D_l_up_m, D_l_up_p)            ! equation 32 Jarosh 2013
  !
  D_l_dn_m = gamma * H_l_dn_m**(n+2_i4b) * f_l_m  ! equation 30 Jarosh 2013
  D_l_dn_p = gamma * H_l_dn_p**(n+2_i4b) * f_l_m  ! equation 30 Jarosh 2013
  D_l_dn_min = min(D_l_dn_m, D_l_dn_p)            ! equation 31 Jarosh 2013
  D_l_dn_max = max(D_l_dn_m, D_l_dn_p)            ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_l_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Sklp(j,i) <= Skl(j,i) .and. H_l_up_m(j,i) <= H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      else if(Sklp(j,i) <= Skl(j,i) .and. H_l_up_m(j,i) > H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_up_m(j,i) <= H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_up_m(j,i) > H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      endif
    enddo
  enddo
  D_l_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Sklm(j,i) .and. H_l_dn_m(j,i) <= H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      else if(Skl(j,i) <= Sklm(j,i) .and. H_l_dn_m(j,i) > H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_dn_m(j,i) <= H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_dn_m(j,i) > H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      endif
    enddo
  enddo
  ! enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k+1/2 index
  H_k_up_m = minus(Hkml, Hkl, Hkpl)
  H_k_up_p = pluss(Hkl, Hkpl, Hkppl)

  ! calculate k-1/2 index
  H_k_dn_m = minus(Hkmml, Hkml, Hkl)
  H_k_dn_p = pluss(Hkml, Hkl, Hkpl)

  ! calculate k flux
  f_k_p = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_m = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up_m = gamma * H_k_up_m**(n+2_i4b) * f_k_p  ! equation 30 Jarosh 2013
  D_k_up_p = gamma * H_k_up_p**(n+2_i4b) * f_k_p  ! equation 30 Jarosh 2013
  D_k_up_min = min(D_k_up_m, D_k_up_p)            ! equation 31 Jarosh 2013
  D_k_up_max = max(D_k_up_m, D_k_up_p)            ! equation 32 Jarosh 2013
  !
  D_k_dn_m = gamma * H_k_dn_m**(n+2_i4b) * f_k_m  ! equation 30 Jarosh 2013
  D_k_dn_p = gamma * H_k_dn_p**(n+2_i4b) * f_k_m  ! equation 30 Jarosh 2013
  D_k_dn_min = min(D_k_dn_m, D_k_dn_p)            ! equation 31 Jarosh 2013
  D_k_dn_max = max(D_k_dn_m, D_k_dn_p)            ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_k_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skpl(j,i) <= Skl(j,i) .and. H_k_up_m(j,i) <= H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      else if(Skpl(j,i) <= Skl(j,i) .and. H_k_up_m(j,i) > H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_up_m(j,i) <= H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_up_m(j,i) > H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      endif
    enddo
  enddo
  D_k_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Skml(j,i) .and. H_k_dn_m(j,i) <= H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      else if(Skl(j,i) <= Skml(j,i) .and. H_k_dn_m(j,i) > H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_dn_m(j,i) <= H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_dn_m(j,i) > H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      endif
    enddo
  enddo
  ! enforce zero diffusion outside the mask
  D_k_up = merge(0._rkind, D_k_up, mask==0)
  D_k_dn = merge(0._rkind, D_k_dn, mask==0)

  ! calculate delta t and t
  divisor = max(maxval(abs(D_k_up)), maxval(abs(D_k_dn)), maxval(abs(D_l_up)), maxval(abs(D_l_dn)))
  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dy**2, dx**2) / divisor
  endif

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy) ! equation 36 Jarosh 2013
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx) ! equation 36 Jarosh 2013
  div_q = div_k + div_l

end subroutine diffusion_MUSCL


! ************************************************************************************************
! private function mass balance: rate distribution after suface height changes
! ************************************************************************************************
function massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), debris(nx,ny), t_total
  integer(i4b), intent(in) :: glacierMask(nx,ny), maxCount, nx, ny, validCount(2)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
  ! Returns
  real(rkind) :: massBalance(nx,ny)
  ! Local variables
  integer(i4b) :: ind, dbr, i, j, k
  
  ! distribute mass balance over surface, using all points in GRU
  do k = 1, nx
    do j = 1, ny
      ! Set mass balance to zero if not on glacier
      if(glacierMask(k,j)==0_i4b)then
        massBalance(k,j) = 0._rkind
      else
        ! set debris
        dbr = 0
        if(debris(k,j) > 0._rkind) dbr = 1
        ! find the index of the elevation in the validElev array
        ind = 0
        do i = 1, validCount(dbr+1)
          if(S(k,j) <= validElev(i,dbr+1)) ind = i
        enddo
        if(ind == 0) ind = validCount(dbr+1)+1 ! elevation is below the lowest valid elevation, extrapolate down
        massBalance(k,j) = (slope(ind,dbr+1) * S(k,j) + intercept(ind,dbr+1))/ t_total
      endif
    enddo
  enddo

  ! convert mass balance to m s-1 from kg m-2 s-1
  massBalance = massBalance / iden_water

end function massBalance


! ************************************************************************************************
! private subroutine run_debrisModel: set up englacial debris advection transport model and run
!   for one sub-time step on a refined grid. Debris input is from rockfall and from englacial debris 
!   emergence, and output from terminus slumping and ice velocity.
!   Follows the implementation of Mayer and Licciulli (2021), from Anderson and Anderson (2016).
! ************************************************************************************************
subroutine run_debrisModel(S, B, debris, dbr_crit, lat_moraine_wid, gamma, n, t_total, m_dot, distance, &
                           dbr_conc, theta_sat, iden_soil, ELA, nx, ny, dx, dy)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dbr_crit, distance(nx,ny), lat_moraine_wid
  real(rkind), intent(in) :: gamma, t_total, m_dot(nx,ny), dbr_conc, theta_sat, iden_soil
  real(rkind), intent(in) :: ELA, dx, dy
  integer(i4b), intent(in) :: n, nx, ny
  real(rkind), intent(inout) :: debris(nx,ny)
  ! Local variables
  real(rkind),parameter :: cfl= 0.5 ! Courant-Friedrichs-Lewy condition
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skl(nx,ny), Skpl(nx,ny), Skml(nx,ny)
  real(rkind) :: S_l_up(nx,ny), S_l_dn(nx,ny), S_k_up(nx,ny), S_k_dn(nx,ny)
  real(rkind) :: slope_l(nx,ny), slope_k(nx,ny), slope(nx,ny), u(nx,ny), v(nx,ny)
  real(rkind) :: mean_slope, ELA_use, topElev, botElev
  real(rkind) :: emergenceElev, englacial_emerg(nx, ny), above_ELA_rockfall_area
  real(rkind) :: rockfall, lat_rockfall(nx, ny), drivingStress(nx,ny)
  real(rkind) :: D(nx,ny), H(nx,ny), dx2, dy2, min_dt
  real(rkind), allocatable:: u_sub(:,:), v_sub(:,:),D_sub(:,:), H_sub(:,:), div_uD(:,:)
  integer(i4b), allocatable:: mask(:,:)
  real(rkind) :: t, max_dt, dt, dt_cfl, deltat
  integer(i4b) :: l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)
  integer(i4b) :: i, j, nx2, ny2, refine=1_i4b ! refine grid for subgrid advection, make odd to have center cell

  H = S - B ! H = ice thickness
  max_dt = 7._rkind * secprday ! maximum time step of 1 week
  min_dt = 0._rkind ! min timestep in seconds 
  t = 0._rkind

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lm = [1, (i, i=1, ny-1)]
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  km = [1, (i, i=1, nx-1)]

  ! calculate glacier ice surface slope with a finite difference
  Sklp  = S(k ,lp)
  Sklm  = S(k ,lm) 
  Skl   = S(k ,l )
  S_l_up = midpt(Sklp, Skl)
  S_l_dn = midpt(Skl, Sklm)
  slope_l= (S_l_up - S_l_dn) / dy ! (m m-1)

  Skpl  = S(kp,l )
  Skml  = S(km,l )
  Skl   = S(k ,l )
  S_k_up = midpt(Skpl, Skl)
  S_k_dn = midpt(Skl, Skml)
  slope_k= (S_k_up - S_k_dn) / dx ! (m m-1)

  ! Mean slope in m/m, set to zero where no glacier (above will include some non-glacier areas on sides)
  mean_slope = sum(merge(sqrt(slope_k**2_i4b + slope_l**2_i4b), 0._rkind, S>B))/count(S>B) 
  
  ! Emergence elevation is below ELA by distance from ELA to lateral moraine width
  !  i.e. emergenceElev = ELA - (topElev - ELA - lat_moraine_wid * mean_slope) 
  ELA_use = ELA ! from mass balance input
  topElev = maxval(merge(S, 0._rkind, S>B)) ! highest elevation on glacier
  botElev = minval(merge(S, 1.e6_rkind, S>B)) ! lowest elevation on glacier
  if(ELA_use > topElev - lat_moraine_wid * mean_slope) ELA_use = topElev - lat_moraine_wid * mean_slope
  emergenceElev = 2._rkind*ELA_use - topElev + lat_moraine_wid * mean_slope ! will not be above ELA

  ! Calculate near-surface debris concentration below emergence elevation
  englacial_emerg = -merge(m_dot, 0._rkind, S<emergenceElev .and. S>B) * dbr_conc * deltat ! (kg m-2)

  ! Add rockfall along the sides of the glacier below the ELA
  ! NOTE: to make rockfall agree with constant englacial debris concentration, influx = outflux, or:
  !   rockfall(kg m-2) * above_ELA_rockfall_area(m2) =  englacial_emerg(kg m-2) * emergence_area(m2)
  above_ELA_rockfall_area = sum(merge(1_i4b, 0_i4b, distance<=lat_moraine_wid .and. S>=ELA_use .and. S>B))*dx*dy
  rockfall = sum(englacial_emerg)*dx*dy/above_ELA_rockfall_area ! (kg m-2)
  lat_rockfall = merge(rockfall, 0._rkind, distance<=lat_moraine_wid .and. S<ELA_use .and. S>=emergenceElev .and. S>B) ! kg m-2
  lat_rockfall = 0._rkind ! turn off rockfall for testing

  ! Sum for total debris influx in m
  D = (lat_rockfall + englacial_emerg)/(iden_soil*(1._rkind - theta_sat)) + debris

  ! Calculate surface velocity using SIA, will be zero where no glacier
  u = gamma * (n + 2_i4b)/(n + 1_i4b) * slope_l**n *(S-B)**(n + 1_i4b) ! (m s-1)
  v = gamma * (n + 2_i4b)/(n + 1_i4b) * slope_k**n *(S-B)**(n + 1_i4b) ! (m s-1)

  ! Make subgrid for velocity advection for stability
  nx2 =(nx-1)*refine + nx
  ny2 =(ny-1)*refine + ny
  dx2 = dx / real(refine+1, rkind)
  dy2 = dy / real(refine+1, rkind)
  allocate(u_sub(nx2, ny2), D_sub(nx2, ny2), v_sub(nx2, ny2), H_sub(nx2, ny2), mask(nx2, ny2), div_uD(nx2, ny2))
  call make_subgrid(u, v, D, H, nx, ny, refine, nx2, ny2, u_sub, v_sub, D_sub, H_sub)
  mask = merge(1_i4b, 0_i4b, H_sub>0._rkind) ! mask for glacier on subgrid

  ! ** movement with subgrid **
  do while (t < t_total)
    dt = t_total - t

    ! Select advection method and call a step
    select case (method)
      ! MUSCL scheme, mass conserving, but more expensive
      case ("MUSCL");   call advection_MUSCL(u_sub, v_sub, D_sub, mask, cfl, max_dt, nx2, ny2, dx2, dy2, div_uD, dt_cfl)
      ! Upwinding scheme, not mass conserving, but stable and quicker
      case ("upwind"); call advection_upwind(u_sub, v_sub, D_sub, mask, cfl, max_dt, nx2, ny2, dx2, dy2, div_uD, dt_cfl)
    end select

    ! update time step
    deltat = min(dt_cfl, dt)
    if(deltat > max_dt) deltat = max_dt
    if(deltat < min_dt) deltat = min_dt
    t = t + deltat

    ! update debris thickness
    D_sub = D_sub + div_uD * deltat
    D_sub = merge(D_sub, 0._rkind, D_sub > 0._rkind ) ! prevent negative debris
  enddo
 
  ! return to original grid
  call average_subgrid(D_sub, nx, ny, refine, debris)
  deallocate(u_sub, v_sub, D_sub, H_sub, mask, div_uD)

  ! remove debris downslope of where slope is too steep to hold debris, > yield stress of 80kPa following Mayer et al. (2021)
  drivingStress = iden_ice * gravity * debris * slope/1000._rkind ! driving stress in kPa
  do j = 1, ny
    do i = 1, nx
      if(drivingStress(i,j) > 80._rkind)then ! 80 kPa yield stress
        if(slope_k(i,j) > 0._rkind .and. slope_l(i,j) > 0._rkind)then
          debris(i:nx,j:ny) = 0._rkind
          exit
        elseif(slope_k(i,j) > 0._rkind .and. slope_l(i,j) < 0._rkind)then
          debris(i:nx, 1:j) = 0._rkind
          exit
        elseif(slope_l(i,j) < 0._rkind .and. slope_k(i,j) > 0._rkind)then
          debris(1:i ,j:ny) = 0._rkind
          exit
        elseif(slope_k(i,j) < 0._rkind .and. slope_k(i,j) < 0._rkind)then
          debris(1:i , 1:j) = 0._rkind
          exit
        endif
      endif
    enddo
  enddo

  ! remove debris where debris thickness is over critical value
  debris = merge(debris, 0._rkind, debris<=dbr_crit) ! Check that debris thickness is not over critical value, effectively making terminal clean ice wedge 

end subroutine run_debrisModel

! ************************************************************************************************
! private function advection_upwind: upwinding advection scheme for debris movement with ice velocity
!  produces numerical diffusion, but is stable
! ************************************************************************************************
subroutine advection_upwind(u, v, D, mask, cfl, max_dt, nx, ny, dx, dy, div_uD, dt_cfl)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: u(nx,ny), v(nx,ny), D(nx,ny), cfl, max_dt, dx, dy
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny
  real(rkind), intent(out) :: div_uD(nx,ny), dt_cfl
  ! Local variables
  real(rkind) :: uklp(nx,ny), uklm(nx,ny), ukl(nx,ny)
  real(rkind) :: vkpl(nx,ny), vkml(nx,ny), vkl(nx,ny)
  real(rkind) :: Dklp(nx,ny), Dklm(nx,ny), Dkl(nx,ny),Dkpl(nx,ny), Dkml(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny)
  real(rkind) :: divisor
  integer(i4b) :: i, j, l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lm = [1, (i, i=1, ny-1)]
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  km = [1, (i, i=1, nx-1)]

  ! indices
  uklp  = u(k ,lp)
  uklm  = u(k ,lm)
  ukl   = u(k ,l )
  vkpl  = v(kp,l )
  vkml  = v(km,l )
  vkl   = v(k ,l )
  !
  Dklp  = D(k ,lp)
  Dklm  = D(k ,lm)
  Dkl   = D(k ,l )
  Dkpl  = D(kp,l )
  Dkml  = D(km,l )

  ! Solving dD/dt + div(D*u) = 0
  !   which is u*dD/dt + u*div(D*u) = 0 assuming du/dt = 0, and then divide by u to get dD/dt = - div(D*u)
  !  Thus, we need to check on sign of u and v to determine upwind direction
  div_l = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(ukl(j,i)>0._rkind)then
        div_l(j,i) = (ukl(j,i)*Dkl(j,i) - uklm(j,i)*Dklm(j,i))/dy
      else
        div_l(j,i) = (uklp(j,i)*Dklp(j,i) - ukl(j,i)*Dkl(j,i))/dy
      endif
    enddo
  enddo
  div_k = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(vkl(j,i)>0._rkind)then
        div_k(j,i) = (vkl(j,i)*Dkl(j,i) - vkml(j,i)*Dkml(j,i))/dx
      else
        div_k(j,i) = (vkpl(j,i)*Dkpl(j,i) - vkl(j,i)*Dkl(j,i))/dx
      endif
    enddo
  enddo
  ! enforce zero advection outside the mask
  div_l = merge(0._rkind, div_l, mask==0_i4b)
  div_k = merge(0._rkind, div_k, mask==0_i4b)
  div_uD = -(div_k + div_l) ! change in debris thickness with time

  ! calculate delta t and t
  divisor = max(maxval(abs(ukl)), maxval(abs(vkl)))
  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dy, dx) / divisor
  endif

end subroutine advection_upwind


! ************************************************************************************************
! private function advection_MUSCL: MUSCL advection scheme for debris movement with ice velocity
!   less numerical diffusion than upwind scheme, but can be unstable if CFL condition is not met
! ************************************************************************************************
subroutine advection_MUSCL(u, v, D, mask, cfl, max_dt, nx, ny, dx, dy, div_uD, dt_cfl)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: u(nx,ny), v(nx,ny), D(nx,ny), cfl, max_dt, dx, dy
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny
  real(rkind), intent(out) :: div_uD(nx,ny), dt_cfl
  ! Local variables
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny), D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: u_l_up(nx,ny), u_l_dn(nx,ny)
  real(rkind) :: v_k_up(nx,ny), v_k_dn(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny)
  real(rkind) :: divisor
  integer(i4b) :: i, j, l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lm = [1, (i, i=1, ny-1)]
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  km = [1, (i, i=1, nx-1)]

  ! MUSCL slope reconstruction for D, u, and v
  D_l_up = pluss(D(k,lm), D(k,l), D(k,lp)) ! D at l+1/2
  D_l_dn = minus(D(k,lm), D(k,l), D(k,lp)) ! D at l-1/2
  D_k_up = pluss(D(km,l), D(k,l), D(kp,l)) ! D at k+1/2
  D_k_dn = minus(D(km,l), D(k,l), D(kp,l)) ! D at k-1/2
  !
  u_l_up = pluss(u(k,lm), u(k,l), u(k,lp)) ! u at l+1/2
  u_l_dn = minus(u(k,lm), u(k,l), u(k,lp)) ! u at l-1/2
  !
  v_k_up = pluss(v(km,l), v(k,l), v(kp,l)) ! v at k+1/2
  v_k_dn = minus(v(km,l), v(k,l), v(kp,l)) ! v at k-1/2

  ! Advection fluxes (upwind, but with MUSCL slopes)
  div_l = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(u_l_up(j,i) > 0._rkind) then
        div_l(j,i) = (u_l_dn(j,i)*D_l_dn(j,i) - u_l_up(j,i)*D_l_up(j,i))/dy
      else
        div_l(j,i) = (u_l_up(j,i)*D_l_up(j,i) - u_l_dn(j,i)*D_l_dn(j,i))/dy
      endif
    enddo
  enddo
  div_k = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(v(j,i) > 0._rkind) then
        div_k(j,i) = (v_k_dn(j,i)*D_k_dn(j,i) - v_k_up(j,i)*D_k_up(j,i))/dx
      else
        div_k(j,i) = (v_k_up(j,i)*D_k_up(j,i) - v_k_dn(j,i)*D_k_dn(j,i))/dx
      endif
    enddo
  enddo

  ! enforce zero advection outside the mask
  div_l = merge(0._rkind, div_l, mask==0_i4b)
  div_k = merge(0._rkind, div_k, mask==0_i4b)
  div_uD = -(div_k + div_l) ! change in debris thickness with time

  ! calculate delta t and t
  divisor = max(maxval(abs(u)), maxval(abs(v)))
  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dy, dx) / divisor
  endif

end subroutine advection_MUSCL

! ************************************************************************************************
! private functions to assist in calculations
! ************************************************************************************************
function minmod(a, b)
  implicit none
  real(rkind), intent(in) :: a(:,:), b(:,:)
  real(rkind), allocatable :: minmod(:,:)
  real(rkind), allocatable :: sign_ab(:,:)

  allocate(minmod(size(a,1), size(a,2)))
  allocate(sign_ab(size(a,1), size(a,2)))
  sign_ab = sign(1.0_rkind, a) + sign(1.0_rkind, b)
  minmod = sign_ab / 2._rkind * min(abs(a), abs(b))
end function minmod

function superbee(r)
  implicit none
  real(rkind), intent(in) :: r(:,:)
  real(rkind), allocatable :: superbee(:,:)

  allocate(superbee(size(r,1), size(r,2)))
  superbee = max(0._rkind, min(2._rkind * r, 1._rkind), min(r, 2._rkind))
end function superbee

function flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dx, dy, n)
  implicit none
  real(rkind), intent(in) :: Skpl(:,:), Skml(:,:), Skplp(:,:), Skmlp(:,:), Sklp(:,:), Skl(:,:), dx, dy
  integer(i4b), intent(in) :: n
  real(rkind), allocatable :: flux(:,:)

  allocate(flux(size(Skpl,1), size(Skpl,2)))
  flux = ( (Skpl - Skml + Skplp - Skmlp)**2_i4b / (4._rkind * dx)**2_i4b & 
                 + (Sklp - Skl)**2_i4b / dy**2_i4b )**((n - 1.0) / 2._rkind)
end function flux

function SIA(Dup, Sup, S, Ddn, Sdn, d)
  implicit none
  real(rkind), intent(in) :: Dup(:,:), Sup(:,:), S(:,:), Ddn(:,:), Sdn(:,:), d
  real(rkind), allocatable :: SIA(:,:)

  allocate(SIA(size(Dup,1), size(Dup,2)))
  SIA = (Dup * (Sup - S) / d - Ddn * (S - Sdn) / d) / d
end function SIA

function midpt(h1, h2)
  implicit none
  real(rkind), intent(in) :: h1(:,:), h2(:,:)
  real(rkind), allocatable :: midpt(:,:)

  allocate(midpt(size(h1,1), size(h1,2)))
  midpt = 0.5_rkind * (h1 + h2)
end function midpt

function pluss(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: pluss(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(pluss(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if(limiter == "minmod")then
    pluss = H - 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if(limiter == "superbee")then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(Hp - H, ones, mask)
    pluss = merge(H - 0.5_rkind * superbee(abs((H - Hm) / divisor)) * (Hp - H), H, mask)
  endif
  deallocate(mask,divisor,ones)
end function pluss

function minus(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: minus(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(minus(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if(limiter == "minmod")then
    minus = H + 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if(limiter == "superbee")then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(H - Hm, ones, mask)
    minus = merge(H + 0.5_rkind * superbee(abs((Hp - H) / divisor)) * (H - Hm), H, mask)
  endif
  deallocate(mask,divisor,ones)
end function minus

subroutine swap_real(a, b)
  real(rkind), intent(inout) :: a, b
  real(rkind) :: temp
  temp = a
  a = b
  b = temp
end subroutine swap_real

subroutine swap_integer(a, b)
  integer(i4b), intent(inout) :: a, b
  integer(i4b) :: temp
  temp = a
  a = b
  b = temp
end subroutine swap_integer

! **************************************************************************************************
! private subroutines to assist in debris advection calculations
! **************************************************************************************************
subroutine make_subgrid(u, v, D, H, nx, ny, refine, nx2, ny2, u_sub, v_sub, D_sub, H_sub)
  implicit none
  ! Arguments
  integer(i4b), intent(in) :: nx, ny, nx2, ny2, refine
  real(rkind), intent(in) :: u(nx,ny), v(nx,ny), D(nx,ny), H(nx,ny)
  real(rkind), intent(inout) :: u_sub(nx2,ny2), v_sub(nx2,ny2), D_sub(nx2,ny2), H_sub(nx2,ny2)
  ! Local variables
  integer :: i, j, k, ii, jj
  real(rkind) :: frac

  ! Fill coarse grid points
  u_sub(1:refine+1:nx2, 1:refine+1:ny2) = u(1:nx, 1:ny)
  v_sub(1:refine+1:nx2, 1:refine+1:ny2) = v(1:nx, 1:ny)
  D_sub(1:refine+1:nx2, 1:refine+1:ny2) = D(1:nx, 1:ny)
  H_sub(1:refine+1:nx2, 1:refine+1:ny2) = H(1:nx, 1:ny)
  
  ! Interpolate between coarse grid points in x-direction
  do i = 1, nx-1
    do j = 1, ny
      ii = (i-1)*(refine+1) + 1
      jj = (j-1)*(refine+1) + 1
      do k = 1, refine
        frac = real(k, rkind) / real(refine+1, rkind)
        u_sub(ii+k, jj) = (1.0_rkind-frac)*u(i,j) + frac*u(i+1,j)
        v_sub(ii+k, jj) = (1.0_rkind-frac)*v(i,j) + frac*v(i+1,j)
        D_sub(ii+k, jj) = (1.0_rkind-frac)*D(i,j) + frac*D(i+1,j)
        H_sub(ii+k, jj) = (1.0_rkind-frac)*H(i,j) + frac*H(i+1,j)
      end do
    end do
  end do
  
  ! Interpolate between coarse grid points in y-direction
  do i = 1, nx
    do j = 1, ny-1
      ii = (i-1)*(refine+1) + 1
      jj = (j-1)*(refine+1) + 1
      do k = 1, refine
        frac = real(k, rkind) / real(refine+1, rkind)
        u_sub(ii, jj+k) = (1.0_rkind-frac)*u(i,j) + frac*u(i,j+1)
        v_sub(ii, jj+k) = (1.0_rkind-frac)*v(i,j) + frac*v(i,j+1)
        D_sub(ii, jj+k) = (1.0_rkind-frac)*D(i,j) + frac*D(i,j+1)
        H_sub(ii, jj+k) = (1.0_rkind-frac)*H(i,j) + frac*H(i,j+1)
      end do
    end do
  end do
  
  ! Interpolate diagonally between coarse grid points
  do i = 1, nx-1
    do j = 1, ny-1
      ii = (i-1)*(refine+1) + 1
      jj = (j-1)*(refine+1) + 1
      do k = 1, refine
        frac = real(k, rkind) / real(refine+1, rkind)
        u_sub(ii+k, jj+k) = (1.0_rkind-frac)*u(i,j) + frac*u(i+1,j+1)
        v_sub(ii+k, jj+k) = (1.0_rkind-frac)*v(i,j) + frac*v(i+1,j+1)
        D_sub(ii+k, jj+k) = (1.0_rkind-frac)*D(i,j) + frac*D(i+1,j+1)
        H_sub(ii+k, jj+k) = (1.0_rkind-frac)*H(i,j) + frac*H(i+1,j+1)
      end do
    end do
  end do
  
  ! Fill boundaries
  u_sub(nx2, :) = u(nx, :)
  u_sub(:, ny2) = u(:, ny)
  v_sub(nx2, :) = v(nx, :)
  v_sub(:, ny2) = v(:, ny)
  D_sub(nx2, :) = D(nx, :)
  D_sub(:, ny2) = D(:, ny)
  H_sub(nx2, :) = H(nx, :)
  H_sub(:, ny2) = H(:, ny)
  
end subroutine make_subgrid


subroutine average_subgrid(D_sub, nx, ny, refine, D)
  implicit none
  ! Arguments
  integer(i4b), intent(in) :: nx, ny, refine
  real(rkind), intent(in) :: D_sub((nx-1)*refine + nx, (ny-1)*refine + ny)
  real(rkind), intent(inout) :: D(nx, ny)
  ! Local variables
  integer :: i, j, ii, jj, k, l
  real(rkind) :: sumD

  ! Average subgrid values back to coarse grid
  do i = 1, nx
    do j = 1, ny
      ii = (i-1)*(refine+1) + 1
      jj = (j-1)*(refine+1) + 1
      sumD = 0._rkind
      do k = 0, refine
        do l = 0, refine
          sumD = sumD + D_sub(ii+k, jj+l)
        end do
      end do
      D(i,j) = sumD / real((refine+1)**2_i4b, rkind)
    end do
  end do

end subroutine average_subgrid  


! ************************************************************************************************
! public subroutine time_updateGlacArea: find date to update glacier domain area, elevation, and layering 
! ************************************************************************************************
subroutine time_updateGlacArea(&
                         ! input
                         now_iyyy, now_im, now_id, now_ih, now_imin, & ! intent(in): current time
                         latitude,                   & ! intent(in): latitude of HRU (degrees, +N, -S)
                         ! output
                         updateJulDay,               & ! intent(inout): julian day of last glacier area update (fraction of day)
                         updateJulDayNext,           & ! intent(inout): julian day of next glacier area update (fraction of day)
                         updateGlacArea,             & ! intent(inout): flag to update glacier area this time step
                         sec_since_last_update,      & ! intent(out):   seconds since last glacier area update
                         ! error control
                         err, message)                 ! intent(out):   error control
  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE time_utils_module,only:compjulday                       ! convert calendar date to julian day
  USE time_utils_module,only:compcalday                       ! convert julian day to calendar date
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  integer(i4b), intent(in)        :: now_iyyy, now_im, now_id, now_ih, now_imin ! current time
  real(rkind), intent(in)         :: latitude                 ! latitude of HRU (degrees, +N, -S)
  real(rkind), intent(inout)      :: updateJulDay             ! julian day of last area update (fraction of day)
  real(rkind), intent(inout)      :: updateJulDayNext         ! julian day of next glacier area update (fraction of day)
  logical, intent(inout)          :: updateGlacArea           ! flag to update glacier area this time step
  real(rkind), intent(out)        :: sec_since_last_update    ! seconds since last glacier area update
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message 
  ! ----- local variables -------------------------------------------------------------------------------------------------
  real(rkind)                     :: currentJulDay            ! current julian day
  integer(i4b)                    :: lowMonth                 ! lowest mass balance month for updating glacier area
  integer(i4b)                    :: nextLowMonthYr           ! next year of lowMonth
  real(rkind)                     :: dsec                     ! seconds fraction of day
  integer(i4b)                    :: iyyy, im, id, ih, imin   ! year, month, day, hour, minute
  integer(i4b),parameter          :: nYears=1                 ! number of years in between glacier area updates (on first of lowMonth)
  character(LEN=256)              :: cmessage                 ! error message of downwind routine
  ! -----------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='time_updateGlacArea/'

  call compjulday(now_iyyy, now_im, now_id, now_ih, now_imin,0._rkind,  & ! input  = year, month, day, hour, minute, second 
                  currentJulDay,err,cmessage)                            ! output = julian day (fraction of day) + error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  currentJulDay = currentJulDay + data_step/secprday ! julian day at end of time step

  ! determine lowMonth based on hemisphere, based on data from Dussaillant et al 2024
  if(latitude > 25._rkind)then
    lowMonth = 10_i4b ! October for Northern Hemisphere
  elseif(latitude < -25._rkind)then
    lowMonth = 4_i4b  ! April for Southern Hemisphere
  else
    lowMonth = 1_i4b  ! January for low latitudes (between 25N and 25S)
  end if

  ! get julian day of next update  (on first of lowMonth)
  if(updateJulDay==realMissing)then ! will only be true at the start of the simulation
    updateJulDay = dJulianStart
    if(now_im>=lowMonth)then ! start of simulation so now is dJulianStart
      nextLowMonthYr = now_iyyy + nYears
    else
      nextLowMonthYr = now_iyyy + nYears - 1_i4b
    endif
    call compjulday(nextLowMonthYr, lowMonth, 1_i4b, 0_i4b, 0_i4b, 0._rkind, updateJulDayNext,err,cmessage) 
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  elseif(updateJuldayNext==realMissing)then ! will only be true at the start of the simulation
    call compcalday(updateJulDay, iyyy, im, id, ih, imin, dsec, err, cmessage) ! get calendar date of last update
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    iyyy = iyyy + nYears
    call compjulday(iyyy, im, id, ih, imin, dsec, updateJuldayNext, err, cmessage) ! next update is nYears after last update
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! determine if glacier area needs to be updated this time step
  sec_since_last_update = (currentJulDay - updateJulDay)*secprday ! seconds since last update
  if(currentJulDay>=updateJuldayNext)then ! update glacier area, reset updateJulDay
    updateGlacArea = .true.
    updateJulDay = updateJulDayNext
    nextLowMonthYr = now_iyyy + nYears ! now_im is always lowMonth here
    call compjulday(nextLowMonthYr, lowMonth, 1_i4b, 0_i4b, 0_i4b, 0._rkind, updateJuldayNext, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

end subroutine time_updateGlacArea


! ************************************************************************************************
! public subroutine updateGlacDomain: update glacier domain area, elevation, and layering 
! ************************************************************************************************
subroutine updateGlacDomain(&
                  ! input
                  iglac,               & ! intent(inout): glacier domain index
                  glac_elev,           & ! intent(in):    elevation of each glacier domain (m) per HRU
                  glac_area,           & ! intent(in):    area of each glacier domain (m2)
                  glac_ablFrac,        & ! intent(in):    fraction of glacier domain area that is ablating
                  glac_debris_thick,   & ! intent(in):    debris thickness of each glacier domain (m) per HRU
                  dom_type,            & ! intent(in):    domain type
                  nSnow,               & ! intent(in):    number of snow layers
                  nLake,               & ! intent(in):    number of lake layers
                  nSoil,               & ! intent(in):    number of soil layers
                  nGlce,               & ! intent(inout): number of glacier ice layers
                  ! output
                  glacMass4AreaChange, & ! intent(inout): mass change (kg m-2)
                  mLayerDepth,         & ! intent(inout): layer thickness (m)
                  mLayerHeight,        & ! intent(inout): layer mid-point height (m)
                  iLayerHeight,        & ! intent(inout): layer interface height (m)
                  DOMarea,             & ! intent(inout): area of each domain (m2)
                  ablFrac,             & ! intent(inout): fraction of glacier domain area that is ablating
                  DOMelev,             & ! intent(inout): elevation of each glacier domain (m) per HRU
                  ! error control
                  err, message)         ! intent(out):   error control
 ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  integer(i4b), intent(inout)     :: iglac                   ! glacier domain index
  real(rkind), intent(in)         :: glac_elev(:)             ! elevation of each glacier domain (m) per HRU
  real(rkind), intent(in)         :: glac_area(:)             ! area of each glacier domain (m2)
  real(rkind), intent(in)         :: glac_ablFrac(:)          ! fraction of glacier domain area that is ablating
  real(rkind), intent(in)         :: glac_debris_thick(:)     ! debris thickness of each glacier domain (m) per HRU
  integer(i4b), intent(in)        :: dom_type                 ! domain type
  integer(i4b), intent(in)        :: nSnow                    ! number of snow layers
  integer(i4b), intent(in)        :: nLake                    ! number of lake layers　(should be 0)
  integer(i4b), intent(in)        :: nSoil                    ! number of soil layers
  integer(i4b), intent(inout)     :: nGlce                    ! number of glacier ice layers
  real(rkind), intent(inout)      :: glacMass4AreaChange      ! mass change (kg m-2)
  real(rkind), intent(inout)      :: mLayerDepth(:)           ! layer thickness (m)
  real(rkind), intent(inout)      :: mLayerHeight(:)          ! layer mid-point height (m)
  real(rkind), intent(inout)      :: iLayerHeight(0:)          ! layer interface height (m)
  real(rkind), intent(inout)      :: DOMarea                  ! area of each domain (m2)
  real(rkind), intent(inout)      :: ablFrac                  ! fraction of glacier domain area that is ablating
  real(rkind), intent(inout)      :: DOMelev                  ! elevation of each glacier domain (m) per HRU
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message 
   ! ----- define local variables ------------------------------------------------------------------------------------------
  integer(i4b)                    :: i                        ! loop index
  real(rkind)                     :: layers_thick             ! depth of layers modifying
  real(rkind)                     :: thick_ratio              ! ratio of new layers thickness to previous thickness
  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='updateGlacDomain/'

  ! update glacier domain elevation, area, and ablating fraction and reset mass change
  DOMelev = glac_elev(iglac) ! realMissing if no area
  DOMarea = glac_area(iglac) ! may be 0
  ablFrac = glac_ablFrac(iglac) ! may be 0
  glacMass4AreaChange = 0._rkind ! reset

  ! update glacier layering
  if(dom_type==glacDbr .and. DOMarea>0._rkind)then ! thickness of average debris cover in HRU changes with debris advection
    ! for zero area, glacDbr needs to keep previous non-zero soil layering so can adjust if debris-covered glacier re-appears
    ! scale soil layer thickness with debris thickness change, keep the same number of layers 
    layers_thick = sum(mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil))
    thick_ratio = glac_debris_thick(iglac)/layers_thick ! glac_debris_thick will be >0 since there is a positive area of debris-covered glacier
    mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)  = mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio
    mLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil) = mLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio
    iLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil) = iLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio            

    ! recalculate the layer heights below soil
    do i=nSnow+nLake+nSoil+1,nSnow+nLake+nSoil+nGlce
      mLayerHeight(i) = mLayerHeight(i) + glac_debris_thick(iglac) - layers_thick
      iLayerHeight(i) = iLayerHeight(i) + glac_debris_thick(iglac) - layers_thick
    enddo

  elseif(DOMarea==0._rkind .and. nSnow>0)then ! no glacier area but still has snow layers, so make snow very thin
    ! scale snow layers so comes back with tiny snow (would prefer 0 thickness but then would have to relayer here)
    layers_thick = sum(mLayerDepth(1:nSnow))
    if(layers_thick>verySmall)then
      thick_ratio = verySmall/layers_thick
      mLayerDepth(1:nSnow) = mLayerDepth(1:nSnow)*thick_ratio
      mLayerHeight(1:nSnow) = mLayerHeight(1:nSnow)*thick_ratio
      iLayerHeight(1:nSnow) = iLayerHeight(1:nSnow)*thick_ratio          

      ! recalculate the layer heights below snow
      do i=nSnow+1,nSnow+nLake+nSoil+nGlce
        mLayerHeight(i) = mLayerHeight(i) + verySmall - layers_thick
        iLayerHeight(i) = iLayerHeight(i) + verySmall - layers_thick
      enddo
    endif

  endif ! ( if changing soil layers or snow layers)

end subroutine updateGlacDomain


end module glacAreaChange_module