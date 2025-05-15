! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C_mask) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module glacierFlowAreaChange_module

! data types
USE nr_type

! access missing values
USE globalData,only:integerMissing     ! missing integer number
USE globalData,only:realMissing        ! missing double precision number

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
USE globalData,only:glacieret          ! domain type for glaciers considered too small for flow                    

implicit none
! define solver settings
character(len=32),parameter :: limiter='superbee' !'minmod'
character(len=32),parameter :: bed_shape='parabolic' !'trapezoid'
character(len=32),parameter :: method='MUSCL' !'upstream'

! privacy
private::run_model,diffusion_MUSCL,diffusion_upstream,minmod,superbee,flux,SIA,H_index,H_plus,H_min
public ::glacierFlowAreaChange
public ::updateGlacierDomain
contains
! ************************************************************************************************
! public subroutine glacierFlowAreaChange: flow of glacier to get new glacier area and elevation
! NOTE: This will eventually run in parallel as a program, but for now it is serial
! ************************************************************************************************
subroutine glacierFlowAreaChange(&
                    ! model control
                    t_total,                 & ! intent(in):    total time to run (s), time since last update
                    nHRU,                    & ! intent(in):    number of HRUs that have a glacier domain
                    nDOM,                    & ! intent(in):    number of glacier domains
                    ndebris,                 & ! intent(in):    number of debris domains in each HRU
                    nclean,                  & ! intent(in):    number of clean domains in each HRU
                    hruInd,                  & ! intent(in):    hruInd of each glacier domain
                    ! glacier topography      
                    nGlacier,                & ! intent(in):    number of glaciers
                    glacInfo,                & ! intent(in):    information for each glacier
                    gridInfo,                & ! intent(in):    information for each grid
                    gridData,                & ! intent(inout): grid data for each glacier
                    ! mass balance
                    massChange,              & ! intent(in):    rchange in glacier water equivalent (kg m-2) in each glacier domain over the nYears
                    elev,                    & ! intent(inout): elevation in each glacier domain (m)
                    ! debris
                    debris_thick_dom,        & ! intent(inout): debris thickness in each glacier domain (m)
                    iden_soil_mean,          & ! intent(in):    mean intrinsic density of soil in each glacier domain (kg m-3)
                    theta_sat_mean,          & ! intent(in):    mean soil porosity in each glacier domain (-)
                    C_constant,              & ! intent(in):    non-spatial concentration for debris advection (kg m-3)
                    dbr_crit,                & ! intent(in):    critical debris thickness start debris-free terminal wedge (m)
                    lat_moraine_wid,         & ! intent(inout): lateral moraine width (m)
                    ! area
                    glacAblArea,             & ! intent(out):   per glacier ablation area (m2)
                    glacAccArea,             & ! intent(out):   per glacier accumulation area (m2)
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
  integer(i4b), intent(in)           :: nGlacier                        ! number of glaciers
  type(glac_info), intent(in)        :: glacInfo(:)                     ! information for each glacier
  type(grid_info), intent(in)        :: gridInfo(:)                     ! information for each grid
  type(grid_double), intent(inout)   :: gridData                        ! data for each grid
  ! mass balance, realMissing value if domain is missing (i.e. glacier does not have one of ablation or accumulation)
  real(rkind), intent(in)            :: massChange(:)                   ! change in glacier water equivalent (kg m-2) in each glacier domain
  real(rkind), intent(inout)         :: elev(:)                         ! elevation of each glacier domain (m)
  ! debris
  real(rkind), intent(inout)         :: debris_thick_dom(:)             ! debris thickness in glacier domain (m)
  real(rkind), intent(in)            :: iden_soil_mean(:)               ! mean intrinsic density of soil (kg m-3)
  real(rkind), intent(in)            :: theta_sat_mean(:)               ! mean soil porosity (-)
  real(rkind), intent(in)            :: C_constant                      ! non-spatial concentration for debris advection (kg m-3)
  real(rkind), intent(in)            :: dbr_crit                        ! critical debris thickness start debris-free terminal wedge (m)
  real(rkind), intent(in)            :: lat_moraine_wid                 ! lateral moraine width (m) 
  ! area 
  real(rkind), intent(out)           :: glacAblArea(nGlacier)           ! per glacier ablation area (m2)
  real(rkind), intent(out)           :: glacAccArea(nGlacier)           ! per glacier accumulation area (m2)
  real(rkind), intent(inout)         :: area(:)                         ! area of each glacier domain (m2)
  real(rkind), intent(out)           :: ablFrac(nDOM)                   ! per domain ablation fraction (-)
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message 
  ! locals
  real(rkind)                        :: elev0(nDOM)                     ! initial elevation of each glacier domain (m)
  real(rkind)                        :: area0(nDOM)                     ! initial area of each glacier domain (m2)
  real(rkind), allocatable           :: surface(:,:)                    ! surface elevation of each glacier domain (m)
  real(rkind), allocatable           :: bed(:,:),bed0(:,:)              ! bed elevation of each glacier domain (m)
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
  real(rkind),parameter              :: thick4area=0.1_rkind            ! an arbitrary small threshold for glacier thickness to be considered as glacier area
  real(rkind),parameter              :: min_thickness=0.01_rkind        ! minimum thickness of debris cover to be considered as debris cover (m)

  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='glacierFlowAreaChange/'
  glacAblArea = 0._rkind
  glacAccArea = 0._rkind
  totVolume = 0._rkind
  ELA_elev = realMissing
  ELA_use = realMissing
  validCount = 0_i4b
  iden_soil = 0._rkind
  theta_sat = 0._rkind

  ! Make one mass balance regression for debris 0 and less than 0 (accumulation and clean ablation)
  !  and another mass balance regression for debris 0 and greater than 0 (accumulation and debris ablation)
  do dbr = 0,1
    if (dbr==0)then
      if (sum(nclean)>0)  validCount(dbr+1) = count(elev/=realMissing .and. debris_thick_dom==0._rkind)
    else 
      if (sum(ndebris)>0) validCount(dbr+1) = count(elev/=realMissing .and. debris_thick_dom >0._rkind)
    end if
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
    if (validCount(dbr+1) == 0) cycle ! no valid points
    if (sum(nclean)==0 .and. dbr==0) cycle ! no clean, skip 
    if (sum(ndebris)==0 .and. dbr==1) cycle ! no debris, skip 
    do iDOM = 1, nDOM
      if (dbr==0)then ! no debris, dbr = 0, index = 1
        if (elev(iDOM)/=realMissing .and. debris_thick_dom(iDOM)==0._rkind) then
          validElev(j,dbr+1) = elev(iDOM)
          validMassChange(j,dbr+1) = massChange(iDOM)
          j = j + 1
        end if
      else ! debris > 0, dbr = 1, index = 2
        if (elev(iDOM)/=realMissing .and. debris_thick_dom(iDOM) >0._rkind) then
          validElev(j,dbr+1) = elev(iDOM)
          validMassChange(j,dbr+1) = massChange(iDOM)
          if (debris_thick_dom(iDOM) > 0._rkind) then ! for now, just take mean of means since is usually depth independent
            iden_soil = iden_soil + iden_soil_mean(iDOM)/validCount(dbr+1)
            theta_sat = theta_sat + theta_sat_mean(iDOM)/validCount(dbr+1)
          end if
          j = j + 1
        end if
      end if
    end do
  
    ! Sort the valid elevations in descending order
    do i = 1, validCount(dbr+1)-1
      do j = 1, validCount(dbr+1)-i
        if (validElev(j,dbr+1) < validElev(j+1,dbr+1)) then
          call swap_real(validElev(j,dbr+1), validElev(j+1,dbr+1))
          call swap_real(validMassChange(j,dbr+1), validMassChange(j+1,dbr+1))
        end if
      end do
    end do

    ! Combine points with the same elevation into a single point that is the average of those points
    i = 1
    do while (i <= validCount(dbr+1))
      n = 1
      sum_mass = validMassChange(i,dbr+1)
      do j = i+1, validCount(dbr+1)
        if (validElev(i,dbr+1) == validElev(j,dbr+1)) then
          sum_mass = sum_mass + validMassChange(j,dbr+1)
          n = n + 1
        else
          exit
        end if
      end do
      if (n > 1) then
        validMassChange(i,dbr+1) = sum_mass / n
        ! Shift remaining elements to the left
        do k = i+1, validCount(dbr+1)-n
          validElev(k,dbr+1) = validElev(k+n,dbr+1)
          validMassChange(k,dbr+1) = validMassChange(k+n,dbr+1)
        end do
        validCount(dbr+1) = validCount(dbr+1) - n + 1
      end if
      i = i + 1
    end do
    
    ! Calculate piecewise linear regression for mass balance as a function of elevation 
    !   also find ELA elevation, assuming monotonically decreasing mass balance with elevation decrease
    !   NOTE: function should probably be capped above and below, but for now just extrapolate
    if (validCount(dbr+1) == 1) then
      slope(1,dbr+1) = 0._rkind
      intercept(1,dbr+1) = validMassChange(1,dbr+1)
      if (validMassChange(1,dbr+1) > 0.0) then
        ELA_elev(dbr+1) = -1.e6_rkind ! all domains are accumulation
      else
        ELA_elev(dbr+1) = 1.e6_rkind ! all domains are ablation
      end if
    else
      ind = 0
      do i = 1, validCount(dbr+1)-1
        slope(i+1,dbr+1)= (validMassChange(i+1,dbr+1) - validMassChange(i,dbr+1)) / (validElev(i+1,dbr+1) - validElev(i,dbr+1))
        intercept(i+1,dbr+1) = validMassChange(i,dbr+1) - slope(i+1,dbr+1) * validElev(i,dbr+1)
        if (validMassChange(i,dbr+1) >= 0._rkind) ind = i
      end do
      if (ind == validCount(dbr+1)-1) ind = validCount(dbr+1)+1 ! all domains accumulation, extrapolate below lowest point
      if (ind == 0) ind = 1 ! all domains ablation, extrapolate above highest point
      slope(1,dbr+1) = slope(2,dbr+1)
      intercept(1,dbr+1) = intercept(2,dbr+1)
      if (slope(1,dbr+1)<0._rkind) then ! don't propogate mass balance inversions
        slope(1,dbr+1) = 1.e-6_rkind
        intercept(1,dbr+1) = validMassChange(1,dbr+1)
      end if
      slope(validCount(dbr+1)+1,dbr+1) = slope(validCount(dbr+1),dbr+1)
      intercept(validCount(dbr+1)+1,dbr+1) = intercept(validCount(dbr+1),dbr+1)
      if (slope(validCount(dbr+1)+1,dbr+1)<0._rkind) then ! don't propogate mass balance inversions
        slope(validCount(dbr+1)+1,dbr+1) = 1.e-6_rkind
        intercept(validCount(dbr+1)+1,dbr+1) = validMassChange(validCount(dbr+1),dbr+1)
      end if
      ELA_elev(dbr+1) = -intercept(ind,dbr+1) / slope(ind,dbr+1)
    end if

  enddo ! end of loop for debris elevation relationships 
  ELA_use = ELA_elev(1) ! default to clean ablation ELA
  if (ELA_use<0._rkind) then ! no clean ablation
    ELA_use = ELA_elev(2)
  end if
  
  print*, "elev0", elev0
  print*, "area0", area0
  print*, "debris_thick_dom0", debris_thick_dom
  print*, "massChange", massChange/1000/t_total !(kg m-2 s-1 to m/s)
  print*, "valid cln",validCount(1), validMassChange(:,1),'xx',validElev(:,1)
  print*, "valid deb",validCount(2), validMassChange(:,2),'xx',validElev(:,2)
 
  ! Initialize new domain areas and elevations
  do i = 1,nDOM
    area(i) = 0._rkind
    ablFrac(i) = 0._rkind
    elev(i) = 0._rkind
    debris_thick_dom(i) = 0._rkind
  end do

  ! Allocate the mapping array from glacier id to index in gridInfo
  allocate(glacid_to_index(nGlacier))

  ! Populate the mapping array
  do iGlac = 1, nGlacier
    glacid_to_index(iGlac) = -1  ! Initialize with an invalid index
    do iGrid = 1, size(gridInfo(:)%grid_id)
      if (gridInfo(iGrid)%grid_id==glacInfo(iGlac)%glac_id) then
        glacid_to_index(iGlac) = iGrid
        exit
      endif
    end do
  end do

  ! run flow for each glacier
  do iGlac = 1, nGlacier
    iGrid = glacid_to_index(iGlac)

    ! set up glacier grid
    ny = gridInfo(iGrid)%ny
    nx = gridInfo(iGrid)%nx
    dx = gridInfo(iGrid)%dx
    dy = gridInfo(iGrid)%dy
    
    ! set height arrays and masks
    allocate(hgt(nx,ny), glacClnMask(nx,ny), glacHiMask(nx,ny), glacLoMask(nx,ny), glacDbrMask(nx,ny), glacAblMask(nx,ny))

    ! set up grid data
    allocate(surface(nx,ny), bed(nx,ny), cell2hru(nx,ny), glacierMask(nx,ny), debris(nx,ny), bed0(nx,ny))
    surface = gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny)
    bed = gridData%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)
    debris = gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny)
    cell2hru = int(gridData%grid(iGrid)%var(iLookGRID%cell2hru)%dat2(1:nx,1:ny))
    glacierMask = int(gridData%grid(iGrid)%var(iLookGRID%glacierMask)%dat2(1:nx,1:ny))

    ! make non-glacier bed have high elevation so glacier does not grow there (this is somewhat arbitrary and unstable, maybe should remove)
    !bed = merge(bed,bed+100._rkind,glacierMask==1)

    ! compute flow
    call run_model(t_total, debris, surface, bed, glacierMask, slope, intercept, validElev, validCount, &
                  maxCount, C_constant, dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA_use, nx, ny, dx, dy, volume)

    totVolume = totVolume+volume ! add volume to total volume, includes debris, not currently used
    
    ! Initialize variables
    hgt = surface - bed

    ! Calculate glacier accumulation and ablation areas
    glacAccArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface>= ELA_use))*dx*dy
    glacAblArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface<  ELA_use))*dx*dy
    ! debugging print
    print*, 'glacAblArea = ', glacAblArea(iGlac), ' glacAccArea = ', glacAccArea(iGlac)," ela = ",ELA_use

    ! Loop through HRUs and calculate domain areas and elevations for each HRU
    ! Order of domains will go HRU 1: clean1, clean2, debris; HRU 2: clean1, clean2, debris etc.
    ! It is possible one of the domains is missing if there is not two clean or debris cover 
    n = 0 ! initialize domain counter
    do iHRU = 1, nHRU
      ! start at last index + 1
      glacClnMask = merge(glacierMask, 0_i4b, cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris< min_thickness)
      glacDbrMask = merge(glacierMask, 0_i4b, cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris>=min_thickness)
      ! if no domain to go into, then give all area to the other domain
      if(nclean(iHRU)==0 .and. sum(glacClnMask)>0) glacDbrMask = glacDbrMask + glacClnMask ! give all area to debris domain (will reduce average domain debris thickness)
      if(ndebris(iHRU)==0 .and. sum(glacDbrMask)>0) glacClnMask = glacClnMask + glacDbrMask ! give all area to clean domain (will make debris thickness zero)

      if (nclean(iHRU)>0)then
        if (nclean(iHRU)>1) then 
          n = n+2 ! two clean domains, now n is last index of clean domains
          ! give higher cells to higher elevation domain
          hiInd = n
          loInd = n-1
          if (elev0(n-1)>elev0(n)) then
            hiInd = n-1
            loInd = n
          end if
          ! keep area ratio the same between the two clean domains
          areaAdd = area0(hiInd)/(area0(n-1)+area0(n)) * sum(glacClnMask) *dx*dy

          if (areaAdd > 0._rkind) then
            ! Flatten the surface array and sort it in descending order
            numCells = nx * ny
            allocate(sortedElev(numCells), sortedIndices(numCells))
            sortedElev = reshape(surface, [numCells])
            sortedIndices = [(i, i=1, numCells)]

            ! Sort the elevations in descending order
            do i = 1, numCells - 1
              do j = i + 1, numCells
                if (sortedElev(i) < sortedElev(j)) then
                  ! Swap elevations
                  call swap_real(sortedElev(i), sortedElev(j))
                  ! Swap indices
                call swap_integer(sortedIndices(i), sortedIndices(j))
                end if
              end do
            end do

            ! Add cells to the mask until the cumulative area matches areaAdd
            cumulativeArea = 0._rkind
            do i = 1, numCells
              ! Get the 2D indices from the 1D index
              j = mod(sortedIndices(i) - 1, nx) + 1_i4b
              k = (sortedIndices(i) - 1) / nx + 1_i4b
            
              ! Add the cell to the mask if in clean area
              if (glacClnMask(j, k) == 1_i4b) then
                glacHiMask(j, k) = 1_i4b
                cumulativeArea = cumulativeArea + dx * dy
              end if
              if (cumulativeArea >= areaAdd) exit
            end do
            deallocate(sortedElev, sortedIndices)
          end if

          ! Set the mask for the other domain
          glacLoMask = merge(glacClnMask, 0_i4b, glacHiMask == 0_i4b)

          ! Calculate areas, elevations, and debris thickness for 2 clean domains
          area(hiInd) = area(hiInd) + sum(glacHiMask)*dx*dy
          elev(hiInd) = elev(hiInd) + sum(surface * glacHiMask) *dx*dy
          glacAblMask = merge(glacHiMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use)
          ablFrac(hiInd) = ablFrac(hiInd)+ sum(glacAblMask) *dx*dy
          debris_thick_dom(hiInd) = 0._rkind

          area(loInd) = area(loInd) + sum(glacLoMask)*dx*dy
          elev(loInd) = elev(loInd) + sum(surface * glacLoMask) *dx*dy
          glacAblMask = merge(glacLoMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use)
          ablFrac(loInd) = ablFrac(loInd)+ sum(glacAblMask) *dx*dy
          debris_thick_dom(loInd) = 0._rkind
        else ! only one clean domain
          n = n+1 ! now n is last index of clean domains
          area(n) = area(n) + sum(glacClnMask)*dx*dy
          elev(n) = elev(n) + sum(surface * glacClnMask) *dx*dy
          debris_thick_dom(n) = 0._rkind
          glacAblMask = merge(glacClnMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use)
          ablFrac(n) = ablFrac(n)+ sum(glacAblMask) *dx*dy
        end if
      end if

      if (ndebris(iHRU)>0)then
        n = n+1 ! currently only one debris domain possible, now n is last index of debris domains
        area(n) = area(n) + sum(glacDbrMask)*dx*dy
        elev(n) = elev(n) + sum(surface * glacDbrMask) *dx*dy
        debris_thick_dom(n) = debris_thick_dom(n) + sum(debris * glacDbrMask) *dx*dy
        glacAblMask = merge(glacDbrMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use)
        ablFrac(n) = ablFrac(n)+ sum(glacAblMask) *dx*dy ! should be 1.0
      end if
    end do

    ! update gridData and deallocate
    gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) = surface
    gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny) = debris
    deallocate(hgt, surface, bed, cell2hru, glacierMask, debris, glacClnMask, glacHiMask, glacLoMask, glacDbrMask, glacAblMask)

  enddo ! end of glacier loop

  ! Set elevations to realMissing if no area in domain
  do iDOM = 1,nDOM
    if (elev(iDOM)==0._rkind) elev(iDOM)=realMissing
    if (area(iDOM)>0._rkind)then 
      elev(iDOM) = elev(iDOM) / area(iDOM)
      ablFrac(iDOM) = ablFrac(iDOM) / area(iDOM)
      debris_thick_dom(iDOM) = debris_thick_dom(iDOM) / area(iDOM)
    end if
  end do
  ! debugging print
  print*, "elev", elev 
  print*, "area", area
  print*, "debris_thick_dom", debris_thick_dom

  deallocate(glacid_to_index, validElev, validMassChange, slope, intercept)

end subroutine glacierFlowAreaChange


! ************************************************************************************************
! private subroutine run_model for setting up flow model and running for each glacier for time period
! ************************************************************************************************
  subroutine run_model(t_total, debris, S, B, glacierMask, slope, intercept, validElev, validCount, &
                      maxCount, C_constant, dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA, nx, ny, dx, dy, volume)
    ! Arguments
    real(rkind), intent(in) :: t_total, dx, dy
    real(rkind), intent(inout) :: debris(nx,ny), S(nx,ny), B(nx,ny)
    real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
    real(rkind), intent(in) ::  C_constant, dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA
    integer(i4b), intent(in) :: validCount(2), maxCount
    integer(i4b), intent(in) :: ny, nx, glacierMask(nx,ny)
    real(rkind), intent(out) :: volume

    ! Local variables
    real(rkind) :: dt, max_dt, min_dt, deltat, div_q(nx,ny), grad_debris(nx,ny), dt_cfl, meanS
    real(rkind) :: gamma, t, m_dot(nx,ny)
    integer(i4b),parameter :: n=3 ! Glen's flow law exponent
    real(rkind),parameter :: A=2.4e-24 ! Modern Glen parameter
    real(rkind),parameter :: cfl= 0.124 ! Courant-Friedrichs-Lewy condition
    integer(i4b) :: i, l(ny), lp(ny), lpp(ny), lm(ny), lmm(ny)
    integer(i4b) :: k(nx), kp(nx), kpp(nx), km(nx), kmm(nx), C_mask(nx,ny)

    ! y direction indices
    l = [(i, i=1, ny)]
    lp = [(i, i=2, ny), ny]
    lpp = [(i, i=3, ny), ny, ny]
    lm = [1, (i, i=1, ny-1)]
    lmm = [1, 1, (i, i=1, ny-2)]
    if (ny == 1) then
      lpp = [1]
      lmm = [1]
    end if
    
    ! x direction indices
    k = [(i, i=1, nx)]
    kp = [(i, i=2, nx), nx]
    kpp = [(i, i=3, nx), nx, nx]
    km = [1, (i, i=1, nx-1)]
    kmm = [1, 1, (i, i=1, nx-2)]
    if (nx == 1) then
      kpp = [1]
      kmm = [1]
    end if

    gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
    max_dt = 31._rkind * secprday! max timestep in seconds, a month
    min_dt = 0._rkind ! min timestep in seconds 
    t = 0._rkind
    !meanS = sum(merge(S - debris,0._rkind,glacierMask==1)) / count(glacierMask==1)
    !print*, "meanS", meanS, maxval(merge(S - debris,0._rkind,glacierMask==1)), minval(merge(S - debris,0._rkind,glacierMask==1))
    !meanS = sum(merge(B,0._rkind,glacierMask==1)) / count(glacierMask==1)
    !print*, "meanB", meanS, maxval(merge(B,0._rkind,glacierMask==1)), minval(merge(B,0._rkind,glacierMask==1))
    !meanS = sum(merge(S - debris-B,0._rkind,glacierMask==1)) / count(glacierMask==1)
    !print*, "meanH", meanS, maxval(merge(S - debris-B,0._rkind,glacierMask==1)), minval(merge(S - debris-B,0._rkind,glacierMask==1)), t

    do while (t < t_total)
      dt = t_total - t

      ! get mass balance rate over surface
      call get_m_dot(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny, m_dot)

      S = S - debris ! remove debris from glacier surface for flow calculation
      ! Select diffusion method and call step
      if (method == "MUSCL") then
        call diffusion_MUSCL(debris, S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, grad_debris, dt_cfl)
      else if (method == "upstream") then
        call diffusion_upstream(debris, S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, l, lp, lm, dx, dy, div_q, grad_debris, dt_cfl)
      end if
  
      ! update time step
      deltat = min(dt_cfl, dt)
      if (deltat > max_dt) deltat = max_dt
      if (deltat < min_dt) deltat = min_dt
      t = t + deltat
  
      ! Update S
      S = S + (m_dot + div_q) * deltat
      S = merge(S, B, S > B)

      meanS = sum(merge(S-B,0._rkind,glacierMask==1)) / count(glacierMask==1)
      print*, "meanH", meanS, maxval(merge(S-B,0._rkind,glacierMask==1)), minval(merge(S-B,0._rkind,glacierMask==1)), t

      ! Check that the glacier is in boundaries, fix small violations
      if (any((S - B) > 0._rkind .and. glacierMask==0)) then
        if (any((S - B) > 10.0 .and. glacierMask==0)) stop 'Glacier exceeds boundaries'
        S = merge(B, S, (S - B) > 0._rkind .and. glacierMask==0)
      end if
      ! Check that glacier surface is not infinite (unstable), bring down to mean glacier height
      ! This is a temporary fix, should be replaced with a more sophisticated method
      if (any((S - B) > 1.e6_rkind .and. glacierMask==1)) then
        meanS = sum(merge(S,0._rkind,glacierMask==1 .and. (S - B) < 1.e6_rkind)) / count(glacierMask==1 .and. (S - B)<=1.e6_rkind)
        S = merge(S,meanS,((S - B)<=1.e6_rkind))
      endif
      
      ! Update debris thickness if there is debris
      if (sum(debris)>0._rkind)then 
        !   NOTE: this uses the englacial debris advection transport model of Anderson and Anderson (2016)
        !         Englacial debris diffusion transport is ignored. The use of advection is valid for glacial
        !         debris transport because advection is defined as the transport of materials due to the bulk
        !         motion of a fluid, and glacier ice is a form of a viscoelastic fluid.
        !   NOTE2: debris thickness < min_thickness is considered as clean ice, so no debris advection
        call get_C_mask(S, B, debris, lat_moraine_wid, ELA, nx, ny, dx, dy, C_mask) ! Calculate near-surface concentration of debris
        debris = debris + ( C_constant*C_mask * m_dot/(1._rkind - theta_sat)/iden_soil + grad_debris ) * deltat
        debris = merge(debris, 0._rkind, S > B) ! remove debris from bedrock
        debris = merge(debris, 0._rkind, debris>=min_thickness) ! Check that debris thickness is not below minimum thickness
        debris = merge(debris, 0._rkind, debris<=dbr_crit) ! Check that debris thickness is not over critical value, effectively making terminal clean ice wedge 
      endif
      ! add debris back to glacier surface
      S = S + debris
    end do ! end of time loop

    ! Calculate volume of glacier (includes debris)
    volume = sum(merge(S-B,0._rkind,glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3

  end subroutine run_model


! ************************************************************************************************
! private subroutines for calculating a time step of glacier flow
! ************************************************************************************************
! ************************************************************************************************
! upwind diffusion scheme
! ************************************************************************************************
subroutine diffusion_upstream(debris, S, B, mask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, l, lp, lm, dx, dy, div_q, grad_debris, dt_cfl)
  implicit none
  real(rkind), intent(in) :: debris(nx,ny), S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, k(nx), kp(nx), km(nx), l(ny), lp(ny), lm(ny), n
  real(rkind), intent(out) :: div_q(nx,ny), grad_debris(nx,ny), dt_cfl
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny), Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkml(nx,ny), Hkl(nx,ny), Hklp(nx,ny), Hklm(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_dn(nx,ny)
  real(rkind) :: H_l_upstream_up(nx,ny), H_l_upstream_dn(nx,ny)
  real(rkind) :: f_l_plus(nx,ny), f_l_min(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_up(nx,ny), H_k_dn(nx,ny)
  real(rkind) :: H_k_upstream_up(nx,ny), H_k_upstream_dn(nx,ny)
  real(rkind) :: f_k_plus(nx,ny), f_k_min(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny), grad_k(nx,ny), grad_l(nx,ny)
  real(rkind) :: divisor
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography  
  
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
  H_l_up = H_index(Hklp, Hkl)
  H_l_dn = H_index(Hkl, Hklm)
  H_l_upstream_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Sklp(j,i) > Skl(j,i)) then
        H_l_upstream_up(j,i) = Hklp(j,i)
      else
        H_l_upstream_up(j,i) = Hkl(j,i)
      end if
    end do
  end do
  H_l_upstream_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) > Sklm(j,i)) then
        H_l_upstream_dn(j,i) = Hkl(j,i)
      else
        H_l_upstream_dn(j,i) = Hklm(j,i)
      end if
    end do
  end do

  ! calculate l flux
  f_l_plus = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_min  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up = gamma * H_l_up**(n+1_i4b) * H_l_upstream_up * f_l_plus
  D_l_dn = gamma * H_l_dn**(n+1_i4b) * H_l_upstream_dn * f_l_min
  ! Enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k upstream 
  H_k_up = H_index(Hkpl, Hkl)
  H_k_dn = H_index(Hkl, Hkml)
  H_k_upstream_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skpl(j,i) > Skl(j,i)) then
        H_k_upstream_up(j,i) = Hkpl(j,i)
      else
        H_k_upstream_up(j,i) = Hkl(j,i)
      end if
    end do
  end do
  H_k_upstream_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) > Skml(j,i)) then
        H_k_upstream_dn(j,i) = Hkl(j,i)
      else
        H_k_upstream_dn(j,i) = Hkml(j,i)
      end if
    end do
  end do

  ! calculate k flux
  f_k_plus = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_min  = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up = gamma * H_k_up**(n+1_i4b) * H_k_upstream_up * f_k_plus
  D_k_dn = gamma * H_k_dn**(n+1_i4b) * H_k_upstream_dn * f_k_min
  ! Enforce zero diffusion outside the mask
  D_k_up = merge(0._rkind, D_k_up, mask==0)
  D_k_dn = merge(0._rkind, D_k_dn, mask==0)

  ! calculate delta t and t
  divisor = max(maxval(abs(D_k_up)), maxval(abs(D_k_dn)), maxval(abs(D_l_up)), maxval(abs(D_l_dn)))

  if (divisor == 0._rkind) then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dx**2_i4b, dy**2_i4b) / divisor
  end if

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy)
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx)
  div_q = div_k + div_l

  ! calculate debris velocity term with a simple finite difference
  if (sum(debris)>0._rkind) then
    Hklp  = debris(k ,lp)*div_l(k ,lp)
    Hklm  = debris(k ,lm)*div_l(k ,lm)
    Hkl   = debris(k ,l )*div_l(k ,l )
    H_l_up = H_index(Hklp, Hkl)
    H_l_dn = H_index(Hkl, Hklm)
    grad_l = (H_l_up - H_l_dn) / dy

    Hkpl  = debris(kp,l )*div_k(kp,l )
    Hkml  = debris(km,l )*div_k(km,l )
    Hkl   = debris(k ,l )*div_k(k ,l )
    H_k_up   = H_index(Hkpl, Hkl)
    H_k_dn = H_index(Hkl, Hkml)
    grad_k = (H_k_up - H_k_dn) / dx

    grad_debris = grad_k + grad_l
  else
    grad_debris = 0._rkind
  end if

end subroutine diffusion_upstream


! ************************************************************************************************
! MUSCL scheme
! ************************************************************************************************
subroutine diffusion_MUSCL(debris, S, B, mask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, grad_debris, dt_cfl)
  implicit none
  real(rkind), intent(in) :: debris(nx,ny), S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, k(nx), kp(nx), km(nx), kpp(nx), kmm(nx), l(ny), lp(ny), lm(ny), lpp(ny), lmm(ny), n
  real(rkind), intent(out) :: div_q(nx,ny), grad_debris(nx,ny), dt_cfl
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny), Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkppl(nx,ny), Hkml(nx,ny), Hkmml(nx,ny), Hkl(nx,ny), Hklp(nx,ny), Hklpp(nx,ny), Hklm(nx,ny), Hklmm(nx,ny)
  real(rkind) :: H_l_min_up(nx,ny), H_l_plus_up(nx,ny)
  real(rkind) :: H_l_min_dn(nx,ny), H_l_plus_dn(nx,ny)
  real(rkind) :: f_l_plus(nx,ny), f_l_min(nx,ny)
  real(rkind) :: D_l_up_m(nx,ny), D_l_up_p(nx,ny), D_l_up_min(nx,ny), D_l_up_max(nx,ny)
  real(rkind) :: D_l_dn_m(nx,ny), D_l_dn_p(nx,ny), D_l_dn_min(nx,ny), D_l_dn_max(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_min_up(nx,ny), H_k_plus_up(nx,ny)
  real(rkind) :: H_k_min_dn(nx,ny), H_k_plus_dn(nx,ny)
  real(rkind) :: f_k_plus(nx,ny), f_k_min(nx,ny)
  real(rkind) :: D_k_up_m(nx,ny), D_k_up_p(nx,ny), D_k_up_min(nx,ny), D_k_up_max(nx,ny)
  real(rkind) :: D_k_dn_m(nx,ny), D_k_dn_p(nx,ny), D_k_dn_min(nx,ny), D_k_dn_max(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_dn(nx,ny)
  real(rkind) :: H_k_up(nx,ny), H_k_dn(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny), grad_k(nx,ny), grad_l(nx,ny)
  real(rkind) :: divisor
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography

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
  H_l_min_up  = H_min(Hklm, Hkl, Hklp)
  H_l_plus_up = H_plus(Hkl, Hklp, Hklpp)

  ! calculate l-1/2 index
  H_l_min_dn  = H_min(Hklmm, Hklm, Hkl)
  H_l_plus_dn = H_plus(Hklm, Hkl, Hklp)

  ! calculate l flux
  f_l_plus = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_min  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up_m = gamma * H_l_min_up**(n+2_i4b) * f_l_plus   ! equation 30 Jarosh 2013
  D_l_up_p = gamma * H_l_plus_up**(n+2_i4b) * f_l_plus  ! equation 30 Jarosh 2013
  D_l_up_min = min(D_l_up_m, D_l_up_p)                  ! equation 31 Jarosh 2013
  D_l_up_max = max(D_l_up_m, D_l_up_p)                  ! equation 32 Jarosh 2013
  !
  D_l_dn_m = gamma * H_l_min_dn**(n+2_i4b) * f_l_min  ! equation 30 Jarosh 2013
  D_l_dn_p = gamma * H_l_plus_dn**(n+2_i4b) * f_l_min ! equation 30 Jarosh 2013
  D_l_dn_min = min(D_l_dn_m, D_l_dn_p)                  ! equation 31 Jarosh 2013
  D_l_dn_max = max(D_l_dn_m, D_l_dn_p)                  ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_l_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Sklp(j,i) <= Skl(j,i) .and. H_l_min_up(j,i) <= H_l_plus_up(j,i)) then
        D_l_up(j,i) = D_l_up_min(j,i)
      else if (Sklp(j,i) <= Skl(j,i) .and. H_l_min_up(j,i) > H_l_plus_up(j,i)) then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if (Sklp(j,i) > Skl(j,i) .and. H_l_min_up(j,i) <= H_l_plus_up(j,i)) then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if (Sklp(j,i) > Skl(j,i) .and. H_l_min_up(j,i) > H_l_plus_up(j,i)) then
        D_l_up(j,i) = D_l_up_min(j,i)
      end if
    end do
  end do
  D_l_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) <= Sklm(j,i) .and. H_l_min_dn(j,i) <= H_l_plus_dn(j,i)) then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      else if (Skl(j,i) <= Sklm(j,i) .and. H_l_min_dn(j,i) > H_l_plus_dn(j,i)) then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if (Skl(j,i) > Sklm(j,i) .and. H_l_min_dn(j,i) <= H_l_plus_dn(j,i)) then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if (Skl(j,i) > Sklm(j,i) .and. H_l_min_dn(j,i) > H_l_plus_dn(j,i)) then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      end if
    end do
  end do
  ! enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k+1/2 index
  H_k_min_up  = H_min(Hkml, Hkl, Hkpl)
  H_k_plus_up = H_plus(Hkl, Hkpl, Hkppl)

  ! calculate k-1/2 index
  H_k_min_dn  = H_min(Hkmml, Hkml, Hkl)
  H_k_plus_dn = H_plus(Hkml, Hkl, Hkpl)

  ! calculate k flux
  f_k_plus = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_min  = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up_m = gamma * H_k_min_up**(n+2_i4b) * f_k_plus   ! equation 30 Jarosh 2013
  D_k_up_p = gamma * H_k_plus_up**(n+2_i4b) * f_k_plus  ! equation 30 Jarosh 2013
  D_k_up_min = min(D_k_up_m, D_k_up_p)                  ! equation 31 Jarosh 2013
  D_k_up_max = max(D_k_up_m, D_k_up_p)                  ! equation 32 Jarosh 2013
  !
  D_k_dn_m = gamma * H_k_min_dn**(n+2_i4b) * f_k_min  ! equation 30 Jarosh 2013
  D_k_dn_p = gamma * H_k_plus_dn**(n+2_i4b) * f_k_min ! equation 30 Jarosh 2013
  D_k_dn_min = min(D_k_dn_m, D_k_dn_p)                  ! equation 31 Jarosh 2013
  D_k_dn_max = max(D_k_dn_m, D_k_dn_p)                  ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_k_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skpl(j,i) <= Skl(j,i) .and. H_k_min_up(j,i) <= H_k_plus_up(j,i)) then
        D_k_up(j,i) = D_k_up_min(j,i)
      else if (Skpl(j,i) <= Skl(j,i) .and. H_k_min_up(j,i) > H_k_plus_up(j,i)) then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if (Skpl(j,i) > Skl(j,i) .and. H_k_min_up(j,i) <= H_k_plus_up(j,i)) then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if (Skpl(j,i) > Skl(j,i) .and. H_k_min_up(j,i) > H_k_plus_up(j,i)) then
        D_k_up(j,i) = D_k_up_min(j,i)
      end if
    end do
  end do
  D_k_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) <= Skml(j,i) .and. H_k_min_dn(j,i) <= H_k_plus_dn(j,i)) then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      else if (Skl(j,i) <= Skml(j,i) .and. H_k_min_dn(j,i) > H_k_plus_dn(j,i)) then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if (Skl(j,i) > Skml(j,i) .and. H_k_min_dn(j,i) <= H_k_plus_dn(j,i)) then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if (Skl(j,i) > Skml(j,i) .and. H_k_min_dn(j,i) > H_k_plus_dn(j,i)) then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      end if
    end do
  end do
  ! enforce zero diffusion outside the mask
  D_k_up = merge(0._rkind, D_k_up, mask==0)
  D_k_dn = merge(0._rkind, D_k_dn, mask==0)

  ! calculate delta t and t
  divisor = max(maxval(abs(D_k_up)), maxval(abs(D_k_dn)), maxval(abs(D_l_up)), maxval(abs(D_l_dn)))
  if (divisor == 0._rkind) then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dy**2, dx**2) / divisor
  end if

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy) ! equation 36 Jarosh 2013
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx) ! equation 36 Jarosh 2013
  div_q = div_k + div_l

  ! calculate debris velocity term with a simple finite difference 
  ! NOTE: maybe this should also use a MUSCL scheme
  if (sum(debris)>0._rkind) then
    Hklp  = debris(k ,lp)*div_l(k ,lp)
    Hklm  = debris(k ,lm)*div_l(k ,lm)
    Hkl   = debris(k ,l )*div_l(k ,l )
    H_l_up = H_index(Hklp, Hkl)
    H_l_dn = H_index(Hkl, Hklm)
    grad_l = (H_l_up - H_l_dn) / dy

    Hkpl  = debris(kp,l )*div_k(kp,l )
    Hkml  = debris(km,l )*div_k(km,l )
    Hkl   = debris(k ,l )*div_k(k ,l )
    H_k_up   = H_index(Hkpl, Hkl)
    H_k_dn = H_index(Hkl, Hkml)
    grad_k = (H_k_up - H_k_dn) / dx

    grad_debris = grad_k + grad_l
  else
    grad_debris = 0._rkind
  end if

end subroutine diffusion_MUSCL

! ************************************************************************************************
! private function for mass balance rate distribution after suface height changes
! ************************************************************************************************
subroutine get_m_dot(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny, m_dot)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), debris(nx,ny), t_total
  integer(i4b), intent(in) :: glacierMask(nx,ny), maxCount, nx, ny, validCount(2)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
  real(rkind), intent(out) :: m_dot(nx,ny)
  integer(i4b) :: ind, dbr, i, j, k
  
  ! distribute mass balance over surface, using all points in GRU
  do k = 1, nx
    do j = 1, ny
      ! Set mass balance to zero if not on glacier
      if (glacierMask(k,j)==0) then
        m_dot(k,j) = 0._rkind
      else
        ! set debris
        dbr = 0
        if (debris(k,j) > 0) dbr = 1
        ! find the index of the elevation in the validElev array
        ind = 0
        do i = 1, validCount(dbr+1)
          if (S(k,j) <= validElev(i,dbr+1)) ind = i
        end do
        if (ind == 0) ind = validCount(dbr+1)+1 ! elevation is below the lowest valid elevation, extrapolate down
        m_dot(k,j) = (slope(ind,dbr+1) * S(k,j) + intercept(ind,dbr+1))/ t_total
      end if
    end do
  end do

  ! convert mass balance to m s-1 from kg m-2 s-1
  m_dot = m_dot / iden_water

end subroutine get_m_dot

! ************************************************************************************************
! private functions for calculate of spatial mask of near-surface concentration of debris
!   Currently adding debris concentration along the sides of the glacier below the ELA, 
!   and where debris currently is, somewhat of a stub
! ************************************************************************************************
subroutine get_C_mask(S, B, debris, lat_moraine_wid0, ELA, nx, ny, dx, dy, C_mask)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), debris(nx,ny), lat_moraine_wid0, ELA, dx, dy
  integer(i4b), intent(in) :: nx, ny
  integer(i4b), intent(out) :: C_mask(nx,ny)
  integer(i4b) :: i, j, k, l, zeroMask(nx, ny)
  real(rkind) :: distance(nx, ny), dist, lat_moraine_wid
  
  ! Compute the Euclidean distance to the nearest zero cell
  zeroMask = merge(0, 1, S - B == 0._rkind) ! mask for cells where S = B
  distance = 1.e6_rkind  ! Initialize with a large value
  do i = 1, nx
    do j = 1, ny
        if (zeroMask(i, j) == 1) then
            distance(i, j) = 0._rkind
        else
            do k = 1, nx
                do l = 1, ny
                    if (zeroMask(k, l) == 1) then
                        dist = sqrt(((i - k) * dx)**2._i4b + ((j - l) * dy)**2._i4b)
                        if (dist < distance(i, j)) then
                            distance(i, j) = dist
                        end if
                    end if
                end do
            end do
        end if
    end do
  end do

  ! calulate the maximum distance to consider for debris concentration, or the width of lateral moraines
  !  based on width of debris concentration near sides of glacier at high elevations
  lat_moraine_wid = lat_moraine_wid0
  if (lat_moraine_wid <= 0._rkind)then
    lat_moraine_wid = 0._rkind
    do i = 1, nx
      do j = 1, ny
        if (S(i,j) > ELA .and. S(i,j) > B(i,j) .and. S(i,j) > 0._rkind) then
          if (distance(i,j) > lat_moraine_wid) then
            lat_moraine_wid = distance(i,j)
          end if
        end if
      end do
    end do
  end if

  ! Create the side mask based on the maximum distance or currently having debris 
  !   and where S>B (has glacier) and S<ELA (in ablation zone)
  C_mask = merge(1, 0, (distance<=lat_moraine_wid .or. debris>0._rkind) .and. S>B .and. S<ELA)

end subroutine get_C_mask

! ************************************************************************************************
! private functions to assist in the glacier flow calculations
! ************************************************************************************************
function minmod(a, b) result(minmod_result)
  implicit none
  real(rkind), intent(in) :: a(:,:), b(:,:)
  real(rkind), allocatable :: minmod_result(:,:)
  real(rkind), allocatable :: sign_ab(:,:)

  allocate(minmod_result(size(a,1), size(a,2)))
  allocate(sign_ab(size(a,1), size(a,2)))
  sign_ab = sign(1.0_rkind, a) + sign(1.0_rkind, b)
  minmod_result = sign_ab / 2._rkind * min(abs(a), abs(b))
end function minmod

function superbee(r) result(superbee_result)
  implicit none
  real(rkind), intent(in) :: r(:,:)
  real(rkind), allocatable :: superbee_result(:,:)

  allocate(superbee_result(size(r,1), size(r,2)))
  superbee_result = max(0._rkind, min(2._rkind * r, 1._rkind), min(r, 2._rkind))
end function superbee

function flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dx, dy, n) result(flux_result)
  implicit none
  real(rkind), intent(in) :: Skpl(:,:), Skml(:,:), Skplp(:,:), Skmlp(:,:), Sklp(:,:), Skl(:,:), dx, dy
  integer(i4b), intent(in) :: n
  real(rkind), allocatable :: flux_result(:,:)

  allocate(flux_result(size(Skpl,1), size(Skpl,2)))
  flux_result = ( (Skpl - Skml + Skplp - Skmlp)**2_i4b / (4._rkind * dx)**2_i4b & 
                 + (Sklp - Skl)**2_i4b / dy**2_i4b )**((n - 1.0) / 2._rkind)
end function flux

function SIA(Dup, Sup, S, Ddn, Sdn, d) result(SIA_result)
  implicit none
  real(rkind), intent(in) :: Dup(:,:), Sup(:,:), S(:,:), Ddn(:,:), Sdn(:,:), d
  real(rkind), allocatable :: SIA_result(:,:)

  allocate(SIA_result(size(Dup,1), size(Dup,2)))
  SIA_result = (Dup * (Sup - S) / d - Ddn * (S - Sdn) / d) / d
end function SIA

function H_index(h1, h2) result(H_index_result)
  implicit none
  real(rkind), intent(in) :: h1(:,:), h2(:,:)
  real(rkind), allocatable :: H_index_result(:,:)

  allocate(H_index_result(size(h1,1), size(h1,2)))
  H_index_result = 0.5_rkind * (h1 + h2)
end function H_index

function H_plus(Hm, H, Hp) result(H_plus_result)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: H_plus_result(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(H_plus_result(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if (limiter == "minmod") then
    H_plus_result = H - 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if (limiter == "superbee") then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(Hp - H, ones, mask)
    H_plus_result = merge(H - 0.5_rkind * superbee(abs((H - Hm) / divisor)) * (Hp - H), H, mask)
  end if
  deallocate(mask,divisor,ones)
end function H_plus

function H_min(Hm, H, Hp) result(H_min_result)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: H_min_result(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(H_min_result(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if (limiter == "minmod") then
    H_min_result = H + 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if (limiter == "superbee") then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(H - Hm, ones, mask)
    H_min_result = merge(H + 0.5_rkind * superbee(abs((Hp - H) / divisor)) * (H - Hm), H, mask)
  end if
  deallocate(mask,divisor,ones)
end function H_min

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

! ************************************************************************************************
! public function to update glacier domain area, elevation, and layering 
! ************************************************************************************************
subroutine updateGlacierDomain(&
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
  real(rkind), intent(inout)      :: iLayerHeight(:)          ! layer interface height (m)
  real(rkind), intent(inout)      :: DOMarea                  ! area of each domain (m2)
  real(rkind), intent(inout)      :: ablFrac                  ! fraction of glacier domain area that is ablating
  real(rkind), intent(inout)      :: DOMelev                  ! elevation of each glacier domain (m) per HRU
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message 
   ! ----- define local variables ------------------------------------------------------------------------------------------
  integer(i4b)                    :: i                        ! loop index
  real(rkind)                     :: soil_thick               ! depth of soil== debris in debris domain of glacier HRU
  real(rkind)                     :: thick_ratio              ! ratio of new debris thickness to previous debris thickness
  real(rkind)                     :: delVol                   ! volume change of glacieret domain
  character(len=256)              :: cmessage                 ! error message
  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='updateGlacierDomain/'

  ! clean and debris ablation domains use area change from glacier flow subroutine
  if (dom_type==glacCln1 .or. dom_type==glacCln2 .or. dom_type==glacDbr)then
    iglac = iglac + 1
    DOMelev = glac_elev(iglac) ! realMissing if no area
    DOMarea = glac_area(iglac) ! may be 0
    ablFrac = glac_ablFrac(iglac) ! may be 0
    glacMass4AreaChange = 0._rkind ! reset
    if (dom_type==glacDbr)then
      ! thickness of average debris cover in HRU changes with debris advection
      ! scale soil layer thickness with debris thickness change, keep the same number of layers 
      soil_thick = sum(mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil))
      thick_ratio = glac_debris_thick(iglac)/soil_thick
      do i = nSnow+nLake+1,nSnow+nLake+nSoil
        mLayerDepth(i)  = mLayerDepth(i)*thick_ratio
        mLayerHeight(i) = mLayerHeight(i)*thick_ratio
        iLayerHeight(i) = iLayerHeight(i)*thick_ratio            
      end do
      do i = nSnow+nLake+nSoil+1,nSnow+nLake+nSoil+nGlce
        mLayerHeight(i) = mLayerHeight(i) + glac_debris_thick(iglac) - soil_thick
        iLayerHeight(i) = iLayerHeight(i) + glac_debris_thick(iglac) - soil_thick
      end do
    end if ! (if debris domain)

  ! glacieret domain uses volume area scaling following Macheret et al., 1988; Chen and Ohmura, 1990; Bahr, 1997
  elseif (dom_type==glacieret)then ! if glacieret grows too big, it will be fixed outside of this module in run_oneGRU
    delVol = glacMass4AreaChange * DOMarea/iden_ice/1.e9_rkind ! km3
    DOMarea = ( delVol/0.033_rkind + (DOMarea/1.e6_rkind)**1.36_rkind )**(1._rkind/1.36_rkind)
    ablFrac = 0._rkind ! glacieret does not keep track of ablating area
    DOMelev = DOMelev + glacMass4AreaChange
    glacMass4AreaChange = 0._rkind ! reset
  endif

  end subroutine updateGlacierDomain


end module glacierFlowAreaChange_module