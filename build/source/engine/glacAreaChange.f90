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
character(len=32),parameter :: bed_shape='parabolic' !'trapezoid'
character(len=32),parameter :: method='MUSCL' !'upstream'

! privacy
private::run_flowModel,volAreaScaling,diffusion_MUSCL,diffusion_upstream,minmod,superbee,flux,SIA,H_index,H_plus,H_min
public ::glacAreaChange
public ::updateGlacDomain
contains
! ************************************************************************************************
! public subroutine glacAreaChange: get new glacier area and elevation based on ice flow model and glacieret area based on vol area change
! NOTE: This will eventually run in parallel as a program, but for now it is serial
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
                    C_constant,              & ! intent(in):    non-spatial concentration for debris advection (kg m-3)
                    dbr_crit,                & ! intent(in):    critical debris thickness start debris-free terminal wedge (m)
                    lat_moraine_wid,         & ! intent(inout): lateral moraine width (m)
                    ! area
                    glacierAreaThresh,       & ! intent(in):    minimum glacier area to be considered a glacier (m2)
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
  real(rkind), intent(in)            :: C_constant                      ! non-spatial concentration for debris advection (kg m-3)
  real(rkind), intent(in)            :: dbr_crit                        ! critical debris thickness start debris-free terminal wedge (m)
  real(rkind), intent(in)            :: lat_moraine_wid                 ! lateral moraine width (m) 
  ! area 
  real(rkind), intent(in)            :: glacierAreaThresh               ! minimum glacier area to be considered a glacier (m2)
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
  real(rkind),parameter              :: min_thickness=0.01_rkind        ! minimum thickness of debris cover to be considered as debris cover (m)
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

  ! Initialize new domain areas and elevations
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
    ! debugging print
    if(printFlag)then 
      volume = sum(merge(glacierMask*(surface - bed), 0._rkind, glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3
      write(*,'(a,i2,a,4(1x,f8.3))') "<GLACIER ",iGlac, " START: ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume
    endif

    if(glacierAccArea(iGlac)+glacierAblArea(iGlac) < glacierAreaThresh)then ! glacier has shrunk below threshold, area may be zero
     ! a glacieret changes area with volume area scaling
     call volAreaScaling(t_total, debris, surface, bed, glacierMask, slope, intercept, validElev, validCount, maxCount, C_constant, &
                         dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA_use, nx, ny, dx, dy, volume, printFlag)
    else ! glacier, run flow model
      ! make non-glacier bed have high elevation so glacier does not grow there (this is somewhat arbitrary and unstable, maybe should remove)
      !bed = merge(bed,bed+100._rkind,glacierMask==1)
      ! compute flow
      call run_flowModel(t_total, debris, surface, bed, glacierMask, slope, intercept, validElev, validCount, maxCount, C_constant, &
                         dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA_use, nx, ny, dx, dy, volume, printFlag)
    endif

    totVolume = totVolume+volume ! add volume to total volume, includes debris, not currently used
    
    ! Initialize variables
    hgt = surface - bed

    ! Calculate glacier accumulation and ablation areas
    glacierAccArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface>= ELA_use))*dx*dy
    glacierAblArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface<  ELA_use))*dx*dy
    ! debugging print
    if(printFlag) write(*,'(a,i2,a,4(1x,f8.3))') ">GLACIER ",iGlac, " END  : ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume

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
        endif
      endif

      if(ndebris(iHRU)>0)then
        n = n+1 ! currently only one debris domain possible, now n is last index of debris domains
        area(n) = area(n) + sum(glacDbrMask)*dx*dy
        elev(n) = elev(n) + sum(surface * glacDbrMask) *dx*dy
        debris_thick_dom(n) = debris_thick_dom(n) + sum(debris * glacDbrMask) *dx*dy
        glacAblMask = merge(glacDbrMask, 0_i4b, cell2hru==hruInd(n) .and. hgt>thick4area .and. surface< ELA_use)
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
! private subroutine volAreaScaling: glacieret area change and debris transport with VA scaling
!   VA scaling is following Macheret et al., 1988; Chen and Ohmura, 1990; Bahr, 1997.
!   After VA scaling for new area, average height of adjoining lowest cells is added to new area increase, or to removed to the
!   bedrock on the lowest cells for new area decrease. Remaining volume change is added/removed as a new height to the entire new area.
!   Debris emerges below the new ELA as in the flow model. This height change is a kluge, but is needed to keep track of volume and
!   for the unlikely case of a glacieret growing into a glacier.
! ************************************************************************************************
  subroutine volAreaScaling(t_total, debris, S, B, glacierMask, slope, intercept, validElev, validCount, maxCount, C_constant, &
                            dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA, nx, ny, dx, dy, volume, printFlag)
    implicit none
    ! Arguments
    real(rkind), intent(in) :: t_total, dx, dy
    real(rkind), intent(inout) :: debris(nx,ny), S(nx,ny), B(nx,ny)
    real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
    real(rkind), intent(in) :: C_constant, dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA
    integer(i4b), intent(in) :: validCount(2), maxCount
    integer(i4b), intent(in) :: ny, nx, glacierMask(nx,ny)
    logical(lgt), intent(in) :: printFlag
    real(rkind), intent(out) :: volume

    ! Local variables
    real(rkind), parameter :: va_exp=1.36_rkind
    real(rkind) :: m_dot(nx,ny), delVol, delArea, area, va_constant, excessVol, H(nx,ny)
    real(rkind) :: around(8), grad_debris(nx,ny), S0(nx,ny)
    integer(i4b) :: C_mask(nx,ny) 

    S = S - debris ! remove debris from glacier surface before volume-area scaling
    S0 = S ! save initial surface for debris transport calculation

    ! initial S and B, calculate coefficients
    volume = sum(merge(S-B,0._rkind,glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3
    area = sum(merge(glacierMask,0_i4b,(S-B)>thick4area))* dx * dy * 1.e-6_rkind ! km2
    if(area <= 0._rkind)then ! no glacier detected coming in
      va_constant = 0.033_rkind ! use a typical value
    else
      va_constant = volume/area**va_exp ! might want to cap this
    endif
    ! debugging print
    if(printFlag)then
      H = merge(S+debris-B,0._rkind,glacierMask==1)
      write(*,'(a,f5.3)') " use VA scaling, va_constant = ", va_constant
      write(*,'(a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = 0.00 mean max H (m) =", sum(H)/count(glacierMask==1), maxval(H),&
            ", mean max debris (m) =", sum(debris)/count(glacierMask==1), maxval(debris)
    endif

    ! get mass balance rate over surface
    m_dot = massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
    delVol = sum(m_dot)*dx*dy*t_total*iden_ice/iden_water * 1.e-9_rkind ! convert from m s-1 to km3 in volume change

    if(volume + delVol <= 0._rkind)then
      delVol = -volume ! glacieret gone
      volume = 0._rkind
      S = B ! glacier surface is bedrock
      debris = 0._rkind ! no debris
      return
    elseif(delVol==0._rkind)then
      return ! no change in volume, so no change surface or debris
    endif
    delArea = (((volume + delVol)/va_constant)**(1._rkind/va_exp) - area)*1.e6_rkind ! change in area from volume-area scaling, m2
    delVol = delVol*1.e9_rkind ! back to m3

    ! Update S and debris velocity, estimates that are a bit klugey
    call updateS_volAreaScaling(S, S0, B, debris, glacierMask, delVol, delArea, thick4area, t_total, m_dot, nx, ny, dx, dy, excessVol, grad_debris)
    
    ! Update debris thickness if there is debris, using englacial debris advection transport model
    if(sum(debris)>0._rkind)then 
      C_mask = emergingDebris(S, B, debris, lat_moraine_wid, ELA, nx, ny, dx, dy) ! Calculate near-surface concentration of debris
      debris = debris + ( C_constant*C_mask * m_dot/(1._rkind - theta_sat)/iden_soil + grad_debris ) * t_total
      debris = merge(debris, 0._rkind, S > B) ! remove debris from bedrock
      debris = merge(debris, 0._rkind, debris>=min_thickness) ! Check that debris thickness is not below minimum thickness
      debris = merge(debris, 0._rkind, debris<=dbr_crit) ! Check that debris thickness is not over critical value, effectively making terminal clean ice wedge 
    endif
    ! add debris back to glacier surface
    S = S + debris    

    ! Recalculate volume of glacier (includes debris)
    volume = sum(merge(S-B,0._rkind,glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3
    ! debugging print
    if(printFlag)then
      H = merge(S-B,0._rkind,glacierMask==1)
      write(*,'(a,f4.2,a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = ", t_total/secprday/365.25_rkind, " mean max H (m) =", &
            sum(H)/count(glacierMask==1), maxval(H), ", mean max debris (m) =", sum(debris)/count(glacierMask==1), maxval(debris)
      write(*,'(a,2(1x,f8.4),a,f6.2)') "  change in area (km2), change in volume (km3) =", &
            delArea*1.e-6_rkind, delVol*1.e-9_rkind, ", excess height (m) added to make volume = ", excessVol/(sum(H) * dx * dy)
    endif

  end subroutine volAreaScaling


! ************************************************************************************************
! private subroutine run_flowModel: set up flow model and debris advection and run for time period
! ************************************************************************************************
  subroutine run_flowModel(t_total, debris, S, B, glacierMask, slope, intercept, validElev, validCount, maxCount, C_constant, &
                           dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA, nx, ny, dx, dy, volume, printFlag)
   implicit none
  ! Arguments
    real(rkind), intent(in) :: t_total, dx, dy
    real(rkind), intent(inout) :: debris(nx,ny), S(nx,ny), B(nx,ny)
    real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
    real(rkind), intent(in) :: C_constant, dbr_crit, min_thickness, lat_moraine_wid, iden_soil, theta_sat, ELA
    integer(i4b), intent(in) :: validCount(2), maxCount
    integer(i4b), intent(in) :: ny, nx, glacierMask(nx,ny)
    logical(lgt), intent(in) :: printFlag
    real(rkind), intent(out) :: volume

    ! Local variables
    real(rkind) :: dt, max_dt, min_dt, deltat, div_q(nx,ny), grad_debris(nx,ny), dt_cfl, meanS
    real(rkind) :: gamma, t, m_dot(nx,ny), H(nx,ny)
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

    gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
    max_dt = 31._rkind * secprday! max timestep in seconds, a month
    min_dt = 0._rkind ! min timestep in seconds 
    t = 0._rkind
    ! debugging print
    if(printFlag)then 
      H = merge(S-B,0._rkind,glacierMask==1)
      write(*,'(a)') " use flow model"
      write(*,'(a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = 0.00 mean max H (m) =", sum(H)/count(glacierMask==1), maxval(H),&
            ", mean max debris (m) =", sum(debris)/count(glacierMask==1), maxval(debris)
      i = 1 ! counter for debugging print, only print once a max_dt
    endif

    do while (t < t_total)
      dt = t_total - t

      ! get mass balance rate over surface
      m_dot = massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)

      S = S - debris ! remove debris from glacier surface for flow calculation
      ! Select diffusion method and call step
      if(method == "MUSCL")then
        call diffusion_MUSCL(debris, S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, grad_debris, dt_cfl)
      else if(method == "upstream")then
        call diffusion_upstream(debris, S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, l, lp, lm, dx, dy, div_q, grad_debris, dt_cfl)
      endif
  
      ! update time step
      deltat = min(dt_cfl, dt)
      if(deltat > max_dt) deltat = max_dt
      if(deltat < min_dt) deltat = min_dt
      t = t + deltat
  
      ! Update S
      S = S + (m_dot + div_q) * deltat
      S = merge(S, B, S > B)

      ! Check that the glacier is in boundaries, fix small violations
      if(any((S - B) > 0._rkind .and. glacierMask==0))then
        if(any((S - B) > 10.0 .and. glacierMask==0)) stop 'Glacier exceeds boundaries'
        S = merge(B, S, (S - B) > 0._rkind .and. glacierMask==0)
      endif
      ! Check that glacier surface is not infinite (unstable), bring down to mean glacier height
      ! This is a temporary fix, should be replaced with a more sophisticated method
      if(any((S - B) > 1.e6_rkind .and. glacierMask==1))then
        meanS = sum(merge(S,0._rkind,glacierMask==1 .and. (S - B) < 1.e6_rkind)) / count(glacierMask==1 .and. (S - B)<=1.e6_rkind)
        S = merge(S,meanS,((S - B)<=1.e6_rkind))
      endif
      
      ! Update debris thickness if there is debris, using englacial debris advection transport model
      if(sum(debris)>0._rkind)then
        C_mask = emergingDebris(S, B, debris, lat_moraine_wid, ELA, nx, ny, dx, dy) ! Calculate near-surface concentration of debris
        debris = debris + ( C_constant*C_mask * m_dot/(1._rkind - theta_sat)/iden_soil + grad_debris ) * deltat
        debris = merge(debris, 0._rkind, S > B) ! remove debris from bedrock
        debris = merge(debris, 0._rkind, debris>=min_thickness) ! Check that debris thickness is not below minimum thickness
        debris = merge(debris, 0._rkind, debris<=dbr_crit) ! Check that debris thickness is not over critical value, effectively making terminal clean ice wedge 
      endif
      ! add debris back to glacier surface
      S = S + debris
      ! debugging print
      if(printFlag)then 
        ! only print once a max_dt
        if(t > i*max_dt)then
          i = i + 1
          H = merge(S-B,0._rkind,glacierMask==1)
          write(*,'(a,f4.2,a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = ", t/secprday/365.25_rkind, " mean max H (m) =", &
                sum(H)/count(glacierMask==1), maxval(H), ", mean max debris (m) =", sum(debris)/count(glacierMask==1), maxval(debris)
        endif
      endif
    enddo ! end of time loop

    ! Calculate volume of glacier (includes debris)
    volume = sum(merge(S-B,0._rkind,glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3

  end subroutine run_flowModel


! ************************************************************************************************
! private subroutines for calculating a time step of glacier flow
! ************************************************************************************************
! ************************************************************************************************
! private subroutine diffusion_upstream: upwind diffusion scheme, not mass conserving, but stable and quicker
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
  real(rkind) :: div_k(nx,ny), div_l(nx,ny)
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
      if(Sklp(j,i) > Skl(j,i))then
        H_l_upstream_up(j,i) = Hklp(j,i)
      else
        H_l_upstream_up(j,i) = Hkl(j,i)
      endif
    enddo
  enddo
  H_l_upstream_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) > Sklm(j,i))then
        H_l_upstream_dn(j,i) = Hkl(j,i)
      else
        H_l_upstream_dn(j,i) = Hklm(j,i)
      endif
    enddo
  enddo

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
      if(Skpl(j,i) > Skl(j,i))then
        H_k_upstream_up(j,i) = Hkpl(j,i)
      else
        H_k_upstream_up(j,i) = Hkl(j,i)
      endif
    enddo
  enddo
  H_k_upstream_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) > Skml(j,i))then
        H_k_upstream_dn(j,i) = Hkl(j,i)
      else
        H_k_upstream_dn(j,i) = Hkml(j,i)
      endif
    enddo
  enddo

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

  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dx**2_i4b, dy**2_i4b) / divisor
  endif

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy)
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx)
  div_q = div_k + div_l

  ! calculate debris velocity term with a simple finite difference
  grad_debris = debris_velocity(debris, div_k, div_l, k, kp, km, l, lp, lm, dx, dy, nx, ny)

end subroutine diffusion_upstream


! ************************************************************************************************
! subroutine diffusion_MUSCL: MUSCL scheme, mass conserving, more accurate, but less stable and slower
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
  real(rkind) :: H_k_min_up(nx,ny), H_k_plus_up(nx,ny), H_k_min_dn(nx,ny), H_k_plus_dn(nx,ny)
  real(rkind) :: f_k_plus(nx,ny), f_k_min(nx,ny)
  real(rkind) :: D_k_up_m(nx,ny), D_k_up_p(nx,ny), D_k_up_min(nx,ny), D_k_up_max(nx,ny)
  real(rkind) :: D_k_dn_m(nx,ny), D_k_dn_p(nx,ny), D_k_dn_min(nx,ny), D_k_dn_max(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_dn(nx,ny), H_k_up(nx,ny), H_k_dn(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny)
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
      if(Sklp(j,i) <= Skl(j,i) .and. H_l_min_up(j,i) <= H_l_plus_up(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      else if(Sklp(j,i) <= Skl(j,i) .and. H_l_min_up(j,i) > H_l_plus_up(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_min_up(j,i) <= H_l_plus_up(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_min_up(j,i) > H_l_plus_up(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      endif
    enddo
  enddo
  D_l_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Sklm(j,i) .and. H_l_min_dn(j,i) <= H_l_plus_dn(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      else if(Skl(j,i) <= Sklm(j,i) .and. H_l_min_dn(j,i) > H_l_plus_dn(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_min_dn(j,i) <= H_l_plus_dn(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_min_dn(j,i) > H_l_plus_dn(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      endif
    enddo
  enddo
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
      if(Skpl(j,i) <= Skl(j,i) .and. H_k_min_up(j,i) <= H_k_plus_up(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      else if(Skpl(j,i) <= Skl(j,i) .and. H_k_min_up(j,i) > H_k_plus_up(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_min_up(j,i) <= H_k_plus_up(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_min_up(j,i) > H_k_plus_up(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      endif
    enddo
  enddo
  D_k_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Skml(j,i) .and. H_k_min_dn(j,i) <= H_k_plus_dn(j,i))then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      else if(Skl(j,i) <= Skml(j,i) .and. H_k_min_dn(j,i) > H_k_plus_dn(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_min_dn(j,i) <= H_k_plus_dn(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_min_dn(j,i) > H_k_plus_dn(j,i))then
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

  ! calculate debris velocity term with a simple finite difference 
  ! NOTE: maybe this should also use a MUSCL scheme
  grad_debris = debris_velocity(debris, div_k, div_l, k, kp, km, l, lp, lm, dx, dy, nx, ny)

end subroutine diffusion_MUSCL

! ************************************************************************************************
! private function  mass balance: rate distribution after suface height changes
! ************************************************************************************************
function massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), debris(nx,ny), t_total
  integer(i4b), intent(in) :: glacierMask(nx,ny), maxCount, nx, ny, validCount(2)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
  real(rkind) :: massBalance(nx,ny)
  integer(i4b) :: ind, dbr, i, j, k
  
  ! distribute mass balance over surface, using all points in GRU
  do k = 1, nx
    do j = 1, ny
      ! Set mass balance to zero if not on glacier
      if(glacierMask(k,j)==0)then
        massBalance(k,j) = 0._rkind
      else
        ! set debris
        dbr = 0
        if(debris(k,j) > 0) dbr = 1
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
! private function debris_velocity: calculation of debris velocity term for debris advection transport
! ************************************************************************************************
function debris_velocity(debris, div_k, div_l, k, kp, km, l, lp, lm, dx, dy, nx, ny)
  implicit none
  real(rkind), intent(in) :: debris(nx,ny), div_k(nx,ny), div_l(nx,ny), dx, dy
  integer(i4b), intent(in) :: k(nx), kp(nx), km(nx), l(ny), lp(ny), lm(ny), nx, ny
  real(rkind) :: debris_velocity(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkml(nx,ny), Hkl(nx,ny)
  real(rkind) :: Hklp(nx,ny), Hklm(nx,ny)
  real(rkind) :: H_k_up(nx,ny), H_k_dn(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_dn(nx,ny)
  integer(i4b) :: i, j 
  real(rkind) :: grad_k(nx,ny), grad_l(nx,ny)

  ! calculate debris velocity term with a simple finite difference
  if(sum(debris)>0._rkind)then
    Hklp  = debris(k ,lp)*div_l(k ,lp)
    Hklm  = debris(k ,lm)*div_l(k ,lm)
    Hkl   = debris(k ,l )*div_l(k ,l )
    H_l_up = H_index(Hklp, Hkl)
    H_l_dn = H_index(Hkl, Hklm)
    grad_l = (H_l_up - H_l_dn) / dy

    Hkpl  = debris(kp,l )*div_k(kp,l )
    Hkml  = debris(km,l )*div_k(km,l )
    Hkl   = debris(k ,l )*div_k(k ,l )
    H_k_up = H_index(Hkpl, Hkl)
    H_k_dn = H_index(Hkl, Hkml)
    grad_k = (H_k_up - H_k_dn) / dx

    debris_velocity = grad_k + grad_l
  else
    debris_velocity = 0._rkind
  endif

end function debris_velocity

! ************************************************************************************************
! private function emergingDebris: calculation of spatial mask of near-surface concentration of debris
!   Currently adding debris concentration along the sides of the glacier below the ELA, 
!   and where debris currently is, somewhat of a stub
!   NOTE: this uses the englacial debris advection transport model of Anderson and Anderson (2016)
!         Englacial debris diffusion transport is ignored. The use of advection is valid for glacial
!         debris transport because advection is defined as the transport of materials due to the bulk
!         motion of a fluid, and glacier ice is a form of a viscoelastic fluid.
!   NOTE2: debris thickness < min_thickness is considered as clean ice, so no debris advection
! ************************************************************************************************
function emergingDebris(S, B, debris, lat_moraine_wid0, ELA, nx, ny, dx, dy)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), debris(nx,ny), lat_moraine_wid0, ELA, dx, dy
  integer(i4b), intent(in) :: nx, ny
  integer(i4b) :: emergingDebris(nx,ny)
  integer(i4b) :: i, j, k, l, zeroMask(nx, ny)
  real(rkind) :: distance(nx, ny), dist, lat_moraine_wid
  
  ! Compute the Euclidean distance to the nearest zero cell
  zeroMask = merge(0, 1, S - B == 0._rkind) ! mask for cells where S = B
  distance = 1.e6_rkind  ! Initialize with a large value
  do i = 1, nx
    do j = 1, ny
        if(zeroMask(i, j) == 1)then
            distance(i, j) = 0._rkind
        else
            do k = 1, nx
                do l = 1, ny
                    if(zeroMask(k, l) == 1)then
                        dist = sqrt(((i - k) * dx)**2._i4b + ((j - l) * dy)**2._i4b)
                        if(dist < distance(i, j))then
                            distance(i, j) = dist
                        endif
                    endif
                enddo
            enddo
        endif
    enddo
  enddo

  ! calulate the maximum distance to consider for debris concentration, or the width of lateral moraines
  !  based on width of debris concentration near sides of glacier at high elevations
  lat_moraine_wid = lat_moraine_wid0
  if(lat_moraine_wid <= 0._rkind)then
    lat_moraine_wid = 0._rkind
    do i = 1, nx
      do j = 1, ny
        if(S(i,j) > ELA .and. S(i,j) > B(i,j) .and. S(i,j) > 0._rkind)then
          if(distance(i,j) > lat_moraine_wid)then
            lat_moraine_wid = distance(i,j)
          endif
        endif
      enddo
    enddo
  endif

  ! Create the side mask based on the maximum distance or currently having debris 
  !   and where S>B (has glacier) and S<ELA (in ablation zone)
  emergingDebris = merge(1, 0, (distance<=lat_moraine_wid .or. debris>0._rkind) .and. S>B .and. S<ELA)

end function emergingDebris

! ************************************************************************************************
! private functions to assist in the glacier flow calculations
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

function H_index(h1, h2)
  implicit none
  real(rkind), intent(in) :: h1(:,:), h2(:,:)
  real(rkind), allocatable :: H_index(:,:)

  allocate(H_index(size(h1,1), size(h1,2)))
  H_index = 0.5_rkind * (h1 + h2)
end function H_index

function H_plus(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: H_plus(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(H_plus(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if(limiter == "minmod")then
    H_plus = H - 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if(limiter == "superbee")then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(Hp - H, ones, mask)
    H_plus = merge(H - 0.5_rkind * superbee(abs((H - Hm) / divisor)) * (Hp - H), H, mask)
  endif
  deallocate(mask,divisor,ones)
end function H_plus

function H_min(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: H_min(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(H_min(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  if(limiter == "minmod")then
    H_min = H + 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if(limiter == "superbee")then
    mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
    divisor = merge(H - Hm, ones, mask)
    H_min = merge(H + 0.5_rkind * superbee(abs((Hp - H) / divisor)) * (H - Hm), H, mask)
  endif
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
! private subroutine updateS_volAreaScaling: update S and debris velocity
!   This is a kluge but needed to keep track of volume and move the debris.
!   Also preserves the ability to turn back into a glacier from a glacieret.
! ************************************************************************************************
subroutine updateS_volAreaScaling(S, S0, B, debris, glacierMask, delVol, delArea, thick4area, t_total, m_dot, nx, ny, dx, dy, excessVol, grad_debris)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S0(nx,ny), B(nx,ny), delVol, delArea, thick4area, t_total, m_dot(nx,ny), dx, dy
  integer(i4b), intent(in) :: glacierMask(nx,ny), nx, ny
  real(rkind), intent(inout) :: S(nx,ny), debris(nx,ny)
  real(rkind), intent(out) :: grad_debris(nx,ny), excessVol

  ! Local variables
  real(rkind) :: cumulativeArea, cumulativeVol, around(8), area
  real(rkind) :: div_k(nx,ny), div_l(nx,ny), div_klp(nx,ny), div_klm(nx,ny), div_kpl(nx,ny), div_kml(nx,ny)
  real(rkind), allocatable :: sortedElev(:)
  integer(i4b), allocatable :: sortedIndices(:)
  integer(i4b) :: i,j,k,l,kp,km,lp,lm, numCells
  integer(i4b) :: k2(nx), l2(ny), kp2(nx), km2(nx), lp2(ny), lm2(ny)

  numCells = nx * ny
  allocate(sortedElev(numCells), sortedIndices(numCells))
  ! Flatten the S array and sort it in descending order
  sortedElev = reshape(S, [numCells])
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

  ! Add or remove cells to/from the mask until the cumulative area matches delArea
  cumulativeArea = 0._rkind
  cumulativeVol = 0._rkind
  if(delArea > 0._rkind)then ! area gain, add to high elevation areas
    do i = 1, numCells
      ! Get the 2D indices from the 1D index
      k = mod(sortedIndices(i) - 1, nx) + 1_i4b
      l = (sortedIndices(i) - 1) / nx + 1_i4b
      kp = min(k+1,nx)
      km = max(k-1,1_i4b)
      lp = min(l+1,ny)
      lm = max(l-1,1_i4b)
      ! Add the cell to the mask if in glacier area and touching an S-B>0 cell
      if(glacierMask(k,l) == 1_i4b .and. S(k,l)-B(k,l)==0._rkind)then
        around = (/S(kp,l)-B(kp,l), S(km,l)-B(km,l), S(k,lp)-B(k,lp), S(k,lm)-B(k,lm),  &
                   S(kp,lp)-B(kp,lp), S(km,lp)-B(km,lp), S(kp,lm)-B(kp,lm), S(km,lm)-B(km,lm)/)
        if(all(around<=0._rkind)) cycle ! not touching any S-B>0 cells, skip
        cumulativeArea = cumulativeArea + dx * dy
        if(cumulativeArea <= delArea)then
          S(k,l) = S(k,l) + sum(around)/count(around>0._rkind) ! add height from average of touching cells
          cumulativeVol = cumulativeVol + (S(k,l)-B(k,l)) * dx * dy
        endif
      endif
      if(cumulativeArea >= delArea) exit
    enddo
    excessVol = delVol - cumulativeVol
  else ! area loss, remove from low elevation areas
    do i = numCells, 1, -1
      ! Get the 2D indices from the 1D index
      k = mod(sortedIndices(i) - 1, nx) + 1_i4b
      l = (sortedIndices(i) - 1) / nx + 1_i4b
      ! Remove the cell from the mask if in glacier area
      if(S(k,l)-B(k,l)>0)then
        cumulativeArea = cumulativeArea + dx * dy
        if(cumulativeArea <= -delArea)then
          cumulativeVol = cumulativeVol + (S(k,l)-B(k,l)) * dx * dy
          S(k,l) = B(k,l) ! remove volume down to bedrock
        endif
      endif
      if(cumulativeArea >= -delArea) exit
    enddo
    excessVol = delVol + cumulativeVol
  endif
  deallocate(sortedElev, sortedIndices)
  
  ! distribute excess (+/-) volume evenly over all S-B>0 cells
  area = sum(merge(glacierMask,0_i4b,(S-B)>thick4area)) * dx * dy ! new area, m2
  if(area > 0._rkind)then
    S = merge( max(B + thick4area, S + excessVol/area), S, S - B > 0._rkind)
  else
    S = B ! no area left, set to bedrock
    debris = 0._rkind ! no debris
    grad_debris = 0._rkind
    excessVol = 0._rkind
    return
  endif

  ! calculate debris velocity term with a simple finite difference
  if(sum(debris)>0._rkind)then
    ! y direction indices
    l2 = [(i, i=1, ny)]
    lp2 = [(i, i=2, ny), ny]
    lm2 = [1, (i, i=1, ny-1)]

    ! x direction indices
    k2 = [(i, i=1, nx)]
    kp2 = [(i, i=2, nx), nx]
    km2 = [1, (i, i=1, nx-1)]

    ! flux divergence estimate, m s-1
    div_klp = (S(k2 ,lp2) - S0(k2 ,lp2))/t_total - m_dot(k2 ,lp2)
    div_klm = (S(k2 ,lm2) - S0(k2 ,lm2))/t_total - m_dot(k2 ,lm2)
    div_kpl = (S(kp2,l2 ) - S0(kp2,l2 ))/t_total - m_dot(kp2,l2 )
    div_kml = (S(km2,l2 ) - S0(km2,l2 ))/t_total - m_dot(km2,l2 )
    div_l = (div_klp - div_klm) / dy
    div_k = (div_kpl - div_kml) / dx
  
    grad_debris = debris_velocity(debris, div_k, div_l, k2, kp2, km2, l2, lp2, lm2, dx, dy, nx, ny)
  else
    grad_debris = 0._rkind
  endif

end subroutine updateS_volAreaScaling

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