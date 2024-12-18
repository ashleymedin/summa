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

module glacFlow_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing      ! missing integer number
USE globalData,only:realMissing         ! missing double precision number

! define data types
USE var_lookup,only:iLookGRID           ! named variables for the glacier grid information
USE data_types,only:&
                    glac_info,       & ! glacier information data structure
                    grid_info,       & ! glacier grid info data structure
                    grid_double        ! x%gru(:)%grid(:)%var(:)%dat2(:,:) (dp)

USE multiconst,only:&
                    secprday,      & ! seconds per day
                    gravity,       & ! gravitational acceleration    (m s-2)
                    iden_water,    & ! intrinsic density of water    (kg m-3)
                    iden_ice         ! intrinsic density of ice      (kg m-3)

implicit none
! define solver settings
character(len=32),parameter :: limiter='superbee' !'minmod'
character(len=32),parameter :: bed_shape='parabolic' !'trapezoid'
character(len=32),parameter :: method='MUSCL' !'upstream'

! privacy
private::run_year,diffusion_MUSCL,diffusion_upstream,minmod,superbee,flux,SIA,H_index,H_plus,H_min
public ::glacFlow
contains
! ************************************************************************************************
! public subroutine glacFlow: flow of glacier to get new glacier area and elevation
! NOTE: This will eventually run in parallel as a program, but for now it is serial
! ************************************************************************************************
subroutine glacFlow(&
                    ! model control
                    nDOM,                    & ! intent(in):    number of glacier domains
                    has_clean,               & ! intent(in):    flag for clean ice
                    has_debris,              & ! intent(in):    flag for debris cover
                    hruInd,                  & ! intent(in):    hruInd of each glacier domain
                    ! glacier topography      
                    nGlacier,                & ! intent(in):    number of glaciers
                    glacInfo,                & ! intent(in):    information for each glacier
                    gridInfo,                & ! intent(in):    information for each grid
                    gridData,                & ! intent(inout): grid data for each glacier
                    ! mass balance
                    nYears,                  & ! intent(in):    number of years to run
                    massChange,              & ! intent(in):    rate of change in glacier water equivalent (kg m-2 s-1) in each glacier domain over the nYears
                    elev,                    & ! intent(inout): elevation of each glacier domain (m)
                    debris_thick_dom,        & ! intent(out):   debris thickness in glacier domain (m)
                    ! area
                    glacAblArea,             & ! intent(out):   per glacier ablation area (m2)
                    glacAccArea,             & ! intent(out):   per glacier accumulation area (m2)
                    area,                    & ! intent(out):   area of each domain (m2)
                    ! error handling
                    err, message)              ! intent(out):   error control
   ! ---------------------------------------------------------------------------------------------
  implicit none
  ! model control
  integer(i4b), intent(in)           :: nDOM                      ! number of glacier domains
  logical, intent(in)                :: has_clean                 ! flag for clean ice
  logical, intent(in)                :: has_debris                ! flag for debris cover
  integer(i8b), intent(in)           :: hruInd(:)                 ! hruInd of each glacier domain
  ! glacier topograpy
  integer(i4b), intent(in)           :: nGlacier                  ! number of glaciers
  type(glac_info), intent(in)        :: glacInfo(:)               ! information for each glacier
  type(grid_info), intent(in)        :: gridInfo(:)               ! information for each grid
  type(grid_double), intent(inout)   :: gridData                  ! data for each grid
   ! mass balance, realMissing value if domain is missing (i.e. a glacier does not have one of ablation or accumulation)
  integer(i4b), intent(in)           :: nYears                    ! number of years to run
  real(rkind), intent(in)            :: massChange(:)             ! rate of change in glacier water equivalent (kg m-2 s-1) in each glacier domain over the nYears
  real(rkind), intent(inout)         :: elev(:)                   ! elevation of each glacier domain (m)
  real(rkind), intent(out)           :: debris_thick_dom(nDOM)    ! debris thickness in glacier domain (m)
  ! area 
  real(rkind), intent(out)           :: glacAblArea(nGlacier)     ! per glacier ablation area (m2)
  real(rkind), intent(out)           :: glacAccArea(nGlacier)     ! per glacier accumulation area (m2)
  real(rkind), intent(out)           :: area(nDOM)                ! area of each glacier domain (m2)
  integer(i4b),intent(out)           :: err                       ! error code
  character(*),intent(out)           :: message                   ! error message 
  ! locals
  real(rkind), allocatable           :: surface(:,:)              ! surface elevation of each glacier domain (m)
  real(rkind), allocatable           :: bed(:,:)                  ! bed elevation of each glacier domain (m)
  integer(i4b), allocatable          :: cell2hru(:,:)             ! map of glacier cell to hru index
  integer(i4b), allocatable          :: glacierMask(:,:)          ! 1-0 mask of glacier domain
  real(rkind)                        :: debris_vol_hru(nDOM)      ! debris volume in each HRU (m3)
  integer(i4b)                       :: i,j,k,n,iGlac,iGrid,iDOM  ! loop indices
  integer(i4b)                       :: i1, i2                    ! indices for mass balance interpolation
  real(rkind)                        :: slope, intercept          ! slope and intercept for linear extrapolation of mass balance
  real(rkind)                        :: sum_mass                  ! sum of mass balance values for averaging
  integer(i4b)                       :: nx, ny                    ! number of grid cells in x and y directions
  real(rkind)                        :: dx, dy                    ! grid cell size in x and y directions
  real(rkind)                        :: stage                     ! fraction of glacier ablation zone covered by debris
  real(rkind)                        :: debris_thick              ! debris thickness for each glacier ablation zone (m)
  real(rkind)                        :: volume                    ! volume of each glacier km3
  real(rkind)                        :: totVolume                 ! total volume of all glaciers km3, might want to send out of routine
  real(rkind), allocatable           :: mb(:,:)                   ! mass balance in each glacier cell over the nYears (m s-1)
  real(rkind)                        :: ELA_elev                  ! ELA in each glacier domain (m s-1)
  real(rkind)                        :: ELA_use                   ! ELA used for domain calculations
  real(rkind), allocatable           :: hgt(:,:)                  ! height of each glacier domain (m)
  integer(i4b)                       :: validCount                ! number of valid points
  real(rkind), allocatable           :: validElev(:)              ! filter out points where equal to realMissing
  real(rkind), allocatable           :: validMassChange(:)        ! filter out points where equal to realMissing
  real(rkind)                        :: temp_elev, temp_mass      ! temporary variables for sorting
  integer(i4b), allocatable          :: glacid_to_index(:)        ! mapping array from glacier id to index in gridInfo
  real(rkind),parameter              :: thick4area=0.1_rkind      ! an arbitrary small threshold for glacier thickness to be considered as glacier area
  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='glacFlow/'
  glacAblArea = 0._rkind
  glacAccArea = 0._rkind
  totVolume = 0._rkind
  validCount = count(elev /= realMissing)
  allocate(validElev(validCount))
  allocate(validMassChange(validCount))

  ! Filter out points where no area and elevation is missing
  j = 1
  do iDOM = 1, nDOM
    if (elev(iDOM) /= realMissing) then
      validElev(j) = elev(iDOM)
      validMassChange(j) = massChange(iDOM)
      j = j + 1
    end if
  end do
  
  ! Simple bubble sort to sort the valid elevations in descending order
  do i = 1, validCount-1
    do j = 1, validCount-i
      if (validElev(j) < validElev(j+1)) then
        ! Swap elevations
        temp_elev = validElev(j)
        validElev(j) = validElev(j+1)
        validElev(j+1) = temp_elev
        ! Swap corresponding massChange values
        temp_mass = validMassChange(j)
        validMassChange(j) = validMassChange(j+1)
        validMassChange(j+1) = temp_mass
      end if
    end do
  end do

  ! Combine points with the same elevation into a single point that is the average of those points
  i = 1
  do while (i <= validCount)
    n = 1
    sum_mass = validMassChange(i)
    do j = i+1, validCount
      if (validElev(i) == validElev(j)) then
        sum_mass = sum_mass + validMassChange(j)
        n = n + 1
      else
        exit
      end if
    end do
    if (n > 1) then
      validMassChange(i) = sum_mass / n
      ! Shift remaining elements to the left
      do k = i+1, validCount-n
        validElev(k) = validElev(k+n)
        validMassChange(k) = validMassChange(k+n)
      end do
      validCount = validCount - n + 1
    end if
    i = i + 1
  end do
  
  ! Calculate the elevation where massChange would be zero, always has 2 points at least
  if (validMassChange(1) <= 0.0) then
    ! Extrapolate above the first point
    i1 = 1
    i2 = 2
  else if (validMassChange(validCount) >= 0.0) then
    ! Extrapolate below the last point
    i1 = validCount - 1
    i2 = validCount
  else
    do i = 1, validCount-1
      if ((validMassChange(i) > 0.0 .and. validMassChange(i+1) < 0.0) .or. &
          (validMassChange(i) < 0.0 .and. validMassChange(i+1) > 0.0)) then
        ! Linear interpolation to find the zero crossing
        i1 = i
        i2 = i + 1
        exit
      end if
    end do
  end if
  ! Calculate slope and intercept for ELA
  slope = (validMassChange(i2) - validMassChange(i1)) / (validElev(i2) - validElev(i1))
  intercept = validMassChange(i1) - slope * validElev(i1)    
  ELA_elev = -intercept / slope
 
  ! Initialize new domain areas and elevations
  do i = 1,nDOM
    area(i) = 0._rkind
    elev(i) = 0._rkind
    debris_vol_hru(i) = 0._rkind
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
    stage = gridInfo(iGrid)%stage
    debris_thick = gridInfo(iGrid)%debris_thick

    ! set up mass balance and height arrays
    allocate(mb(nx,ny), hgt(nx,ny))
    mb = 0._rkind
    hgt = 0._rkind

    ! set up grid data
    allocate(surface(nx,ny), bed(nx,ny), cell2hru(nx,ny), glacierMask(nx,ny))
    surface = gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny)
    bed = gridData%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)
    cell2hru = int(gridData%grid(iGrid)%var(iLookGRID%cell2hru)%dat2(1:nx,1:ny))
    glacierMask = int(gridData%grid(iGrid)%var(iLookGRID%glacierMask)%dat2(1:nx,1:ny))

    ! distribute mass balance over surface, using all points in GRU
    do k = 1, nx
      do j = 1, ny
        ! Set mass balance to zero if not on glacier
        if (glacierMask(k,j)==0) then
          mb(k,j) = 0._rkind
        else
          if (validCount==1) then
            mb(k,j) = validMassChange(1)
          else ! find the interval in which surface(k,j) lies
            if (surface(k,j) >= validElev(1)) then
              ! Extrapolate above the first point
              i1 = 1
              i2 = 2
            else if (surface(k,j) <= validElev(validCount)) then
              ! Extrapolate below the last point
              i1 = validCount - 1
              i2 = validCount
            else
              ! Interpolate between points
              do n = 1, validCount-1
                if (surface(k,j) <= validElev(n) .and. surface(k,j) >= validElev(n+1)) then
                  i1 = n
                  i2 = n + 1
                  exit
                end if
              end do
            end if
            ! Calculate slope and intercept for mass balance
            slope = (validMassChange(i2) - validMassChange(i1)) / (validElev(i2) - validElev(i1))
            intercept = validMassChange(i1) - slope * validElev(i1)
            
            mb(k,j) = slope * surface(k,j) + intercept 
          end if
          !mb(k,j) = 3.0_rkind/365.0_rkind/secprday/900.0_rkind ! 3 kg/m2 ice /year to m s-1 to match the test case
        end if ! end of glacier mask check
      end do ! end of ny loop
    end do ! end of nx loop

    ! convert mass balance to m s-1 from kg m-2 s-1
    mb = mb / iden_water

    ! compute flow
    call run_year(nYears,surface, bed, mb, glacierMask, nx, ny, dx, dy, volume)
    totVolume = totVolume+volume
    
    ! Initialize variables
    hgt = surface - bed

    ! Force glacier to have an accumulation area and ablation area, even if they are tiny, so that the domain is always defined
    ELA_use = ELA_elev
    if (ELA_elev > maxval(merge(surface,0._rkind,hgt>thick4area))) then
      ELA_use = maxval(merge(surface,0._rkind, hgt>thick4area))-0.0001_rkind
    else if (ELA_elev < minval(merge(surface,1.e6_rkind, hgt>thick4area))) then
      ELA_use = minval(merge(surface,1.e6_rkind, hgt>thick4area))+0.0001_rkind
    end if
    
    ! Calculate glacier accumulation and ablation areas
    glacAccArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface>= ELA_use))*dx*dy
    glacAblArea(iGlac) = sum(merge(glacierMask, 0_i4b, hgt>thick4area .and. surface<  ELA_use))*dx*dy
    print*, 'glacAblArea = ', glacAblArea(iGlac), ' glacAccArea = ', glacAccArea(iGlac),ELA_use

    ! Loop through HRUs and calculate domain areas and elevations for each HRU
    ! Order of domains will go HRU 1: Acc, Abl clean, Abl debris HRU 2: Acc, Abl clean, Abl debris etc.
    ! It is possible one of the ablation domains is missing if there is no clean or debris cover
    n = 1
    if (has_clean)  n=n+1
    if (has_debris) n=n+1
    do k = 1, nDOM / n
      area(n*(k-1)+1) = area(n*(k-1)+1) + sum(merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface>=ELA_use))*dx*dy
      if (has_clean) area(n*(k-1)+2) = area(n*(k-1)+2) + sum(merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use))*dx*dy*(1.0_rkind-stage)
      if (has_debris) area(n*(k-1)+n) = area(n*(k-1)+n) + sum(merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use))*dx*dy*stage
      
      ! Calculate mean elevation for accumulation and ablation areas
      if (count(cell2hru==hruInd(k) .and. hgt>thick4area .and. surface>=ELA_use) > 0) then
        elev(n*(k-1)+1) = elev(n*(k-1)+1) + sum(surface * merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface>=ELA_use)) / &
                          count(cell2hru==hruInd(k) .and. hgt>thick4area .and. surface>=ELA_use) /nGlacier
      end if
      if (count(cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use) > 0) then
        elev(n*(k-1)+n)  = elev(n*(k-1)+n) + sum(surface * merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use)) / &
                          count(cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use) /nGlacier
        if (has_clean) elev(n*(k-1)+2) =  elev(n*(k-1)+n)  ! only one of debris or clean, then already set        
      end if
      ! thickness of average debris cover in HRU changes with area change if there is more than one glacier in an HRU and they have different debris thickness 
      debris_vol_hru(n*(k-1)+n) = debris_vol_hru(n*(k-1)+n) + sum(merge(glacierMask, 0_i4b, cell2hru==hruInd(k) .and. hgt>thick4area .and. surface< ELA_use))*dx*dy*stage*debris_thick
    end do
    ! update gridData and deallocate
    gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) = surface
    deallocate(mb, hgt, surface, bed, cell2hru, glacierMask)

  enddo ! end of glacier loop

  ! Set elevations to realMissing if no area in domain
  do iDOM = 1,nDOM
    if (elev(iDOM)==0._rkind) elev(iDOM)=realMissing
    if (has_debris) then ! debris thickness is calculated for each HRU
      if (area(iDOM)>0._rkind)then 
        debris_thick_dom(iDOM) = debris_vol_hru(iDOM) / area(iDOM) ! if ablation clean or accumulation, will be 0
      else
        debris_thick_dom(iDOM) = 0._rkind ! will skip this domain
      end if
    else
      debris_thick_dom(iDOM) = 0._rkind
    end if
  end do
  print*,"elev ", elev, "area", area, "massChange", massChange, "debris_thick_dom", debris_thick_dom

  deallocate(glacid_to_index, validElev, validMassChange)

end subroutine glacFlow


! ************************************************************************************************
! private subroutine run_year for setting up flow model and running for each glacier for time period
! ************************************************************************************************
  subroutine run_year(y_end, S, B, m_dot, glacierMask, nx, ny, dx, dy, volume)
    ! Arguments
    real(rkind), intent(in) :: dx, dy
    real(rkind), intent(inout) :: S(nx,ny), B(nx,ny), m_dot(nx,ny)
    integer(i4b), intent(in) :: y_end, ny, nx, glacierMask(nx,ny)
    real(rkind), intent(out) :: volume

    ! Local variables
    real(rkind) :: t_total, dt, max_dt, min_dt, deltat, div_q(nx,ny), dt_cfl, meanS
    real(rkind) :: gamma, t
    integer(i4b) :: i
    integer(i4b),parameter :: n=3 ! Glen's flow law exponent,
    real(rkind),parameter :: A=2.4e-24 ! Modern Glen parameter
    real(rkind),parameter :: cfl= 0.124 ! Courant-Friedrichs-Lewy condition
    integer(i4b) :: l(ny), lp(ny), lpp(ny), lm(ny), lmm(ny)
    integer(i4b) :: k(nx), kp(nx), kpp(nx), km(nx), kmm(nx)

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

    t_total = y_end * secprday* 365._rkind
    gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
    max_dt = 31._rkind * secprday! max timestep in seconds, a month
    min_dt = secprday !0._rkind ! min timestep in seconds 
    t = 0._rkind

    do while (t < t_total)
      dt = t_total - t
      ! Select diffusion method and call step
      if (method == "MUSCL") then
        call diffusion_MUSCL(S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, dt_cfl)
      else if (method == "upstream") then
        call diffusion_upstream(S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, l, lp, lm, dx, dy, div_q, dt_cfl)
      end if
  
      ! update time step
      deltat = min(dt_cfl, dt)
      if (deltat > max_dt) deltat = max_dt
      if (deltat < min_dt) deltat = min_dt
      t = t + deltat
  
      ! Update S
      S = S + (m_dot + div_q) * deltat
      S = merge(S, B, S > B)

      ! Check that the glacier is in boundaries, fix small violations
      if (any((S - B) > 0._rkind .and. glacierMask==0)) then
        if (any((S - B) > 10.0 .and. glacierMask==0)) stop 'Glacier exceeds boundaries'
        S = merge(B, S, (S - B) > 0._rkind .and. glacierMask==0)
      end if
      ! Check that glacier surface is not infinite (unstable), bring down to mean glacier height
      ! This is a temporary fix, should be replaced with a more sophisticated method
      if (any((S - B) > 1.e6_rkind .and. glacierMask==1)) then
        meanS = sum(merge(S,0._rkind,(glacierMask==1 .and. (S - B) < 1.e6_rkind))) / count(glacierMask==1 .and. (S - B)<=1.e6_rkind)
        S = merge(S,meanS,((S - B)<=1.e6_rkind))
      endif
    end do

    ! Calculate volume, not currently used but could be useful
    volume = sum(merge(S-B,0._rkind,glacierMask==1)) * dx * dy * 1.e-9_rkind ! km3

  end subroutine run_year


! ************************************************************************************************
! private subroutines for calculating a time step of glacier flow
! ************************************************************************************************
! ************************************************************************************************
! upwind diffusion scheme
! ************************************************************************************************
subroutine diffusion_upstream(S, B, mask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, l, lp, lm, dx, dy, div_q, dt_cfl)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, k(nx), kp(nx), km(nx), l(ny), lp(ny), lm(ny), n
  real(rkind), intent(out) :: div_q(nx,ny), dt_cfl
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny), Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkml(nx,ny), Hkl(nx,ny), Hklp(nx,ny), Hklm(nx,ny)
  real(rkind) :: H_l_up(nx,ny), H_l_down(nx,ny)
  real(rkind) :: H_l_upstream_up(nx,ny), H_l_upstream_down(nx,ny)
  real(rkind) :: f_l_plus(nx,ny), f_l_min(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_up(nx,ny), H_k_down(nx,ny)
  real(rkind) :: H_k_upstream_up(nx,ny), H_k_upstream_down(nx,ny)
  real(rkind) :: f_k_plus(nx,ny), f_k_min(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: div_l(nx,ny), div_k(nx,ny)
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
  H_l_up   = H_index(Hklp, Hkl)
  H_l_down = H_index(Hkl, Hklm)
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
  H_l_upstream_down = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) > Sklm(j,i)) then
        H_l_upstream_down(j,i) = Hkl(j,i)
      else
        H_l_upstream_down(j,i) = Hklm(j,i)
      end if
    end do
  end do

  ! calculate l flux
  f_l_plus = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_min  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up = gamma * H_l_up**(n+1_i4b) * H_l_upstream_up * f_l_plus
  D_l_dn = gamma * H_l_down**(n+1_i4b) * H_l_upstream_down * f_l_min
  ! Enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k upstream 
  H_k_up   = H_index(Hkpl, Hkl)
  H_k_down = H_index(Hkl, Hkml)
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
  H_k_upstream_down = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if (Skl(j,i) > Skml(j,i)) then
        H_k_upstream_down(j,i) = Hkl(j,i)
      else
        H_k_upstream_down(j,i) = Hkml(j,i)
      end if
    end do
  end do

  ! calculate k flux
  f_k_plus = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_min  = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up = gamma * H_k_up**(n+1_i4b) * H_k_upstream_up * f_k_plus
  D_k_dn = gamma * H_k_down**(n+1_i4b) * H_k_upstream_down * f_k_min
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

end subroutine diffusion_upstream

! ************************************************************************************************
! MUSCL scheme
! ************************************************************************************************
subroutine diffusion_MUSCL(S, B, mask, gamma, n, cfl, max_dt, nx, ny, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, dt_cfl)
  implicit none
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, k(nx), kp(nx), km(nx), kpp(nx), kmm(nx), l(ny), lp(ny), lm(ny), lpp(ny), lmm(ny), n
  real(rkind), intent(out) :: div_q(nx,ny), dt_cfl
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny), Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkppl(nx,ny), Hkml(nx,ny), Hkmml(nx,ny), Hkl(nx,ny), Hklp(nx,ny), Hklpp(nx,ny), Hklm(nx,ny), Hklmm(nx,ny)
  real(rkind) :: H_l_min_up(nx,ny), H_l_plus_up(nx,ny)
  real(rkind) :: H_l_min_down(nx,ny), H_l_plus_down(nx,ny)
  real(rkind) :: f_l_plus(nx,ny), f_l_min(nx,ny)
  real(rkind) :: D_l_up_m(nx,ny), D_l_up_p(nx,ny), D_l_up_min(nx,ny), D_l_up_max(nx,ny)
  real(rkind) :: D_l_dn_m(nx,ny), D_l_dn_p(nx,ny), D_l_dn_min(nx,ny), D_l_dn_max(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_min_up(nx,ny), H_k_plus_up(nx,ny)
  real(rkind) :: H_k_min_down(nx,ny), H_k_plus_down(nx,ny)
  real(rkind) :: f_k_plus(nx,ny), f_k_min(nx,ny)
  real(rkind) :: D_k_up_m(nx,ny), D_k_up_p(nx,ny), D_k_up_min(nx,ny), D_k_up_max(nx,ny)
  real(rkind) :: D_k_dn_m(nx,ny), D_k_dn_p(nx,ny), D_k_dn_min(nx,ny), D_k_dn_max(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: div_l(nx,ny), div_k(nx,ny)
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
  H_l_min_down  = H_min(Hklmm, Hklm, Hkl)
  H_l_plus_down = H_plus(Hklm, Hkl, Hklp)

  ! calculate l flux
  f_l_plus = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_min  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up_m = gamma * H_l_min_up**(n+2_i4b) * f_l_plus   ! equation 30 Jarosh 2013
  D_l_up_p = gamma * H_l_plus_up**(n+2_i4b) * f_l_plus  ! equation 30 Jarosh 2013
  D_l_up_min = min(D_l_up_m, D_l_up_p)                  ! equation 31 Jarosh 2013
  D_l_up_max = max(D_l_up_m, D_l_up_p)                  ! equation 32 Jarosh 2013
  !
  D_l_dn_m = gamma * H_l_min_down**(n+2_i4b) * f_l_min  ! equation 30 Jarosh 2013
  D_l_dn_p = gamma * H_l_plus_down**(n+2_i4b) * f_l_min ! equation 30 Jarosh 2013
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
      if (Skl(j,i) <= Sklm(j,i) .and. H_l_min_down(j,i) <= H_l_plus_down(j,i)) then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      else if (Skl(j,i) <= Sklm(j,i) .and. H_l_min_down(j,i) > H_l_plus_down(j,i)) then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if (Skl(j,i) > Sklm(j,i) .and. H_l_min_down(j,i) <= H_l_plus_down(j,i)) then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if (Skl(j,i) > Sklm(j,i) .and. H_l_min_down(j,i) > H_l_plus_down(j,i)) then
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
  H_k_min_down  = H_min(Hkmml, Hkml, Hkl)
  H_k_plus_down = H_plus(Hkml, Hkl, Hkpl)

  ! calculate k flux
  f_k_plus = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_min  = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up_m = gamma * H_k_min_up**(n+2_i4b) * f_k_plus   ! equation 30 Jarosh 2013
  D_k_up_p = gamma * H_k_plus_up**(n+2_i4b) * f_k_plus  ! equation 30 Jarosh 2013
  D_k_up_min = min(D_k_up_m, D_k_up_p)                  ! equation 31 Jarosh 2013
  D_k_up_max = max(D_k_up_m, D_k_up_p)                  ! equation 32 Jarosh 2013
  !
  D_k_dn_m = gamma * H_k_min_down**(n+2_i4b) * f_k_min  ! equation 30 Jarosh 2013
  D_k_dn_p = gamma * H_k_plus_down**(n+2_i4b) * f_k_min ! equation 30 Jarosh 2013
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
      if (Skl(j,i) <= Skml(j,i) .and. H_k_min_down(j,i) <= H_k_plus_down(j,i)) then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      else if (Skl(j,i) <= Skml(j,i) .and. H_k_min_down(j,i) > H_k_plus_down(j,i)) then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if (Skl(j,i) > Skml(j,i) .and. H_k_min_down(j,i) <= H_k_plus_down(j,i)) then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if (Skl(j,i) > Skml(j,i) .and. H_k_min_down(j,i) > H_k_plus_down(j,i)) then
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

end subroutine diffusion_MUSCL


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

function flux(sjpl, sjml, sjplp, sjmlp, sjlp, sjl, dj, dl, n) result(flux_result)
  implicit none
  real(rkind), intent(in) :: sjpl(:,:), sjml(:,:), sjplp(:,:), sjmlp(:,:), sjlp(:,:), sjl(:,:), dj, dl
  integer(i4b), intent(in) :: n
  real(rkind), allocatable :: flux_result(:,:)

  allocate(flux_result(size(sjpl,1), size(sjpl,2)))
  flux_result = ( (sjpl - sjml + sjplp - sjmlp)**2_i4b / (4._rkind * dj)**2_i4b & 
                 + (sjlp - sjl)**2_i4b / dl**2_i4b )**((n - 1.0) / 2._rkind)
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


end module glacFlow_module