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

USE multiconst,only:&
                    secperday,    & ! seconds per day
                    gravity,      & ! gravitational acceleration    (m s-2)
                    iden_ice        ! intrinsic density of ice      (kg m-3)

character(len=32),parameter :: limiter='superbee' !'minmod'
character(len=32),parameter :: bed_shape='parabolic' !'trapezoid'
character(len=32),parameter :: method='MUSCL' !'upstream'

! privacy
implicit none
private::run_year,diffusion_MUSCL,diffusion_upstream,minmod,superbee,flux,SIA,H_index,H_plus,H_min
public ::glacFlow
contains
! ************************************************************************************************
! public subroutine glacFlow: flow of glacier to get new glacier area and elevation
! NOTE: This will eventually run in parallel as a program, but for now it is serial
! ************************************************************************************************
subroutine glacFlow(&
                    ! model control
                    gruId,              & ! intent(in):    gruId
                    nDOM,               & ! intent(in):    number of glacier domains
                    hruId,              & ! intent(in):    hruId of each glacier domain
                    ! glacier topography
                    nGlacier,           & ! intent(in):    number of glaciers
                    surface,            & ! intent(in):    surface elevation of each glacier domain (m)
                    dx, dy,             & ! intent(in):    grid spacing (m) by glacier
                    Nx0, Ny0,           & ! intent(in):    number of grid cells in x and y directions by glacier
                    nxgrid, nygrid,     & ! intent(in):    max number of grid cells in x and y directions by glacier
                    bed,                & ! intent(in):    bed elevation of each glacier domain (m)
                    cell2hru,           & ! intent(in):    map of glacier cell to hru
                    glacierMask,        & ! intent(in):    mask of glacier domain
                    ! mass balance
                    GWE_delta,          & ! intent(in):    change in glacier water equivalent (m s-1) in each glacier domain
                    elev,               & ! intent(inout): elevation of each glacier domain (m)
                    ! area
                    glacAblArea,        & ! intent(out):   per glacier ablation area (m2)
                    glacAccArea,        & ! intent(out):   per glacier accumulation area (m2)
                    area,           & ! intent(out):   area of each domain (m2)
                    elev,           & ! intent(out):   elev of each domain (m2)
                    ! error handling
                    err, message)         ! intent(out):   error control
   ! ---------------------------------------------------------------------------------------------
  implicit none
  ! model control
  integer(i8b), intent(in)    :: gruId(:)                  ! gruId
  integer(i4b), intent(in)    :: nDOM                      ! number of glacier domains
  integer(i8b), intent(in)    :: hruId(:)                  ! hruId of each glacier domain
  ! glacier topography
  integer(i4b), intent(in)    :: nGlacier                  ! number of glaciers
  real(rkind), intent(in)     :: surface(:,:,:)            ! surface elevation of each glacier domain (m)
  real(rkind), intent(in)     :: dx(:), dy(:)              ! grid spacing (m) by glacier
  integer(i4b), intent(in)    :: Nx0(:), Ny0(:)            ! number of grid cells in x and y directions by glacier
  real(rkind), intent(in)     :: nxgrid,nygrid             ! max number of grid cells in x and y directions by glacier
  real(rkind), intent(in)     :: bed(:,:,:)                ! bed elevation of each glacier domain (m)
  real(rkind), intent(in)     :: cell2hru(:,:,:)           ! map of glacier cell to hru
  real(rkind), intent(in)     :: glacierMask(:,:,:)        ! mask of glacier domain
  ! mass balance, realMissing value if domain is missing (i.e. a glacier does not have one of ablation or accumulation)
  real(rkind), intent(in)     :: GWE_delta(:)              ! change in glacier water equivalent (m s-1) in each glacier domain
  real(rkind), intent(inout)  :: elev(:)                   ! elevation of each glacier domain (m)
  ! area 
  real(rkind), intent(out)    :: glacAblArea(nGlacier)     ! per glacier ablation area (m2)
  real(rkind), intent(out)    :: glacAccArea(nGlacier)     ! per glacier accumulation area (m2)
  real(rkind), intent(out)    :: area(nDOM)                ! area of each glacier domain (m2)
  integer(i4b),intent(out)    :: err                       ! error code
  character(*),intent(out)    :: message                   ! error message 
  ! locals
  integer(i4b)                :: i,j,k                     ! loop indices
  integer(i4b),parameter      :: nYr = 1                   ! number of years to run
  real(rkind)                 :: Ny, Nx                    ! number of grid cells in x and y directions
  real(rkind)                 :: volume                    ! volume of each glacier km3
  real(rkind)                 :: totVolume                 ! total volume of all glaciers km3, might want to send out of routine
  real(rkind)                 :: mb(nygrid,nxgrid)         ! mass balance in each glacier domain (m s-1)
  real(rkind)                 :: ELA_elev                  ! ELA in each glacier domain (m s-1)
  real(rkind)                 :: hgt(nygrid,nxgrid)          ! height of each glacier domain (m)
  real(rkind)                 :: dly                       ! area of each grid cell (m2)
  integer(i4b)                :: validCount                ! number of valid points
  real(rkind), allocatable    :: validElev(:)              ! filter out points where equal to realMissing
  real(rkind), allocatable    :: validGWE_delta(:)         ! filter out points where equal to realMissing
  real(rkind)                 :: temp_elev, temp_GWE       ! temporary variables for sorting
  real(rkind)                 :: slope, intercept          ! slope and intercept for linear extrapolation of GWE
  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='glacFlow/'
  glacAblArea = 0._rkind
  glacAccArea = 0._rkind
  totVolume = 0._rkind
  validCount = count(elev /= realMissing)
  allocate(validElev(validCount))
  allocate(validGWE_delta(validCount))
  
  ! Filter out points where no area and elevation is missing
  j = 1
  do i = 1, nDOM
    if (elev(i) /= realMissing) then
      validElev(j) = elev(i)
      validGWE_delta(j) = GWE_delta(i)
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
        ! Swap corresponding GWE_delta values
        temp_GWE = validGWE_delta(j)
        validGWE_delta(j) = validGWE_delta(j+1)
        validGWE_delta(j+1) = temp_GWE
      end if
    end do
  end do
  
  ! Calculate the elevation where GWE_delta would be zero
  if (validGWE_delta(1) <= 0.0) then
    ! Extrapolate above the first point
    slope = (validGWE_delta(2) - validGWE_delta(1)) / (validElev(2) - validElev(1))
    intercept = validGWE_delta(1) - slope * validElev(1)
    ELA_elev = -intercept / slope
  else if (validGWE_delta(validCount) >= 0.0) then
    ! Extrapolate below the last point
    slope = (validGWE_delta(validCount) - validGWE_delta(validCount-1)) / (validElev(validCount) - validElev(validCount-1))
    intercept = validGWE_delta(validCount) - slope * validElev(validCount)
    ELA_elev = -intercept / slope
  else
    do i = 1, validCount-1
      if ((validGWE_delta(i) > 0.0 .and. validGWE_delta(i+1) < 0.0) .or. (validGWE_delta(i) < 0.0 .and. validGWE_delta(i+1) > 0.0)) then
        ! Linear interpolation to find the zero crossing
        slope = (validGWE_delta(i+1) - validGWE_delta(i)) / (validElev(i+1) - validElev(i))
        intercept = validGWE_delta(i) - slope * validElev(i)
        ELA_elev = -intercept / slope
        exit
      end if
    end do
  end if
  deallocate(validElev, validGWE_delta)

  ! Initialize domain areas and elevations
  do i = 1,nDOM
    area(i) = 0._rkind
    elev(i) = 0._rkind
  end do

  ! run flow for each glacier
  do i = 1,nGlacier
    ! set up glacier
    Nx = Nx0(i)
    Ny = Ny0(i)
    mb = 0._rkind
    hgt = 0._rkind

    ! distribute mass balance over surface, using all points in GRU
    do j = 1, Nx
      do k = 1, Ny
        ! Find the interval in which surface(k,j,i) lies
        if (surface(k,j,i) >= elev(1)) then
            ! Extrapolate above the first point
          slope = (GWE_delta(2) - GWE_delta(1)) / (elev(2) - elev(1))
          intercept = GWE_delta(1) - slope * elev(1)
          mb(k,j) = slope * surface(k,j,i) + intercept
        else if (surface(j) <= elev(nDOM)) then
          ! Extrapolate below the last point
          slope = (GWE_delta(nDOM) - GWE_delta(nDOM-1)) / (elev(nDOM) - elev(nDOM-1))
          intercept = GWE_delta(nDOM) - slope * elev(nDOM)
          mb(k,j) = slope * surface(k,j,i) + intercept
        else
          ! Interpolate between points
          do n = 1, nDOM-1
            if (surface(k,j,i) <= elev(n) .and. surface(k,j,i) >= elev(n+1)) then
              slope = (GWE_delta(n+1) - GWE_delta(n)) / (elev(n+1) - elev(n))
              intercept = GWE_delta(n) - slope * elev(n)
              mb(k,j) = slope * surface(k,j,i) + intercept
              exit
            end if
          end do
        end if
        ! Set mass balance to zero if not on glacier
        if (glacierMask(k,j,i)==0) then
          mb(k,j) = 0._rkind
        end if
      end do
    end do

    ! compute flow
    call run_year(nYr,surface(1:Ny,1:Nx,i), bed(1:Ny,1:Nx,i), mb(1:Ny,1:Nx), glacierMask(1:Ny,1:Nx,i), Ny, Nx, dx(i), dy(i), volume)
    totVolume = totVolume+volume
    
    ! Initialize variables
    hgt(1:Ny,1:Nx) = surface(1:Ny,1:Nx,i) - bed(1:Ny,1:Nx,i)
    dly = dx(i) * dy(i)
    
    ! Calculate glacier ablation and accumulation areas
    glacAblArea(i) = sum(merge(dly, 0.0, (hgt(1:Ny,1:Nx)>0) .and. (surface(1:Ny,1:Nx,i)<  ELA_elev)))
    glacAccArea(i) = sum(merge(dly, 0.0, (hgt(1:Ny,1:Nx)>0) .and. (surface(1:Ny,1:Nx,i)>= ELA_elev)))
    
    ! Loop through HRUs and calculate domain areas and elevations
    ! Order of domains will go HRU 1: Acc, Abl, HRU 2: Acc, Abl, etc.
    do k = 1, nDOM / 2
      area(2*k-1) = area(2*k-1) + sum(merge(dly, 0.0, (cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (hgt(1:Ny,1:Nx)>0) .and. (surface(1:Ny,1:Nx,i)>= ELA_elev)))
      area(2*k)   = area(2*k)   + sum(merge(dly, 0.0, (cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (hgt(1:Ny,1:Nx)>0) .and. (surface(1:Ny,1:Nx,i)<  ELA_elev)))
      
      ! Calculate mean elevation for accumulation and ablation areas
      if (count((cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i)>=ELA_elev)) > 0) then
        elev(2*k-1) = elev(2*k-1) + sum(surface(1:Ny,1:Nx,i) * merge(1.0, 0.0, (cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i) >= ELA_elev))) / &
                          count((cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i)>=ELA_elev))
      end if
      if (count((cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i) <ELA_elev)) > 0) then
        elev(2*k)   = elev(2*k)   + sum(surface(1:Ny,1:Nx,i) * merge(1.0, 0.0, (cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i) < ELA_elev))) / &
                          count((cell2hru(1:Ny,1:Nx,i)==hruId(k)) .and. (surface(1:Ny,1:Nx,i)< ELA_elev))
      end if
    end do
  enddo ! end of glacier loop

  ! Set elevations to realMissing if no area
  do i = 1,nDOM
    if (elev(i)==0._rkind) elev(i)=realMissing
  end do

end subroutine glacFlow


! ************************************************************************************************
! private subroutine run_year for setting up flow model and running for each glacier for time period
! ************************************************************************************************
  subroutine run_year(y_end, S, B, m_dot, glacierMask, Ny, Nx, dx, dy, volume)
    ! Arguments
    real(rkind), intent(in) :: y_end, dx, dy
    real(rkind), intent(inout) :: S(:,:), B(:,:), m_dot(:,:)
    integer(i4b), intent(in) :: Nx, Ny
    logical(lgt), intent(in) :: glacierMask(:)
    real(rkind), intent(out) :: volume

    ! Local variables
    real(rkind) :: t_total, dt, max_dt, min_dt, deltat
    integer(i4b) :: i, j
    integer(i4b),parameter :: n=3 ! Glen's flow law exponent,
    real(rkind),parameter :: A=2.4e-24 ! Modern Glen parameter
    real(rkind),parameter :: cfl=0.124 ! Courant-Friedrichs-Lewy condition
    integer(i4b) :: l(Nx), lp(Nx), lpp(Nx), lm(Nx), lmm(Nx)
    integer(i4b) :: k(Ny), kp(Ny), kpp(Ny), km(Ny), kmm(Ny)

  ! x direction indices
    l = [(i, i=0, Nx-1)]
    lp = [(i, i=1, Nx-1), Nx-1]
    lpp = [(i, i=2, Nx-1), Nx-1, Nx-1]
    lm = [0, (i, i=0, Nx-2)]
    lmm = [0, 0, (i, i=0, Nx-3)]
    if (Nx == 1) then
      lpp = [0]
      lmm = [0]
    end if
    
    ! y direction indices
    k = [(i, i=0, Ny-1)]
    kp = [(i, i=1, Ny-1), Ny-1]
    kpp = [(i, i=2, Ny-1), Ny-1, Ny-1]
    km = [0, (i, i=0, Ny-2)]
    kmm = [0, 0, (i, i=0, Ny-3)]
    if (Ny == 1) then
      kpp = [0]
      kmm = [0]
    end if

    t_total = y_end * secprday* 365._rkind
    gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
    max_dt = 31._rkind * secprday! max timestep in seconds, a month
    min_dt = 0._rkind !min timestep in seconds 
    t = 0._rkind
    do while (t < t_total)
      dt = t_total - t

      ! Select diffusion method and call step
      if (method == "MUSCL") then
        call diffusion_MUSCL(S, B, gamma, n, cfl, max_dt, Ny, Nx, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, dt_cfl)
      else if (method == "upstream") then
        call diffusion_upstream(S, B, gamma, n, cfl, max_dt, Ny, Nx, k, kp, km, l, lp, lm, dx, dy, div_q, dt_cfl)
      end if
  
      ! update time step
      deltat = min(dt_cfl, dt)
      if (deltat > max_dt) deltat = max_dt
      if (deltat < min_dt) deltat = min_dt
      t = t + deltat
  
      ! Update S
      S = S + (m_dot + div_q) * deltat
      S = max(S, B)
  
      ! Check that the glacier is in boundaries, fix small violations
      if (any((S - B)(.not. glacierMask) > 0._rkind)) then
        if (any((S - B)(.not. glacierMask) > 10.0)) stop 'Glacier exceeds boundaries'
        S(.not. glacierMask) = B(.not. glacierMask)
      end if
    end do

    volume = sum(S(glacierMask)-B(glacierMask)) * dx * dy * 1.e-9_rkind ! km3

  end subroutine run_year


! ************************************************************************************************
! private subroutines for calculating a time step of glacier flow
! ************************************************************************************************
! ************************************************************************************************
! upwind diffusion scheme
! ************************************************************************************************
subroutine diffusion_upstream(S, B, gamma, n, cfl, max_dt, Ny, Nx, k, kp, km, l, lp, lm, dx, dy, div_q, dt_cfl)
  implicit none
  real(rkind), intent(in) :: S(:,:), B(:,:), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: Ny, Nx, k(:), kp(:), km(:), l(:), lp(:), lm(:), n
  real(rkind), intent(out) :: div_q(Ny, Nx), dt_cfl
  real(rkind) :: H(Ny, Nx)
  real(rkind) :: Sklp(Ny), Sklm(Ny), Skplp(Ny), Skplm(Ny), Skpl(Ny), Skl(Ny), Skmlp(Ny), Skmlm(Ny), Skml(Ny)
  real(rkind) :: Hkpl(Ny), Hkml(Ny), Hkl(Ny), Hklp(Ny), Hklm(Ny)
  real(rkind) :: H_l_up(Ny, Nx), H_l_down(Ny, Nx)
  real(rkind) :: H_l_upstream_up(Ny, Nx), H_l_upstream_down(Ny, Nx)
  real(rkind) :: f_l_plus(Ny, Nx), f_l_min(Ny, Nx)
  real(rkind) :: D_l_up(Ny, Nx), D_l_dn(Ny, Nx)
  real(rkind) :: H_k_up(Ny, Nx), H_k_down(Ny, Nx)
  real(rkind) :: H_k_upstream_up(Ny, Nx), H_k_upstream_down(Ny, Nx)
  real(rkind) :: f_k_plus(Ny, Nx), f_k_min(Ny, Nx)
  real(rkind) :: D_k_up(Ny, Nx), D_k_dn(Ny, Nx)
  real(rkind) :: div_l(Ny, Nx), div_k(Ny, Nx)
  real(rkind) :: divisor
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography  
  
  ! indices
  Sklp = S(k, :, lp)
  Sklm = S(k, :, lm)
  Skplp = S(kp, :, lp)
  Skplm = S(kp, :, lm)
  Skpl = S(kp, :, l)
  Skl = S(k, :, l)
  Skmlp = S(km, :, lp)
  Skmlm = S(km, :, lm)
  Skml = S(km, :, l)
  !
  Hkpl = H(kp, :, l)
  Hkml = H(km, :, l)
  Hkl = H(k, :, l)
  Hklp = H(k, :, lp)
  Hklm = H(k, :, lm)

  ! calculate l upstream
  call H_index(Hklp, Hkl, H_l_up)
  call H_index(Hkl, Hklm, H_l_down)
  H_l_upstream_up = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Sklp(j, i) > Skl(j, i)) then
        H_l_upstream_up(j, i) = Hklp(j, i)
      else
        H_l_upstream_up(j, i) = Hkl(j, i)
      end if
    end do
  end do
  H_l_upstream_down = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skl(j, i) > Sklm(j, i)) then
        H_l_upstream_down(j, i) = Hkl(j, i)
      else
        H_l_upstream_down(j, i) = Hklm(j, i)
      end if
    end do
  end do

  ! calculate l flux
  call flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dx, dy, f_l_plus)
  call flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dx, dy, f_l_min)

  ! calculate l Diffusivity
  D_l_up = gamma * H_l_up**(n+1_i4b) * H_l_upstream_up * f_l_plus
  D_l_dn = gamma * H_l_down**(n+1_i4b) * H_l_upstream_down * f_l_min

  ! calculate k upstream 
  call H_index(Hkpl, Hkl, H_k_up)
  call H_index(Hkl, Hkml, H_k_down)
  H_k_upstream_up = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skpl(j, i) > Skl(j, i)) then
        H_k_upstream_up(j, i) = Hkpl(j, i)
      else
        H_k_upstream_up(j, i) = Hkl(j, i)
      end if
    end do
  end do
  H_k_upstream_down = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skl(j, i) > Skml(j, i)) then
        H_k_upstream_down(j, i) = Hkl(j, i)
      else
        H_k_upstream_down(j, i) = Hkml(j, i)
      end if
    end do
  end do

  ! calculate k flux
  call flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dy, dx, f_k_plus)
  call flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dy, dx, f_k_min)

  ! calculate k Diffusivity
  D_k_up = gamma * H_k_up**(n+1_i4b) * H_k_upstream_up * f_k_plus
  D_k_dn = gamma * H_k_down**(n+1_i4b) * H_k_upstream_down * f_k_min

  ! calculate delta t and t
  divisor = maxval(abs(D_k_up), abs(D_k_dn), abs(D_l_up), abs(D_l_dn))
  if (divisor == 0._rkind) then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dx**2_i4b, dy**2_i4b) / divisor
  end if

  ! Calculate the time step values
  call SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dx, div_l)
  call SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dy, div_k)
  div_q = div_k + div_l

end subroutine diffusion_upstream

! ************************************************************************************************
! MUSCL scheme
! ************************************************************************************************
subroutine diffusion_MUSCL(S, B, gamma, n, cfl, max_dt, Ny, Nx, k, kp, km, kpp, kmm, l, lp, lm, lpp, lmm, dx, dy, div_q, dt_cfl)
  implicit none
  real(rkind), intent(in) :: S(:,:), B(:,:), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: Ny, Nx, k(:), kp(:), km(:), kpp(:), kmm(:), l(:), lp(:), lm(:), lpp(:), lmm(:), n
  real(rkind), intent(out) :: div_q(Ny, Nx), dt_cfl
  real(rkind) :: H(Ny, Nx)
  real(rkind) :: Sklp(Ny), Sklm(Ny), Skplp(Ny), Skplm(Ny), Skpl(Ny), Skl(Ny), Skmlp(Ny), Skmlm(Ny), Skml(Ny)
  real(rkind) :: Hkpl(Ny), Hkppl(Ny), Hkml(Ny), Hkmml(Ny), Hkl(Ny), Hklp(Ny), Hklpp(Ny), Hklm(Ny), Hklmm(Ny)
  real(rkind) :: H_l_min_up(Ny, Nx), H_l_plus_up(Ny, Nx)
  real(rkind) :: H_l_min_down(Ny, Nx), H_l_plus_down(Ny, Nx)
  real(rkind) :: f_l_plus(Ny, Nx), f_l_min(Ny, Nx)
  real(rkind) :: D_l_up_m(Ny, Nx), D_l_up_p(Ny, Nx), D_l_up_min(Ny, Nx), D_l_up_max(Ny, Nx)
  real(rkind) :: D_l_dn_m(Ny, Nx), D_l_dn_p(Ny, Nx), D_l_dn_min(Ny, Nx), D_l_dn_max(Ny, Nx)
  real(rkind) :: D_l_up(Ny, Nx), D_l_dn(Ny, Nx)
  real(rkind) :: H_k_min_up(Ny, Nx), H_k_plus_up(Ny, Nx)
  real(rkind) :: H_k_min_down(Ny, Nx), H_k_plus_down(Ny, Nx)
  real(rkind) :: f_k_plus(Ny, Nx), f_k_min(Ny, Nx)
  real(rkind) :: D_k_up_m(Ny, Nx), D_k_up_p(Ny, Nx), D_k_up_min(Ny, Nx), D_k_up_max(Ny, Nx)
  real(rkind) :: D_k_dn_m(Ny, Nx), D_k_dn_p(Ny, Nx), D_k_dn_min(Ny, Nx), D_k_dn_max(Ny, Nx)
  real(rkind) :: D_k_up(Ny, Nx), D_k_dn(Ny, Nx)
  real(rkind) :: div_l(Ny, Nx), div_k(Ny, Nx)
  real(rkind) :: divisor
  integer(i4b) :: i, j

  H = S - B ! H = ice thickness, S = Surface height, B = bed topography

  ! indices
  Sklp = S(k, :, lp)
  Sklm = S(k, :, lm)
  Skplp = S(kp, :, lp)
  Skplm = S(kp, :, lm)
  Skpl = S(kp, :, l)
  Skl = S(k, :, l)
  Skmlp = S(km, :, lp)
  Skmlm = S(km, :, lm)
  Skml = S(km, :, l)
  !
  Hkpl = H(kp, :, l)
  Hkppl = H(kpp, :, l)
  Hkml = H(km, :, l)
  Hkmml = H(kmm, :, l)
  Hkl = H(k, :, l)
  Hklp = H(k, :, lp)
  Hklpp = H(k, :, lpp)
  Hklm = H(k, :, lm)
  Hklmm = H(k, :, lmm)

  ! calculate l+1/2 index
  call H_min(Hklm, Hkl, Hklp, H_l_min_up)
  call H_plus(Hkl, Hklp, Hklpp, H_l_plus_up)

  ! calculate l-1/2 index
  call H_min(Hklmm, Hklm, Hkl, H_l_min_down)
  call H_plus(Hklm, Hkl, Hklp, H_l_plus_down)

  ! calculate l flux
  call flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dx, dy, f_l_plus)
  call flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dx, dy, f_l_min)

  ! calculate l Diffusivity
  D_l_up_m = gamma * H_l_min_up**(n+2_i4b) * f_l_plus   ! equation 30 Jarosh 2013
  D_l_up_p = gamma * H_l_plus_up**(n+2_i4b) * f_l_plus  ! equation 30 Jarosh 2013
  D_l_up_min = min(D_l_up_m, D_l_up_p)                ! equation 31 Jarosh 2013
  D_l_up_max = max(D_l_up_m, D_l_up_p)                ! equation 32 Jarosh 2013
  !
  D_l_dn_m = gamma * H_l_min_down**(n+2_i4b) * f_l_min  ! equation 30 Jarosh 2013
  D_l_dn_p = gamma * H_l_plus_down**(n+2_i4b) * f_l_min ! equation 30 Jarosh 2013
  D_l_dn_min = min(D_l_dn_m, D_l_dn_p)                ! equation 31 Jarosh 2013
  D_l_dn_max = max(D_l_dn_m, D_l_dn_p)                ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_l_up = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Sklp(j, i) <= Skl(j, i) .and. H_l_min_up(j, i) <= H_l_plus_up(j, i)) then
        D_l_up(j, i) = D_l_up_min(j, i)
      else if (Sklp(j, i) <= Skl(j, i) .and. H_l_min_up(j, i) > H_l_plus_up(j, i)) then
        D_l_up(j, i) = D_l_up_max(j, i)
      else if (Sklp(j, i) > Skl(j, i) .and. H_l_min_up(j, i) <= H_l_plus_up(j, i)) then
        D_l_up(j, i) = D_l_up_max(j, i)
      else if (Sklp(j, i) > Skl(j, i) .and. H_l_min_up(j, i) > H_l_plus_up(j, i)) then
        D_l_up(j, i) = D_l_up_min(j, i)
      end if
    end do
  end do
  D_l_dn = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skl(j, i) <= Sklm(j, i) .and. H_l_min_down(j, i) <= H_l_plus_down(j, i)) then
        D_l_dn(j, i) = D_l_dn_min(j, i)
      else if (Skl(j, i) <= Sklm(j, i) .and. H_l_min_down(j, i) > H_l_plus_down(j, i)) then
        D_l_dn(j, i) = D_l_dn_max(j, i)
      else if (Skl(j, i) > Sklm(j, i) .and. H_l_min_down(j, i) <= H_l_plus_down(j, i)) then
        D_l_dn(j, i) = D_l_dn_max(j, i)
      else if (Skl(j, i) > Sklm(j, i) .and. H_l_min_down(j, i) > H_l_plus_down(j, i)) then
        D_l_dn(j, i) = D_l_dn_min(j, i)
      end if
    end do
  end do

  ! calculate k+1/2 index
  call H_min(Hkml, Hkl, Hkpl, H_k_min_up)
  call H_plus(Hkl, Hkpl, Hkppl, H_k_plus_up)

  ! calculate k-1/2 index
  call H_min(Hkmml, Hkml, Hkl, H_k_min_down)
  call H_plus(Hkml, Hkl, Hkpl, H_k_plus_down)

  ! calculate k flux
  call flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dy, dx, f_k_plus)
  call flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dy, dx, f_k_min)

  ! calculate k Diffusivity
  D_k_up_m = gamma * H_k_min_up**(n+2_i4b) * f_k_plus   ! equation 30 Jarosh 2013
  D_k_up_p = gamma * H_k_plus_up**(n+2_i4b) * f_k_plus  ! equation 30 Jarosh 2013
  D_k_up_min = min(D_k_up_m, D_k_up_p)                ! equation 31 Jarosh 2013
  D_k_up_max = max(D_k_up_m, D_k_up_p)                ! equation 32 Jarosh 2013
  !
  D_k_dn_m = gamma * H_k_min_down**(n+2_i4b) * f_k_min  ! equation 30 Jarosh 2013
  D_k_dn_p = gamma * H_k_plus_down**(n+2_i4b) * f_k_min ! equation 30 Jarosh 2013
  D_k_dn_min = min(D_k_dn_m, D_k_dn_p)                ! equation 31 Jarosh 2013
  D_k_dn_max = max(D_k_dn_m, D_k_dn_p)                ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_k_up = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skpl(j, i) <= Skl(j, i) .and. H_k_min_up(j, i) <= H_k_plus_up(j, i)) then
        D_k_up(j, i) = D_k_up_min(j, i)
      else if (Skpl(j, i) <= Skl(j, i) .and. H_k_min_up(j, i) > H_k_plus_up(j, i)) then
        D_k_up(j, i) = D_k_up_max(j, i)
      else if (Skpl(j, i) > Skl(j, i) .and. H_k_min_up(j, i) <= H_k_plus_up(j, i)) then
        D_k_up(j, i) = D_k_up_max(j, i)
      else if (Skpl(j, i) > Skl(j, i) .and. H_k_min_up(j, i) > H_k_plus_up(j, i)) then
        D_k_up(j, i) = D_k_up_min(j, i)
      end if
    end do
  end do
  D_k_dn = 0._rkind
  do j = 1, Ny
    do i = 1, Nx
      if (Skl(j, i) <= Skml(j, i) .and. H_k_min_down(j, i) <= H_k_plus_down(j, i)) then
        D_k_dn(j, i) = D_k_dn_min(j, i)
      else if (Skl(j, i) <= Skml(j, i) .and. H_k_min_down(j, i) > H_k_plus_down(j, i)) then
        D_k_dn(j, i) = D_k_dn_max(j, i)
      else if (Skl(j, i) > Skml(j, i) .and. H_k_min_down(j, i) <= H_k_plus_down(j, i)) then
        D_k_dn(j, i) = D_k_dn_max(j, i)
      else if (Skl(j, i) > Skml(j, i) .and. H_k_min_down(j, i) > H_k_plus_down(j, i)) then
        D_k_dn(j, i) = D_k_dn_min(j, i)
      end if
    end do
  end do

  ! calculate delta t and t
  divisor = maxval(abs(D_k_up), abs(D_k_dn), abs(D_l_dn), abs(D_l_dn))
  if (divisor == 0._rkind) then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dx**2, dy**2) / divisor
  end if

  ! Calculate the time step values
  call SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dx, div_l) ! equation 36 Jarosh 2013
  call SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dy, div_k) ! equation 36 Jarosh 2013
  div_q = div_k + div_l

end subroutine diffusion_MUSCL


! ************************************************************************************************
! private functions to assist in the glacier flow calculations
! ************************************************************************************************
function minmod(a, b) result(minmod_result)
  implicit none
  real(rkind), intent(in) :: a, b
  real(rkind) :: minmod_result
  real(rkind) :: sign

  sign = sign(a) + sign(b)
  minmod_result = sign / 2._rkind * min(abs(a), abs(b))
end function minmod

function superbee(r) result(superbee_result)
  implicit none
  real(rkind), intent(in) :: r
  real(rkind) :: superbee_result

  superbee_result = max(0._rkind, min(2._rkind * r, 1._rkind), min(r, 2._rkind))
end function superbee

function flux(sjpl, sjml, sjplp, sjmlp, sjlp, sjl, dj, dl, n) result(flux_result)
  implicit none
  real(rkind), intent(in) :: sjpl, sjml, sjplp, sjmlp, sjlp, sjl, dj, dl, n
  real(rkind) :: flux_result

  flux_result = ((sjpl - sjml + sjplp - sjmlp)**2_i4b / (4._rkind * dj)**2_i4b + (sjlp - sjl)**2_i4b / dl**2_i4b)**((n - 1.0) / 2._rkind)
end function flux

function SIA(Dup, Sup, S, Ddn, Sdn, d) result(SIA_result)
  implicit none
  real(rkind), intent(in) :: Dup, Sup, S, Ddn, Sdn, d
  real(rkind) :: SIA_result

  SIA_result = (Dup * (Sup - S) / d - Ddn * (S - Sdn) / d) / d
end function SIA

function H_index(h1, h2) result(H_index_result)
  implicit none
  real(rkind), intent(in) :: h1, h2
  real(rkind) :: H_index_result

  H_index_result = 0.5_rkind * (h1 + h2)
end function H_index

function H_plus(Hm, H, Hp) result(H_plus_result)
  implicit none
  real(rkind), intent(in) :: Hm, H, Hp
  real(rkind) :: H_plus_result

  if (limiter == "minmod") then
    H_plus_result = H - 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if (limiter == "superbee") then
    if (Hp /= H .and. Hp /= Hm .and. H /= Hm) then
      H_plus_result = H - 0.5_rkind * superbee(abs((H - Hm) / (Hp - H))) * (Hp - H)
    else
      H_plus_result = H
    end if
  end if
end function H_plus

function H_min(Hm, H, Hp, limiter) result(H_min_result)
  implicit none
  real(rkind), intent(in) :: Hm, H, Hp
  real(rkind) :: H_min_result

  if (limiter == "minmod") then
    H_min_result = H + 0.5_rkind * minmod(H - Hm, Hp - H) * (Hp - H)
  else if (limiter == "superbee") then
    if (Hp /= H .and. Hp /= Hm .and. H /= Hm) then
      H_min_result = H + 0.5_rkind * superbee(abs((Hp - H) / (H - Hm))) * (H - Hm)
    else
      H_min_result = H
    end if
  end if
end function H_min

end module glacFlow_module