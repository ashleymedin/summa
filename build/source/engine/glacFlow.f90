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
subroutine glacFlow()
  implicit none
  ! model control
  integer(i8b), intent(in)    :: gruId(:)             ! gruId
  integer(i4b), intent(in)    :: nDOM                 ! number of domains
  integer(i8b), intent(in)    :: hruId(:)             ! hruId
  ! mass balance
  real(rkind), intent(in)     :: GWE_deltaYr(:)       ! change in glacier water equivalent (m s-1) in each glacier domain
  real(rkind), intent(inout)  :: elev(:)              ! elevation of each glacier domain (m)
  ! area 
  real(rkind), intent(inout)  :: glacAblArea(:)       ! per glacier ablation area (m2)
  real(rkind), intent(inout)  :: glacAccArea(:)       ! per glacier accumulation area (m2)
  real(rkind), intent(out)    :: area(nDOM)           ! area of each glacier domain (m2)
  ! locals
  integer(i4b)                :: i,j
  integer(i4b)                :: nGlacier
  integer(i4b),parameter      :: Ny==1, Nx==100, nYr = 1
  real(rkind)                 :: volume ! volume of glacier km3
  real                        :: surface(Ny,Nx),base(Ny,Nx), mb(Ny,Nx)
  
  ! get starting surface and base from file and number of glaciers
  ! This logic makes sense if assuming multiple glaciers in each HRU and one HRU per GRU, or one glacier in each GRU with multiple HRUs
  ! If some glaciers are not in a particular glacier HRU, this logic will not capture that

  nGlacier = size(glacAblArea)
  ! get by gruID()
  area_all(1:nGlacier) = 
  surface_all
  base_all
  hruCount_all
  ! starting values can get from shapefile
  ! maybe call this with a flag to initialize progressive? Or read into coldState_glac
  abl_elev = (zmin_m + zmed_m)/2._rkind
  acc_elev = (zmax_m + zmed_m)/2._rkind
  abl_area = area_km2/2._rkind * 1.e6_rkind
  acc_area = abl_area


  ! run flow for each glacier
  do i = 1,nGlacier
    hruCount = 0
    if(glacAblArea(i)+glacAccArea(i)==0._rkind) cycle ! not growing glaciers

    do k = 1,nGlacier
      if (area_all(k).ne. glacAblArea(i)+glacAccArea(i))cycle
      surface = surface_all(:,:,k)
      base = base_all(:,:,k)
      hruCount = hruCount_all(k) ! number of HRUs in glacier
      exit
    enddo
    if(hruCount==0)then
      message = trim(message)//"cannot find glacier topography info in file"
      err=20; return
    endif
    
    allocate(hruGlac(hruCount),mb0(hruCount*2),elev0(hruCount*2))
    hruGlac =  ! get these
    do j = 1,nDOM
      if (any(hruGlac)==hruId(j))then
        if(elev(j).ne.realMissing)then
          !use to extrapolate

        endif
      endif
    enddo

    ! distribute mass balance over surface
    ! Build mass balance gradient from elev, mb_in, 
    ! If use elevation of domains outside of HRU in glacier, then storage change won't map correctly
    ! Map to surface, extrapolate linearly
    ! Order will go HRU 1: Acc, Abl, HRU 2: Acc, Abl unless have zero area then skip
    mb =
    m_dot = ! convert to m/s
    do j = 1, Ny
      do k = 1, Nx
        if (.not. glacierMask(j, k)) then
          m_dot(j, k) = 0._rkind
        end if
      end do
    end do

    call run_year(nYr,surface, base, mb, glacierMask, x, Ny, Nx, dx, dy, volume)

    ! based on mass balance, use width to get accumulation area and ablation area and elevations


    ! HOW MAP THIS BACK TO HRUS?


  enddo ! end of glacier loop

end subroutine glacFlow


! ************************************************************************************************
! private subroutine run_year for setting up flow model and running for each glacier for time period
! ************************************************************************************************
  subroutine run_year(y_end, S, B, m_dot, glacierMask, x, Ny, Nx, dx, dy, volume)
    ! Arguments
    real(rkind), intent(in) :: y_end, dx, dy
    real(rkind), intent(inout) :: S(:,:), B(:,:), m_dot(:,:)
    integer(i4b), intent(in) :: Nx, Ny
    logical(lgt), intent(in) :: glacierMask(:)
    real(rkind), intent(in) :: x(:)
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