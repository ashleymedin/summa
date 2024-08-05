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

module qTimeDelay_module

! data types
USE nrtype

! constants, time information
USE multiconst,only:secprhour  ! number of seconds in an hour
USE globalData,only:data_step  ! length of the data step (s)

! look-up values for the sub-grid routing method
USE mDecisions_module,only:      &
 timeDelay,&  ! time-delay histogram
 qInstant     ! instantaneous routing

implicit none
private
public::qOverland
public::qGlacier
contains


 ! *************************************************************************************************************
 ! public subroutine qOverland: compute the time delay in runoff in a basin (places runoff in future time steps)
 ! *************************************************************************************************************
 subroutine qOverland(&
                      ! input
                      ixRouting,             &  ! index for routing method
                      averageTotalRunoff,    &  ! total runoff to the channel from all active components (m s-1)
                      fracFuture,            &  ! fraction of runoff in future time steps (m s-1)
                      qFuture,               &  ! runoff in future time steps (m s-1)
                      ! output
                      averageInstantRunoff,  &  ! instantaneous runoff (m s-1)
                      averageRoutedRunoff,   &  ! routed runoff (m s-1)
                      err,message)              ! error control
 implicit none
 ! input
 integer(i4b),intent(in)    :: ixRouting                 ! index for routing method
 real(rkind),intent(in)     :: averageTotalRunoff        ! total runoff to the channel from all active components (m s-1)
 real(rkind),intent(in)     :: fracFuture(:)             ! fraction of runoff in future time steps (m s-1)
 real(rkind),intent(inout)  :: qFuture(:)                ! runoff in future time steps (m s-1)
 ! output
 real(rkind),intent(out)    :: averageInstantRunoff      ! instantaneous runoff (m s-1)
 real(rkind),intent(out)    :: averageRoutedRunoff       ! routed runoff (m s-1)
 integer(i4b),intent(out)   :: err                       ! error code
 character(*),intent(out)   :: message                   ! error message
 ! internal
 real(rkind),parameter      :: valueMissing=-9999._rkind ! missing value
 integer(i4b)               :: nTDH                      ! number of points in the time-delay histogram
 integer(i4b)               :: iFuture                   ! index in time delay histogram
 ! initialize error control
 err=0; message='qOverland/'

 ! assign instantaneous runoff (m s-1)  (Note: this variable is redundant with averageTotalRunoff, could remove)
 averageInstantRunoff = averageTotalRunoff

 ! compute routed runoff (m s-1)
 select case(ixRouting)  ! (select option for sub-grid routing)
  ! ** instantaneous routing
  case(qInstant)
   averageRoutedRunoff = averageInstantRunoff

  ! ** time delay histogram
  case(timeDelay)
   ! identify number of points in the time-delay histogram
   nTDH = size(qFuture)
   ! place a fraction of runoff in future steps
   qFuture(1:nTDH) = qFuture(1:nTDH) + averageInstantRunoff*fracFuture(1:nTDH)
   ! save the routed runoff
   averageRoutedRunoff = qFuture(1)
   ! move array back
   do iFuture=2,nTDH
    qFuture(iFuture-1) = qFuture(iFuture)
   end do
   qFuture(nTDH) = 0._rkind

  ! ** error checking
  case default; err=20; message=trim(message)//'cannot find option for sub-grid routing'; return

 end select ! (select option for sub-grid routing)
 ! For open water SUMMA doesn't run any calculations
 !  the values for any output variables in the netCDF will stay at the value at which they were initialized, which may be a large negative
 ! Coast may be similarly large and negative
 !if (averageRoutedRunoff < 0._rkind) averageRoutedRunoff = valueMissing

 end subroutine qOverland


 ! *************************************************************************************************************
 ! public subroutine qGlacier: compute the time delay in glacier melt to route to stream in a basin (places runoff in future time steps)
 !      Note: might want to divide ablation reservoir into the standard two parts: one for snow melt and one for ice melt
 !            and then route each separately.  Also should make k_abl and k_acc parameters in the input file.
 ! *************************************************************************************************************
 subroutine qGlacier(&
                     ! input
                     glacAblMelt,           &  ! total melt into ablation reservoir (m s-1)
                     glacAccMelt,           &  ! total melt into accumulation reservoir (m s-1)
                     glacAblArea,           &  ! per glacier acumulation area (m2)
                     glacAccArea,           &  ! per glacier ablation area (m2)
                     qAblFuture,            &  ! per glacier ablation reservoir runoff in future time steps (m s-1)
                     qAccFuture,            &  ! per glacier accumlation reservoir runoff in future time steps (m s-1)
                     ! output
                     glacierRoutedRunoff,   &  ! routed glacier runoff (m s-1)
                     err,message)              ! error control
 implicit none
 ! input
 real(rkind),intent(in)     :: data_step                 ! time step for the data (s)
 real(rkind),intent(in)     :: glacAblMelt               ! total melt into ablation reservoir (m s-1)
 real(rkind),intent(in)     :: glacAccMelt               ! total melt into accumulation reservoir (m s-1)
 real(rkind),intent(in)     :: glacAblArea(:)            ! per glacier acumulation area (m2)
 real(rkind),intent(in)     :: glacAccArea(:)            ! per glacier ablation area (m2)
 real(rkind),intent(inout)  :: qAblFuture(:)             ! per glacier ablation reservoir runoff in future time steps (m s-1)
 real(rkind),intent(inout)  :: qAccFuture(:)             ! per glacier accumlation reservoir runoff in future time steps (m s-1)
 ! output
 real(rkind),intent(out)    :: glacierRoutedRunoff       ! routed glacier runoff (m s-1)
 integer(i4b),intent(out)   :: err                       ! error code
 character(*),intent(out)   :: message                   ! error message
 ! internal
 real(rkind)                :: qAbl                      ! hourly ablation reservoir runoff (m s-1)
 real(rkind)                :: qAcc                      ! hourly accumulation reservoir runoff (m s-1)
 real(rkind)                :: frac                      ! fraction of glacier area
 real(rkind)                :: glacAblTotal              ! total ablation area (m2)
 real(rkind)                :: glacAccTotal              ! total accumulation area (m2)
 integer(i4b)               :: nGlacier                  ! number of glaciers in the basin
 integer(i4b)               :: iGlacier                  ! index for glaciers
 real(rkind),parameter      :: k_abl=10                  ! storage coefficient ablation reservoir (hours)
 real(rkind),parameter      :: k_acc=400                 ! storage coefficient accumulation reservoir (hours)
 ! initialize error control
 err=0; message='qGlacier/' 

 nGlacier = size(qAblFuture) ! number of glaciers in the basin
 glacierRoutedRunoff = 0._rkind

 glacAblTotal = sum(glacAblArea)
 glacAccTotal = sum(glacAccArea)

 do iGlacier=1,nGlacier
   ! ablation reservoir runoff (m s-1)
   frac = glacAblArea(iGlacier)/glacAblTotal
   qAbl = qAblFuture(iGlacier) + glacAblMelt*frac - glacAblMelt*frac*exp(-1._rkind/k_abl)
   qAblFuture(iGlacier) = qAbl*exp(-data_step/secprhour/k_abl) ! place runoff in future time steps 

   ! accumulation reservoir runoff (m s-1)
   frac = glacAccArea(iGlacier)/glacAccTotal
   qAcc = qAccFuture(iGlacier) + glacAccMelt(iGlacier)*frac - glacAccMelt(iGlacier)*frac*exp(-1._rkind/k_acc)
   qAccFuture(iGlacier) = qAcc*exp(-data_step/secprhour/k_acc) ! place runoff in future time steps 

   ! routed glacier runoff (m s-1)
   glacierRoutedRunoff = glacierRoutedRunoff + qAbl + qAcc 
 end do

 end subroutine qGlacier


end module qTimeDelay_module
