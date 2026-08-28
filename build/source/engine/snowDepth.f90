module snowDepth_module

! data types
USE nr_type
USE data_types,only:&
                    var_ilength,   & ! x%var(:)%dat            (i4b)
                    var_dlength,   & ! x%var(:)%dat            (rkind)
                    zLookup          ! x%z(:)%var(:)%lookup(:) (rkind)

! constants
USE multiconst,only:&
                    Tfreeze,      & ! freezing temperature (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_air,     & ! intrinsic density of air (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
USE globalData,only:verySmall       ! a small number

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG        ! named variables for structure elements
USE var_lookup,only:iLookDIAG        ! named variables for structure elements
USE var_lookup,only:iLookFLUX        ! named variables for structure elements
USE var_lookup,only:iLookPARAM       ! named variables for structure elements
USE var_lookup,only:iLookINDEX       ! named variables for structure elements

! privacy
implicit none
private
public::snowDepth

contains

! ************************************************************************************************
! public subroutine snowDepth: compute snow depth for one sub timestep
! ************************************************************************************************
subroutine snowDepth(&
                        dt_sub,                 & ! intent(in):    time step (s)
                        nSnow,                  & ! intent(in):    number of snow layers
                        scalarSnowSublimation,  & ! intent(in):    scalar sublimation of snow (kg m-2)
                        mLayerVolFracLiq,       & ! intent(inout): volumetric fraction of liquid water in each layer (-)
                        mLayerVolFracIce,       & ! intent(inout): volumetric fraction of ice in each layer (-)
                        mLayerTemp,             & ! intent(in):    temperature of each layer (K)
                        mLayerMeltFreeze,       & ! intent(in):    volumetric melt in each layer (kg m-3)
                        mpar_data,              & ! intent(in):    model parameters
                        ! output
                        tooMuchSublim,          & ! intent(out):   flag to denote that there was too much sublimation in a given time step
                        mLayerDepth,            & ! intent(inout): depth of each layer (m)
                        ! error control
                        err,message)              ! intent(out):   error control

  ! -----------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(qp),intent(in)                  :: dt_sub                 ! time step (s)
  integer(i4b),intent(in)              :: nSnow                  ! number of snow layers
  real(rkind),intent(in)               :: scalarSnowSublimation  ! scalar sublimation of snow (kg m-2)
  real(rkind),intent(inout)            :: mLayerVolFracLiq(:)    ! volumetric fraction of liquid water in each layer (-)
  real(rkind),intent(inout)            :: mLayerVolFracIce(:)    ! volumetric fraction of ice in each layer (-)
  real(rkind),intent(in)               :: mLayerTemp(:)          ! temperature of each layer (K)
  real(rkind),intent(in)               :: mLayerMeltFreeze(:)    ! volumetric melt in each layer (kg m-3)
  type(var_dlength),intent(in)         :: mpar_data              ! model parameters
  logical(lgt)                         :: tooMuchSublim          ! flag to denote that there was too much sublimation in a given time step
  real(rkind),intent(inout)            :: mLayerDepth(:)         ! depth of each layer (m)
  integer(i4b),intent(out)             :: err                    ! error code
  character(*),intent(out)             :: message                ! error message
  ! local variables
  character(len=256)                   :: cmessage               ! error message
  integer(i4b)                         :: iSnow                  ! index of snow layers
  real(rkind)                          :: massLiquid             ! mass liquid water (kg m-2)

  ! * compute change in ice content of the top snow layer due to sublimation...
  ! ---------------------------------------------------------------------------
  ! initialize the flags
  tooMuchSublim=.false.  ! too much sublimation (merge snow layers)
  ! NOTE: this is done BEFORE densification
  if(nSnow>0)then ! snow layers exist

    ! try to remove ice from the top layer
    iSnow=1

    ! save the mass of liquid water (kg m-2)
    massLiquid = mLayerDepth(iSnow)*mLayerVolFracLiq(iSnow)*iden_water

    ! add/remove the depth of snow gained/lost by frost/sublimation (m)
    ! NOTE: assume constant density
    mLayerDepth(iSnow) = mLayerDepth(iSnow) + dt_sub*scalarSnowSublimation/(mLayerVolFracIce(iSnow)*iden_ice)

    ! check that we did not remove the entire layer
    if(mLayerDepth(iSnow) < verySmall)then
      tooMuchSublim=.true.
      return
    endif

    ! update the volumetric fraction of liquid water
    mLayerVolFracLiq(iSnow) = massLiquid / (mLayerDepth(iSnow)*iden_water)

  ! no snow
  else

    ! no snow: check that sublimation is zero
    if(abs(scalarSnowSublimation) > verySmall)then
      message=trim(message)//'sublimation of snow has been computed when no snow exists'
      err=20; return
    end if

  end if  ! (if snow layers exist)


  ! *** account for compaction and cavitation in the snowpack...
  ! ------------------------------------------------------------
  if(nSnow>0)then
    call snowDensify(&
                    ! intent(in): variables
                    dt_sub,                                            & ! intent(in): time step (s)
                    nSnow,                                             & ! intent(in): number of snow layers
                    mLayerTemp(1:nSnow),                               & ! intent(in): temperature of each layer (K)
                    mLayerMeltFreeze(1:nSnow),                         & ! intent(in): volumetric melt in each layer (kg m-3)
                    ! intent(in): parameters
                    mpar_data%var(iLookPARAM%densScalGrowth)%dat(1),   & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1),   & ! intent(in): temperature scaling factor for grain growth (K-1)
                    mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1),  & ! intent(in): rate of grain growth (s-1)
                    mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1),   & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1),   & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                    mpar_data%var(iLookPARAM%baseViscosity)%dat(1),    & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                    ! intent(inout): state variables
                    mLayerDepth(1:nSnow),                              & ! intent(inout): depth of each layer (m)
                    mLayerVolFracLiq(1:nSnow),                         & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                    mLayerVolFracIce(1:nSnow),                         & ! intent(inout):  volumetric fraction of ice after itertations (-)
                    ! output: error control
                    err,cmessage)                     ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
  end if  ! if snow layers exist

end subroutine snowDepth


! ************************************************************************************************
 ! public subroutine snowDensify: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine snowDensify(&
                       ! intent(in): variables
                       dt,                             & ! intent(in): time step (s)
                       nSnow,                          & ! intent(in): number of snow layers
                       mLayerTemp,                     & ! intent(in): temperature of each layer (K)
                       mLayerMeltFreeze,               & ! intent(in): volumnetric melt in each layer (kg m-3)
                       ! intent(in): parameters
                       densScalGrowth,                 & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                       tempScalGrowth,                 & ! intent(in): temperature scaling factor for grain growth (K-1)
                       grainGrowthRate,                & ! intent(in): rate of grain growth (s-1)
                       densScalOvrbdn,                 & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                       tempScalOvrbdn,                 & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                       baseViscosity,                  & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                       ! intent(inout): state variables
                       mLayerDepth,                    & ! intent(inout): depth of each layer (m)
                       mLayerVolFracLiqNew,            & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceNew,            & ! intent(inout):  volumetric fraction of ice after itertations (-)
                       ! output: error control
                       err,message)                      ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute change in snow density over the time step
 implicit none
 ! intent(in): variables
 real(rkind),intent(in)              :: dt                       ! time step (seconds)
 integer(i4b),intent(in)             :: nSnow                    ! number of snow layers
 real(rkind),intent(in)              :: mLayerTemp(:)            ! temperature of each snow layer after iterations (K)
 real(rkind),intent(in)              :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
 ! intent(in): parameters
 real(rkind),intent(in)              :: densScalGrowth           ! density scaling factor for grain growth (kg-1 m3)
 real(rkind),intent(in)              :: tempScalGrowth           ! temperature scaling factor for grain growth (K-1)
 real(rkind),intent(in)              :: grainGrowthRate          ! rate of grain growth (s-1)
 real(rkind),intent(in)              :: densScalOvrbdn           ! density scaling factor for overburden pressure (kg-1 m3)
 real(rkind),intent(in)              :: tempScalOvrbdn           ! temperature scaling factor for overburden pressure (K-1)
 real(rkind),intent(in)              :: baseViscosity            ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
 ! intent(inout): state variables
 real(rkind),intent(inout)           :: mLayerDepth(:)           ! depth of each layer (m)
 real(rkind),intent(inout)           :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water in each snow layer after iterations (-)
 real(rkind),intent(inout)           :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice in each snow layer after iterations (-)
 ! intent(out): error control
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 integer(i4b)                        :: iSnow                       ! index of snow layers
 real(rkind)                         :: chi1,chi2,chi3,chi4,chi5    ! multipliers in the densification algorithm (-)
 real(rkind)                         :: halfWeight                  ! half of the weight of the current snow layer (kg m-2)
 real(rkind)                         :: weightSnow                  ! total weight of snow above the current snow layer (kg m-2)
 real(rkind)                         :: CR_grainGrowth              ! compaction rate for grain growth (s-1)
 real(rkind)                         :: CR_ovrbdnPress              ! compaction rate associated with over-burden pressure (s-1)
 real(rkind)                         :: CR_metamorph                ! compaction rate for metamorphism (s-1)
 real(rkind)                         :: massIceOld                  ! mass of ice in the snow layer (kg m-2)
 real(rkind)                         :: massLiqOld                  ! mass of liquid water in the snow layer (kg m-2)
 real(rkind)                         :: scalarDepthNew              ! updated layer depth (m)
 real(rkind)                         :: scalarDepthMin              ! minimum layer depth (m)
 real(rkind)                         :: volFracIceLoss              ! volumetric fraction of ice lost due to melt (-)
 real(rkind), dimension(nSnow)       :: mLayerVolFracAirNew         ! volumetric fraction of air in each layer after compaction (-)
 real(rkind),parameter               :: snwden_min=100._rkind       ! minimum snow density for reducing metamorphism rate (kg m-3)
 real(rkind),parameter               :: snwDensityMax=550._rkind    ! maximum snow density for collapse under melt (kg m-3)
 real(rkind),parameter               :: wetSnowThresh=0.01_rkind    ! threshold to discriminate between "wet" and "dry" snow
 real(rkind),parameter               :: minLayerDensity=40._rkind   ! minimum snow density allowed for any layer (kg m-3)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="snowDensify/"

 ! NOTE: still need to process the case of "snow without a layer"
 if(nSnow==0)return

 ! initialize the weight of snow above each layer (kg m-2)
 weightSnow = 0._rkind

 ! loop through snow layers
 do iSnow=1,nSnow

  ! save mass of liquid water and ice (mass does not change)
  massIceOld = iden_ice*mLayerVolFracIceNew(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
  massLiqOld = iden_water*mLayerVolFracLiqNew(iSnow)*mLayerDepth(iSnow) ! (kg m-2)

  ! *** compute the compaction associated with grain growth (s-1)
  ! compute the base rate of grain growth (-)
  if(mLayerVolFracIceNew(iSnow)*iden_ice <snwden_min) chi1=1._rkind
  if(mLayerVolFracIceNew(iSnow)*iden_ice>=snwden_min) chi1=exp(-densScalGrowth*(mLayerVolFracIceNew(iSnow)*iden_ice - snwden_min))
  ! compute the reduction of grain growth under colder snow temperatures (-)
  chi2 = exp(-tempScalGrowth*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the acceleration of grain growth in the presence of liquid water (-)
  if(mLayerVolFracLiqNew(iSnow) > wetSnowThresh)then; chi3=2._rkind  ! snow is "wet"
  else; chi3=1._rkind; end if                                         ! snow is "dry"
  ! compute the compaction associated with grain growth (s-1)
  CR_grainGrowth = grainGrowthRate*chi1*chi2*chi3

  ! **** compute the compaction associated with over-burden pressure (s-1)
  halfWeight = (massIceOld + massLiqOld)/2._rkind  ! half of the weight from the current layer
  weightSnow = weightSnow + halfweight             ! there is some over-burden pressure from the layer itself, add half
  chi4 = exp(-tempScalOvrbdn*(Tfreeze - mLayerTemp(iSnow)))       ! increase in compaction under colder snow temperatures (-)
  chi5 = exp(-densScalOvrbdn*mLayerVolFracIceNew(iSnow)*iden_ice) ! increase in compaction under low density snow (-)
  CR_ovrbdnPress = (weightSnow/baseViscosity)*chi4*chi5 ! compaction associated with over-burden pressure (s-1)
  weightSnow = weightSnow + halfweight                  ! update the snow weight, add half of the weight from the current layer

  ! *** compute the compaction rate associated with snowmelt (s-1)
  ! NOTE: loss of ice due to snowmelt is implicit, so can be updated directly
  if(iden_ice*mLayerVolFracIceNew(iSnow) < snwDensityMax)then ! only collapse layers if below a critical density
    volFracIceLoss = max(0._rkind,mLayerMeltFreeze(iSnow)/iden_ice)  ! volumetric fraction of ice lost due to melt (-)
    scalarDepthNew = mLayerDepth(iSnow) * mLayerVolFracIceNew(iSnow)/(mLayerVolFracIceNew(iSnow) + volFracIceLoss) ! adjust snow depth to account for cavitation
  else ! do not allow collapse of snow layers
    scalarDepthNew = mLayerDepth(iSnow)
  end if

  ! **** compute the total compaction rate associated with metamorphism
  CR_metamorph = CR_grainGrowth + CR_ovrbdnPress
  ! update depth due to metamorphism (implicit solution)
  scalarDepthNew = scalarDepthNew/(1._rkind + CR_metamorph*dt)
  ! ensure that the new depth is in line with the maximum amount of compaction that can occur given the masses of ice and liquid in the layer
  scalarDepthMin = (massIceOld / iden_ice) + (massLiqOld / iden_water)
  mLayerDepth(iSnow) = max(scalarDepthMin, scalarDepthNew)

  ! update volumetric ice and liquid water content
  mLayerVolFracIceNew(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
  mLayerVolFracLiqNew(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)
  mLayerVolFracAirNew(iSnow) = 1.0_rkind - mLayerVolFracIceNew(iSnow) - mLayerVolFracLiqNew(iSnow)

 end do  ! looping through snow layers

 ! check depth
 if(any(mLayerDepth(1:nSnow) < 0._rkind))then
  do iSnow=1,nSnow
   write(*,'(a,1x,i4,1x,4(f12.5,1x))') 'iSnow, mLayerDepth(iSnow)', iSnow, mLayerDepth(iSnow)
  end do
  message=trim(message)//'unreasonable value for snow depth'
  err=20; return
 end if

 ! check for low/high snow density
 if(any(mLayerVolFracIceNew(1:nSnow)*iden_ice + mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracAirNew(1:nSnow)*iden_air < minLayerDensity) .or. &
    any(mLayerVolFracIceNew(1:nSnow) + mLayerVolFracLiqNew(1:nSnow) + mLayerVolFracAirNew(1:nSnow) > 1._rkind))then
  do iSnow=1,nSnow
   write(*,*) 'iSnow, volFracIce, density = ', iSnow, mLayerVolFracIceNew(iSnow),  mLayerVolFracIceNew(iSnow)*iden_ice
  end do
  message=trim(message)//'unreasonable value for snow density'
  err=20; return
 end if

 end subroutine snowDensify

end module snowDepth_module

