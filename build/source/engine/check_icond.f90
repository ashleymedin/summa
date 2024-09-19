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

module check_icond_module
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! access domain types
USE globalData,only:upland          ! domain type for upland areas
USE globalData,only:glacAcc         ! domain type for glacier accumulation areas
USE globalData,only:glacAbl         ! domain type for glacier ablation areas
USE globalData,only:wetland         ! domain type for wetland areas

implicit none
private
public::check_icond
contains

 ! ************************************************************************************************
 ! public subroutine check_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine check_icond(nGRU,                          & ! number of GRUs
                        bvarData,                      & ! model basin variables
                        progData,                      & ! model prognostic (state) variables
                        diagData,                      & ! model diagnostic variables
                        mparData,                      & ! model parameters
                        indxData,                      & ! layer index data
                        lookupData,                    & ! lookup table data
                        attrData,                      & ! model attributes
                        checkEnthalpy,                 & ! flag if to check enthalpy for consistency
                        use_lookup,                    & ! flag to use the lookup table for soil enthalpy                             
                        err,message)                     ! error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookBVAR                           ! variable lookup structure
 USE var_lookup,only:iLookPARAM                          ! variable lookup structure
 USE var_lookup,only:iLookPROG                           ! variable lookup structure
 USE var_lookup,only:iLookDIAG                           ! variable lookup structure
 USE var_lookup,only:iLookINDEX                          ! variable lookup structure
 USE var_lookup,only:iLookATTR                           ! variable lookup structure
 USE globalData,only:gru_struc                           ! gru-hru mapping structures
 USE data_types,only:gru_doubleVec                       ! gru double precision structure no hru no domain
 USE data_types,only:gru_hru_doubleVec                   ! hru double precision structure no domain
 USE data_types,only:gru_hru_dom_doubleVec               ! full double precision structure
 USE data_types,only:gru_hru_dom_intVec                  ! full integer structure
 USE data_types,only:gru_hru_dom_z_vLookup               ! full lookup structure
 USE data_types,only:gru_hru_double                      ! hru double precision structure no domain no depth
 USE globalData,only:iname_soil,iname_snow,iname_ice,iname_lake ! named variables to describe the type of layer
 USE multiconst,only:&
                       LH_fus,    &                      ! latent heat of fusion                (J kg-1)
                       iden_ice,  &                      ! i,trinsic density of ice             (kg m-3)
                       iden_water,&                      ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &                      ! gravitational acceleration           (m s-2)
                       Tfreeze                           ! freezing point of pure water         (K)
 USE snow_utils_module,only:fracliquid                   ! compute volumetric fraction of liquid water in snow based on temperature
 USE updatState_module,only:updateSnLaIc                   ! update snow states
 USE updatState_module,only:updateSoil                   ! update soil states
 USE enthalpyTemp_module,only:T2enthTemp_cas             ! convert temperature to enthalpy for canopy air space
 USE enthalpyTemp_module,only:T2enthTemp_veg             ! convert temperature to enthalpy for vegetation
 USE enthalpyTemp_module,only:T2enthTemp_snLaIc          ! convert temperature to enthalpy for snow, lake, and ice
 USE enthalpyTemp_module,only:T2enthTemp_soil            ! convert temperature to enthalpy for soil
 
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 integer(i4b),intent(in)                   :: nGRU           ! number of grouped response units
 type(gru_doubleVec),intent(inout)         :: bvarData       ! basin variables
 type(gru_hru_dom_doubleVec),intent(inout) :: progData       ! prognostic vars
 type(gru_hru_dom_doubleVec),intent(inout) :: diagData       ! diagnostic vars
 type(gru_hru_dom_doubleVec),intent(in)    :: mparData       ! parameters
 type(gru_hru_dom_intVec),intent(in)       :: indxData       ! layer indexes
 type(gru_hru_dom_z_vLookup),intent(in)    :: lookupData     ! lookup table data
 type(gru_hru_double),intent(in)           :: attrData       ! attributes
 logical(lgt),intent(in)                   :: checkEnthalpy  ! if true either need enthTemp as starting residual value, or for state variable initialization
 logical(lgt),intent(in)                   :: use_lookup     ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
 integer(i4b),intent(out)                  :: err            ! error code
 character(*),intent(out)                  :: message        ! returned error message
 ! locals
 character(len=256)             :: cmessage              ! downstream error message
 integer(i4b)                   :: i,iGRU,iHRU,iDOM      ! loop index
 ! temporary variables for realism checks
 integer(i4b)                   :: iLayer                ! index of model layer
 integer(i4b)                   :: iSoil                 ! index of soil layer
 real(rkind)                    :: fLiq                  ! fraction of liquid water on the vegetation canopy (-)
 real(rkind)                    :: vGn_m                 ! van Genutchen "m" parameter (-)
 real(rkind)                    :: tWat                  ! total water on the vegetation canopy (kg m-2)
 real(rkind)                    :: scalarTheta           ! liquid water equivalent of total water [liquid water + ice] (-)
 real(rkind)                    :: h1,h2                 ! used to check depth and height are consistent
 real(rkind)                    :: kappa                 ! constant in the freezing curve function (m K-1)
 integer(i4b)                   :: nSnow                 ! number of snow layers
 integer(i4b)                   :: nLake                 ! number of lake layers
 integer(i4b)                   :: nSoil                 ! number of soil layers
 integer(i4b)                   :: nIce                  ! number of ice layers
 integer(i4b)                   :: nLayers               ! total number of layers
 integer(i4b)                   :: nState                ! total number of states
 real(rkind),parameter          :: xTol=1.e-10_rkind     ! small tolerance to address precision issues
 real(rkind),parameter          :: canIceTol=1.e-3_rkind ! small tolerance to allow existence of canopy ice for above-freezing temperatures (kg m-2)
 real(rkind)                    :: remaining_area        ! remaining area of the HRU
 real(rkind)                    :: remaining_elev        ! remaining elevation of the HRU
 real(rkind)                    :: glacAblAreaTot       ! total basin glacier ablation area from bvarData (m2)
 real(rkind)                    :: glacAccAreaTot       ! total basin glacier accumulation area from bvarData (m2)
 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="check_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! Check that the initial conditions do not conflict with parameters, structure, etc.
 ! --------------------------------------------------------------------------------------------------------

 ! check and correct domain area and elevation, and ensure that the area is positive, and make backwards compatible
 do iGRU = 1,nGRU
   glacAblAreaTot = 0.0_rkind
   glacAccAreaTot = 0.0_rkind
   do iHRU=1,gru_struc(iGRU)%hruCount
     ! update the HRU area and elevation
     remaining_area = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
     remaining_elev = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)*attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%elevation)
     do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
       if (gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type.ne.upland) then
         remaining_area = remaining_area - progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
         remaining_elev = remaining_elev - progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) * progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1)
         if (gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type==glacAbl) then
           glacAblAreaTot = glacAblAreaTot + progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)
         else if (gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type==glacAcc) then
           glacAccAreaTot = glacAccAreaTot + progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) 
         end if
       end if
     end do
     do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
       if (gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type==upland) then
         progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) = remaining_area
         if(remaining_area>0.0_rkind) then 
           progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1) = remaining_elev/remaining_area
         else
           if (remaining_area<-xTol) write(*,'(A,E22.16,A)') 'Warning: area of upland HRU (=', remaining_area, ') < 0. Resetting to 0.0'
           progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1) = 0.0_rkind
           progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) = attrData%gru(iGRU)%hru(iHRU)%var(iLookATTR%elevation)
         end if
       end if
     end do
   end do
   if (abs(glacAblAreaTot-sum(bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat))>xTol) then
     write(*,'(A,E22.16,A,E22.16,A)') 'Warning: glacier ablation domain area (=', glacAblAreaTot, ') does not match the sum of the basin areas (=', sum(bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat), '). Resetting basin areas to be equal.'
     bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat(:) = glacAblAreaTot/size(bvarData%gru(iGRU)%var(iLookBVAR%glacAblArea)%dat)
   end if
   if (abs(glacAccAreaTot-sum(bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat))>xTol) then
     write(*,'(A,E22.16,A,E22.16,A)') 'Warning: glacier accumulation domain area (=', glacAccAreaTot, ') does not match the sum of the basin areas (=', sum(bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat), '). Resetting basin areas to be equal.'
     bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat(:) = glacAccAreaTot/size(bvarData%gru(iGRU)%var(iLookBVAR%glacAccArea)%dat)
   end if

 enddo

 ! check for realistic values of albedo
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
    ! ensure the spectral average albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinWinter)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinWinter)%dat(1)
    ! ensure the visible albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxVisible)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxVisible)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinVisible)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinVisible)%dat(1)
    ! ensure the nearIR albedo is realistic
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) > mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxNearIR)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMaxNearIR)%dat(1)
    if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) < mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinNearIR)%dat(1)) &
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMinNearIR)%dat(1)
   end do
  end do
 end do

 ! ensure the initial conditions are consistent with the constitutive functions
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1,gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! associate local variables with variables in the data structures
    associate(&
    ! state variables in the canopy air space
    scalarCanairTemp     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanairTemp)%dat(1)     ,& ! canopy air temperature (K)
    scalarCanairEnthalpy => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarCanairEnthalpy)%dat(1) ,& ! canopy air enthalpy (J m-3)
    ! state variables in the vegetation canopy
    scalarCanopyTemp     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyTemp)%dat(1)     ,& ! canopy temperature (K)
    scalarCanopyEnthTemp => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) ,& ! canopy temperature component of enthalpy (J m-3)
    scalarCanopyEnthalpy => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1) ,& ! canopy enthalpy (J m-3)
    scalarCanopyLiq      => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyLiq)%dat(1)      ,& ! mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyIce      => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyIce)%dat(1)      ,& ! mass of ice on the vegetation canopy (kg m-2)
    heightCanopyTop      => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%heightCanopyTop)%dat(1)     ,& ! height of the top of the canopy layer (m)
    heightCanopyBottom   => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%heightCanopyBottom)%dat(1)  ,& ! height of the bottom of the canopy layer (m)
    specificHeatVeg      => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%specificHeatVeg)%dat(1)     ,& ! specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation    => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%maxMassVegetation)%dat(1)   ,& ! maximum mass of vegetation (kg m-2)
    ! state variables in the layer domains
    mLayerTemp           => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerTemp)%dat              ,& ! temperature (K)
    mLayerEnthTemp       => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%mLayerEnthTemp)%dat          ,& ! temperature component of enthalpy (J m-3)
    mLayerEnthalpy       => diagData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookDIAG%mLayerEnthalpy)%dat          ,& ! enthalpy (J m-3)
    mLayerVolFracLiq     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq)%dat        ,& ! volumetric fraction of liquid water in each snow layer (-)
    mLayerVolFracIce     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce)%dat        ,& ! volumetric fraction of ice in each snow layer (-)
    mLayerMatricHead     => progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerMatricHead)%dat        ,& ! matric head (m)
    mLayerLayerType      => indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat              ,& ! type of layer (ix_soil or ix_snow)
    ! depth varying soil properties
    soil_dens_intr       => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%soil_dens_intr)%dat         ,& ! intrinsic soil density             (kg m-3)
    vGn_alpha            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_alpha)%dat              ,& ! van Genutchen "alpha" parameter (m-1)
    vGn_n                => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n)%dat                  ,& ! van Genutchen "n" parameter (-)
    theta_sat            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_sat)%dat              ,& ! soil porosity (-)
    theta_res            => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_res)%dat              ,& ! soil residual volumetric water content (-)
    ! snow parameters
    snowfrz_scale        => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%snowfrz_scale)%dat(1)       ,& ! scaling parameter for the snow freezing curve (K-1)
    FCapil               => mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%FCapil)%dat(1)               & ! fraction of pore space in tension storage (-)
    )  ! (associate local variables with model parameters)

    ! compute the constant in the freezing curve function (m K-1)
    kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

    ! check canopy ice content for unrealistic situations
    if(scalarCanopyIce > canIceTol .and. scalarCanopyTemp > Tfreeze)then
     ! ice content > threshold, terminate run
     write(message,'(A,E22.16,A,E9.3,A,F7.3,A,F7.3,A)') trim(message)//'canopy ice (=',scalarCanopyIce,') > canIceTol (=',canIceTol,') when canopy temperature (=',scalarCanopyTemp,') > Tfreeze (=',Tfreeze,')'
     err=20; return
    else if(scalarCanopyIce > 0._rkind .and. scalarCanopyTemp > Tfreeze)then
     ! if here, ice content < threshold. Could be sublimation on previous timestep or simply wrong input. Print a warning
	  write(*,'(A,E22.16,A,F7.3,A,F7.3,A)') 'Warning: canopy ice content in restart file (=',scalarCanopyIce,') > 0 when canopy temperature (=',scalarCanopyTemp,') > Tfreeze (=',Tfreeze,'). Continuing.',NEW_LINE('a')
    end if
    scalarTheta = scalarCanopyIce + scalarCanopyLiq

   if(checkEnthalpy)then ! enthalpy as state variable (cold start often only has temperature)
      call T2enthTemp_cas(&
                  scalarCanairTemp,       & ! intent(in): canopy air temperature (K)
                  scalarCanairEnthalpy)     ! intent(out): enthalpy of the canopy air space (J m-3)
 
      call T2enthTemp_veg(&
                  (heightCanopyTop-heightCanopyBottom), & ! intent(in): canopy depth (m)
                  specificHeatVeg,        & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                  maxMassVegetation,      & ! intent(in): maximum mass of vegetation (kg m-2)
                  snowfrz_scale,          & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                  scalarCanopyTemp,       & ! intent(in): canopy temperature (K)
                  scalarTheta,            & ! intent(in): canopy water content (kg m-2)
                  scalarCanopyEnthTemp)     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
     scalarCanopyEnthalpy = scalarCanopyEnthTemp - LH_fus * scalarCanopyIce/ (heightCanopyTop-heightCanopyBottom)
    end if

    ! number of layers
    nSnow   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
    nLake   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
    nSoil   = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
    nIce    = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nIce
    nLayers = nSnow + nLake + nSoil + nIce

    ! loop through all layers
    do iLayer=1,nLayers

     ! *****
     ! * check that the initial volumetric fraction of liquid water and ice is reasonable...
     ! *************************************************************************************
     select case(mlayerLayerType(iLayer))

      ! ***** snow, ice, lake, volume expansion allowed
      case(iname_snow, iname_lake, iname_ice)
       iSoil       = integerMissing
       vGn_m       = realMissing
       scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
       ! (check liquid water)
       if(mLayerVolFracLiq(iLayer) < 0._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < 0: layer = ',iLayer; err=20; return; end if
       if(mLayerVolFracLiq(iLayer) > 1._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > 1: layer = ',iLayer; err=20; return; end if
       ! (check ice)
       if (mlayerLayerType(iLayer)==iname_snow) then
         if(mLayerVolFracIce(iLayer) > 0.80_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > 0.80: layer = ',iLayer; err=20; return; end if
         if(scalarTheta > 0.80_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > 0.80: layer = '    ,iLayer; err=20; return; end if
       else ! glacier ice or lake (could be all ice)
         if(mLayerVolFracIce(iLayer) > 1._rkind  )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > 1: layer = ',iLayer; err=20; return; end if
         if(scalarTheta > 1._rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > 1: layer = '      ,iLayer; err=20; return; end if
       end if
       if (mlayerLayerType(iLayer)==iname_lake) then ! lake could be all liquid
         if(mLayerVolFracIce(iLayer) < 0._rkind  )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0: layer = '   ,iLayer; err=20; return; end if
       else if (mlayerLayerType(iLayer)==iname_ice) then ! glacier ice should be mostly ice
         if(mLayerVolFracIce(iLayer) < 0.80_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.80: layer = ',iLayer; err=20; return; end if
       else if (mlayerLayerType(iLayer)==iname_snow) then ! 
         if(mLayerVolFracIce(iLayer) < 0.05_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.05: layer = ',iLayer; err=20; return; end if
       end if
       ! check total water
       if(scalarTheta < 0.05_rkind)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < 0.05: layer = ',iLayer; err=20; return; end if

      ! ***** soil, no volume expansion
      case(iname_soil)
       iSoil       = iLayer - nSnow - nLake
       vGn_m       = 1._rkind - 1._rkind/vGn_n(iSoil)
       scalarTheta = mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)
       ! (check liquid water)
       if(mLayerVolFracLiq(iLayer) < theta_res(iSoil)-xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < theta_res: layer = ',iLayer; err=20; return; end if
       if(mLayerVolFracLiq(iLayer) > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > theta_sat: layer = ',iLayer; err=20; return; end if
       ! (check ice)
       if(mLayerVolFracIce(iLayer) < 0._rkind             )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0: layer = '        ,iLayer; err=20; return; end if
       if(mLayerVolFracIce(iLayer) > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > theta_sat: layer = ',iLayer; err=20; return; end if
       ! check total water
       if(scalarTheta < theta_res(iSoil)-xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < theta_res: layer = ',iLayer; err=20; return; end if
       if(scalarTheta > theta_sat(iSoil)+xTol)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > theta_sat: layer = ',iLayer; err=20; return; end if

      case default
       write(*,*) 'Cannot recognize case in initial vol water/ice check: type=', mlayerLayerType(iLayer)
       err=20; message=trim(message)//'cannot identify layer type'; return
     end select

     ! *****
     ! * check that the initial conditions are consistent with the constitutive functions...
     ! *************************************************************************************
     select case(mLayerLayerType(iLayer))

      ! ** snow, lake, ice
      case(iname_snow, iname_lake, iname_ice)

       ! check that snow temperature is less than freezing
       if(mLayerTemp(iLayer) > Tfreeze)then
        message=trim(message)//'initial snow temperature is greater than freezing'
        err=20; return
       end if

       ! ensure consistency among state variables
       call updateSnLaIc(&
                       ! input
                       mLayerTemp(iLayer),             & ! intent(in): temperature (K)
                       scalarTheta,                    & ! intent(in): volumetric fraction of total water (-)
                       snowfrz_scale,                  & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                       ! output
                       mLayerVolFracLiq(iLayer),       & ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce(iLayer),       & ! intent(out): volumetric fraction of ice (-)
                       fLiq,                           & ! intent(out): fraction of liquid water (-)
                       err,cmessage)                     ! intent(out): error control
       if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

       if(checkEnthalpy)then ! enthalpy as state variable (cold start often only has temperature)
          call T2enthTemp_snLaIc(&
                       ! input
                       snowfrz_scale,                  & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                       mLayerTemp(iLayer),             & ! intent(in):  layer temperature (K)
                       scalarTheta,                    & ! intent(in):  volumetric total water content (-)
                       ! output
                       mLayerEnthTemp(iLayer),         & ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
                       err,cmessage)                     ! intent(out): error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
          mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
       endif

      ! ** soil
      case(iname_soil)

      ! ensure consistency among state variables
<<<<<<< HEAD
       call updateSoil(&
                      ! input
=======
      call updateSnow(&
                      mLayerTemp(iLayer),             & ! intent(in): temperature (K)
                      scalarTheta,                    & ! intent(in): volumetric fraction of total water (-)
                      snowfrz_scale,                  & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                      mLayerVolFracLiq(iLayer),       & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),       & ! intent(out): volumetric fraction of ice (-)
                      fLiq,                           & ! intent(out): fraction of liquid water (-)
                      err,cmessage)                     ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

      if(checkEnthalpy)then ! enthalpy as state variable (cold start often only has temperature)
         call T2enthTemp_snow(&
                      snowfrz_scale,                  & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                      mLayerTemp(iLayer),             & ! intent(in):  layer temperature (K)
                      scalarTheta,                    & ! intent(in):  volumetric total water content (-)
                      mLayerEnthTemp(iLayer))           ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
         mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
      endif

     ! ** soil
     case(iname_soil)

      ! ensure consistency among state variables
      call updateSoil(&
>>>>>>> develop
                      mLayerTemp(iLayer),              & ! intent(in): layer temperature (K)
                      mLayerMatricHead(iLayer-nSnow-nLake),  & ! intent(in): matric head (m)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      scalarTheta,                     & ! intent(out): volumetric fraction of total water (-)
                      mLayerVolFracLiq(iLayer),        & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),        & ! intent(out): volumetric fraction of ice (-)
                      err,cmessage)                      ! intent(out): error control
       if(err/=0)then; message=trim(cmessage)//trim(cmessage); return; end if  ! (check for errors)

<<<<<<< HEAD
       if(checkEnthalpy)then ! enthalpy as state variable (cold start often only has temperature)
        call T2enthTemp_soil(&
                      ! input
=======
      if(checkEnthalpy)then ! enthalpy as state variable (cold start often only has temperature)
         call T2enthTemp_soil(&
>>>>>>> develop
                      use_lookup,                      & ! intent(in):  flag to use the lookup table for soil enthalpy
                      soil_dens_intr(iSoil),           & ! intent(in):  intrinsic soil density (kg m-3)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      iSoil,                           & ! intent(in):  index of the control volume within the domain
                      lookupData%gru(iGRU)%hru(iHRU)%dom(iDOM),  & ! intent(in):  lookup table data structure
                      realMissing,                     & ! intent(in):  lower value of integral (not computed)
                      mLayerTemp(iLayer),              & ! intent(in):  layer temperature (K)
<<<<<<< HEAD
                      mLayerMatricHead(iLayer-nSnow-nLake),  & ! intent(in):  matric head (m)
                     ! output
                      mLayerEnthTemp(iLayer),          & ! intent(out): temperature component of enthalpy soil layer (J m-3)
                      err,cmessage)                      ! intent(out): error control      
        if(err/=0)then; message=trim(cmessage)//trim(cmessage); return; end if  ! (check for errors)
        mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
       endif
=======
                      mLayerMatricHead(iLayer-nSnow),  & ! intent(in):  matric head (m)
                      mLayerEnthTemp(iLayer))            ! intent(out): temperature component of enthalpy soil layer (J m-3)
         mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
      endif
>>>>>>> develop

      case default; err=10; message=trim(message)//'unknown case for model layer'; return
     end select

    end do  ! (looping through layers)
    
    ! end association to variables in the data structures
    end associate

    ! if snow layers exist, compute snow depth and SWE
    if(nSnow > 0)then
     progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSWE)%dat(1) = sum( (progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                                                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice)  * &
                                                                                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
    end if  ! if snow layers exist

    ! check that the layering is consistent
    do iLayer=1,nLayers
     h1 = sum(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerDepth)%dat(1:iLayer)) ! sum of the depths up to the current layer
     h2 = progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%iLayerHeight)%dat(iLayer) - progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%iLayerHeight)%dat(0)  ! difference between snow-atm interface and bottom of layer
     if(abs(h1 - h2) > 1.e-6_rkind)then
      write(message,'(a,1x,i0,a,f5.3,a,f5.3)') trim(message)//'mis-match between layer depth and layer height; layer = ', iLayer, '; sum depths = ',h1,'; height = ',h2
      err=20; return
     end if
    end do

   end do ! iDOM
  end do ! iHRU
 end do ! iGRU

 end subroutine check_icond

end module check_icond_module
