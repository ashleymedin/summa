

module computResidWithPrime_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil
USE globalData,only:iname_ice       ! named variables for ice
USE globalData,only:iname_lake      ! named variables for lake

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! privacy
implicit none
public::computResidWithPrime
contains

! **********************************************************************************************************
! public subroutine computResidWithPrime: compute the residual vector
! **********************************************************************************************************
subroutine computResidWithPrime(&
                      ! input: model control
                      dt,                        & ! intent(in):  length of the time step (seconds)
                      nSnow,                     & ! intent(in):  number of snow layers
                      nLake,                     & ! intent(in):  number of lake layers
                      nSoil,                     & ! intent(in):  number of soil layers
                      nLayers,                   & ! intent(in):  total number of layers
                      enthalpyStateVec,               & ! intent(in):  flag if enthalpy is state variable
                      ! input: flux vectors
                      sMul,                      & ! intent(in):  state vector multiplier (used in the residual calculations)
                      fVec,                      & ! intent(in):  flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanairTempPrime,     & ! intent(in):  prime value for the temperature of the canopy air space (K s-1)
                      scalarCanopyTempPrime,     & ! intent(in):  prime value for the temperature of the vegetation canopy (K s-1)
                      scalarCanopyWatPrime,      & ! intent(in):  prime value for the water on the vegetation canopy (kg m-2 s-1)
                      mLayerTempPrime,           & ! intent(in):  prime vector of the temperature of each layer (K s-1)
                      scalarAquiferStoragePrime, & ! intent(in):  prime value for storage of water in the aquifer (m s-1)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIcePrime,      & ! intent(in):  prime value for the ice on the vegetation canopy (kg m-2 s-1)
                      scalarCanopyLiqPrime,      & ! intent(in):  prime value for the liq on the vegetation canopy (kg m-2 s-1)
                      mLayerVolFracIcePrime,     & ! intent(in):  prime vector of the volumetric ice in each layer (s-1)
                      mLayerVolFracWatPrime,     & ! intent(in):  prime vector of the volumetric water in each layer (s-1)
                      mLayerVolFracLiqPrime,     & ! intent(in):  prime vector of the volumetric liq in each layer (s-1)
                      ! input: enthalpy terms
                      scalarCanopyCmTrial,       & ! intent(in):  Cm of vegetation canopy (J kg K-1)
                      mLayerCmTrial,             & ! intent(in):  Cm of each layer (J kg K-1)
                      scalarCanairEnthalpyPrime, & ! intent(in):  prime value for the enthalpy of the canopy air space (W m-3)
                      scalarCanopyEnthalpyPrime, & ! intent(in):  prime value for the of enthalpy of the vegetation canopy (W m-3)
                      mLayerEnthalpyPrime,       & ! intent(in):  prime vector of the of enthalpy of each layer (W m-3)
                      ! input: data structures
                      prog_data,                 & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,                 & ! intent(in):  model diagnostic variables for a local HRU
                      flux_data,                 & ! intent(in):  model fluxes for a local HRU
                      indx_data,                 & ! intent(in):  index data
                      ! output
                      rAdd,                      & ! intent(out): additional (sink) terms on the RHS of the state equation
                      rVec,                      & ! intent(out): residual vector
                      err,message)                 ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  real(rkind),intent(in)          :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
  integer(i4b),intent(in)         :: nLake                     ! number of lake layers
  integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the layer domains
  logical(lgt),intent(in)         :: enthalpyStateVec               ! flag if enthalpy is state variable
  ! input: flux vectors
  real(qp),intent(in)             :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
  real(rkind),intent(in)          :: fVec(:)                   ! flux vector
  ! input: state variables (already disaggregated into scalars and vectors)
  real(rkind),intent(in)          :: scalarCanairTempPrime     ! prime value for temperature of the canopy air space (K s-1)
  real(rkind),intent(in)          :: scalarCanopyTempPrime     ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),intent(in)          :: scalarCanopyWatPrime      ! prime value for canopy total water content (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerTempPrime(:)        ! prime vector of temperature of each snow/soil layer (K s-1) content
  real(rkind),intent(in)          :: scalarAquiferStoragePrime ! prime value of aquifer storage (m s-1)
  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
  real(rkind),intent(in)          :: scalarCanopyIcePrime      ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: scalarCanopyLiqPrime      ! prime value for the liq on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerVolFracIcePrime(:)  ! prime vector of volumetric fraction of ice (s-1)
  real(rkind),intent(in)          :: mLayerVolFracWatPrime(:)  ! prime vector of the volumetric water in each layer (s-1)
  real(rkind),intent(in)          :: mLayerVolFracLiqPrime(:)  ! prime vector of the volumetric water in each layer (s-1)
  ! input: enthalpy terms
  real(qp),intent(in)             :: scalarCanopyCmTrial       ! Cm of vegetation canopy (-)
  real(qp),intent(in)             :: mLayerCmTrial(:)          ! Cm of each layer (-)
  real(rkind),intent(in)          :: scalarCanairEnthalpyPrime ! prime value for enthalpy of the canopy air space (W m-3)
  real(rkind),intent(in)          :: scalarCanopyEnthalpyPrime ! prime value for enthalpy of the vegetation canopy (W m-3)
  real(rkind),intent(in)          :: mLayerEnthalpyPrime(:)    ! prime vector of enthalpy of each layer (W m-3)
  ! input: data structures
  type(var_dlength),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
  type(var_dlength),intent(in)    :: diag_data                 ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)    :: flux_data                 ! model fluxes for a local HRU
  type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
  ! output
  real(rkind),intent(out)         :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
  real(qp),intent(out)            :: rVec(:)   ! NOTE: qp      ! residual vector
  integer(i4b),intent(out)        :: err                       ! error code
  character(*),intent(out)        :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  integer(i4b)                     :: iLayer                   ! index of layer within the layer domains
  integer(i4b),parameter           :: ixVegVolume=1            ! index of the desired vegetation control volumne (currently only one veg layer)
  real(rkind)                      :: scalarCanopyHydPrime     ! trial value for canopy water (kg m-2), either liquid water content or total water content
  real(rkind),dimension(nLayers)   :: mLayerVolFracHydPrime    ! vector of volumetric water content (-), either liquid water content or total water content
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! link to the necessary variables for the residual computations
  associate(&
    ! canopy and layer depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp]      canopy depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model fluxes (sink terms in the soil domain)
    mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)
    mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
    ! number of state variables of a specific type
    nSlicSoilNrg            => indx_data%var(iLookINDEX%nSlicSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the layer domains
    nSlicSoilHyd            => indx_data%var(iLookINDEX%nSlicSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the layer domains
    nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! model indices
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ixSlicSoilNrg           => indx_data%var(iLookINDEX%ixSlicSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixSlicSoilHyd           => indx_data%var(iLookINDEX%ixSlicSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in layer domains
    layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in layer domains
    ) ! association to necessary variables for the residual computations
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computResidWithPrime/"

    ! ---
    ! * compute sink terms...
    ! -----------------------

    ! intialize additional terms on the RHS as zero
    rAdd(:) = 0._rkind

    ! add melt freeze terms only if not using enthalpy terms 
    ! NOTE: would need to use these if were using enthTemp terms
    if(.not.enthalpyStateVec)then
      ! compute energy associated with melt freeze for the vegetation canopy
      if(ixVegNrg/=integerMissing) rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*scalarCanopyIcePrime/canopyDepth   ! energy associated with melt/freeze (J m-3)
 
      ! compute energy associated with melt/freeze for snow
      ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
      if(nSlicSoilNrg>0)then
        do concurrent (iLayer=1:nLayers,ixSlicSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the layer domains)
          select case( layerType(iLayer) )
            case(iname_snow, iname_lake, iname_ice); rAdd( ixSlicSoilNrg(iLayer) ) = rAdd( ixSlicSoilNrg(iLayer) ) + LH_fus*iden_ice * mLayerVolFracIcePrime(iLayer)
            case(iname_soil);                        rAdd( ixSlicSoilNrg(iLayer) ) = rAdd( ixSlicSoilNrg(iLayer) ) + LH_fus*iden_water * mLayerVolFracIcePrime(iLayer)
          end select
        end do  ! looping through non-missing energy state variables in the layer domains
      endif

    endif

    ! sink terms soil hydrology (-)
    ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
    ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
    ! NOTE 3: rAdd for all other Hyd is =0, and is defined in the initialization above
    ! NOTE 4: same sink terms for matric head and liquid matric potential
    if(nSoilOnlyHyd>0)then
      do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the layer domains)
       rAdd( ixSoilOnlyHyd(iLayer) ) = rAdd( ixSoilOnlyHyd(iLayer) ) + ( ( mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow+nLake) - mLayerCompress(iLayer) )*dt
      end do  ! looping through non-missing energy state variables in the layer domains
    endif

    ! ---
    ! * compute the residual vector...
    ! --------------------------------

    ! compute the residual vector for the vegetation canopy
    ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
    ! --> energy balance
    if(enthalpyStateVec)then
      if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = scalarCanairEnthalpyPrime - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
      if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = scalarCanopyEnthalpyPrime - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
    else
      if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = sMul(ixCasNrg) * scalarCanairTempPrime - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
      if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = sMul(ixVegNrg) * scalarCanopyTempPrime + scalarCanopyCmTrial * scalarCanopyWatPrime/canopyDepth &
                                                   - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
    endif                                               
    ! --> mass balance
    if(ixVegHyd/=integerMissing)then
      scalarCanopyHydPrime = merge(scalarCanopyWatPrime, scalarCanopyLiqPrime, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
      rVec(ixVegHyd) = sMul(ixVegHyd)*scalarCanopyHydPrime - ( fVec(ixVegHyd)*dt + rAdd(ixVegHyd) )
    endif

    ! compute the residual vector for the snow and soil sub-domains for energy
    if(nSlicSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSlicSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the layer domains)
        if(enthalpyStateVec)then
          rVec( ixSlicSoilNrg(iLayer) ) = mLayerEnthalpyPrime(iLayer) - ( fVec( ixSlicSoilNrg(iLayer) )*dt + rAdd( ixSlicSoilNrg(iLayer) ) )
        else
          rVec( ixSlicSoilNrg(iLayer) ) = sMul( ixSlicSoilNrg(iLayer) ) * mLayerTempPrime(iLayer) + mLayerCmTrial(iLayer) * mLayerVolFracWatPrime(iLayer) &
                                         - ( fVec( ixSlicSoilNrg(iLayer) )*dt + rAdd( ixSlicSoilNrg(iLayer) ) )
        endif
      end do  ! looping through non-missing energy state variables in the layer domains
    endif

    ! compute the residual vector for the snow and soil sub-domains for hydrology
    ! NOTE: residual depends on choice of state variable
    if(nSlicSoilHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSlicSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the layer domains)
        ! (get the correct state variable)
        mLayerVolFracHydPrime(iLayer) = merge(mLayerVolFracWatPrime(iLayer), mLayerVolFracLiqPrime(iLayer), (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
        ! (compute the residual)
        rVec( ixSlicSoilHyd(iLayer) ) = mLayerVolFracHydPrime(iLayer) - ( fVec( ixSlicSoilHyd(iLayer) )*dt + rAdd( ixSlicSoilHyd(iLayer) ) )
      end do  ! looping through non-missing energy state variables in the layer domains
    endif

    ! compute the residual vector for the aquifer
    if(ixAqWat/=integerMissing)  rVec(ixAqWat) = sMul(ixAqWat)*scalarAquiferStoragePrime - ( fVec(ixAqWat)*dt + rAdd(ixAqWat) )

    ! print result
    if(globalPrintFlag)then
      write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
    endif

    ! check
    if(any(isNan(rVec)))then
      message=trim(message)//'vector of residuals contains NaN value(s) ' ! formerly known as the Indian bread error
      write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      err=20; return
    endif
    
  end associate

end subroutine computResidWithPrime

end module computResidWithPrime_module
