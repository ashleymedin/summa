


module summaSolveSundialsIDA_module


!======= Inclusions ===========
USE, intrinsic :: iso_c_binding
USE nrtype
USE type4IDA

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:model_decisions ! model decision structure

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes
USE globalData,only:diag_meta       ! metadata on the model diagnostic variables
USE globalData,only:prog_meta       ! metadata on the model prognostic variables
USE globalData,only:deriv_meta      ! metadata on the model derivatives

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J K-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookDERIV     ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup         ! data vector

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  qbaseTopmodel ! TOPMODEL-ish baseflow parameterization

! privacy
 implicit none
 private::setInitialCondition
 private::setSolverParams
 private::find_rootdir
 public::layerDisCont4IDA
 public::summaSolveSundialsIDA

contains

!-------------------
! * public subroutine summaSolveSundialsIDA: solve F(y,y') = 0 by IDA (y is the state vector)
! ------------------
subroutine summaSolveSundialsIDA(                         &
                      dt,                      & ! intent(in):    data time step
                      atol,                    & ! intent(in):    absolute telerance
                      rtol,                    & ! intent(in):    relative tolerance
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nStat,                   & ! intent(in):    total number of state variables
                      ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors
                      stateVecInit,            & ! intent(in):    initial state vector
                      sMul,                    & ! intent(inout): state vector multiplier (USEd in the residual calculations)
                      dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                      ! input: data structures
                      lookup_data,             & ! intent(in):    lookup tables
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(in):    index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                      ! output
                      ixSaturation,            & ! intent(inout) index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      idaSucceeds,             & ! intent(out):   flag to indicate if ida successfully solved the problem in current data step
                      tooMuchMelt,             & ! intent(inout):   flag to denote that there was too much melt
                      dt_out,                  & ! intent(out):   time step sum for entire data window at termination of sundials
                      stateVec,                & ! intent(out):   model state vector
                      stateVecPrime,           & ! intent(out):   derivative of model state vector
                      err,message              & ! intent(out):   error control
                    )

  !======= Inclusions ===========
  USE fida_mod                                    ! Fortran interface to IDA
  USE fsundials_context_mod                       ! Fortran interface to SUNContext
  USE fnvector_serial_mod                         ! Fortran interface to serial N_Vector
  USE fsunmatrix_dense_mod                        ! Fortran interface to dense SUNMatrix
  USE fsunlinsol_dense_mod                        ! Fortran interface to dense SUNLinearSolver
  USE fsunmatrix_band_mod                         ! Fortran interface to banded SUNMatrix
  USE fsunlinsol_band_mod                         ! Fortran interface to banded SUNLinearSolver
  USE fsunnonlinsol_newton_mod                    ! Fortran interface to Newton SUNNonlinearSolver
  USE fsundials_matrix_mod                        ! Fortran interface to generic SUNMatrix
  USE fsundials_nvector_mod                       ! Fortran interface to generic N_Vector
  USE fsundials_linearsolver_mod                  ! Fortran interface to generic SUNLinearSolver
  USE fsundials_nonlinearsolver_mod               ! Fortran interface to generic SUNNonlinearSolver
  USE allocspace_module,only:allocLocal           ! allocate local data structures
  USE eval8summaSundials_module,only:eval8summa4IDA         ! DAE/ODE functions
  USE eval8summaSundials_module,only:eval8summaSundials     ! residual of DAE
  USE computJacobSundials_module,only:computJacob4IDA       ! system Jacobian
  USE tol4IDA_module,only:computWeight4IDA        ! weigth required for tolerances
  USE var_derive_module,only:calcHeight           ! height at layer interfaces and layer mid-point

  !======= Declarations =========
  implicit none

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! calling variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(qp),intent(in)             :: dt                     ! data time step
  real(qp),intent(inout)          :: atol(:)                ! vector of absolute tolerances
  real(qp),intent(inout)          :: rtol(:)                ! vector of relative tolerances
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nStat                  ! total number of state variables
  integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  ! input: state vectors
  real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
  real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
  real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
  ! input: data structures
  type(zLookup),intent(in)        :: lookup_data            ! lookup tables
  type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  type(var_d),        intent(in)  :: forc_data              ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
  type(var_ilength),  intent(in)  :: indx_data              ! indices defining model states and layers
  ! input-output: data structures
  type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
  ! output: state vectors
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
  real(rkind),intent(inout)       :: stateVecPrime(:)       ! model state vector (y')
  logical(lgt),intent(out)        :: idaSucceeds            ! flag to indicate if IDA is successful
  logical(lgt),intent(inout)      :: tooMuchMelt            ! flag to denote that there was too much melt
  real(qp),intent(out)            :: dt_out                 ! time step sum for entire data window at termination of sundials
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  type(N_Vector),           pointer :: sunvec_y             ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_yp            ! sundials derivative vector
  type(N_Vector),           pointer :: sunvec_av            ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A             ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS         ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS        ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem              ! IDA memory
  type(c_ptr)                       :: sunctx               ! SUNDIALS simulation context
  type(eqnsData),           target  :: eqns_data            ! IDA type
  integer(i4b)                      :: retval, retvalr      ! return value
  logical(lgt)                      :: feasible             ! feasibility flag
  real(qp)                          :: t0                   ! staring time
  real(qp)                          :: dt_last(1)           ! last time step
  real(qp)                          :: dt_diff              ! difference from previous timestep
  integer(kind = 8)                 :: mu, lu               ! in banded matrix mode in Sundials type
  integer(c_long)                   :: nState               ! total number of state variables in Sundials type
  real(rkind)                       :: rVec(nStat)          ! residual vector
  integer(i4b)                      :: iVar, i              ! indices
  integer(i4b)                      :: nRoot                ! total number of roots (events) to find
  real(qp)                          :: tret(1)              ! time in data window
  real(qp)                          :: tretPrev             ! previous time in data window
  integer(i4b),allocatable          :: rootsfound(:)        ! crossing direction of discontinuities
  integer(i4b),allocatable          :: rootdir(:)           ! forced crossing direction of discontinuities
  logical(lgt)                      :: tinystep             ! if step goes below small size
  logical(lgt),parameter            :: offErrWarnMessage = .true. ! flag to turn IDA warnings off, default true
  logical(lgt),parameter            :: detect_events = .true.     ! flag to do event detection and restarting, default true
  logical(lgt),parameter            :: use_fdJac = .false.        ! flag to use finite difference Jacobian, default false

  ! -----------------------------------------------------------------------------------------------------

  ! initialize error control
  err=0; message="summaSolveSundialsIDA/"

  nState = nStat
  idaSucceeds = .true.
  ! fill eqns_data which will be required later to call eval8summaSundials
  eqns_data%dt                      = dt
  eqns_data%nSnow                   = nSnow
  eqns_data%nSoil                   = nSoil
  eqns_data%nLayers                 = nLayers
  eqns_data%nState                  = nState
  eqns_data%ixMatrix                = ixMatrix
  eqns_data%firstSubStep            = firstSubStep
  eqns_data%computeVegFlux          = computeVegFlux
  eqns_data%scalarSolution          = scalarSolution

  allocate( eqns_data%atol(nState) )
  eqns_data%atol = atol

  allocate( eqns_data%rtol(nState) )
  eqns_data%rtol = rtol

  allocate( eqns_data%sMul(nState) )
  eqns_data%sMul                    = sMul

  allocate( eqns_data%dMat(nState) )
  eqns_data%dMat                    = dMat

  ! allocate space for the temporary prognostic variable structure
  call allocLocal(prog_meta(:),eqns_data%prog_data,nSnow,nSoil,err,message)
  if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
  eqns_data%prog_data               = prog_data

  ! allocate space for the temporary diagnostic variable structure
  call allocLocal(diag_meta(:),eqns_data%diag_data,nSnow,nSoil,err,message)
  if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
  eqns_data%diag_data               = diag_data

  ! allocate space for the temporary flux variable structure
  call allocLocal(flux_meta(:),eqns_data%flux_data,nSnow,nSoil,err,message)
  if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
  eqns_data%flux_data               = flux_data

  ! allocate space for the derivative structure
  call allocLocal(deriv_meta(:),eqns_data%deriv_data,nSnow,nSoil,err,message)
  if(err/=0)then; err=20; message=trim(message)//trim(message); return; end if
  eqns_data%deriv_data              = deriv_data

  eqns_data%lookup_data             = lookup_data
  eqns_data%type_data               = type_data
  eqns_data%attr_data               = attr_data
  eqns_data%mpar_data               = mpar_data
  eqns_data%forc_data               = forc_data
  eqns_data%bvar_data               = bvar_data
  eqns_data%indx_data               = indx_data

  ! allocate space
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel)then
    allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
  else
    allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
  end if
  allocate( eqns_data%mLayerMatricHeadLiqTrial(nSoil) )
  allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
  allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
  allocate( eqns_data%mLayerVolFracWatTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracWatPrev(nLayers) )
  allocate( eqns_data%mLayerTempTrial(nLayers) )
  allocate( eqns_data%mLayerTempPrev(nLayers) )
  allocate( eqns_data%mLayerVolFracIceTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracIcePrev(nLayers) )
  allocate( eqns_data%mLayerVolFracLiqTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracLiqPrev(nLayers) )
  allocate( eqns_data%mLayerEnthalpyTrial(nLayers) )
  allocate( eqns_data%mLayerEnthalpyPrev(nLayers) )
  allocate( eqns_data%fluxVec(nState) )
  allocate( eqns_data%resSink(nState) )

  retval = FSUNContext_Create(c_null_ptr, sunctx)

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(nState, stateVec, sunctx)
  if (.not. associated(sunvec_y)) then; err=20; message='summaSolveSundialsIDA: sunvec = NULL'; return; endif

  sunvec_yp => FN_VMake_Serial(nState, stateVecPrime, sunctx)
  if (.not. associated(sunvec_yp)) then; err=20; message='summaSolveSundialsIDA: sunvec = NULL'; return; endif

  ! Initialize solution vectors
  call setInitialCondition(nState, stateVecInit, sunvec_y, sunvec_yp)

  ! Create memory
  ida_mem = FIDACreate(sunctx)
  if (.not. c_associated(ida_mem)) then; err=20; message='summaSolveSundialsIDA: ida_mem = NULL'; return; endif

  ! Attach user data to memory
  eqns_data%ida_mem = ida_mem
  retval = FIDASetUserData(ida_mem, c_loc(eqns_data))
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDASetUserData'; return; endif

  ! Initialize memory
  t0 = 0._rkind
  retval = FIDAInit(ida_mem, c_funloc(eval8summa4IDA), t0, sunvec_y, sunvec_yp)
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDAInit'; return; endif

  ! set tolerances
  retval = FIDAWFtolerances(ida_mem, c_funloc(computWeight4IDA))
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDAWFtolerances'; return; endif

  ! initialize rootfinding problem and allocate space, counting roots
  if(detect_events)then
    nRoot = 0
    if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing) nRoot = nRoot+1
    if(nSnow>0)then
      do i = 1,nSnow
        if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing) nRoot = nRoot+1
      enddo
    endif
    if(nSoil>0)then
      do i = 1,nSoil
        if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing) nRoot = nRoot+1
        if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing) nRoot = nRoot+1
      enddo
    endif
    allocate( rootsfound(nRoot) )
    allocate( rootdir(nRoot) )
    rootdir = 0
    retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4IDA))
    if (retval /= 0) then; err=20; message='solveByIDA: error in FIDARootInit'; return; endif
  else ! will not use, allocate at something
    nRoot = 1
    allocate( rootsfound(nRoot) )
    allocate( rootdir(nRoot) )
  endif

  ! define the form of the matrix
  select case(ixMatrix)
    case(ixBandMatrix)
      mu = ku; lu = kl;
      ! Create banded SUNMatrix for use in linear solves
      sunmat_A => FSUNBandMatrix(nState, mu, lu, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolveSundialsIDA: sunmat = NULL'; return; endif

      ! Create banded SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Band(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolveSundialsIDA: sunlinsol = NULL'; return; endif

    case(ixFullMatrix)
      ! Create dense SUNMatrix for use in linear solves
      sunmat_A => FSUNDenseMatrix(nState, nState, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolveSundialsIDA: sunmat = NULL'; return; endif

      ! Create dense SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolveSundialsIDA: sunlinsol = NULL'; return; endif

      ! check
    case default;  err=20; message='summaSolveSundialsIDA: error in type of matrix'; return

  end select  ! form of matrix

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDASetLinearSolver'; return; endif

  ! Set the user-supplied Jacobian routine
  if(.not.use_fdJac)then
    retval = FIDASetJacFn(ida_mem, c_funloc(computJacob4IDA))
    if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDASetJacFn'; return; endif
  endif

  ! Create Newton SUNNonlinearSolver object
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
  if (.not. associated(sunnonlin_NLS)) then; err=20; message='summaSolveSundialsIDA: sunnonlinsol = NULL'; return; endif

  ! Attach the nonlinear solver
  retval = FIDASetNonlinearSolver(ida_mem, sunnonlin_NLS)
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDASetNonlinearSolver'; return; endif

  ! Enforce the solver to stop at end of the time step
  retval = FIDASetStopTime(ida_mem, dt)
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in FIDASetStopTime'; return; endif

  ! Set solver parameters such as maximum order, number of iterations, ...
  call setSolverParams(dt, nint(mpar_data%var(iLookPARAM%maxiter)%dat(1)), ida_mem, retval)
  if (retval /= 0) then; err=20; message='summaSolveSundialsIDA: error in setSolverParams'; return; endif

  ! Disable error messages and warnings
  if(offErrWarnMessage) then
    retval = FIDASetErrFile(ida_mem, c_null_ptr)
    retval = FIDASetNoInactiveRootWarn(ida_mem)
  endif

  ! need the following values for the first substep
  eqns_data%scalarCanopyTempPrev     = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)
  eqns_data%scalarCanopyIcePrev      = prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)
  eqns_data%scalarCanopyLiqPrev      = prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)
  eqns_data%scalarCanopyEnthalpyPrev = diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)
  eqns_data%mLayerTempPrev(:)        = prog_data%var(iLookPROG%mLayerTemp)%dat(:)
  eqns_data%mLayerMatricHeadPrev(:)  = prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
  eqns_data%mLayerVolFracWatPrev(:)  = prog_data%var(iLookPROG%mLayerVolFracWat)%dat(:)
  eqns_data%mLayerVolFracIcePrev(:)  = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(:)
  eqns_data%mLayerVolFracLiqPrev(:)  = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(:)
  eqns_data%mLayerEnthalpyPrev(:)    = diag_data%var(iLookDIAG%mLayerEnthalpy)%dat(:)
  eqns_data%scalarAquiferStoragePrev = prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)
  eqns_data%ixSaturation             = ixSaturation

  tinystep = .false.

  !**********************************************************************************
  !****************************** Main Solver ***************************************
  !************************* loop on one_step mode **********************************
  !**********************************************************************************

  tret(1) = t0           ! initial time
  do while(tret(1) < dt)

    tretPrev = tret(1)
    ! call this at beginning of step to reduce root bouncing (only looking in one direction)
    if(detect_events .and. .not.tinystep)then
      call find_rootdir(eqns_data, rootdir)
      retval = FIDASetRootDirection(ida_mem, rootdir)
      if (retval /= 0) then; err=20; message='solveByIDA: error in FIDASetRootDirection'; return; endif
    endif

    eqns_data%firstFluxCall = .false.
    eqns_data%firstSplitOper = .true.
    ! call IDASolve, advance solver just one internal step
    retvalr = FIDASolve(ida_mem, dt, tret, sunvec_y, sunvec_yp, IDA_ONE_STEP)
    if( retvalr < 0 )then
      idaSucceeds = .false.
      exit
    endif

    tooMuchMelt = .false.
    feasible = .true.
    ! loop through non-missing energy state variables in the snow domain to see if need to merge
    do concurrent (i=1:nSnow,indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)
      if (stateVec(indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)) > Tfreeze) tooMuchMelt = .true. !need to merge
    enddo
    if(tooMuchMelt)exit

    ! get the last stepsize and difference from previous end time, not necessarily the same
    retval = FIDAGetLastStep(ida_mem, dt_last)
    dt_diff = tret(1) - tretPrev

    ! compute the flux and the residual vector for a given state vector
    call eval8summaSundials(&
                  ! input: model control
                  dt_diff,                            & ! intent(in):    current stepsize
                  eqns_data%dt,                       & ! intent(in):    total data step
                  eqns_data%nSnow,                    & ! intent(in):    number of snow layers
                  eqns_data%nSoil,                    & ! intent(in):    number of soil layers
                  eqns_data%nLayers,                  & ! intent(in):    number of layers
                  eqns_data%nState,                   & ! intent(in):    number of state variables in the current subset
                  .false.,                            & ! intent(in):    check for feasibility once outside Sundials loop
                  eqns_data%firstSubStep,             & ! intent(in):    flag to indicate if we are processing the first sub-step
                  eqns_data%firstFluxCall,            & ! intent(inout): flag to indicate if we are processing the first flux call
                  eqns_data%firstSplitOper,           & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                  eqns_data%computeVegFlux,           & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  eqns_data%scalarSolution,           & ! intent(in):    flag to indicate the scalar solution
                  ! input: state vectors
                  stateVec,                           & ! intent(in):    model state vector
                  stateVecPrime,                      & ! intent(in):    model state vector
                  eqns_data%sMul,                     & ! intent(inout): state vector multiplier (used in the residual calculations)
                  ! input: data structures
                  model_decisions,                    & ! intent(in):    model decisions
                  eqns_data%lookup_data,              & ! intent(in):    lookup data
                  eqns_data%type_data,                & ! intent(in):    type of vegetation and soil
                  eqns_data%attr_data,                & ! intent(in):    spatial attributes
                  eqns_data%mpar_data,                & ! intent(in):    model parameters
                  eqns_data%forc_data,                & ! intent(in):    model forcing data
                  eqns_data%bvar_data,                & ! intent(in):    average model variables for the entire basin
                  eqns_data%prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                  ! input-output: data structures
                  eqns_data%indx_data,                & ! intent(inout): index data
                  eqns_data%diag_data,                & ! intent(inout): model diagnostic variables for a local HRU
                  eqns_data%flux_data,                & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                  eqns_data%deriv_data,               & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: here we need to pass some extra variables that do not get updated in in the Sundials loops
                  eqns_data%scalarCanopyTempTrial,    & ! intent(in):    trial value of canopy temperature (K)
                  eqns_data%scalarCanopyTempPrev,     & ! intent(in):    previous value of canopy temperature (K)
                  eqns_data%scalarCanopyIceTrial,     & ! intent(out):   trial value for mass of ice on the vegetation canopy (kg m-2)
                  eqns_data%scalarCanopyIcePrev,      & ! intent(in):    value for mass of ice on the vegetation canopy (kg m-2)
                  eqns_data%scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                  eqns_data%scalarCanopyLiqPrev,      & ! intent(in):    value of canopy liquid water (kg m-2)
                  eqns_data%scalarCanopyEnthalpyTrial,& ! intent(out):   trial value for enthalpy of the vegetation canopy (J m-3)
                  eqns_data%scalarCanopyEnthalpyPrev, & ! intent(in):    value for enthalpy of the vegetation canopy (J m-3)
                  eqns_data%mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                  eqns_data%mLayerTempPrev,           & ! intent(in):    vector of layer temperature (K)
                  eqns_data%mLayerMatricHeadLiqTrial, & ! intent(out):   trial value for liquid water matric potential (m)
                  eqns_data%mLayerMatricHeadTrial,    & ! intent(out):   trial value for total water matric potential (m)
                  eqns_data%mLayerMatricHeadPrev,     & ! intent(in):    value for total water matric potential (m)
                  eqns_data%mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                  eqns_data%mLayerVolFracWatPrev,     & ! intent(in):    vector of volumetric total water content (-)
                  eqns_data%mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                  eqns_data%mLayerVolFracIcePrev,     & ! intent(in):    vector of volumetric ice water content (-)
                  eqns_data%mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                  eqns_data%mLayerVolFracLiqPrev,     & ! intent(in):    vector of volumetric liquid water content (-)
                  eqns_data%scalarAquiferStorageTrial,& ! intent(out):   trial value of storage of water in the aquifer (m)
                  eqns_data%scalarAquiferStoragePrev, & ! intent(in):    value of storage of water in the aquifer (m)
                  eqns_data%mLayerEnthalpyPrev,       & ! intent(in):    vector of enthalpy for snow+soil layers (J m-3)
                  eqns_data%mLayerEnthalpyTrial,      & ! intent(out):   trial vector of enthalpy for snow+soil layers (J m-3)
                  ! input-output: baseflow
                  eqns_data%ixSaturation,             & ! intent(inout): index of the lowest saturated layer
                  eqns_data%dBaseflow_dMatric,        & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                  ! output: flux and residual vectors
                  feasible,                           & ! intent(out):   flag to denote the feasibility of the solution
                  eqns_data%fluxVec,                  & ! intent(out):   flux vector
                  eqns_data%resSink,                  & ! intent(out):   additional (sink) terms on the RHS of the state equation
                  rVec,                               & ! intent(out):   residual vector
                  eqns_data%err,eqns_data%message)      ! intent(out):   error control

    ! sum of fluxes smoothed over timestep
    do iVar=1,size(flux_meta)
      flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + eqns_data%flux_data%var(iVar)%dat(:) * dt_diff
    end do

    ! sum of mLayerCmpress over the time step, since using prev value
    ! Note: correct for difference between dt_diff and dt_last(1), e.g. if dt_last is longer than dt_diff should be adding less than the whole
    mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:) &
                                    * ( eqns_data%mLayerMatricHeadTrial(:) - eqns_data%mLayerMatricHeadPrev(:) )*dt_diff/dt_last(1)
print*,dt_diff,dt_last(1),tret(1),tretPrev
    ! save required quantities for next step
    eqns_data%scalarCanopyTempPrev     = eqns_data%scalarCanopyTempTrial
    eqns_data%scalarCanopyIcePrev      = eqns_data%scalarCanopyIceTrial
    eqns_data%scalarCanopyLiqPrev      = eqns_data%scalarCanopyLiqTrial
    eqns_data%scalarCanopyEnthalpyPrev = eqns_data%scalarCanopyEnthalpyTrial
    eqns_data%mLayerTempPrev(:)        = eqns_data%mLayerTempTrial(:)
    eqns_data%mLayerMatricHeadPrev(:)  = eqns_data%mLayerMatricHeadTrial(:)
    eqns_data%mLayerVolFracWatPrev(:)  = eqns_data%mLayerVolFracWatTrial(:)
    eqns_data%mLayerVolFracIcePrev(:)  = eqns_data%mLayerVolFracIceTrial(:)
    eqns_data%mLayerVolFracLiqPrev(:)  = eqns_data%mLayerVolFracLiqTrial(:)
    eqns_data%mLayerEnthalpyPrev(:)    = eqns_data%mLayerEnthalpyTrial(:)
    eqns_data%scalarAquiferStoragePrev = eqns_data%scalarAquiferStorageTrial

    ! Restart for where vegetation and layers cross freezing point
    if(detect_events)then
      if (retvalr .eq. IDA_ROOT_RETURN) then !IDASolve succeeded and found one or more roots at tret(1)
        ! rootsfound[i]= +1 indicates that gi is increasing, -1 g[i] decreasing, 0 no root
        !retval = FIDAGetRootInfo(ida_mem, rootsfound)
        !if (retval < 0) then; err=20; message='solveByIDA: error in FIDAGetRootInfo'; return; endif
        !print '(a,f15.7,2x,17(i2,2x))', "time, rootsfound[] = ", tret(1), rootsfound
        ! Reininitialize solver for running after discontinuity and restart
        retval = FIDAReInit(ida_mem, tret(1), sunvec_y, sunvec_yp)
        if (retval /= 0) then; err=20; message='solveByIDA: error in FIDAReInit'; return; endif
        if(dt_last(1) < 0.01_rkind)then ! don't keep calling if step is small
          retval = FIDARootInit(ida_mem, 0, c_funloc(layerDisCont4IDA))
          tinystep = .true.
        else
          retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4IDA))
          tinystep = .false.
        endif
        if (retval /= 0) then; err=20; message='solveByIDA: error in FIDARootInit'; return; endif
      endif
    endif

  enddo ! while loop on one_step mode until time dt

  !****************************** End of Main Solver ***************************************

  err               = eqns_data%err
  message           = eqns_data%message
  if( .not. feasible) idaSucceeds = .false.

  if(idaSucceeds)then
    ! copy to output data
    diag_data     = eqns_data%diag_data
    flux_data     = eqns_data%flux_data
    deriv_data    = eqns_data%deriv_data
    ixSaturation  = eqns_data%ixSaturation
    dt_out        = tret(1) ! should be dt, probably do not need to keep this
  endif

  ! free memory
  deallocate( eqns_data%sMul )
  deallocate( eqns_data%dMat )
  deallocate( eqns_data%dBaseflow_dMatric )
  deallocate( eqns_data%mLayerMatricHeadLiqTrial )
  deallocate( eqns_data%mLayerMatricHeadTrial )
  deallocate( eqns_data%mLayerMatricHeadPrev )
  deallocate( eqns_data%mLayerVolFracWatTrial )
  deallocate( eqns_data%mLayerVolFracWatPrev )
  deallocate( eqns_data%mLayerVolFracIceTrial )
  deallocate( eqns_data%mLayerTempPrev )
  deallocate( eqns_data%mLayerTempTrial )
  deallocate( eqns_data%mLayerVolFracIcePrev )
  deallocate( eqns_data%mLayerVolFracLiqPrev )
  deallocate( eqns_data%mLayerEnthalpyTrial )
  deallocate( eqns_data%mLayerEnthalpyPrev )
  deallocate( eqns_data%fluxVec )
  deallocate( eqns_data%resSink )
  deallocate( rootsfound )
  deallocate( rootdir )

  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)

end subroutine summaSolveSundialsIDA

! ----------------------------------------------------------------
! SetInitialCondition: routine to initialize u and up vectors.
! ----------------------------------------------------------------
subroutine setInitialCondition(neq, y, sunvec_u, sunvec_up)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fsundials_nvector_mod
  USE fnvector_serial_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)  :: sunvec_u  ! solution N_Vector
  type(N_Vector)  :: sunvec_up ! derivative N_Vector
  integer(c_long) :: neq
  real(rkind)        :: y(neq)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  real(c_double), pointer :: up(:)

  ! get data arrays from SUNDIALS vectors
  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
  up(1:neq) => FN_VGetArrayPointer(sunvec_up)

  uu = y
  up = 0._rkind

end subroutine setInitialCondition

! ----------------------------------------------------------------
! setSolverParams: private routine to set parameters in ida solver
! ----------------------------------------------------------------
subroutine setSolverParams(dt,nonlin_iter,ida_mem,retval)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fida_mod   ! Fortran interface to IDA

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind),intent(in)      :: dt                 ! time step
  integer,intent(in)          :: nonlin_iter        ! maximum number of nonlinear iterations, default = 4, set in parameters
  type(c_ptr),intent(inout)   :: ida_mem            ! IDA memory
  integer(i4b),intent(out)    :: retval             ! return value


  !======= Internals ============
  real(qp),parameter          :: coef_nonlin = 0.33 ! Coeff. in the nonlinear convergence test, default = 0.33
  integer,parameter           :: max_order = 5      ! maximum BDF order,  default = 5
  integer,parameter           :: acurtest_fail = 50 ! maximum number of error test failures, default = 10
  integer,parameter           :: convtest_fail = 50 ! maximum number of convergence test failures, default = 10
  integer(kind = 8),parameter :: max_step = 999999  ! maximum number of steps,  default = 500
  real(qp),parameter          :: h_init = 0         ! initial stepsize
  real(qp)                    :: h_max              ! maximum stepsize,  default = infinity

  ! Set the maximum BDF order
  retval = FIDASetMaxOrd(ida_mem, max_order)
  if (retval /= 0) return

  ! Set Coeff. in the nonlinear convergence test
  retval = FIDASetNonlinConvCoef(ida_mem, coef_nonlin)
  if (retval /= 0) return

  ! Set maximun number of nonliear iterations
  retval = FIDASetMaxNonlinIters(ida_mem, nonlin_iter)
  if (retval /= 0) return

  !  Set maximum number of convergence test failures
  retval = FIDASetMaxConvFails(ida_mem, convtest_fail)
  if (retval /= 0) return

  !  Set maximum number of error test failures
  retval = FIDASetMaxErrTestFails(ida_mem, acurtest_fail)
  if (retval /= 0) return

  ! Set maximum number of steps
  retval = FIDASetMaxNumSteps(ida_mem, max_step)
  if (retval /= 0) return

  ! Set maximum stepsize
  h_max = dt
  retval = FIDASetMaxStep(ida_mem, h_max)
  if (retval /= 0) return

  ! Set initial stepsize
  retval = FIDASetInitStep(ida_mem, h_init)
  if (retval /= 0) return

end subroutine setSolverParams


! ----------------------------------------------------------------------------------------
! find_rootdir: private routine to determine which direction to look for the root, by
!  determining if the variable is greater or less than the root. Need to do this to prevent
!  bouncing around solution
! ----------------------------------------------------------------------------------------
 subroutine find_rootdir(eqns_data,rootdir)

 !======= Inclusions ===========
 use, intrinsic :: iso_c_binding
 use fsundials_nvector_mod
 use fnvector_serial_mod
 use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
 use globalData,only:integerMissing     ! missing integer
 use var_lookup,only:iLookINDEX         ! named variables for structure elements
 use multiconst,only:Tfreeze            ! freezing point of pure water (K)

 !======= Declarations =========
 implicit none

 ! calling variables
 type(eqnsData),intent(in)  :: eqns_data  ! equations data
 integer(i4b),intent(inout) :: rootdir(:) ! root function directions to search

 ! local variables
 integer(i4b)               :: i,ind     ! indices
 integer(i4b)               :: nState    ! number of states
 integer(i4b)               :: nSnow     ! number of snow layers
 integer(i4b)               :: nSoil     ! number of soil layers
 real(rkind)                :: xPsi      ! matric head at layer (m)
 real(rkind)                :: TcSoil    ! critical point when soil begins to freeze (K)

 ! get equations data variables
 nState = eqns_data%nState
 nSnow = eqns_data%nSnow
 nSoil = eqns_data%nSoil

 ! initialize
 ind = 0

 ! identify the critical point when vegetation begins to freeze
 if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
   ind = ind+1
   rootdir(ind) = 1
   if(eqns_data%scalarCanopyTempPrev > Tfreeze) rootdir(ind) = -1
 endif

 if(nSnow>0)then
   do i = 1,nSnow
     ! identify the critical point when the snow layer begins to freeze
     if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
       ind = ind+1
       rootdir(ind) = 1
       if(eqns_data%mLayerTempPrev(i) > Tfreeze) rootdir(ind) = -1
     endif
   end do
 endif

 if(nSoil>0)then
   do i = 1,nSoil
     xPsi = eqns_data%mLayerMatricHeadPrev(i)
     ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
     if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
       ind = ind+1
       rootdir(ind) = 1
       if(xPsi > 0._rkind ) rootdir(ind) = -1
     endif
     ! identify the critical point when the soil layer begins to freeze
     if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
       ind = ind+1
       TcSoil = crit_soilT(xPsi)
       rootdir(ind) = 1
       if(eqns_data%mLayerTempPrev(i+nSnow) > TcSoil) rootdir(ind) = -1
     endif
   end do
 endif

 end subroutine find_rootdir


! ----------------------------------------------------------------------------------------
! layerDisCont4IDA: The root function routine to find soil matrix potential = 0,
!  soil temp = critical frozen point, and snow and veg temp = Tfreeze
! ----------------------------------------------------------------------------------------
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------------------------------
 integer(c_int) function layerDisCont4IDA(t, sunvec_u, sunvec_up, gout, user_data) &
      result(ierr) bind(C,name='layerDisCont4IDA')

 !======= Inclusions ===========
 use, intrinsic :: iso_c_binding
 use fsundials_nvector_mod
 use fnvector_serial_mod
 use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
 use globalData,only:integerMissing     ! missing integer
 use var_lookup,only:iLookINDEX         ! named variables for structure elements
 use multiconst,only:Tfreeze            ! freezing point of pure water (K)

 !======= Declarations =========
 implicit none

 ! calling variables
 real(c_double), value      :: t         ! current time
 type(N_Vector)             :: sunvec_u  ! solution N_Vector
 type(N_Vector)             :: sunvec_up ! derivative N_Vector
 real(c_double)             :: gout(999) ! root function values, if (nVeg + nSnow + 2*nSoil)>999, problem
 type(c_ptr),    value      :: user_data ! user-defined data

 ! local variables
 integer(i4b)               :: i,ind     ! indices
 integer(i4b)               :: nState    ! number of states
 integer(i4b)               :: nSnow     ! number of snow layers
 integer(i4b)               :: nSoil     ! number of soil layers
 real(rkind)                :: xPsi      ! matric head at layer (m)
 real(rkind)                :: TcSoil    ! critical point when soil begins to freeze (K)

 ! pointers to data in SUNDIALS vectors
 real(c_double), pointer :: uu(:)
 type(eqnsData), pointer :: eqns_data      ! equations data

 !======= Internals ============
 ! get equations data from user-defined data
 call c_f_pointer(user_data, eqns_data)
 nState = eqns_data%nState
 nSnow = eqns_data%nSnow
 nSoil = eqns_data%nSoil

 ! get data array from SUNDIALS vector
 uu(1:nState) => FN_VGetArrayPointer(sunvec_u)

 ! initialize
 ind = 0

 ! identify the critical point when vegetation begins to freeze
 if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
   ind = ind+1
   gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)) - Tfreeze
 endif

 if(nSnow>0)then
   do i = 1,nSnow
     ! identify the critical point when the snow layer begins to freeze
     if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
       ind = ind+1
       gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)) - Tfreeze
     endif
   end do
 endif

 if(nSoil>0)then
   do i = 1,nSoil
     ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
     if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
       ind = ind+1
       xPsi = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
       gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
     else
       xPsi = eqns_data%prog_data%var(iLookPROG%mLayerMatricHead)%dat(i)
     endif
     ! identify the critical point when the soil layer begins to freeze
     if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
       ind = ind+1
       TcSoil = crit_soilT(xPsi)
       gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)) - TcSoil
     endif
   end do
 endif

 ! return success
 ierr = 0
 return

 end function layerDisCont4IDA

end module summaSolveSundialsIDA_module
