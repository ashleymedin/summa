module type4kinsol

! data types
USE nrtype
USE, intrinsic :: iso_c_binding      

USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions
implicit none

type eqnsData
  real(rkind)                 :: dt                           ! data step
  integer(i4b)                :: nSnow                        ! number of snow layers
  integer(i4b)                :: nSoil                        ! number of soil layers
  integer(i4b)                :: nLayers                      ! total number of layers
  integer(i4b)                :: nState                       ! total number of state variables
  integer(i4b)                :: ixMatrix                     ! form of matrix (dense or banded)
  logical(lgt)                :: firstSubStep                 ! flag to indicate if we are processing the first sub-step
  logical(lgt)                :: firstFluxCall                ! flag to indicate if we are processing the first flux call
  logical(lgt)                :: firstSplitOper               ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                :: computeVegFlux               ! flag to indicate if computing fluxes over vegetation
  logical(lgt)                :: scalarSolution               ! flag to denote if implementing the scalar solution
  type(zLookup)               :: lookup_data                  ! lookup tables
  type(var_i)                 :: type_data                    ! type of vegetation and soil
  type(var_d)                 :: attr_data                    ! spatial attributes
  type(var_dlength)           :: mpar_data                    ! model parameters
  type(var_d)                 :: forc_data                    ! model forcing data
  type(var_dlength)           :: bvar_data                    ! model variables for the local basin
  type(var_dlength)           :: prog_data                    ! prognostic variables for a local HRU
  type(var_ilength)           :: indx_data                    ! indices defining model states and layers
  type(var_dlength)           :: diag_data                    ! diagnostic variables for a local HRU
  type(var_dlength)           :: flux_data                    ! model fluxes for a local HRU
  type(var_dlength)           :: deriv_data                   ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind), allocatable             :: stateVecPrev(:)        ! state vector from the previous iteration to help with infeasibility
  real(qp), allocatable       :: sMul(:)                      ! state vector multiplier (used in the residual calculations)
  real(rkind), allocatable    :: dMat(:)                      ! diagonal of the Jacobian matrix
  real(rkind), allocatable    :: fluxVec(:)                   ! flux vector
  real(qp), allocatable       :: resSink(:)                   ! additional (sink) terms on the RHS of the state equation
  real(rkind), allocatable    :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)
  integer(i4b)                :: ixSaturation                 ! index of the lowest saturated layer
  integer(i4b)                :: err                          ! error code
  character(len = 50)         :: message                      ! error message
  logical(lgt)                         :: feasible               ! flag to denote the feasibility of the solution
  real(rkind)                          :: fEval                  ! function evaluation
  real(rkind),allocatable              :: fScale(:)              ! function scaling vector
  logical(lgt)                         :: firstStateiteration   ! flag to denote if we computed an iteration so we know to save the state
end type eqnsData


end module type4kinsol