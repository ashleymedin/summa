module mizuroute_coupling

  USE nr_type, only: i4b, rkind
  USE summa_type, only:summa1_type_dec

  ! mizuRoute public interface
  use public_var, only: integerMissing

  implicit none
  private

  ! SUMMA physical constants, restated here because this module must not see the SUMMA module
  ! tree (multiconst sits next to globalData and var_lookup, whose names mizuRoute shares)
  real(rkind), parameter :: iden_water = 1000._rkind   ! intrinsic density of liquid water (kg m-3), as multiconst
  real(rkind), parameter :: Cp_water   = 4181._rkind   ! specific heat of liquid water (J kg-1 K-1), as multiconst
  real(rkind), parameter :: Tfreeze    = 273.16_rkind  ! freezing point of pure water (K), as multiconst

  public :: init_mizuroute_from_summa
  public :: route_mizuroute_from_summa
  public :: define_mizuroute_output_from_summa
  public :: write_mizuroute_output_from_summa
  public :: get_mizuroute_streamflow
  public :: init_stream_network_from_summa
  public :: get_mizuroute_reach_hydraulics
  public :: remap_lateral_energy

  ! *****************************************************************************
  ! SUMMA--mizuRoute coupling interface
  ! *****************************************************************************
  !
  ! This module provides the thin interface between SUMMA and mizuRoute.
  ! It translates between the SUMMA host-model interface and explicitly named
  ! mizuRoute interfaces, but should not contain implementation logic belonging
  ! to either model.
  !
  ! IMPORTANT:
  ! Do not USE modules with names shared by SUMMA and mizuRoute in this module
  ! (e.g., globalData, var_lookup, read_param_module, popMetadat_module).
  ! Both models define modules with these names, and exposing them in the
  ! coupling layer can result in ambiguous or incorrect module resolution.
  !
  ! Native SUMMA logic should remain in SUMMA routines, and native mizuRoute
  ! logic should remain in routines compiled as part of the MIZUROUTE target.
  ! This coupling module should only:
  !
  !   - access the SUMMA structure exposed at the coupling interface;
  !   - transfer data between SUMMA and mizuRoute structures; and
  !   - call uniquely named mizuRoute initialization, routing, and I/O routines.
  !
  ! Keeping this interface thin isolates the SUMMA and mizuRoute module
  ! namespaces and avoids the need to rename shared module names in either model.
  ! *****************************************************************************

  ! *****************************************************************************
  ! Spatial organization of runoff passed from SUMMA to mizuRoute
  ! *****************************************************************************
  !
  ! SUMMA uses a hierarchical spatial structure in which one or more HRUs are
  ! contained within each GRU:
  !
  !       SUMMA HRU_1 ----\
  !       SUMMA HRU_2 -----+--> SUMMA GRU
  !             ...       /
  !       SUMMA HRU_n ----/
  !
  ! SUMMA can optionally simulate lateral flow among HRUs within a GRU. Runoff
  ! from the HRUs is then aggregated to the GRU level. SUMMA subsequently
  ! applies its runoff-routing calculation, which represents the time delay
  ! associated with routing through the unresolved river network.
  !
  ! The current coupling therefore supplies one routed runoff value per
  ! SUMMA GRU:
  !
  !       SUMMA HRUs
  !            |
  !            | optional lateral flow among HRUs
  !            v
  !       aggregate runoff to SUMMA GRU
  !            |
  !            | SUMMA time-delay routing through
  !            | the unresolved river network
  !            v
  !       SUMMA GRU runoff
  !            |
  !            | ------------------
  !            | coupling interface
  !            | ------------------
  !            v
  !       mizuRoute host-model runoff
  !            |
  !            | optional spatial remapping
  !            v
  !       mizuRoute river-network HRU runoff
  !            |
  !            | aggregate runoff to reaches
  !            v
  !       lateral inflow to river reaches
  !            |
  !            | route through explicit river network
  !            v
  !       routed streamflow
  !
  !
  ! Coupling interface
  ! ------------------
  !
  ! The coupling interface contains the runoff values and corresponding SUMMA
  ! GRU identifiers:
  !
  !       coupling(:)%id      SUMMA GRU IDs
  !       coupling(:)%qsim    routed SUMMA GRU runoff [m s-1]
  !
  ! These are transferred to the native mizuRoute runoff structure:
  !
  !       runoff%hru_id(:)    IDs of the host-model runoff elements
  !       runoff%sim(:)       runoff on those elements [m s-1]
  !
  ! Thus:
  !
  !       coupling(:)%id   --> runoff%hru_id(:)
  !       coupling(:)%qsim --> runoff%sim(:)
  !
  ! The name runoff%hru_id can be confusing in the coupled configuration.
  ! It does not imply that these IDs are the mizuRoute river-network HRU IDs.
  ! Rather, runoff%hru_id identifies the spatial elements on which runoff is
  ! supplied by the host land model. For the current SUMMA coupling these
  ! elements are GRUs, so runoff%hru_id contains SUMMA GRU IDs.
  !
  !
  ! Spatial remapping
  ! -----------------
  !
  ! mizuRoute may receive runoff on a spatial domain that differs from the
  ! HRUs associated with its river network. When this occurs, runoff is
  ! spatially remapped before river-network routing:
  !
  !       runoff%hru_id / runoff%sim
  !       host-model runoff elements
  !                    |
  !                    | spatial remapping
  !                    v
  !       runoff%basinRunoff
  !       runoff on river-network HRUs
  !
  ! Within the remapping machinery, the host-model runoff elements are called
  ! qHRUs ("runoff HRUs") to distinguish them from the destination
  ! river-network HRUs. The remapping structure therefore contains:
  !
  !       remap%qhru_id(:)    IDs of source runoff elements (qHRUs)
  !       remap%hru_id(:)     IDs of destination river-network HRUs
  !
  ! The qHRU terminology is specific to the remapping relationship.
  ! remap%qhru_id refers back to the host-model IDs stored in runoff%hru_id.
  ! The IDs in remap%qhru_id are matched against runoff%hru_id to locate the
  ! corresponding values in runoff%sim. Those values are then remapped onto
  ! the river-network HRUs identified by remap%hru_id.
  !
  ! For example, if SUMMA supplies runoff from a single GRU with ID 101 that
  ! contributes to multiple mizuRoute river-network HRUs:
  !
  !       runoff%hru_id       = [101]
  !       runoff%sim          = [runoff from SUMMA GRU 101]
  !
  !       remap%qhru_id       = [101, 101, 101, ...]
  !       remap%hru_id        = [river-network HRU IDs ...]
  !
  ! Each occurrence of qhru_id=101 therefore refers to the same source runoff
  ! value stored in runoff%sim(1). The remapping weights distribute that runoff
  ! onto the corresponding river-network HRUs, producing runoff%basinRunoff.
  !
  ! Spatial remapping is optional. If the host land model already supplies
  ! runoff on the mizuRoute river-network HRUs, the remapping step is skipped
  ! and runoff%sim is used directly as runoff%basinRunoff.
  !
  !
  ! River-network routing
  ! ---------------------
  !
  ! After the optional spatial-remapping step, the coupled workflow is:
  !
  !       runoff%sim
  !       runoff supplied by SUMMA
  !              |
  !              | optional spatial remapping
  !              v
  !       runoff%basinRunoff
  !       runoff on river-network HRUs
  !              |
  !              | basin2reach:
  !              | aggregate river-network HRUs to reaches and
  !              | convert runoff depth to lateral reach inflow
  !              v
  !       lateral inflow to river reaches
  !              |
  !              | route_network
  !              v
  !       routed streamflow
  !
  ! Spatial remapping and river-network routing are therefore distinct
  ! operations. Spatial remapping only reconciles the spatial discretization
  ! of runoff supplied by the host land model with the HRUs associated with
  ! the mizuRoute river network. It can be omitted when those spatial
  ! discretizations already coincide.
  !
  ! The coupled implementation does not include the mizuRoute runoff-routing
  ! routine for the unresolved river network. That component is deliberately
  ! excluded from the set of mizuRoute source files compiled and linked into
  ! SUMMA, because the corresponding time-delay routing is always performed by
  ! SUMMA before runoff crosses the coupling interface.
  !
  ! In the coupled configuration, mizuRoute therefore begins with the runoff
  ! supplied by SUMMA and is responsible only for spatial remapping (when
  ! required), aggregation of runoff to river reaches, and routing through the
  ! explicit river network.
  !
  ! Stream temperature exchange
  ! ---------------------------
  !
  ! When SUMMA carries stream domains (the reach water column of a GRU), a
  ! second exchange runs after the routing step. mizuRoute hands back, per
  ! reach, the discharge leaving the reach, the water volume and the lateral
  ! inflow it received (get_mizuroute_reach_hydraulics), and the energy flux
  ! carried by the GRU runoff (coupling(:)%esim) is remapped onto the reaches
  ! exactly like the runoff itself (remap_lateral_energy). SUMMA then walks
  ! the reaches in routing order and solves the water column of each stream
  ! domain, so temperature is routed downstream by SUMMA while water is routed
  ! by mizuRoute. The per-reach arrays live in summaStruct%stream_net; this
  ! module only fills the mizuRoute side of them.
  !
  ! Possible HRU-level coupling
  ! ---------------------------
  !
  ! Although the current implementation passes runoff at SUMMA GRU resolution,
  ! the coupling interface is not fundamentally restricted to GRUs. SUMMA
  ! runoff could instead be supplied directly at HRU resolution, with SUMMA
  ! HRU IDs stored in runoff%hru_id. In that configuration the SUMMA HRUs would
  ! be the source qHRUs and mizuRoute could perform the required spatial
  ! aggregation and routing through the unresolved river network. HRU-level
  ! coupling is not currently implemented.
  !
  ! *****************************************************************************
contains
  
  !-----------------------------------------------------------------------
  ! Initialize mizuRoute within the SUMMA data structures
  !-----------------------------------------------------------------------
  subroutine init_mizuroute_from_summa(summaStruct, ierr, message)
  USE public_var,     only: iulog
  USE nr_utils,       only: match_index
  USE init_mizuRoute, only: init_mizuroute_domain 

  type(summa1_type_dec), intent(inout) :: summaStruct
  integer,               intent(out)   :: ierr
  character(*),          intent(out)   :: message
  real(rkind)                          :: length_conv
  real(rkind)                          :: time_conv
  integer(i4b)                         :: iGRU
  integer(i4b)                         :: nSpace(1:2) = integerMissing
  integer(i4b)                         :: n_write
  character(len=256)                   :: cmessage
  
  ierr = 0
  message = 'init_mizuroute_from_summa/'
  associate(info => summaStruct%config%mizu_info, domain => summaStruct%mizu_domain)

  ! -----------------------------------------------------------------------
  ! Define host-model information required by mizuRoute
  ! -----------------------------------------------------------------------
 
  ! general info
  info%is_print     = .true.
  info%do_mizuroute = .true.
  info%do_remapping = allocated(info%remap%remap_file)
 
  ! logging
  iulog = summaStruct%config%iulog_summa

  ! time information
  n_write           = summaStruct%n_write
  info%dt_landmodel = summaStruct%data_step
  
  ! SUMMA provides runoff on a one-dimensional basin (GRU) domain 
  nSpace(1) = summaStruct%nGRU_local
  nSpace(2) = integerMissing
  info%is_gridded = (nSpace(2) /= integerMissing)
  
  ! ---- initialize unit conversions (multipliers) ----
  length_conv = 1._rkind   ! no conversion needed: summa runoff length = m
  time_conv   = 1._rkind   ! no conversion needed: summa runoff time = s-1
 
  ! -----------------------------------------------------------------------
  ! Initialize the mizuRoute domain
  !
  ! This performs the major mizuRoute initialization operations, including:
  !   - reading routing configuration and parameter information
  !   - reading the river-network topology
  !   - constructing the river-network data structures
  !   - allocating the runoff and river-routing data structures
  !   - populating the host-model runoff IDs
  !   - reading spatial-remapping information, when required
  !   - constructing the indices required for spatial remapping
  ! -----------------------------------------------------------------------
  call init_mizuroute_domain(summaStruct%instance_parallel%rank, &
                             info, domain, nSpace, n_write,      &
                             summaStruct%coupling(:)%id,         &
                             length_conv, time_conv,             &
                             ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end associate
  
  end subroutine init_mizuroute_from_summa
  
  !-----------------------------------------------------------------------
  ! Network routing in mizuRoute
  !-----------------------------------------------------------------------
  subroutine route_mizuroute_from_summa(modelTimeStep, summaStruct, ierr, message)
  USE network_routing_module, only: network_routing
  integer(i4b),          intent(in)    :: modelTimeStep
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  integer(i4b)       :: iGRU
  integer(i4b)       :: idx_buff
  character(len=256) :: cmessage
  
  ierr    = 0
  message = 'route_mizuroute_from_summa/'
  associate(info => summaStruct%config%mizu_info, domain => summaStruct%mizu_domain)
  
    ! Determine the index of the output buffer (if writePerStep n_write=1)
    idx_buff = merge(1, modelTimeStep, summaStruct%n_write == 1)
  
    ! Transfer summa routed runoff into the mizuRoute runoff structure
    domain%river_network%core%runoff%sim(:) = summaStruct%coupling(:)%qsim

    ! Route the complete runoff field
    call network_routing(idx_buff,             &
                         domain%river_network, &
                         domain%remap%routing, &
                         info%do_remapping,    &
                         ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end associate
  
  end subroutine route_mizuroute_from_summa
  
  !-----------------------------------------------------------------------
  ! Define mizuRoute output based on the SUMMA model structure
  !-----------------------------------------------------------------------
  subroutine define_mizuroute_output_from_summa(ncid, summaStruct, ierr, message)
  USE mizuroute_output_module, only: define_mizuroute_output
  integer(i4b),          intent(in)  :: ncid
  type(summa1_type_dec), intent(in)  :: summaStruct
  integer(i4b),          intent(out) :: ierr
  character(*),          intent(out) :: message
  character(len=256) :: cmessage
  
  ierr = 0
  message = 'define_mizuroute_output_from_summa/'
  call define_mizuroute_output(ncid, summaStruct%config%mizu_info, summaStruct%mizu_domain, ierr, cmessage, &
                               write_stream=any(summaStruct%stream_net%ixDOM > 0))
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine define_mizuroute_output_from_summa

  !-----------------------------------------------------------------------
  ! Write mizuRoute output from the SUMMA model structure
  !-----------------------------------------------------------------------
  subroutine write_mizuroute_output_from_summa(ncid, istart, numtim, summaStruct, ierr, message)
  USE mizuroute_output_module, only: write_mizuroute_output
  integer(i4b),          intent(in)    :: ncid
  integer(i4b),          intent(in)    :: istart
  integer(i4b),          intent(in)    :: numtim
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  character(len=256) :: cmessage

  ierr = 0
  message = 'write_mizuroute_output_from_summa/'
  if(any(summaStruct%stream_net%ixDOM > 0))then
    call write_mizuroute_output(ncid, istart, numtim, summaStruct%config%mizu_info, summaStruct%mizu_domain, ierr, cmessage, &
                                tReach=summaStruct%stream_net%tOutHist, vReach=summaStruct%stream_net%velHist)
  else
    call write_mizuroute_output(ncid, istart, numtim, summaStruct%config%mizu_info, summaStruct%mizu_domain, ierr, cmessage)
  endif
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine write_mizuroute_output_from_summa

  !-----------------------------------------------------------------------
  ! Get mizuRoute streamflow
  !-----------------------------------------------------------------------
  subroutine get_mizuroute_streamflow(modelTimeStep, summaStruct, simFlow)
  integer(i4b),          intent(in)  :: modelTimeStep
  type(summa1_type_dec), intent(in)  :: summaStruct
  real(rkind),           intent(out) :: simFlow
  integer(i4b) :: idx_buff
  integer(i4b) :: ixSeg

  idx_buff = merge(1, modelTimeStep, summaStruct%n_write == 1)
  ixSeg    = summaStruct%config%mizu_info%ntopo%ixSegOut
  simFlow = summaStruct%mizu_domain%river_network%driver%method(1)%streamflow(ixSeg,idx_buff)

  end subroutine get_mizuroute_streamflow

  !-----------------------------------------------------------------------
  ! Build the river network as seen by the stream temperature model
  !-----------------------------------------------------------------------
  ! The SUMMA side supplies, per GRU, the reach id its stream HRU stands for
  ! (0 = the reach the GRU drains to) and where that stream domain lives.
  ! Everything else comes from the mizuRoute topology.
  subroutine init_stream_network_from_summa(summaStruct, streamSegId, ixStreamHRU, ixStreamDOM, domArea, ierr, message)
  USE var_lookup, only: ixNTOPO, ixHRU2SEG
  USE public_var, only: realMissing, iulog
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(in)    :: streamSegId(:)   ! per GRU: reach id of the stream HRU (0 = reach mapped from the GRU id)
  integer(i4b),          intent(in)    :: ixStreamHRU(:)   ! per GRU: index of the stream HRU within the GRU (0 = none)
  integer(i4b),          intent(in)    :: ixStreamDOM(:)   ! per GRU: index of the stream domain within that HRU
  real(rkind),           intent(in)    :: domArea(:)       ! per GRU: planform area of the stream domain (m2)
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  integer(i4b)                         :: iGRU, iSeg, iHRU, jSeg, nUps, maxUps
  real(rkind)                          :: reachArea        ! reach planform area from the routing geometry (m2)

  ierr = 0
  message = 'init_stream_network_from_summa/'
  associate(net => summaStruct%stream_net, core => summaStruct%mizu_domain%river_network%core, &
            driver => summaStruct%mizu_domain%river_network%driver)

  net%nSeg = core%topology%n_seg
  maxUps = 0
  do iSeg=1,net%nSeg
    maxUps = max(maxUps, size(core%ntopo(iSeg)%UREACHI))
  end do
  allocate(net%segId(net%nSeg), net%rchOrder(net%nSeg), net%nUps(net%nSeg), net%ixUps(max(maxUps,1),net%nSeg), &
           net%ixGRU(net%nSeg), net%ixHRU(net%nSeg), net%ixDOM(net%nSeg), net%length(net%nSeg),                 &
           net%qUp(net%nSeg), net%qLat(net%nSeg), net%qOut(net%nSeg), net%vol(net%nSeg), net%depth(net%nSeg),     &
           net%velocity(net%nSeg), net%eLat(net%nSeg), net%tLat(net%nSeg), net%tUp(net%nSeg), net%tOut(net%nSeg),   &
           net%tOutHist(net%nSeg,summaStruct%n_write), net%velHist(net%nSeg,summaStruct%n_write), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating the stream network'; return; endif
  net%tOutHist(:,:) = Tfreeze; net%velHist(:,:) = 0._rkind

  ! topology
  net%ixUps(:,:) = 0
  do iSeg=1,net%nSeg
    net%segId(iSeg)    = core%ntopo(iSeg)%REACHID
    net%rchOrder(iSeg) = core%topology%ntopo(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
    net%length(iSeg)   = core%param(iSeg)%RLENGTH
    nUps = size(core%ntopo(iSeg)%UREACHI)
    net%nUps(iSeg) = nUps
    if(nUps>0) net%ixUps(1:nUps,iSeg) = core%ntopo(iSeg)%UREACHI(1:nUps)
  end do

  ! the stream domain standing for each reach
  net%ixGRU(:) = 0; net%ixHRU(:) = 0; net%ixDOM(:) = 0
  do iGRU=1,summaStruct%nGRU_local
    if(ixStreamHRU(iGRU) < 1) cycle
    if(streamSegId(iGRU) > 0)then
      jSeg = findloc(driver%seg_id, streamSegId(iGRU), dim=1)
      if(jSeg < 1)then
        write(message,'(a,i0,a,i0,a)') trim(message)//'streamSegId ',streamSegId(iGRU),' of the stream HRU in GRU ',summaStruct%gru_struc(iGRU)%gru_id,' is not a reach of the river network'
        ierr=20; return
      endif
    else
      ! the reach the GRU drains to: the routing HRU with the GRU id (only meaningful without spatial remapping)
      iHRU = findloc(driver%hru_id, int(summaStruct%gru_struc(iGRU)%gru_id,kind=i4b), dim=1)
      if(iHRU < 1)then
        write(message,'(a,i0,a)') trim(message)//'GRU ',summaStruct%gru_struc(iGRU)%gru_id,' has a stream HRU but no routing HRU with its id; give the reach with streamSegId'
        ierr=20; return
      endif
      jSeg = core%topology%hru2seg(iHRU)%var(ixHRU2SEG%hruSegIndex)%dat(1)
    endif
    if(net%ixGRU(jSeg) > 0)then
      write(message,'(a,i0,a)') trim(message)//'reach ',net%segId(jSeg),' is represented by more than one stream HRU'
      ierr=20; return
    endif
    net%ixGRU(jSeg) = iGRU
    net%ixHRU(jSeg) = ixStreamHRU(iGRU)
    net%ixDOM(jSeg) = ixStreamDOM(iGRU)
    ! the column is the reach: its area should be the reach planform area, or the residence time is off by the ratio
    reachArea = core%param(jSeg)%RLENGTH*core%param(jSeg)%R_WIDTH
    if(abs(domArea(iGRU) - reachArea) > 0.1_rkind*reachArea) &
      write(iulog,'(a,i0,a,es10.3,a,es10.3,a)') ' WARNING: stream domain of GRU ',summaStruct%gru_struc(iGRU)%gru_id, &
        ' has area ',domArea(iGRU),' m2 but its reach is ',reachArea,' m2 (length x width); the reach volume per unit area uses the reach'
  end do

  ! initial values
  net%qUp(:) = 0._rkind; net%qLat(:) = 0._rkind; net%qOut(:) = 0._rkind
  net%vol(:) = realMissing; net%depth(:) = 0._rkind; net%velocity(:) = 0._rkind
  net%eLat(:) = 0._rkind; net%tLat(:) = Tfreeze; net%tUp(:) = Tfreeze; net%tOut(:) = Tfreeze

  end associate
  end subroutine init_stream_network_from_summa

  !-----------------------------------------------------------------------
  ! Reach hydraulics after a routing step: discharge, lateral inflow, volume, depth, velocity
  !-----------------------------------------------------------------------
  ! Uses the first configured routing method. The volume is only carried by
  ! the IRF, KW, MC and DW methods; for the others the depth falls back to
  ! the Manning normal depth. The depth passed to SUMMA is the reach volume
  ! spread over the reach planform area (length x width of the routing
  ! geometry), so the SUMMA column has the mean depth of the reach.
  subroutine get_mizuroute_reach_hydraulics(summaStruct, ierr, message)
  USE public_var, only: realMissing
  USE hydraulic,  only: flow_depth, flow_area
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  integer(i4b)                         :: iSeg
  real(rkind)                          :: yNorm, aFlow, reachArea
  real(rkind), parameter               :: minVol = 1.e-6_rkind

  ierr = 0
  message = 'get_mizuroute_reach_hydraulics/'
  associate(net => summaStruct%stream_net, core => summaStruct%mizu_domain%river_network%core)

  do iSeg=1,net%nSeg
    net%qOut(iSeg) = core%flux(iSeg)%ROUTE(1)%REACH_Q
    net%qUp(iSeg)  = core%flux(iSeg)%ROUTE(1)%REACH_INFLOW
    net%qLat(iSeg) = core%flux(iSeg)%BASIN_QR(1)
    net%vol(iSeg)  = core%flux(iSeg)%ROUTE(1)%REACH_VOL(1)
    reachArea = net%length(iSeg)*core%param(iSeg)%R_WIDTH
    ! Manning normal depth for the discharge: the depth where the routing method carries no volume, and the
    ! floor for the velocity (a reach whose volume has not caught up with its discharge is not moving at Q L/V)
    yNorm = 0._rkind
    if(net%qOut(iSeg) > 0._rkind) &
      yNorm = flow_depth(net%qOut(iSeg), core%param(iSeg)%R_WIDTH, core%param(iSeg)%SIDE_SLOPE, &
                         core%param(iSeg)%R_SLOPE, core%param(iSeg)%R_MAN_N, &
                         zf=core%param(iSeg)%FLDP_SLOPE, bankDepth=core%param(iSeg)%R_DEPTH)
    if(net%vol(iSeg) > minVol .and. reachArea > 0._rkind)then
      net%depth(iSeg) = net%vol(iSeg)/reachArea
    else
      net%vol(iSeg)   = realMissing
      net%depth(iSeg) = yNorm
    endif
    if(net%qOut(iSeg) > 0._rkind)then
      aFlow = flow_area(max(net%depth(iSeg),yNorm), core%param(iSeg)%R_WIDTH, core%param(iSeg)%SIDE_SLOPE, &
                        zf=core%param(iSeg)%FLDP_SLOPE, bankDepth=core%param(iSeg)%R_DEPTH)
      net%velocity(iSeg) = net%qOut(iSeg)/max(aFlow, minVol)
    else
      net%velocity(iSeg) = 0._rkind
    endif
  end do

  end associate
  end subroutine get_mizuroute_reach_hydraulics

  !-----------------------------------------------------------------------
  ! Remap the energy flux carried by the GRU runoff onto the reaches
  !-----------------------------------------------------------------------
  ! The energy flux density (W m-2) goes through the same spatial remapping
  ! and basin-to-reach aggregation as the runoff depth, so dividing by the
  ! reach lateral inflow gives its flow-weighted temperature.
  subroutine remap_lateral_energy(summaStruct, ierr, message)
  USE process_remap_module, only: remap_runoff
  USE process_remap_module, only: basin2reach
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer(i4b),          intent(out)   :: ierr
  character(*),          intent(out)   :: message
  real(rkind), allocatable             :: qSave(:)         ! the runoff, put back after the remapping
  real(rkind), allocatable             :: basinNrg(:)      ! energy flux density on the routing HRUs (W m-2)
  real(rkind), parameter               :: rhoCp = iden_water*Cp_water ! rho_w * Cp_w (J m-3 K-1)
  real(rkind), parameter               :: minFlow = 1.e-12_rkind
  integer(i4b)                         :: iSeg
  character(len=256)                   :: cmessage

  ierr = 0
  message = 'remap_lateral_energy/'
  associate(net => summaStruct%stream_net, info => summaStruct%config%mizu_info, &
            core => summaStruct%mizu_domain%river_network%core, domain => summaStruct%mizu_domain)

  allocate(basinNrg(size(core%runoff%basinRunoff)), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating basinNrg'; return; endif

  ! the energy flux density takes the place of the runoff for the remapping, then the runoff is restored
  qSave = core%runoff%sim
  core%runoff%sim(:) = summaStruct%coupling(:)%esim
  if(info%do_remapping)then
    call remap_runoff(core%runoff, domain%remap%routing, basinNrg, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
    basinNrg = core%runoff%sim
  endif
  core%runoff%sim = qSave

  ! aggregate to the reaches (W), no lower limit since this is not a runoff
  call basin2reach(basinNrg, core%ntopo, core%param, net%eLat, ierr, cmessage, limitRunoff=.false.)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  do iSeg=1,net%nSeg
    if(net%qLat(iSeg) > minFlow)then
      net%tLat(iSeg) = max(net%eLat(iSeg)/(rhoCp*net%qLat(iSeg)), Tfreeze)
    else
      net%tLat(iSeg) = Tfreeze
    endif
  end do

  end associate
  end subroutine remap_lateral_energy

end module mizuroute_coupling
