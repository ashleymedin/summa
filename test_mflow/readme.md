# SUMMA MODFLOW case study
This folder contains a case study to show how a typical SUMMA setup looks coupling with MODFLOW.

  ! ****************************************************************************************
  ! *** Thin BMI coupler: SUMMA land model  <-->  MODFLOW 6 groundwater model          ***
  ! ****************************************************************************************
  !
  ! Each SUMMA data step:
  !   1. (feedback) the MODFLOW 6 water-table head from the previous step is written into
  !      SUMMA as the prescribed-head lower boundary condition of the soil column
  !      (BMI input  "soil_water_sat-zone_top__head", parameter "lowerBoundHead").
  !   2. SUMMA advances one step.
  !   3. the drainage out the base of the SUMMA soil column
  !      (BMI output "soil_water__drainage_volume_flux", flux "scalarSoilDrainage")
  !      is regridded onto the MODFLOW 6 grid and written into the RCH package
  !      RECHARGE array.
  !   4. MODFLOW 6 advances one step (prepare/do/finalize_time_step).
  !   5. the new MODFLOW 6 head field is read back and aggregated per SUMMA HRU,
  !      ready to be applied at step 1 of the next iteration (explicit, one-step lag).
  !
  ! Coupling requires SUMMA to be built with -DUSE_MODFLOW6=ON (sets MODFLOW_ACTIVE) and
  ! the SUMMA model decision  groundwatr = modflow  with  bcLowrSoiH = presHead.
  !
  ! The MODFLOW 6 model (read from mfsim.nam in the working directory) must:
  !   * use length unit metres and TDIS TIME_UNITS SECONDS,
  !   * have exactly one MODFLOW time step per SUMMA forcing data step: MODFLOW delt
  !     must equal SUMMA's data_step (taken from the forcing file's data_step
  !     attribute, constant for the run; the coupler reads it at runtime and checks
  !     it against MODFLOW delt).
  !   * contain an RCH package with READASARRAYS,
  !   * be a single GWF model discretised with DIS.
  !
  ! Configuration (Fortran namelist), default file name "summa_modflow6.config":
  !
  !   &coupler
  !     mf6_model_name     = 'MYMODEL' ! GWF model name, as in mfsim.nam (upper case)
  !     rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
  !     bflow_package_name = 'CHD'     ! head-dependent boundary package (CHD/DRN/RIV/GHB) whose
  !                                    !    simulated flow feeds back per HRU as scalarAquiferBaseflow
  !                                    !    ('' => skip the baseflow feedback)
  !     map_file           = ''        ! optional HRU->cell weight file; if blank a nearest-cell
  !                                    !    map is built from the MODFLOW 6 DIS grid geometry
  !     mf6_epsg           = 0         ! EPSG code of the MODFLOW grid's projected CRS, used only to
  !                                    !    reproject SUMMA's lon/lat HRU centres before the built-in
  !                                    !    nearest-cell map (nHRU>1; see nearest_hru).  Only WGS84 UTM
  !                                    !    is supported (32601-32660 N, 32701-32760 S); 0 = no
  !                                    !    reprojection (fine when nHRU==1, or with an explicit map_file)
  !     feedback           = .true.    ! .false. => one-way (SUMMA drainage -> MODFLOW only)
  !   /
  !
  ! With feedback = .true. the coupler writes two groundwater quantities back into SUMMA
  ! each step so its water balance and routed streamflow include the aquifer:
  !     scalarAquiferStorage  = Sy * (MODFLOW water table - soil-column base)
  !     scalarAquiferBaseflow = MODFLOW <bflow_package_name> outflow over the HRU footprint
  ! (scalarAquiferRecharge is not exchanged - SUMMA sets it to its own soil drainage.)
  !
  ! Usage:  summa_modflow6.exe <fileManager.txt> [summa_modflow6.config]