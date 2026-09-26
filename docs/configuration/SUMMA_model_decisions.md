# Model Decisions in SUMMA

Model decisions select which process parameterizations and numerical methods SUMMA uses for a
run. They are supplied in the [model decisions file](../input_output/SUMMA_input.md#infile_model_decisions),
whose path is given by the `decisionsFile` entry of the [file manager](../input_output/SUMMA_input.md#infile_master_configuration).

* The file has one decision per line as `<keyword> <option>`; lines may appear in any order.
* Keywords and options are **case sensitive**.
* The authoritative list of keywords is `build/source/dshare/get_ixname.f90`
  (`function get_ixdecisions`); the authoritative list of options, defaults and the rules that
  couple decisions together is `build/source/engine/mDecisions.f90`. If this page and the code
  disagree, the code wins.
* Simulation start/end time (`simStartTime`, `simEndTime`) and time-zone handling
  (`tmZoneInfo`) are **not** model decisions — they moved to the [file manager](../input_output/SUMMA_input.md#infile_master_configuration).
* Every decision must be present in the file. Several options are named `notPopulatedYet`;
  where accepted (noted below) they resolve to the default so that inactive decisions can be
  given a placeholder value.

## Summary

There are 44 model decisions. Defaults (used for `notPopulatedYet` where accepted) are in **bold**.

| # | Keyword | Options | Purpose |
|---|---------|---------|---------|
| 1 | [soilCatTbl](#soilcattbl) | STAS, STAS-RUC, ROSETTA | soil-category lookup table |
| 2 | [vegeParTbl](#vegepartbl) | USGS, MODIFIED_IGBP_MODIS_NOAH, plumberCABLE, plumberCHTESSEL, plumberSUMMA | vegetation-category lookup table |
| 3 | [soilStress](#soilstress) | NoahType, CLM_Type, SiB_Type | soil-moisture control on stomatal resistance |
| 4 | [stomResist](#stomresist) | BallBerry, Jarvis, simpleResistance, BallBerryFlex, BallBerryTest | stomatal resistance |
| 5 | [bbTempFunc](#bbtempfunc) | q10Func, Arrhenius | Ball-Berry: leaf-temperature control on photosynthesis |
| 6 | [bbHumdFunc](#bbhumdfunc) | humidLeafSurface, scaledHyperbolic | Ball-Berry: humidity control on stomatal resistance |
| 7 | [bbElecFunc](#bbelecfunc) | linear, linearJmax, quadraticJmax | Ball-Berry: dependence of photosynthesis on PAR |
| 8 | [bbCO2point](#bbco2point) | origBWB, Leuning | Ball-Berry: use of the CO2 compensation point |
| 9 | [bbNumerics](#bbnumerics) | NoahMPsolution, newtonRaphson | Ball-Berry: iterative solution method |
| 10 | [bbAssimFnc](#bbassimfnc) | colimitation, minFunc | Ball-Berry: controls on carbon assimilation |
| 11 | [bbCanIntg8](#bbcanintg8) | constantScaling, laiScaling | Ball-Berry: leaf-to-canopy scaling of photosynthesis |
| 12 | [num_method](#num_method) | **homegrown** (alias `itertive`), kinsol, ida | numerical method for the model equations |
| 13 | [fDerivMeth](#fderivmeth) | numericl, analytic | flux derivatives for the Jacobian |
| 14 | [LAI_method](#lai_method) | monTable, specified | source of LAI/SAI |
| 15 | [cIntercept](#cintercept) | **notPopulatedYet**, sparseCanopy, storageFunc | canopy interception |
| 16 | [f_Richards](#f_richards) | **mixdform** | form of Richards' equation |
| 17 | [groundwatr](#groundwatr) | qTopmodl, bigBuckt, noXplict, modflow, modLatflow | groundwater parameterization |
| 18 | [hc_profile](#hc_profile) | constant, pow_prof, exp_prof | hydraulic-conductivity profile with depth |
| 19 | [bcUpprTdyn](#bcupprtdyn) | presTemp, nrg_flux, zeroFlux | upper boundary condition, thermodynamics |
| 20 | [bcLowrTdyn](#bclowrtdyn) | presTemp, zeroFlux, presFlux | lower boundary condition, thermodynamics |
| 21 | [bcUpprSoiH](#bcupprsoih) | presHead, liq_flux | upper boundary condition, soil hydrology |
| 22 | [bcLowrSoiH](#bclowrsoih) | presHead, bottmPsi, drainage, zeroFlux | lower boundary condition, soil hydrology |
| 23 | [veg_traits](#veg_traits) | Raupach_BLM1994, CM_QJRMS1988, vegTypeTable | vegetation roughness length and displacement height |
| 24 | [rootProfil](#rootprofil) | **powerLaw**, doubleExp | rooting profile |
| 25 | [canopyEmis](#canopyemis) | simplExp, difTrans | canopy emissivity |
| 26 | [snowIncept](#snowincept) | stickySnow, lightSnow | maximum canopy snow interception capacity |
| 27 | [windPrfile](#windprfile) | exponential, logBelowCanopy | canopy wind profile |
| 28 | [astability](#astability) | standard, louisinv, mahrtexp | atmospheric stability function |
| 29 | [compaction](#compaction) | consettl, anderson | snow compaction / densification |
| 30 | [snowLayers](#snowlayers) | jrdn1991, CLM_2010 | snow layer combination / sub-division rules |
| 31 | [thCondSnow](#thcondsnow) | tyen1965, melr1977, jrdn1991, smnv2000 | snow thermal conductivity |
| 32 | [thCondSoil](#thcondsoil) | funcSoilWet, mixConstit, hanssonVZJ | soil thermal conductivity |
| 33 | [canopySrad](#canopysrad) | noah_mp, CLM_2stream, UEB_2stream, NL_scatter, BeersLaw | canopy shortwave radiation |
| 34 | [alb_method](#alb_method) | conDecay, varDecay | snow albedo |
| 35 | [spatial_gw](#spatial_gw) | localColumn, singleBasin | spatial representation of groundwater |
| 36 | [subRouting](#subrouting) | timeDlay, qInstant | within-basin (sub-grid) routing |
| 37 | [snowDenNew](#snowdennew) | **hedAndPom**, anderson, pahaut_76, constDens | new-snow density |
| 38 | [snowUnload](#snowunload) | **meltDripUnload**, windUnload | unloading of snow from the canopy |
| 39 | [nrgConserv](#nrgconserv) | closedForm, enthalpyForm, enthalpyFormAN | state variable / residual form for the energy equations |
| 40 | [aquiferIni](#aquiferini) | **fullStart**, emptyStart | initial fill level of the aquifer |
| 41 | [infRateMax](#infratemax) | **topmodel_GA**, GreenAmpt, noInfExc | maximum infiltration rate |
| 42 | [surfRun_SE](#surfrun_se) | **homegrown_SE**, FUSEPRMS, FUSEAVIC, FUSETOPM, zero_SE | saturation-excess surface runoff |
| 43 | [read_force](#read_force) | **readPerStep**, readFullSeries | how forcing data are read |
| 44 | [write_buff](#write_buff) | **writePerStep**, writeFullSeries | how model output is buffered before writing |
| 45 | [deepTherml](#deeptherml) | **none**, aquiferTemp, airTempGW, bedrockLyrs | deep thermal state below the hydrologically active soil column |
| 46 | [hyporhTdyn](#hyporhtdyn) | **none**, proxy | hyporheic exchange in a stream domain |

---

<a id="soilcattbl"></a>
## 1. soilCatTbl — soil-category dataset

Selects the block of `SOILPARM.TBL` used for soil properties keyed on `soilTypeIndex`.

| Option | Description |
|---|---|
| STAS | STATSGO soil classes (Noah/Noah-MP default) |
| STAS-RUC | STATSGO classes with the RUC land-surface-model parameter set |
| ROSETTA | classes/parameters from the ROSETTA pedotransfer database |

<a id="vegepartbl"></a>
## 2. vegeParTbl — vegetation-category dataset

Selects the block of `VEGPARM.TBL` used for vegetation properties keyed on `vegTypeIndex`.

| Option | Description |
|---|---|
| USGS | 24-class USGS land-use classification |
| MODIFIED_IGBP_MODIS_NOAH | 20-class modified IGBP/MODIS classification |
| plumberCABLE | vegetation parameters tuned for the PLUMBER CABLE configuration |
| plumberCHTESSEL | vegetation parameters tuned for the PLUMBER CH-TESSEL configuration |
| plumberSUMMA | vegetation parameters tuned for the PLUMBER SUMMA configuration |

<a id="soilstress"></a>
## 3. soilStress — soil-moisture control on stomatal resistance

| Option | Description |
|---|---|
| NoahType | thresholded linear function of volumetric liquid water content |
| CLM_Type | thresholded linear function of matric head ([CLM technical note, 2010](https://doi.org/10.1029/2011MS00045)) |
| SiB_Type | exponential of the log of matric head |

<a id="stomresist"></a>
## 4. stomResist — stomatal resistance

| Option | Description |
|---|---|
| BallBerry | Ball-Berry conductance model (as in Noah-MP) |
| Jarvis | Jarvis environmental-stress model |
| simpleResistance | prescribed minimum stomatal resistance scaled by environmental factors |
| BallBerryFlex | flexible Ball-Berry scheme; enables decisions 5–11 |
| BallBerryTest | flexible Ball-Berry scheme, testing configuration; enables decisions 5–11 |

Decisions 5–11 are read only when `stomResist` is `BallBerryFlex` or `BallBerryTest`;
otherwise their values are ignored (use `notPopulatedYet`).

<a id="bbtempfunc"></a>
## 5. bbTempFunc — Ball-Berry leaf-temperature function

| Option | Description |
|---|---|
| q10Func | Q10 function, as in CLM4 and Noah-MP |
| Arrhenius | Arrhenius functions, as in CLM5 and CABLE |

<a id="bbhumdfunc"></a>
## 6. bbHumdFunc — Ball-Berry humidity function

| Option | Description |
|---|---|
| humidLeafSurface | humidity at the leaf surface (Bonan et al., 2011) |
| scaledHyperbolic | scaled hyperbolic function (Leuning et al., 1995) |

<a id="bbelecfunc"></a>
## 7. bbElecFunc — Ball-Berry electron-transport (PAR) function

| Option | Description |
|---|---|
| linear | linear function, as in CLM4 and Noah-MP |
| linearJmax | linear Jmax function, as in CABLE |
| quadraticJmax | quadratic Jmax function, as in SSiB and CLM5 |

<a id="bbco2point"></a>
## 8. bbCO2point — Ball-Berry CO2 compensation point

| Option | Description |
|---|---|
| origBWB | original Ball-Woodrow-Berry formulation (no explicit compensation point) |
| Leuning | Leuning formulation using the CO2 compensation point |

<a id="bbnumerics"></a>
## 9. bbNumerics — Ball-Berry iterative solution method

| Option | Description |
|---|---|
| NoahMPsolution | fixed-point iteration, maximum 3 iterations (as in Noah-MP and CLM4) |
| newtonRaphson | full Newton-Raphson iteration to convergence |

<a id="bbassimfnc"></a>
## 10. bbAssimFnc — Ball-Berry carbon assimilation

| Option | Description |
|---|---|
| colimitation | smooth co-limitation of the three assimilation controls (Collatz et al. 1991; Sellers et al. 1996) |
| minFunc | take the minimum of the three assimilation controls |

<a id="bbcanintg8"></a>
## 11. bbCanIntg8 — Ball-Berry leaf-to-canopy scaling

| Option | Description |
|---|---|
| constantScaling | constant scaling factor from leaf to canopy |
| laiScaling | exponential function of LAI (Leuning et al., 1995, eq. 9) |

<a id="num_method"></a>
## 12. num_method — numerical method

| Option | Description |
|---|---|
| homegrown | SUMMA's built-in backward-Euler solver (concepts from Numerical Recipes); constant sub-step with adaptive retries |
| itertive | accepted as a backward-compatible alias for `homegrown` |
| kinsol | backward Euler with the SUNDIALS KINSOL nonlinear solver; constant step size. Requires SUMMA built with `-DUSE_SUNDIALS=ON`. [KINSOL docs](https://sundials.readthedocs.io/en/latest/kinsol/) |
| ida | adaptive-step implicit differential-algebraic solution with SUNDIALS IDA. Requires SUMMA built with `-DUSE_SUNDIALS=ON`. [IDA docs](https://sundials.readthedocs.io/en/latest/ida/) |

<a id="fderivmeth"></a>
## 13. fDerivMeth — flux derivatives for the Jacobian

| Option | Description |
|---|---|
| numericl | finite-difference (numerical) derivatives |
| analytic | analytical derivatives. Required for `num_method = kinsol` or `ida`; with `homegrown` either option may be used |

<a id="lai_method"></a>
## 14. LAI_method — source of LAI and SAI

| Option | Description |
|---|---|
| monTable | LAI/SAI taken directly from the monthly vegetation-class table |
| specified | LAI/SAI computed from the green-vegetation fraction and the `winterSAI` / `summerLAI` parameters |

<a id="cintercept"></a>
## 15. cIntercept — canopy interception

| Option | Description |
|---|---|
| notPopulatedYet | undefined (backward compatibility) |
| sparseCanopy | a fixed fraction of rainfall reaches the ground as throughfall; canopy drainage above a storage threshold |
| storageFunc | throughfall is a function of relative canopy storage; 100% throughfall at canopy capacity |

<a id="f_richards"></a>
## 16. f_Richards — form of Richards' equation

| Option | Description |
|---|---|
| mixdform | mixed (head/moisture) form of Richards' equation |

The moisture-based form was removed; `mixdform` (or `notPopulatedYet`) is the only accepted value.

<a id="groundwatr"></a>
## 17. groundwatr — groundwater parameterization

| Option | Description |
|---|---|
| qTopmodl | TOPMODEL-style baseflow from a per-column store |
| bigBuckt | lumped "big bucket" aquifer model |
| noXplict | no explicit groundwater; soil drainage leaves the column |
| modflow | couple with MODFLOW 6 (no aquifer modelled in SUMMA) |
| modLatflow | couple with MODFLOW 6 but allow lateral flow in the soil column (no aquifer modelled in SUMMA) |

See also [`spatial_gw`](#spatial_gw), which sets whether the store is per column or per basin.

With `modflow`, each data step SUMMA sends the soil-column drainage to the MODFLOW 6
recharge package and receives the MODFLOW water-table head back as the soil-column lower
boundary condition, so it **requires** [`bcLowrSoiH`](#bclowrsoih) `= presHead`. It is
only available in a build configured with `-DUSE_MODFLOW6=ON` and must be run through the
`summa_modflow6` coupler executable (which also drives MODFLOW 6); a plain `summa` run
with `groundwatr = modflow` is rejected at start-up. MODFLOW then supplies
`scalarAquiferBaseflow` (routed into streamflow) and `scalarAquiferStorage`, and
`scalarAquiferRecharge` is set equal to the soil-column drainage.

`modLatflow` is `modflow` plus lateral flow through the soil column itself, for hillslopes
where water moves downslope through the soil as well as recharging the aquifer. It shares
every requirement of `modflow` above, and additionally **requires**
[`hc_profile`](#hc_profile) `= exp_prof`: the lateral transmissivity is the vertical
integral of the conductivity over the soil column alone, since MODFLOW carries everything
below it, and the exponential profile is the one that integrates to a finite base rather
than assuming a shallow aquifer of its own. That in turn requires
[`infRateMax`](#infratemax) `= topmodel_GA` (or `noInfExc`). The lateral flow appears as
`basin__ColumnOutflow` and is added to total runoff alongside the MODFLOW baseflow.

<a id="hc_profile"></a>
## 18. hc_profile — hydraulic-conductivity profile

| Option | Description |
|---|---|
| constant | saturated hydraulic conductivity constant with depth |
| pow_prof | power-law decrease with depth, normalized to `k_soil` at `compactedDepth` and reaching zero at the base of the soil |
| exp_prof | exponential decrease with depth, `K(z) = k_soil * exp(-f_hydCond * z)`, finite at the base of the soil |

`constant` means a homogeneous column, so it pairs with [`infRateMax`](#infratemax) `GreenAmpt`. `pow_prof` and `exp_prof` vary with
depth and so require `topmodel_GA` (or `noInfExc`), which evaluates the conductivity at the wetting front.

[`groundwatr`](#groundwatr) `qTopmodl` requires `pow_prof`, the only option whose transmissivity represents a shallow aquifer.
`exp_prof` instead stops at an impermeable base, so it is not selectable with `qTopmodl`; it is the profile for lateral flow with
no shallow aquifer of its own, which is glacier debris over ice (used internally whatever this decision is set to) or soil over an
external aquifer, and it is accordingly required by `modLatflow`.

`pow_prof` cannot be used with [`bcLowrSoiH`](#bclowrsoih) `presHead`: its conductivity is exactly zero at the base of the soil, so a
prescribed head can drive no drainage. That rules it out for the MODFLOW-coupled options, which require `presHead`. Use `exp_prof`, which decays with depth but stays finite there.

The decay rate for `exp_prof` is the parameter `f_hydCond` (m-1), default 0.75, a soil-column value leaving ~5% of the surface
conductivity at 4 m. Supraglacial debris is 1-5 m-1 over a 0.3-1 m depth, so glacier runs should set it explicitly.

<a id="bcupprtdyn"></a>
## 19. bcUpprTdyn — upper boundary condition, thermodynamics

| Option | Description |
|---|---|
| presTemp | prescribed surface temperature |
| nrg_flux | energy flux computed from the surface energy balance |
| zeroFlux | zero energy flux at the upper boundary |

<a id="bclowrtdyn"></a>
## 20. bcLowrTdyn — lower boundary condition, thermodynamics

| Option | Description |
|---|---|
| presTemp | prescribed temperature at the bottom of the soil column |
| zeroFlux | zero energy flux at the bottom of the soil column |
| presFlux | prescribed energy flux into the bottom of the soil column, the geothermal heat flux |

`presTemp` holds the base of the column at `lowerBoundTemp`, which is only physical if the
column reaches the depth where the annual temperature signal is damped out. `zeroFlux`
insulates the base instead, so a shallow column has no thermal memory below the seasonal
cycle and its bottom temperature tracks the surface with a lag.

`presFlux` prescribes the energy flux `lowerBoundNrgFlux` (W m-2, default 0.06, the continental
average geothermal heat flux; roughly 0.03 in old cratons and 0.1 in tectonically active
regions) entering the base of the column. Like `presTemp` it is a prescribed constant, not a
computed or time-varying flux. It is the right lower boundary for a deep column in permafrost
or cold-region work: the long-term temperature gradient at the base is
`lowerBoundNrgFlux`/`thCond_soil` rather than zero, so the deep column keeps a realistic
temperature and the water draining from it is not simply a lagged copy of the surface.
The flux is applied only where the bottom layer is soil; under a glacier or lake column the
base stays zero flux, since the impermeable ice at the base of a glacier column has no way to
drain the melt the flux would produce.

<a id="bcupprsoih"></a>
## 21. bcUpprSoiH — upper boundary condition, soil hydrology

| Option | Description |
|---|---|
| presHead | prescribed head (prescribed volumetric liquid water content for the mixed form) |
| liq_flux | prescribed liquid water flux (infiltration) at the soil surface |

<a id="bclowrsoih"></a>
## 22. bcLowrSoiH — lower boundary condition, soil hydrology

| Option | Description |
|---|---|
| presHead | prescribed matric head at the bottom of the soil column |
| bottmPsi | flux computed from the matric head gradient in the lowest layer |
| drainage | free (gravity) drainage |
| zeroFlux | zero liquid water flux at the bottom of the soil column |

[`groundwatr`](#groundwatr) `= modflow` forces `presHead` here (the prescribed head is the
coupled MODFLOW 6 water table); any other value is rejected at start-up.

<a id="veg_traits"></a>
## 23. veg_traits — roughness length and displacement height

| Option | Description |
|---|---|
| Raupach_BLM1994 | [Raupach (1994)](https://doi.org/10.1007/BF00709229) simplified expressions |
| CM_QJRMS1988 | [Choudhury and Monteith (1988)](https://doi.org/10.1002/qj.49711448006) four-layer heat-budget model |
| vegTypeTable | constant values taken from the vegetation-type table |

<a id="rootprofil"></a>
## 24. rootProfil — rooting profile

| Option | Description |
|---|---|
| powerLaw | power-law root density with depth (also selected by `notPopulatedYet`) |
| doubleExp | double-exponential profile (Zeng, 2001) |

<a id="canopyemis"></a>
## 25. canopyEmis — canopy emissivity

| Option | Description |
|---|---|
| simplExp | simple exponential function of LAI+SAI |
| difTrans | function of the diffuse transmissivity of the canopy |

<a id="snowincept"></a>
## 26. snowIncept — canopy snow interception capacity

| Option | Description |
|---|---|
| stickySnow | maximum interception capacity increases with temperature (increased cohesion in warm conditions) ([Andreadis et al., 2009](https://doi.org/10.1029/2008WR007042)) |
| lightSnow | maximum interception capacity an inverse function of new-snow density ([Hedstrom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4)) |

<a id="windprfile"></a>
## 27. windPrfile — canopy wind profile

| Option | Description |
|---|---|
| exponential | exponential wind-speed decay that extends to the ground surface |
| logBelowCanopy | logarithmic wind-speed profile below the canopy |

<a id="astability"></a>
## 28. astability — atmospheric stability function

| Option | Description |
|---|---|
| standard | standard Monin-Obukhov similarity (after Anderson, 1976) |
| louisinv | Louis (1979) inverse-power function |
| mahrtexp | Mahrt (1987) exponential function |

<a id="compaction"></a>
## 29. compaction — snow densification

| Option | Description |
|---|---|
| consettl | constant settlement rate |
| anderson | semi-empirical method of Anderson (1976) (destructive metamorphism + overburden) |

<a id="snowlayers"></a>
## 30. snowLayers — snow layer combination and sub-division

| Option | Description |
|---|---|
| jrdn1991 | SNTHERM rules applied identically to all layers; grows to as many as 100 layers (Jordan, 1991) |
| CLM_2010 | CLM rules; combination/sub-division depend on layer index, giving up to a 5-layer snowpack ([CLM technical note, 2010](https://doi.org/10.1029/2011MS00045)) |

<a id="thcondsnow"></a>
## 31. thCondSnow — snow thermal conductivity

| Option | Description |
|---|---|
| tyen1965 | Yen (1965) |
| melr1977 | Mellor (1977) |
| jrdn1991 | Jordan (1991), as used in SNTHERM |
| smnv2000 | Smirnova et al. (2000) |

<a id="thcondsoil"></a>
## 32. thCondSoil — soil thermal conductivity

| Option | Description |
|---|---|
| funcSoilWet | function of soil wetness (Kersten-number approach) |
| mixConstit | volume-weighted mixture of soil constituents |
| hanssonVZJ | Hansson et al. (VZJ 2004); tuned to the Mizoguchi laboratory freezing experiment |

<a id="canopysrad"></a>
## 33. canopySrad — canopy shortwave radiation

| Option | Description |
|---|---|
| noah_mp | full Noah-MP two-stream implementation (includes its own albedo) |
| CLM_2stream | CLM two-stream model |
| UEB_2stream | UEB two-stream model (Mahat and Tarboton, 2012) |
| NL_scatter | simplified scattering method of Nijssen and Lettenmaier (1999) |
| BeersLaw | Beer's-law extinction (as in VIC) |

<a id="alb_method"></a>
## 34. alb_method — snow albedo

| Option | Description |
|---|---|
| conDecay | constant decay time scale (as in VIC, CLASS) |
| varDecay | variable decay from destructive metamorphism and soot content (BATS-style) |

<a id="spatial_gw"></a>
## 35. spatial_gw — spatial representation of groundwater

| Option | Description |
|---|---|
| localColumn | a separate groundwater store for each soil column |
| singleBasin | a single groundwater store shared across the whole GRU |

<a id="subrouting"></a>
## 36. subRouting — within-basin routing

| Option | Description |
|---|---|
| timeDlay | route runoff through a time-delay histogram (gamma-distribution unit hydrograph) |
| qInstant | no routing; runoff is delivered instantaneously |

<a id="snowdennew"></a>
## 37. snowDenNew — new-snow density

| Option | Description |
|---|---|
| hedAndPom | temperature-dependent exponential relation ([Hedstrom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4)); also selected by `notPopulatedYet` |
| anderson | Anderson (1976) temperature relation |
| pahaut_76 | Pahaut (1976); depends on air temperature and wind speed (Col de Porte) |
| constDens | constant new-snow density, taken directly from the `constSnowDen` parameter |

<a id="snowunload"></a>
## 38. snowUnload — unloading of intercepted snow

| Option | Description |
|---|---|
| meltDripUnload | temperature-driven unloading plus liquid drip; controlled by `snowUnloadingCoeff` and `ratioDrip2Unloading` ([Hedstrom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4); [Storck et al., 2002](https://doi.org/10.1029/2002WR001281)); also selected by `notPopulatedYet` |
| windUnload | temperature- and wind-driven unloading; controlled by `rateTempUnloading`, `rateWindUnloading` and minimum-temperature / minimum-windspeed thresholds ([Roesch et al., 2001](https://doi.org/10.1007/s003820100153)) |

<a id="nrgconserv"></a>
## 39. nrgConserv — form of the energy equations

Chooses the state variable / residual formulation for energy conservation.

| Option | Description |
|---|---|
| closedForm | temperature with a closed-form heat capacity |
| enthalpyForm | enthalpy as the state variable, with a temperature–enthalpy lookup table for soil |
| enthalpyFormAN | enthalpy as the state variable, with an analytical temperature–enthalpy relation for soil |

`enthalpyForm` / `enthalpyFormAN` require `num_method` to be `homegrown`, `kinsol` or `ida`.
With the `itertive` alias, the value is forced to `closedForm` for backward compatibility.

<a id="aquiferini"></a>
## 40. aquiferIni — initial aquifer fill level

| Option | Description |
|---|---|
| fullStart | start from the aquifer value in the initial-conditions file (default; for a cold start this is typically full, since draining to equilibrium is easier than filling) |
| emptyStart | start from an empty aquifer; intended only for comparing solution methods, not for realistic simulation |

`emptyStart` applies to a cold start only. Every restart file SUMMA writes carries the global
attribute `aqEmptyStarted = 1`, and a run that reads such a file keeps the aquifer storage the
file holds, exactly as `fullStart` would, with a warning in the log. An initial-conditions file
written by hand has no such attribute, so the aquifer is emptied there. This keeps a continued
run, or a spun-up calibration sample, from re-emptying the aquifer at every restart.

With `spatial_gw = singleBasin` the basin-average aquifer storage is not part of the restart
file, so that configuration always starts at the `aquiferIni` fill level -- full (1 m) or empty
-- whatever the initial-conditions file holds.

<a id="infratemax"></a>
## 41. infRateMax — maximum infiltration rate

| Option | Description |
|---|---|
| topmodel_GA | Green-Ampt with the [`hc_profile`](#hc_profile) conductivity evaluated at the wetting front (default; also selected by `notPopulatedYet` when `num_method = itertive`) |
| GreenAmpt | Green-Ampt maximum infiltration rate |
| noInfExc | maximum infiltration rate set very high, effectively disabling infiltration-excess runoff (saturation-excess runoff can still occur) |

This must match [`hc_profile`](#hc_profile): a depth-varying profile requires `topmodel_GA` or `noInfExc`, and `constant` with
`topmodel_GA` gives a deprecation warning.

<a id="surfrun_se"></a>
## 42. surfRun_SE — saturation-excess surface runoff

| Option | Description |
|---|---|
| homegrown_SE | SUMMA's built-in saturation-excess procedure (default; also selected by `notPopulatedYet`) |
| FUSEPRMS | PRMS saturation-excess formulation, as implemented in FUSE |
| FUSEAVIC | ARNO/VIC saturation-excess formulation, as implemented in FUSE |
| FUSETOPM | TOPMODEL saturation-excess formulation, as implemented in FUSE |
| zero_SE | no saturation-excess surface runoff |

<a id="read_force"></a>
## 43. read_force — how forcing is read

Renamed from `readForcing` in earlier versions.

| Option | Description |
|---|---|
| readPerStep | read one time step of forcing at a time (default; also selected by `notPopulatedYet`) |
| readFullSeries | read the whole forcing series for a file in one buffered read |

<a id="write_buff"></a>
## 44. write_buff — how output is buffered

Renamed from `writeOutput` in earlier versions.

| Option | Description |
|---|---|
| writePerStep | write model output every time step (default; also selected by `notPopulatedYet`) |
| writeFullSeries | buffer a whole output file in memory and write it once |

<a id="deeptherml"></a>
## 45. deepTherml — deep thermal state below the soil column

The temperature of the water that groundwater hands to the channel. The choices are whichever
mechanisms exist in the code, not a fixed set decided up front.

| Option | Description |
|---|---|
| none | the base of the hydrologically active soil column itself, floored at freezing (default; also selected by `notPopulatedYet`) |
| aquiferTemp | a well-mixed temperature carried by the big-bucket aquifer store |
| airTempGW | groundwater temperature scaled from the air temperature, after Wade et al. (2024) |
| bedrockLyrs | the deepest soil layers carry heat alone, extending the column into bedrock |

`none` is what SUMMA did before this decision existed, so an existing configuration keeps its
behaviour on upgrade. The only change is that the temperature is floored at freezing: the water
leaving the column is liquid even where the layer it left is frozen.

`aquiferTemp` gives the big-bucket aquifer its own temperature, `scalarAquiferTemp` (a new,
optional restart variable). Recharge arrives at the temperature of the base of the soil column
and mixes into the store, while baseflow and transpiration leave at the store's own temperature
and so do not change it:

```latex
S \frac{dT}{dt} = R\,(T_{rech} - T)
```

so the store relaxes towards the recharge temperature with a time constant of storage over
recharge, which for a real aquifer is months to years. That damping is the point. How much you
get depends on how much water the store holds: with `aquiferScaleFactor` and
`aquiferBaseflowRate` set so the bucket is nearly empty, the aquifer temperature simply tracks
the base of the soil column. It requires [`groundwatr`](#groundwatr) `= bigBuckt`, since there
is no store to carry a temperature otherwise, and is rejected at start-up with `qTopmodl` or
`noXplict`.

`airTempGW` is the cheap alternative of Wade et al. (2024, *Environmental Modelling and
Software* 171:105866, eq. 9). The groundwater temperature is bounded below by deep groundwater,
approximated by the mean annual air temperature, and above by the ground surface, approximated
by a smoothed daily air temperature, and one coefficient picks the effective sourcing depth:

```latex
T_{GW} = C_{ATGW}\,(AT_D - AT) + AT
```

`C_ATGW` (0-1) is a **basin** parameter, in `basinParamInfo.txt`, since Wade et al. tune it per
stream order; `AT_D` is the air temperature averaged over `gwTempWindow` days (2-14 in Wade et
al.) and `AT` its annual mean. `C_ATGW = 0` is deep, temporally invariant groundwater and
`C_ATGW = 1` is shallow groundwater that follows the ground surface; the 0.5 default is a
neutral starting point for a calibration, not a recommendation. Both means are kept as running
means of the forcing, so no extra input is needed, but the annual mean starts at the first air
temperature the run sees and needs a year of spin-up before it means what its name says. They
are written to the restart file (`scalarAirTempWindow`, `scalarAirTempAnnual`), so a spun-up run
carries them forward.

In permafrost `airTempGW` needs care: its lower bound is the mean annual air temperature, which
is below freezing in the colder interior and northern zones, while real sub-permafrost or talik
groundwater is at or above 0 C. The result is floored at freezing, so the formula there mostly
returns 0 C rather than anything physical. That is why it is a per-basin opt-in rather than a
default: it is useful in parts of Alaska, not all of it. `aquiferTemp` with a geothermal lower
boundary ([`bcLowrTdyn`](#bclowrtdyn) `presFlux`) is the process-based cold-region path.

`bedrockLyrs` extends the soil column downward into bedrock. The deepest `nBedrock` soil layers
keep their energy state but carry no hydrology, the way the `noThetaChange` layers at the base
of a glacier column do: the Richards solve never reaches them, water leaves the column at the
base of the active layers, and what they add is thermal memory. Their thermal conductivity and
porosity come from `thCond_bedrock` and `theta_sat_bedrock` rather than the usual soil texture
parameters, and `theta_sat_bedrock` also fixes their liquid water content, held fully saturated
for the whole run.

`nBedrock` is an optional variable in the initial conditions file, beside `nSoil` and
`nLakeFrz`, and is fixed for the run. It is per domain because it has to agree with that
domain's `nSoil`, and it is rejected on a glacier column, whose ice already owns the layers
below, and where it would leave no hydrologically active soil.

Pair it with [`bcLowrTdyn`](#bclowrtdyn) `presFlux`. The bedrock layers are then started on the
steady gradient that flux holds,

```latex
T(z) = T_{base} + \frac{q\,(z - z_{base})}{k}
```

anchored at the deepest active layer, so the column begins in equilibrium with its own lower
boundary and needs no spin-up, in the same way glacier ice is not spun up. A 30 m column on
gulkana, eight active layers over seven bedrock layers, drifts 0.09 K at its base over fifteen
months, against 1.36 K for the same column with every layer hydrologically active and started
uniform.

`bedrockLyrs` is mutually exclusive with `aquiferTemp` and `airTempGW`: `deepTherml` is a single
choice, so nothing needs to guard against combining them.

<a id="hyporhtdyn"></a>
## 46. hyporhTdyn — hyporheic exchange in a stream domain

| Option | Description |
|---|---|
| none | no exchange with the bed (default; also selected by `notPopulatedYet`) |
| proxy | a tuned fraction of the reach flow returns at the temperature the reach had `hypLag` hours ago, after Wade et al. (2024) |

SUMMA's stream domain otherwise assumes zero exchange with the bed. `proxy` is the conceptual
representation of Wade et al. (2024, *Environmental Modelling and Software* 171:105866,
eqs. 11-12): a fraction `hypFrac` of the reach flow leaves into the bed and returns at
`scalarHypTemp`, the mean of the reach's own outlet temperature over the previous `hypLag`
hours:

```latex
\Phi_{hyp} = \frac{H_{frac}\,Q}{V}\,(T_{hyp} - T)
```

It is written as one more `(T_in - T)` exchange alongside the upstream, lateral and surface
inflows, so it moves no water: the flow returns what it took, and the mean reach temperature is
unchanged while the diurnal swing is damped. In Wade et al. this was the single largest
accuracy gain of the whole study, cutting headwater daily-maximum error from 3.70 to 1.44 C.

`hypFrac` (0-1) and `hypLag` (1-24 h) are both tuned; the 0.2 and 6 h defaults are a starting
point, not a recommendation. Wade et al. tune them per stream order, which here means one value
per reach: a stream domain is the only domain with area in its HRU, so an HRU-level parameter
already gives one value per reach. The reach's temperature history is kept in the restart
variable `hypTempPast`; a run that starts without one builds the history as it goes and leaves
the exchange out of the energy balance until it has some.

The mechanism is not blocked by permafrost, since streams usually keep an unfrozen talik even
under continuous permafrost, but it has not been validated there. `hypFrac` and `hypLag` are an
empirical fit rather than derived SUMMA physics, which is why this is opt-in and separate from
[`deepTherml`](#deeptherml): the two compose, and either can be used without the other.
