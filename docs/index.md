# Structure for Unifying Multiple Modeling Alternatives: SUMMA

[![Build Status](https://travis-ci.org/NCAR/summa.svg?branch=develop)](https://travis-ci.org/NCAR/summa)
[![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/NCAR/SUMMA/master/COPYING)
[![Documentation Status](https://readthedocs.org/projects/summa/badge/?version=latest)](http://summa.readthedocs.org/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.800772.svg)](https://doi.org/10.5281/zenodo.800772)

SUMMA (Clark et al., [2015a](#clark_2015a);[b](#clark_2015b);[c](#clark_2015c); [2021](#clark_2021) is a hydrologic modeling framework that can be used for the systematic analysis of alternative model conceptualizations with respect to flux parameterizations, spatial configurations, and numerical solution techniques. It can be used to configure a wide range of hydrological model alternatives and we anticipate that systematic model analysis will help researchers and practitioners understand reasons for inter-model differences in model behavior. When applied across a large sample of catchments, SUMMA may provide insights in the dominance of different physical processes and regional variability in the suitability of different modeling approaches. An important application of SUMMA is selecting specific physics options to reproduce the behavior of existing models – these applications of "**model mimicry**" can be used to define reference (benchmark) cases in structured model comparison experiments, and can help diagnose weaknesses of individual models in different hydroclimatic regimes.

SUMMA is built on a common set of conservation equations and a common numerical solver, which together constitute the  “**structural core**” of the model. Different modeling approaches can then be implemented within the structural core, enabling a controlled and systematic analysis of alternative modeling options, and providing insight for future model development.

The important modeling features are:

 1. The formulation of the conservation model equations is cleanly separated from their numerical solution;

 2. Different model representations of physical processes (in particular, different flux parameterizations) can be used within a common set of conservation equations; and

 3. The physical processes can be organized in different spatial configurations, including model elements of different shape and connectivity (e.g., nested multi-scale grids and HRUs).


## Documentation
SUMMA documentation is available [online](http://summa.readthedocs.io/) and remains a work in progress. Additional SUMMA information including publications, test data sets, and sample applications can be found on the [SUMMA web site](http://www.ral.ucar.edu/projects/summa) at NCAR.

## Developer Guidelines for Contributing New Modular Components
New modular components may be added by using similar existing modular components as a template. The following steps (if applicable) may be used as a guideline. This process is illustrated using the addition of a new surface hydrology flux parameterization as an example.
1. Identify a similar model component
    * Identify the appropriate subdirectory within SUMMA's `source` directory such as:
        * `driver`: high-level program operations including the main driver
        * `dshare`: modules related to data storage and access
        * `engine`: low-level operations for physical and numerical processes
            * e.g., applies to flux calculations
    * Identify the appropriate source file and module
        * source files have self-explanatory names
            * e.g., `soilLiqFlx.f90` corresponds to operations for liquid water fluxes in soil
        * each source file generally contains one module
            * e.g., `soilLiqFlx.f90` contains `soilLiqFlx_module`
    * Identify the appropriate procedure
        * isolate the module procedure
            * e.g., within `soilLiqFlx_module`, the `surfaceFlx` module subroutine handles operations for surface hydrology fluxes
        * isolate the internal procedure
            * e.g., within the `contains` block of `surfaceFlx`, we have `update_surfaceFlx_prescribedHead` containing operations for specifying a prescirbed pressure head surface boundary condition
            * `update_surfaceFlx_prescribedHead` may be used as a template for our example contribution
        * note that procedure names in SUMMA are organized using the terms *initialize*, *update*, and *finalize* to categorize operation types
            * *initialize* procedures are used for initial setup steps (initialization of variables, memory allocation, etc.)
            * *update* procedures are used for major computational operations (e.g., flux calculations)
            * *finalize* procedures are for post-processing and final error control checks
2. Determine input and output variables
    * Found by examining dummy variables in argument lists
        * Note that internal procedures inherit the dummy variables from the applicable module procedure by default
        * e.g., for the `update_surfaceFlx_prescribedHead` internal subroutine, the argument list of the `surfaceFlx` module subroutine applies: `subroutine surfaceFlx(io_soilLiqFlx,in_surfaceFlx,io_surfaceFlx,out_surfaceFlx)`
    * Dummy variables may be objects with multiple data and procedure components
        * Such objects are declared using derived types (most commonly defined in `data_types.f90`)
        * Objects may be used to concisely interface data between the procedure and the caller
        * for SUMMA objects, the nomenclature `in_foobar`, `io_foobar`, and `out_foobar` is used for objects that interface input, input-output, and output data between the `foobar` procedure and its caller, respectively
    * The `intent` attribute within dummy variable declarations indicates usage for input, input-output, or output
        * e.g., within `surfaceFlx` we have `type(in_type_surfaceFlx) ,intent(in)    :: in_surfaceFlx`, indicating the `in_surfaceFlx` object is for input data only
        * as noted above, the `in_surfaceFlx` object interfaces input data between the `surfaceFlx` module subroutine and its caller (the `soilLiqFlx` module subroutine)
3. Create a skeleton of the new procedure
    * Choose a self explanatory name for the new procedure
        * e.g. `update_surfaceFlx_example_flux`
    * e.g., at the conclusion of this step, we would have a skeleton within the `contains` block of `surfaceFlx` similar to the following:


    ````fortran
    subroutine update_surfaceFlx_example_flux
    ! main computations for the calculation of an example flux
    
    end subroutine update_surfaceFlx_example_flux
    ````

    * For new module procedures, the argument list from the template routine should be adjusted to match the needs of the new procedure
    * For internal procedures (such as the example above), adding new data to interface objects from the corresponding module procedure may be required (see next step)
4. Update derived type definitions for interface objects
    * It may be desirable to add data components to existing objects related to the template procedure
        * e.g., adding a new numerical constant to be used in calculating a surface hydrology flux would require interfacing that data to the `update_surfaceFlx_example_flux` subroutine, which can be done using the `in_surfaceFlx` object
    * Derived type definitions for interface objects are found in `source/dshare/data_types.f90`
        * e.g., for the `in_surfaceFlx` object we have the `in_type_surfaceFlx` derived type in the `data_types` module:
     
        ```fortran
        type, public :: in_type_surfaceFlx ! intent(in) data
          ! input: model control
          logical(lgt) :: firstSplitOper   ! flag indicating if desire to compute infiltration
          logical(lgt) :: deriv_desired    ! flag to indicate if derivatives are desired
          integer(i4b) :: ixRichards       ! index defining the option for the Richards equation (moisture or mixdform)
          integer(i4b) :: bc_upper         ! index defining the type of boundary conditions
          integer(i4b) :: nRoots           ! number of layers that contain roots
          integer(i4b) :: ixIce            ! index of lowest ice layer
          integer(i4b) :: nSoil            ! number of soil layers
          ! [...] ! additional data components here
         contains
          procedure :: initialize => initialize_in_surfaceFlx
        end type in_type_surfaceFlx
        ```

        * adding a new numerical constant (say `example_flux_constant`) may be done as follows:

        ```fortran
        type, public :: in_type_surfaceFlx ! intent(in) data
          ! input: model control
          logical(lgt) :: firstSplitOper   ! flag indicating if desire to compute infiltration
          logical(lgt) :: deriv_desired    ! flag to indicate if derivatives are desired
          integer(i4b) :: ixRichards       ! index defining the option for the Richards equation (moisture or mixdform)
          integer(i4b) :: bc_upper         ! index defining the type of boundary conditions
          integer(i4b) :: nRoots           ! number of layers that contain roots
          integer(i4b) :: ixIce            ! index of lowest ice layer
          integer(i4b) :: nSoil            ! number of soil layers
          ! input: values for the example flux
          real(rkind)  :: example_flux_constant ! numerical constant for example flux
          ! [...] ! additional data components here
         contains
          procedure :: initialize => initialize_in_surfaceFlx
        end type in_type_surfaceFlx
        ```

        * note that SUMMA uses the following `kind` parameters: `lgt` for logical variables, `i4b` for integer variables, and `rkind` for real variables
    * Additionally, we have procedure components for *initialize* and *finalize* operations for data interfacing
        * e.g., `call in_surfaceFlx % initialize` points to the `initialize_in_surfaceFlx` class procedure (in the `contains` block of the `data_types` module) for initializing data components:

        ```fortran
        subroutine initialize_in_surfaceFlx(in_surfaceFlx,nRoots,ixIce,nSoil,ibeg,iend,in_soilLiqFlx,io_soilLiqFlx,&
                                           &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                           &iLayerHeight,dHydCond_dTemp,iceImpedeFac)
         class(in_type_surfaceFlx),intent(out) :: in_surfaceFlx ! input object for surfaceFlx
         ! [...] ! additional variable declarations here

         associate(&
          ! model control
          firstSplitOper         => in_soilLiqFlx % firstSplitOper,                      & ! flag to compute infiltration
          deriv_desired          => in_soilLiqFlx % deriv_desired,                       & ! flag indicating if derivatives are desired
          ixRichards             => model_decisions(iLookDECISIONS%f_Richards)%iDecision,& ! index of the form of the Richards equation
          ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision & ! index defining the type of boundary conditions
         &)
          ! intent(in): model control
          in_surfaceFlx % firstSplitOper = firstSplitOper          ! flag indicating if desire to compute infiltration
          in_surfaceFlx % deriv_desired  = deriv_desired           ! flag indicating if derivatives are desired
          in_surfaceFlx % ixRichards     = ixRichards              ! index defining the form of the Richards equation (moisture or mixdform)
          in_surfaceFlx % bc_upper       = ixBcUpperSoilHydrology  ! index defining the type of boundary conditions (Neumann or Dirichlet)
          in_surfaceFlx % nRoots         = nRoots                  ! number of layers that contain roots
          in_surfaceFlx % ixIce          = ixIce                   ! index of lowest ice layer
          in_surfaceFlx % nSoil          = nSoil                   ! number of soil layers
         end associate

         ! [...] ! additional associate blocks here
        end subroutine initialize_in_surfaceFlx
        ```
        
        * new data components, such as `example_flux_constant`, must be applied within the procedure components:

        ```fortran
        subroutine initialize_in_surfaceFlx(in_surfaceFlx,nRoots,ixIce,nSoil,ibeg,iend,in_soilLiqFlx,io_soilLiqFlx,&
                                           &model_decisions,prog_data,mpar_data,flux_data,diag_data,&
                                           &iLayerHeight,dHydCond_dTemp,iceImpedeFac,example_flux_constant)
         class(in_type_surfaceFlx),intent(out) :: in_surfaceFlx ! input object for surfaceFlx
         ! [...] ! additional variable declarations here
         real(rkind),intent(in) :: example_flux_constant ! declaration for new constant

         associate(&
          ! model control
          firstSplitOper         => in_soilLiqFlx % firstSplitOper,                      & ! flag to compute infiltration
          deriv_desired          => in_soilLiqFlx % deriv_desired,                       & ! flag indicating if derivatives are desired
          ixRichards             => model_decisions(iLookDECISIONS%f_Richards)%iDecision,& ! index of the form of the Richards equation
          ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision & ! index defining the type of boundary conditions
         &)
          ! intent(in): model control
          in_surfaceFlx % firstSplitOper = firstSplitOper          ! flag indicating if desire to compute infiltration
          in_surfaceFlx % deriv_desired  = deriv_desired           ! flag indicating if derivatives are desired
          in_surfaceFlx % ixRichards     = ixRichards              ! index defining the form of the Richards equation (moisture or mixdform)
          in_surfaceFlx % bc_upper       = ixBcUpperSoilHydrology  ! index defining the type of boundary conditions (Neumann or Dirichlet)
          in_surfaceFlx % nRoots         = nRoots                  ! number of layers that contain roots
          in_surfaceFlx % ixIce          = ixIce                   ! index of lowest ice layer
          in_surfaceFlx % nSoil          = nSoil                   ! number of soil layers
         end associate
         
         ! [...] ! additional associate blocks here
        
         ! assignment statements for the example flux
         in_surfaceFlx % example_flux_constant = example_flux_constant ! numerical constant for example flux
        
        end subroutine initialize_in_surfaceFlx
        ```

        * for the above example, we have added a dummy variable for the new example flux constant and an assignment statement to initialize the new data component `in_surfaceFlx % example_flux_constant`
        * note that the corresponding call to `in_surfaceFlx` within the `soilLiqFlx` subroutine would need to be updated to include the additional argument `example_flux_constant`
5. Add operations to the skeleton procedure
    * add main operations within the skeleton procedure created in the above steps to complete the new modular component
        * for the example flux parameterization (using a toy model of constant infiltration), we have:
     
        ````fortran
        subroutine update_surfaceFlx_example_flux
         ! main computations for the calculation of an example flux
         associate(&
          ! input: flux at the upper boundary
          scalarRainPlusMelt => in_surfaceFlx % scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
          ! input: numerical constants
          example_flux_constant => in_surfaceFlx % example_flux_constant
          ! output: runoff and infiltration
          scalarSurfaceRunoff       => out_surfaceFlx % scalarSurfaceRunoff       , & ! surface runoff (m s-1)
          scalarSurfaceInfiltration => out_surfaceFlx % scalarSurfaceInfiltration   & ! surface infiltration (m s-1)
         &)
          scalarSurfaceInfiltration = example_flux_constant                          ! toy model of constant infiltration
          scalarSurfaceRunoff       = scalarRainPlusMelt - scalarSurfaceInfiltration ! compute surface runoff
         end associate
        end subroutine update_surfaceFlx_example_flux
        ````

## Credits
SUMMA's initial implementation is described in two papers published in [Water Resources Research](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1944-7973). If you use SUMMA, please credit these two publications.

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, 2015a: A unified approach for process-based hydrologic modeling: Part 1. Modeling concept. _Water Resources Research_, [doi:10.1002/2015WR017198](http://dx.doi.org/10.1002/2015WR017198).<a id="clark_2015a"></a>

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015b: A unified approach for process-based hydrologic modeling: Part 2. Model implementation and case studies. _Water Resources Research_, [doi:10.1002/2015WR017200](http://dx.doi.org/10.1002/2015WR017200).<a id="clark_2015b"></a>
 
 * Clark, M. P., Zolfaghari, R., Green, K. R., Trim, S., Knoben, W. J. M., Bennett, A., Nijssen, B., Ireson, A., Spiteri, R. J., 2021: The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests. _Journal of Hydrometeorology_, [doi:10.1175/JHM-D-20-0175.1](http://dx.doi.org/10.1175/JHM-D-20-0175.1).<a id="clark_2021"></a>

In addition, an NCAR technical note describes the SUMMA implementation in detail:

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015c: The structure for unifying multiple modeling alternatives (SUMMA), Version 1.0: Technical Description. _NCAR Technical Note NCAR/TN-514+STR_, 50 pp., [doi:10.5065/D6WQ01TD](http://dx.doi.org/10.5065/D6WQ01TD).<a id="clark_2015c"></a>


## License
SUMMA is distributed under the GNU Public License Version 3. For details see the file `COPYING` in the SUMMA root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
