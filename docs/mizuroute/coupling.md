# Coupling architecture

The SUMMA--mizuRoute coupling implements the design described in
[Coupling design](design.md) through a deliberately thin software boundary
between SUMMA and mizuRoute.

Rather than incorporating mizuRoute logic throughout the SUMMA codebase, the
coupled implementation combines a small host-model interface, a compatibility
layer, and selected native mizuRoute components.

At the highest level, the implementation consists of three layers:

1. **The host-model interface**, which performs the small number of explicit
   exchanges between SUMMA and mizuRoute.
2. **The compatibility layer**, which adapts dependencies expected by the
   native mizuRoute source to the coupled environment.
3. **The native mizuRoute components**, which perform river-network
   initialization and routing.

Model-specific implementation logic remains on the appropriate side of this
boundary: SUMMA logic remains in SUMMA and routing logic remains in mizuRoute.

## Compatibility shims

A key part of the coupling development was separating the mizuRoute routing
components from assumptions associated with the standalone mizuRoute
executable. Some native mizuRoute modules have substantial dependencies on
infrastructure used by the standalone model for global state, numerical types,
parallel I/O, model initialization, and diagnostic calculations. These
dependencies are either unnecessary in a coupled SUMMA simulation or duplicate
functionality already provided by SUMMA.

Rather than modifying the upstream mizuRoute source to create a SUMMA-specific
version, the coupled build currently uses a small set of **compatibility
shims**. These shims provide the module names, data definitions, constants, or
routines expected by the native mizuRoute source while adapting those
dependencies to the coupled environment.

The principal compatibility shims include:

- **`nrtype_shim`**, which provides the numerical type definitions expected by
  the mizuRoute source;
- **`globalData_shim`**, which provides the subset of global model state needed
  by the routing components;
- **`pio_utils_shim`**, which provides the NetCDF-related definitions required
  by mizuRoute without bringing the standalone ParallelIO infrastructure into
  the coupled model;
- **`init_model_data_shim`**, which provides the initialization functionality
  required to construct the mizuRoute model state within the coupled
  environment; and
- **`water_balance_shim`**, which provides the water-balance functionality
  required by the routing code without introducing the complete standalone
  implementation.

The CMake configuration explicitly includes these compatibility modules
alongside the required native mizuRoute components.

The shims are maintained on the SUMMA side of the interface rather than in the
mizuRoute Git submodule. The mizuRoute submodule therefore remains an
unmodified version of the upstream model, while the compatibility requirements
of the coupled application are isolated within SUMMA.

The current coupling uses the upstream **mizuRoute v3.1.1** release. No source
files within the mizuRoute Git submodule have been modified for the SUMMA
coupling. All changes required to compile and use mizuRoute within SUMMA are
contained in the SUMMA codebase, principally in the coupling interface,
compatibility shims, and CMake build configuration.

The shims are not intended to define the long-term coupling interface. Their
current use largely reflects dependencies within mizuRoute modules that make
it difficult to compile individual initialization and routing components
without also introducing infrastructure required by the standalone model.

The longer-term goal is to work with the mizuRoute developers to reduce these
dependencies and establish cleaner interfaces between the routing algorithms
and the surrounding model infrastructure. As those interfaces are improved,
the compatibility shims should be substantially reduced or, where possible,
eliminated.

Importantly, reducing the need for shims should not require moving
SUMMA-specific code into mizuRoute. The goal is instead to make the native
mizuRoute components sufficiently self-contained that they can be called
directly by SUMMA or by other host models.

## Coupling points

With these dependencies isolated, SUMMA interacts with mizuRoute at only five
well-defined points in the simulation workflow:

1. **Provide mizuRoute configuration to SUMMA.**  
   mizuRoute options are supplied through the TOML configuration used by
   SUMMA. The configuration is activated with the `-c` or `--config`
   command-line option. The existing SUMMA control files are unchanged.

2. **Define mizuRoute NetCDF output.**  
   When coupled routing is enabled, the dimensions and variables needed for
   mizuRoute results are defined within the SUMMA NetCDF output file. SUMMA
   retains responsibility for managing the output file.

3. **Initialize mizuRoute.**  
   Before time stepping, the interface provides the information required to
   initialize mizuRoute. The mizuRoute initialization code reads the routing
   information, constructs the river-network topology, initializes the routing
   data structures, and prepares the selected routing method.

4. **Run river-network routing.**  
   During model execution, SUMMA provides the runoff field required by
   mizuRoute and invokes the network-routing operation. The routing
   calculations themselves remain within mizuRoute.

5. **Write mizuRoute NetCDF output.**  
   Routed quantities are made available to the SUMMA output system and written
   into the SUMMA NetCDF output files. This keeps file management and model I/O
   under the control of the host model rather than requiring the standalone
   mizuRoute I/O workflow.

These five coupling points define the practical boundary between SUMMA and
mizuRoute. The interface is intentionally narrow so that the routing
implementation remains localized within mizuRoute and the amount of
SUMMA-specific coupling code remains small.

## Runoff exchange and routing responsibilities

The coupled configuration separates runoff routing through the unresolved
drainage network (handled by SUMMA) from routing through the explicit river
network (handled by mizuRoute).

SUMMA computes runoff within its HRUs, aggregates runoff to the GRU level, and
applies its existing time delay routing (using a parameterized unit hydrograph)
through the unresolved river network. The resulting routed GRU runoff is then
passed across the coupling interface to mizuRoute.

Within mizuRoute, the runoff supplied by SUMMA is stored on the spatial
elements defined by the host model. These elements may differ from the HRUs
associated with the mizuRoute river network. When necessary, mizuRoute
spatially remaps the supplied runoff onto the river-network HRUs before
aggregating the resulting runoff to river reaches and routing flow through the
explicit river network.

The coupled workflow can therefore be summarized as:

```text
SUMMA HRUs
    |
    | optional lateral flow among HRUs
    v
aggregate runoff to SUMMA GRUs
    |
    | SUMMA routing through the unresolved river network
    v
routed SUMMA GRU runoff
    |
    | ------------------
    | coupling interface
    | ------------------
    v
mizuRoute runoff input
    |
    | optional spatial remapping
    v
mizuRoute river-network HRUs
    |
    | aggregate runoff to reaches
    v
lateral reach inflow
    |
    | explicit river-network routing
    v
routed streamflow
```

The mizuRoute version of the parameterized unit hydrograph to represent
routing through the unresolved drainage network is deliberately not included
in the set of mizuRoute source files compiled and linked into SUMMA. The
corresponding time-delay routing is already performed by SUMMA before runoff
crosses the coupling interface. The coupled mizuRoute components therefore
provide spatial remapping (where required), aggregation of runoff to river
reaches, and routing through the explicit river network.

Spatial remapping is optional. If runoff from the host model is already defined
on the mizuRoute river-network HRUs, the remapping step can be skipped.

The current SUMMA coupling supplies runoff at GRU resolution. The interface
could in principle instead supply runoff from individual SUMMA HRUs and rely on
mizuRoute for the subsequent spatial aggregation and unresolved-network
routing, but this configuration is not currently implemented.

## Stream temperature exchange

When SUMMA carries **stream HRUs** (see the [input description](../input_output/SUMMA_input.md#infile_stream)),
a second exchange runs after every routing step. mizuRoute keeps routing the
water; SUMMA solves the temperature of the water column of each reach and
routes heat downstream itself. The formulation follows Wanders et al. (2019,
*WRR*, DynWat) after van Beek et al. (2012, *WRR*): a well-mixed reach whose
energy balance is the surface exchange (SUMMA's usual bare-surface fluxes over
open water or ice), conduction to the bed and the heat carried by the water
entering from upstream and from the local catchment,

```text
rho_w c_p d(h T)/dt = H_surface + H_bed + rho_w c_p [ Q_up (T_up - T) + q_lat (T_lat - T) + q_sfc (T_sfc - T) ] / A
```

with `A` the reach planform area. Ice is treated as in van Beek et al.: when
the water column reaches the freezing point and the energy balance is still
negative, ice forms in the water layers and, once per step, rises into an
**ice cover** at the top of the lake layers (`lakeIceCover.f90`, index
`nLakeFrz` counts the cover layers). The cover is a lake layer handled like
glacier ice: impermeable, sharp freezing curve, melt squeezed to its surface
and into the flow, thinning as it melts; the water beneath stays mixed and
advected by the reach flow. A cover thinner than `lakeIceMinThick` (5 mm,
the breakup thickness of Wanders et al.) returns to the water. Snow builds
on the cover; rain and snowmelt on the cover run off into the reach.

The sequence within a time step is:

```text
land HRUs of every GRU (run_oneGRU, parallel over GRUs)
    |  runoff and the heat it carries, both through the SUMMA unit hydrograph:
    |  surface runoff at the surface-layer temperature, drainage and baseflow
    |  at the bottom-soil temperature, glacier melt at the freezing point
    v
coupling(:)%qsim, coupling(:)%esim          (m s-1 and W m-2 per GRU)
    |
    | ------------------ coupling interface ------------------
    v
mizuRoute routes the water                  (route_mizuroute_from_summa)
    |
    v
reach discharge, upstream inflow, volume,   (get_mizuroute_reach_hydraulics)
lateral inflow; depth = volume / (length x width), velocity = Q length / volume
    |
    v
energy flux remapped and aggregated to      (remap_lateral_energy: the same
reaches, divided by the lateral inflow       remap_runoff + basin2reach as the runoff)
    |
    | ------------------ back to SUMMA ------------------
    v
network pass, reaches in routing order      (run_streamNetwork, serial)
    |  T_up from the outflow of the upstream reaches, then the stream domain of
    |  the reach is run (run_oneHRU with streamPass=.true.), and its column
    |  temperature becomes the outflow temperature of the reach; a reach with
    |  no stream HRU only mixes what flows into it
    v
T_reach, v_reach output; scalarStream* fluxes and scalarStreamTemp per stream HRU
```

The per-reach arrays live in `summa1_type_dec%stream_net` (`data_types.f90`,
type `stream_network`). The coupling module only fills the mizuRoute side of
them; the network pass and the column physics are SUMMA code
(`streamtemp.f90`, the lake layers of `coupled_em`, the advective source in
`snowLakeSoilGlceNrgFlux.f90`, `lakeLiqFlux` in `snowLakeGlceLiqFlux.f90`,
`lakeResize`/`lakePrescribeDepth` in `snowLakeGlceDepth.f90`).

Approximations of this first implementation: one well-mixed column per reach
(no longitudinal sub-reaches); the liquid depth of the column is prescribed
from the reach volume once per step (each lake layer keeping its share of the
column, with a floor so a layer that melted or sublimated away refills), so the
water mass stays entirely in mizuRoute; the ice cover is at most one layer
(`nLakeIceLayers_poss` in `globalData.f90`) and forms only once enough ice
is present to stand (twice `lakeIceMinThick`); no bed seepage
or hyporheic exchange; no shortwave penetration below the top lake layer;
mizuRoute never sees the ice, so winter depth and velocity are open-water
values; and open-water evaporation exceeding all other runoff of a GRU is not
debited from the channel (the network takes no negative lateral inflow). The
network pass is serial and, like the rest of the coupled mizuRoute, needs the
whole domain on one process.

## Build-system separation

The software boundary described above is also reflected in the CMake build
system. mizuRoute support is optional and is enabled with:

```bash
-DUSE_MIZUROUTE=ON
```

The default is `OFF`, so SUMMA applications that do not require river routing
do not need to include the mizuRoute routing components.

CMake does not build the standalone mizuRoute executable and link SUMMA against
it. Instead, the coupled build assembles the routing capability from the
portions of the mizuRoute source required for initialization and routing,
together with the compatibility and coupling layers described above.

The source dependencies are organized into groups for:

- numerical and data types;
- base mizuRoute utilities;
- compatibility shims;
- river-network topology;
- initialization;
- runoff remapping;
- routing methods;
- output; and
- the SUMMA--mizuRoute coupling interface.

The CMake source configuration separates these components into `MIZU_TYPE`,
`MIZU_BASE`, `MIZU_SHIM`, `MIZU_TOPO`, `MIZU_INIT`, `MIZU_RMAP`,
`MIZU_ROUTE`, `MIZU_OUTPUT`, and `MIZU_COUPLING`.

This decomposition allows the coupled build to select the mizuRoute
functionality required by SUMMA without introducing the complete standalone
mizuRoute execution environment. Where a native mizuRoute module has
dependencies that are inappropriate for the coupled application, the
corresponding compatibility shim is compiled instead.

The resulting build therefore combines the native mizuRoute river-network
algorithms with a small amount of host-specific compatibility code and a thin
runtime interface.

## TOML configuration dependency

TOML configuration support is handled somewhat differently from mizuRoute
itself. The `toml-f` parser is included in SUMMA as a Git submodule and is
built as a standard SUMMA dependency rather than only when mizuRoute is
enabled.

The TOML configuration currently provides the additional options required for
capabilities such as coupled mizuRoute execution. The existing SUMMA control
files remain unchanged. The longer-term intention, however, is to migrate
SUMMA control information to TOML, so TOML support is not treated as an
optional mizuRoute-specific dependency.
