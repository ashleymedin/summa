# Documentation still to write

Work that has landed in code but is not yet described in the user-facing docs.
Delete an entry once its docs are written.

## The deep thermal state and stream temperature work (Sep 2026)

Covered already: `docs/configuration/SUMMA_model_decisions.md` has `bcLowrTdyn = presFlux`,
`deepTherml` and `hyporhTdyn`; `docs/input_output/SUMMA_input.md` lists the new restart
variables; `docs/whats-new.md` and `docs/assets/changes_fromV3Summa.txt` carry all of it.

Still to do:

- `docs/input_output/SUMMA_output.md` does not list the new output variables:
  `scalarAquiferTemp`, `scalarAirTempWindow`, `scalarAirTempAnnual`, `scalarFrostTableDepth`,
  `scalarActiveLayerDepth`, `scalarHypTemp`.
- No parameter reference page lists `lowerBoundNrgFlux`, `gwTempWindow`, `hypFrac`, `hypLag`
  (local) or `C_ATGW` (basin). They are in `localParamInfo.txt` / `basinParamInfo.txt` defaults
  and in the decisions page, but there is no single parameter table to add them to.
- `docs/configuration/SUMMA_configuration.md` says nothing about setting up a deep soil column
  for permafrost: column depth, the `parSoil` depth dimension in `trialParams.nc` for bedrock
  properties, and which lower boundary conditions to pair with it.
- `utils/test/README.md` and `utils/test/test_regression/README.md` do not mention that the
  bundled stream-temperature test also exercises `hyporhTdyn`.

## Known limitation, not yet fixed

- A glacier debris column still leans on the clamped lower-boundary flux
  (`max(scalarGlceMelt, min(0, scalarDrainage))`) to bring melt up into the debris. Giese et al.
  (2020, The Cryosphere 14:1555) instead add the melt of the top ice layer to the debris as a
  mass source, with zero Darcy conductivity at the interface, and give every debris layer a
  linear-reservoir runoff sink. See the discussion in `changes_fromV3Summa.txt` entry 76.
