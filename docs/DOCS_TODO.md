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

## Design note worth writing up

- Melt enters a glacier debris column as a bounded lower-boundary flux,
  `max(scalarGlceMelt, min(0, scalarDrainage))`, with a finite conductivity at the base so
  capillary suction draws it up. Uptake then follows the column's own moisture state and stops
  when the column saturates. Giese et al. (2020, The Cryosphere 14:1555) instead inject the melt
  as a mass source with zero conductivity at the interface, which forces water in whether or not
  the debris can take it, and needs a per-layer runoff sink to remove it again.
