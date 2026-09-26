# Documentation still to write

Work that has landed in code but is not yet described in the user-facing docs.
Delete an entry once its docs are written.

## The deep thermal state and stream temperature work (Sep 2026)

Covered already: `docs/configuration/SUMMA_model_decisions.md` has `bcLowrTdyn = presFlux`,
`deepTherml` (including `bedrockLyrs`) and `hyporhTdyn`; `docs/configuration/SUMMA_configuration.md`
has the deep bedrock column setup (`nBedrock`, `bedrockDepth`, `nBedrockLayers`); `docs/input_output/SUMMA_input.md`
lists the new restart variables including `nBedrock`; `docs/whats-new.md` and
`docs/assets/changes_fromV3Summa.txt` carry all of it.

Still to do:

- `docs/input_output/SUMMA_output.md` does not list the new output variables:
  `scalarAquiferTemp`, `scalarAirTempWindow`, `scalarAirTempAnnual`, `scalarFrostTableDepth`,
  `scalarActiveLayerDepth`, `scalarHypTemp`. NOTE: this is not unique to this work — the file has
  never enumerated individual output variables for anything (it documents dimensions and file
  types, not variables); adding one would be a new kind of page, not a gap in existing coverage.
- No parameter reference page lists `lowerBoundNrgFlux`, `gwTempWindow`, `hypFrac`, `hypLag`,
  `thCond_bedrock`, `theta_sat_bedrock` (local) or `C_ATGW` (basin). Same caveat as above: no
  SUMMA parameter, old or new, has a docs table — they live in `localParamInfo.txt` /
  `basinParamInfo.txt` and the decisions page only.
- `utils/test/README.md` and `utils/test/test_regression/README.md` do not mention that the
  bundled stream-temperature test also exercises `hyporhTdyn`. NOTE: checked — it doesn't yet;
  `make_stream_domain.py` copies the Provo settings unmodified, which have `hyporhTdyn = none`.
  Turning it on is a test-config change (and would need its own reference-output check), not a
  docs fix; flagging rather than doing it silently.
- No test case exercises `deepTherml = bedrockLyrs`; it has been run only from a hand-built
  cold state in a scratch directory.

