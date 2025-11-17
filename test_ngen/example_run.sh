#./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/gauge_01073000/settings/example_realization_config_w_summa_bmi_routing.json

./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/gauge_01073000/settings/example_realization_config_w_summa_bmi.json
python -m nwm_routing -V4 -f  ./test/data/routing/ngen_routing_config_unit_test.yaml

#./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./data/gauge_01073000/example_bmi_multi_realization_config_w_routing.json

rm test/data/routing/*.parquet
