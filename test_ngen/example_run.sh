#./cmake_build/ngen ./data/gauge_01073000/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/example_realization_config_w_summa_bmi_routing.json
#rm test/data/routing/*.parquet

./cmake_build/ngen ./data/gauge_01073000/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/example_realization_config_w_summa_bmi.json

#./cmake_build/ngen ./data/gauge_01073000/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./data/gauge_01073000/example_bmi_multi_realization_config_w_routing.json