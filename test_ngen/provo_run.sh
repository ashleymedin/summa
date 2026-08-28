#./cmake_build/ngen ./extern/summa/summa/test_ngen/domain_provo/settings/gage-10154200_subset.gpkg '' ./extern/summa/summa/test_ngen/domain_provo/settings/gage-10154200_subset.gpkg'' ./extern/summa/summa/test_ngen/settings/provo_realization_config_w_summa_bmi_routing.json

./cmake_build/ngen ./extern/summa/summa/test_ngen/domain_provo/settings/gage-10154200_subset.gpkg '' ./extern/summa/summa/test_ngen/domain_provo/settings/gage-10154200_subset.gpkg '' ./extern/summa/summa/test_ngen/domain_provo/settings/provo_realization_config_w_summa_bmi.json
python -m nwm_routing -V4 -f  ./extern/summa/summa/test_ngen/domain_provo/settings/provo_routing.yaml

rm extern/summa/summa/test_ngen/domain_provo/simulations/*.parquet