# SUMMA Recipes

This page contains instructions for certain common actions one might do when working with the SUMMA source code.

## Add a variable

To add a new variable, some infrastructure must be set up before the variable can be used. This all happens inside `./source/dshare/`:

1. Add a new entry to `var_lookup.f90` under an appropriate category for your variable. For example, `scalarSaturatedArea` fits well under "Soil hydrology" in "(7) Diagnostic variables". Make not of the structure you've added the variable under. In this case, this would be `iLook_diag`.
2. Further down in `var_lookup.f90`, increment the correct data structure (here: `iLook_diag`) by 1. For example, if the current structure ends at 110, add the value 111.
3. Add a corresponding entry in `popMetadat.f90`. This defines how the variable will described in SUMMA's output files. Keep the ordering consistent with `var_lookup.f90`.
4. Add a corresponding entry to the correct data structure function in `get_ixname.f90`.

If everything has gone well SUMMA should compile and run as normal. If not, SUMMA will exit the run early and provide diagnostics. Common errors are:

* Data structure in `var_lookup.f90` not incremented (step 2):

```
FATAL ERROR: summa_initialize/summa_defineGlobalData/checkStruc/problem with structure constructor iLookDIAG [element=111]
```

* `popMetadat.f90` not updated (step 3):

```
FATAL ERROR: summa_initialize/summa_defineGlobalData/checkStruc/checkPopulated/diag_meta structure is not populated for named variable # 75 in structure iLookDIAG
```

* `get_ixname.f90` not updated (step 4):

```
FATAL ERROR: summa_initialize/summa_defineGlobalData/checkStruc/checkPopulated/get_ixUnknown/variable scalarSaturatedArea is not found in any structure
```