# SUMMA Output Files

<a id="outfile_file_formats"></a>
## Output file formats
All SUMMA output files are in [NetCDF format](SUMMA_input#infile_format_nc).

<a id="outfile_dimensions"></a>
## Output file dimensions
SUMMA output files can have the following dimensions (as defined in `build/source/netcdf/def_output.f90`). Dimensions may be present even in output files where they are not actually used. Most of these dimensions are pretty self-explanatory, except perhaps the `[mid|ifc][Snow|Soil|Toto]andTime` dimensions, which combine depth and time information as a work-around for the lack of support for variable-length or ragged arrays in earlier versions of NetCDF-4. While SUMMA will likely move to support these variable-length arrays in its output, we currently use a different organization. The dimensions indicated by `ifc` are associated with variables that are specified at the interfaces between layers including the very top and bottom. For example, the flux into or out of a layer would be arranged along an `ifc` dimension. The dimensions indicated by `mid` are associated with variables that are specified at the mid-point of each layer (or layer-average). `Snow`, `Soil`, and `Toto` indicate snow layers, soil layers, and all layers, respectively. This is explained in detail in the [model history file](#outfile_history) section.

| Dimension | long name | notes |
|-----------|-----------|-------|
| gru       | dimension for the GRUs | Variables and parameters that vary by GRU |
| hru       | dimension for the HRUs | Variables and parameters that vary by HRU |
| depth     | dimension for soil depth | Variables and parameters that are defined for a fixed number of layers |
| scalar    | dimension for scalar variables | Scalar variables and parameters (degenerate dimension) |
| spectral  | dimension for the number of spectral bands | Variables and parameters that vary for different spectral regimes |
| time      | dimension for the time step | Time-varying variables and parameters |
| tdh       | dimension for the time delay routing vectors | Variables and parameters that are held in memory as part of routing routines |
| midSnow   | dimension for midSnow | Time-varying variables and parameters at the mid-point of each snow layer |
| midSoil   | dimension for midSoil | Time-varying variables and parameters at the mid-point of each soil layer |
| midToto   | dimension for midToto | Time-varying variables and parameters at the mid-point of each layer in the combined soil and snow profile |
| ifcSnow   | dimension for ifcSnow | Time-varying variables and parameters at the interfaces between snow layers (including top and bottom) |
| ifcSoil   | dimension for ifcSoil | Time-varying variables and parameters at the interfaces between soil layers (including top and bottom) |
| ifcToto   | dimension for ifcToto | Time-varying variables and parameters at the interfaces between all layers in the combined soil and snow profile (including top and bottom) |

<a id="outfile_restart"></a>
## Restart or state file
A SUMMA restart file is in [NetCDF forma](SUMMA_input#infile_format_nc) and is written by `build/source/netcdf/modelwrite.f90:writeRestart()`. This file is also an input file because it specifies the initial conditions at the start of a model simulation. It is described in more detail in the [SUMMA input](SUMMA_input#infile_initial_conditions) documentation. Note that when the file is written, the time for which it is valid is included as part of the model file name.

<a id="outfile_history"></a>
## Model history files
SUMMA history files are in [NetCDF format](SUMMA_input#infile_format_nc) and describe the time evolution of SUMMA variables and parameters. The files are written by the `writeParam`, `writeData`, and `writeTime` subroutines in `build/source/netcdf/modelwrite.f90`. SUMMA output is pretty flexible. You can output many time-varying model variables and parameters, including summary statistics. You can specify what you want to output in the  [output control file](SUMMA_input#infile_output_control), which is one of SUMMA's required input files.
