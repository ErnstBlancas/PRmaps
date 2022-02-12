# PRmaps
## Input files
The input files are:
  - POSCAR with the cell of the compound
  - FORCE_CONSTANTS with the [Phonopy FORCE_CONSTANTS format](https://phonopy.github.io/phonopy/input-files.html#format-of-force-constants) (.hdf5 format is not acepted)
  - BTE.qpoints
  - BTE.w_final
 
## Creating output file with the results
The eigenvectors generated with Phonopy are used to obtain the participation rates and the atomic participation rates. To generate this file the wrapper has the mode write `-m w`. 
 - The size of the supercell is specified with `--supercell_size x y z`  
 - The density of the grid is specified with `--mesh_grid x y z` and is saved as `mesh.yaml`. The mesh grid must be equal or smaller than the BTE mesh use
 - The output file name could be set with `--output_yaml name`

For example to create a [5, 5, 3] supercell with a [12, 12, 12] mesh grid, and save the results in `BiOCuSe_results.yaml` the code is:
```
python main.py -m w --supercell_size 5 5 3 --mesh_grid 12 12 12 --output_yaml BiOCuSe_results
```


## Generating the results
The .yaml file with the results could be read with the mode read `m -r`. With this mode, different plots could be produced. 
 - The different modes could be selected with `-p [m, s or n]`
 - The input file is specified with `--output_yaml name`
 - The size of the heatmap is defined with `--heatmap_xy x y`.
 - The range of the colorbar could be set with:
   - `--robust True` the range is fixed according to the maximum and the minimum
   - `--robust False` in this case the range must be specified with `--heatmap_range min max`
 - The interpolation method use for smooth the heatmap could be specified with `--interpolation_method method_name`. The methods that could be use are from [Matplotlib `.inshow(interpolation)`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html)
 - The colormap used in the plots could be changed with `--cmap name`. The colour schemes available are from [Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)

 
### Multiplot mode
With the multiplot mode `-p m` one single .pdf file with all the heatmaps. The name of the file could be changed with `--multiplot_name name`. For example:
```
python main.py -m r -p m --output_yaml BiOCuSe_results --heatmap_xy 15 15 --robust True --cmap Blues --interpolation_method bicubic --multiplot_name BiOCuSe_multiplot
```
creates `BiOCuSe_multiplot.pdf` with a [15,15] heatmap, with the range fix to the Total maximum and minimum and using the bicubic interpolation method. To obtain the same figure but with a different range (0 to 7000):

```
python main.py -m r -p m --output_yaml BiOCuSe_results --heatmap_xy 15 15 --robust False --heatmap_range 0 7000 --cmap Blues --interpolation_method bicubic --multiplot_name BiOCuSe_multiplot
```

### Singleplot mode
With the single plot mode `-p s` a .pdf file is produced for each element in the compound and with the total distribution. All the previous parameters work except `--multiplot_name`. For example:
```
python main.py -m r -p s --output_yaml BiOCuSe_results --heatmap_xy 15 15 --robust True --cmap Blues --interpolation_method bicubic
```
the command generates 5 .pdf files with the name of the compound and the element (for example `BiOCuSe_Cu.pdf`) and one file with the total distribution (`BiOCuSe_Total.pdf`).

### No plot mode
With the no plot mode `-p n` the data used for the plots is saved as `raw_data`. Only the heatmap size and the input file are requested:
```
python main.py -m r -p n --output_yaml BiOCuSe_results --heatmap_xy 15 15
```
explicar datos que se produce?


