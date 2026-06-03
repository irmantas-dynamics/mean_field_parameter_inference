Needs staticly compiled liblsoda library:
https://github.com/sdwfrost/liblsoda

To compile static library replace Makefile  in `liblsoda-master/src` with the one that is included in `Supplementary_note_2/liblsoda_make` directory.

To compile library `losslib.c` for fast loss function calculation in `Makefile` document:
* change the variable `SRC_DIR` that shows path to the `liblsoda` source code directory;
* change the variable `LSODA_DIR` that shows path to  the statically compiled `liblsoda` library (can be set equal to `SRC_DIR`);


Parameter inference is done using time series in folder `time_series` precalculated by matlab function `QIF_AD/Noninvasive/QIF_AD_micro_time_series.m` for finite networks and `QIF_AD/Noninvasive/QIF_AD_macro_time_series_series.m` for infinite networks case.

Please copy these time series from folder `QIF_AD/Noninvasive/time_series`.

* `QIF_AD_DiffEvolution.py` - python script that estimates parameters using  Differential Evolution algorithm from SciPy package;

* `pymoo_CMAES_C.py` - python script that estimates parameters using  Covariance matrix adaptation evolution strategy (CMA-ES)  algorithm from pymoo package;

* `pymoo_PSO_C.py` - python script that estimates parameters using  Particle swarm optimization (PSO)  algorithm from pymoo package;


* Results are polotted by matlab script `results_plot.m`
