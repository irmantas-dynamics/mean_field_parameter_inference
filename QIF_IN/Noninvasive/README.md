1. Time series generation is performed by matlab script (output to `time_series` folder):
   * `QIF_IN_micro_time_series.m` - from the microscopic equations;
   * `QIF_IN_macro_time_series.m` - from macroscopic equations;

2. Optimization is performed by python script (ouputs to `optimization_results` folder):
    * `QIF_IN_DiffEvolution_noninv.py`

3. Hidden variables prediction (Fig. 3 in manuscript), from infered parameters is performed by matlab script:
    * `QIF_IN_time_series_prediction.m`
