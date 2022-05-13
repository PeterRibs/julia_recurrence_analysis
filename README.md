# Recurrence analysis

In the field, I recorded the dance of a specie named Swallow-tailed Manakin. After tracking the videos using Python and obtaining the time series, I used the JULIA language to perform the Recurrence analysis. Its objective is to perceive patterns in a time series.

Recurrence Quantification Analysis is made upon the recurrence matrix (or recurrence plots) which is a square matrix of binary states (0 or 1). The rows and columns have the length of the time series and each cell (i,j) of this matrix assumes the value 1 if the i-th term of the time series has the same value (with a precision of + r, a parameter of the analysis) of the j-th term, otherwise, it assumes 0.

We extracted six metrics of the RQA:
- Recurrence rate (RR): the percentage of cells in the recurrence matrix that are recurrent, which indicates the number of different heights that are recurrent in a time series. 
- Determinism (DET): the percentage of cells in the recurrence matrix that form diagonal lines, indicating the percentage of the temporal series that has parts of sequential data repeated over the series. Indicates how repetitive the display is. 
- Diagonal entropy (DLEntr): the Shannon entropy of these diagonal lines, indicating the length variation of the repeated sequential data. Detects subtle variations of how this predictability occurs. 
- Laminarity (LAM): the percentage of cells in the recurrence matrix that form vertical lines (sequential data with no variation). Indicates the percentage of time series composed of stationary moments, that is, how long do males remain stationary. 
- Microstates entropy (MCEntr) measures the Shannon entropy of sub-matrices of size 2x2, which are randomly sampled from the recurrence matrix.

## Tecnology:

- Julia
- Jupyter Notebook

## Packages:

* `DelimitedFiles`
* `DifferentialEquations`
* `DynamicalSystems`
* `LaTeXStrings`
* `Random`
* `Statistics`
* `Sundials`




