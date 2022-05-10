# HOI
Retrieving high-order information multiplets from data using the O-information

Work in progress.

It uses the gaussian copula estimator (https://github.com/robince/gcmi). Other estimators are possible of course (at the moment the fast bootstrap works for the gaussian copula/covariance based ones).

The main files are

*hoi_exhaustive_loop_zerolag.m
*hoi_exhaustive_loop_lagged.m

which find the significant redundant and synergistic multiplets up to the desired order.

The python version https://github.com/PranavMahajan25/HOI_toolbox does not contain neither FDR correction, nor sped-up bootstrap. 

Tentative todo list

* update the python repo
* plotting (already started by @renzocom)
* currently we select the first nbest multiplets. I suspect that this number (a couple of dozens) is very small compared to the whole set of significant multiplets, given the nmber of combinations, so we should actually wnantify how many multiplets are in the tail of the distribution, and fix nbest accordingly (the execution time scales linearly). Update 02.12.2021 I have modified the main function to save all the values, and added the plot_hist_nbest.m function to visualize these numbers.

# References
Rosas, F. E., Mediano, P. A., Gastpar, M., & Jensen, H. J. (2019). Quantifying high-order interdependencies via multivariate extensions of the mutual information. Physical Review E, 100(3), 032305. https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.032305

Stramaglia, S., Scagliarini, T., Daniels, B. C., & Marinazzo, D. (2021). Quantifying dynamical high-order interdependencies from the o-information: an application to neural spiking dynamics. Frontiers in Physiology, 11, 1784. https://www.frontiersin.org/articles/10.3389/fphys.2020.595736/full
