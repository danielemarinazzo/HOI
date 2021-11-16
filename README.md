# HOI
Retrieving high-order information multiplets from data using the O-information

Work in progress.

The files starting with hoi_* (basically the zero lag ones) are fully optimized.

The lagged version is fully functional (including the proper FDR correction), but is not optimized for speed yet.

The python version https://github.com/PranavMahajan25/HOI_toolbox does not contain neither FDR correction, nor speed-up with bootstrap. The equivalent matlab files in this repo are the ones starting with old_* which in turn use o_information_boot.m and o_information_lagged_boot.m (this latter is still the current one for the lagged version).

# Tentative todo list

* update the python repo
* optimize the lagged version
* plotting (already started by @renzocom)
* currently we select the first nbest multiplets. I suspect that this number (a couple of dozens) is very small compared to the whole set of significant multiplets, given the nmber of combinations, so we should actually wnantify how many multiplets are in the tail of the distribution, and fix nbest accordingly (the execution time scales linearly).

# References
Rosas, F. E., Mediano, P. A., Gastpar, M., & Jensen, H. J. (2019). Quantifying high-order interdependencies via multivariate extensions of the mutual information. Physical Review E, 100(3), 032305. https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.032305

Stramaglia, S., Scagliarini, T., Daniels, B. C., & Marinazzo, D. (2021). Quantifying dynamical high-order interdependencies from the o-information: an application to neural spiking dynamics. Frontiers in Physiology, 11, 1784. https://www.frontiersin.org/articles/10.3389/fphys.2020.595736/full
