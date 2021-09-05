## mixed-phase clouds analytical framework

The jupyter notebooks collected here implement an framework for anaytical treatment of mixed phase clouds,
as it was presented in a series of publications (co-)authored by A. Korolev.


### Publications

For the analytical treatment of glaciation time
- Korolev, A., Mazin, I.P., 2003. Supersaturation of water vapor in clouds. Journal of the Atmospheric Sciences 60, 2957–2974. https://doi.org/10.1175/1520-0469(2003)060¡2957:SOWVIC¿2.0.CO;2
- Korolev, A., Isaac, G., 2003. Phase transformation of mixed-phase clouds. Q.J.R. Meteorol. Soc. 129, 19–38. https://doi.org/10.1256/qj.01.203
- Korolev, A., 2008. Rates of phase transformations in mixed-phase clouds. Q.J.R. Meteorol. Soc. 134, 595–608. https://doi.org/10.1002/qj.230
- Pinsky, M., Khain, A., Korolev, A., 2014. Analytical Investigation of Glaciation Time in Mixed-Phase Adiabatic Cloud Volumes. Journal of the Atmospheric Sciences 71, 4143–4157. https://doi.org/10.1175/JAS-D-13-0359.1

For the growth model
- Chen, J.-P., Lamb, D., 1994. The theoretical basis for the parameterization of ice crystal habits: Growth by vapor deposition. Journal of Atmospheric Sciences 51, 1206–1222. https://doi.org/10.1175/1520-0469(1994)051<1206:TTBFTP>2.0.CO;2
- Sulia, K.J., Harrington, J.Y., 2011. Ice aspect ratio influences on mixed-phase clouds: Impacts on phase partitioning in parcel models. J. Geophys. Res. 116. https://doi.org/10.1029/2011JD016298
- Zhang, D., Wang, Z., Heymsfield, A., Fan, J., Luo, T., 2014. Ice Concentration Retrieval in Stratiform Mixed-Phase Clouds Using Cloud Radar Reflectivity Measurements and 1D Ice Growth Model Simulations. Journal of the Atmospheric Sciences 71, 3613–3635. https://doi.org/10.1175/JAS-D-13-0354.1


Single plots from the some of the publications are reproduced in the respective notebooks.



### Parameterizations

A couple of temperature-dependent parametrization are implemented for

- latent heat of condensation
- latent heat of sublimation
- water vapor diffusion in air
- thermal conductivity
- heat capacity of dry air
- heat capacity of water vapor
- saturation vapor pressure over water
- saturation vapor pressure over ice
- density of supercooled water
- dynamic viscosity of air

with a focus on the temperature intervall -40 to 0°C.
The can be accessed via the `pstore` variable. An overview with plots is provided in `parametrizations_overview.ipynb`

### Acknowledgements

The implementation of Korolev and Mazin (2003)  is based on an earlier version by Johannes Bühl.
Thanks to Alexei Korolev, sharing his original skript helped a lot identifying bugs and cross-checking this implementation.


### License
See the LICENSE file for more information.
Copyright 2021, Martin Radenz MIT License