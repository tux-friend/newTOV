# newTOV
solver for stellar equations with tabulated EOS, switch between Newtonian or TOV equations

## Description
Creates a 1D profile for a star with given EoS table. This 1D profile then can be used for e.g. SPH simulations.
It can be switched between Newtonian stellar equations and TOV equations. For solving the ODE system the 4th order Runge-Kutta method is applied.

The EOS driver routine of Evan O'Connor et al. is used (https://github.com/evanoconnor/EOSdriver).

The form of the TOV equations used with the baryon density can be found in https://arxiv.org/abs/1306.4034 by Kaplan et al.
EOS tables are available on https://stellarcollapse.org/equationofstate

The code is build upon better understanding the equations and process of building neutron star models. For sure it can be made much faster 
and efficient. Feel free to commit and pull requests...
