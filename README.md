[![DOI](https://zenodo.org/badge/638126735.svg)](https://zenodo.org/badge/latestdoi/638126735)

```
      &     (/
  &&   &&    &
    &&   &&   &&
&&   &&&  &&  &&/
    &&&( &&& &&& &&
      #%%%&&%&&%%%&&
            ((%%%%%%%%&&
        %%# %##%#((####%&%
          (%%##%%&%#########%%
          /%%###&%%%%#######&&%#
            %%####%%%####%##%&&&&%%
              /%%%%%%%#%#%#%#%#&&&&%&&%
                #%&&&%%&%#%#&&&&%#%&%&&&&/
                /%%%#%(##%%%#%&&&%%%##%#&&
                  /%%(##%&&&%%#&&&%%%###%%##
                  (%%%%%%%#%%#%&&&########((
                    %#&%#%%%%%%%&&%%&######((/
                    &&&&#(/%#%%&&%########((///
                      /%&&&%%%%%%&###%####(((#((
                        %&&&&%&%%&&&###(#(((((((((/
                        #&&&&&&%%%%%#((((((((((((//
                            #&&&&&&%%&&####((((//////,,,/ # .#%/&#
                              &&&&&&#(((((((//((#(  . ..#%%,   ##%#
                                ###%%%%%(///// //  ,,,,.,,/ //((/
                                  &&((/    /    ,,,   (((%##((#####(#/
                                /((  /  //( /////////((#(#################//#
                          (/#((((//((///(((((((((((#(#%%#&%&#&%&&%#%&%#%#######%%%%(
                      %%&%%%&&&//((###((((      %&%%%%%%%%%#%%###%%%###%%%%%%%&&&&%&&&
                  %%&&%&%%%#%/                        (#%%%%%#%###%###(##%%%#(%####&%%&&#
                                                            /##/#/(######(((###%%%####%%&%
                                                                        ((%%((####%#########%%
                                                                              #%# ##%(#%%%%%#  #%
                                                                                      &&% &&%    /(
```

Code for generating **O**ptimal **S**tability **P**olynomials in **R**oots for **E**xplicit Time **I**ntegration.

## Dependencies

* [`IpOpt`](https://github.com/coin-or/Ipopt) is the core package (optimizer). Following the [installation instructions](https://coin-or.github.io/Ipopt/INSTALL.html) should suffice, no special installation directory is required.
For the linear solver only [`MUMPS`](https://github.com/coin-or-tools/ThirdParty-Mumps) has been applied both in combination with [`METIS`](https://github.com/KarypisLab/METIS) and without.
* [`NAG dco/c++`](https://www.nag.com/content/downloads-dco-c-versions) is used to compute the necessary derivatives algorithmically. `dco/c++` is proprietary software, but chances are that you can obtain an academic license (`NAG Campus`) if you are working in research.
After obtaining `dco/c++` and licensing it, you need to change the path in the Makefiles (line 7) accordingly, i.e., `DCO_PATH=YOUR/PATH/TO/DCO`.
* _Optional_: If you want to compute also the monomial coefficients of the stability polynomial the usage of higher precision datatypes is necessary, where I resort to the implementation by [`Boost`](https://github.com/boostorg/multiprecision).

## Building
After obtaining and licensing the dependencies, execute in both directories `Feasibility_Problem` and `Optimization_Problem` 
```
make -j NUMTHREADS
```
where you can specify the `NUMTHREADS` according to your machine, e.g. `8`.
This builds object files and binaries in the corresponding directories `obj` and `bin`.

## Usage

Best starting point are the examples.
In order to carry out the optimization of a stability polynomial of degree $S$ and linear order of accuracy $p$ you require the spectrum and a reference timestep $\Delta t_\text{Ref}$ with a corresponding reference stage count $S_\text{Ref}$.
The syntax of calling the feasible point searchers/optimizers is then 
```
./Roots_Real(Imag).exe S p S_ref dt_ref Spectrum
```
If the spectrum itself does not form a convex hull, you need to supply the path to the files containing the real and imaginary part, respectively.
A call would then look like this:
```
./Roots_Real(Imag).exe S p S_ref dt_ref Spectrum PathToHullPoints
```
Again, this is best seen in the examples.

`Roots_Real.exe` looks for the parameter file `Roots_Real.opt` and `Roots_RealImag.exe` accordingly for `Roots_RealImag.opt` in the working directory.
If none of these files is present, default `Ipopt` options are used.

## Affiliation
![rwth_acom_en_rgb_cropped](https://github.com/DanielDoehring/OSPREI/assets/75639095/8de4aa2f-e2b6-463f-b2ce-a2652e546710)

This project was developed at the [Institute for Applied and Computational Mathematics (ACoM) at RWTH Aachen University](https://www.acom.rwth-aachen.de/).

## Acknowledgements

This code is a result of research performed in the research unit "Structure-Preserving Numerical Methods for Bulk- and
Interface Coupling of Heterogeneous Models (SNuBIC)"

![snubic-logo6](https://github.com/DanielDoehring/OSPREI/assets/75639095/a88000ca-d43c-4ac2-94d8-a0636fc7f32c)

This project has benefited from funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)
through the research unit FOR 5409 "Structure-Preserving Numerical Methods for Bulk- and
Interface Coupling of Heterogeneous Models (SNuBIC)" (project number 463312734).

![DFG](https://github.com/DanielDoehring/OSPREI/assets/75639095/a0ebb057-a0c5-4dc9-96df-03ac74894d2a)
