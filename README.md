# OSPREI

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

Code for generation of **O**ptimal **S**tabiliy **P**olynomials for **E**xplicit time **I**ntegration.

## Dependencies

* [`IpOpt`](https://github.com/coin-or/Ipopt) is the core pacakge (optimizer). Following the [installation instructions](https://coin-or.github.io/Ipopt/INSTALL.html) should suffice, no special installation directory is required.
For the linear solver [`MUMPS`](https://github.com/coin-or-tools/ThirdParty-Mumps) has been tested, both in combination with [`METIS`](https://github.com/KarypisLab/METIS) and without.
* [`dco/c++`](https://www.nag.com/content/downloads-dco-c-versions) is used to compute the necessary derivatives algorithmically. `dco/c++` is proprietary, but chances are that you can obtain an academic license.
After obtaining `dco/c++` and licensing it, you need to change the path in the Makefiles accordingly `DCO_PATH=YOUR/PATH/TO/DCO`.

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