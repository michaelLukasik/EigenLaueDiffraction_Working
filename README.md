# EigenLaueDiffraction_Working

A working directory for my xray-diffraction simulation for periodic crystaline structures, ranging from simple minerals such as table salt to biological structures such as DNA. 
The project takes a quantum mechanical based approach to calculating the far (and less precicely, near) field diffraction patterns for vairous materials. Currently, the simulation is suited for crystals <100 units in size (or roughly <5000 atoms) at a resolution scale of ~1nm at the diffraction plane. 

Typical approaches to opical-path tracing employ a fourier-optics approach to solve the wave equation via approximation to first order and using a fourier series to simulate "propogation" of wave over some distance. 

This approach has two glaring flaws that are often overlooked in favour of simplicity of calculations: 
1) The implicit "squashing" of the third spacial directions when approximating the radiating stucture as a 2-Dimensional Appeture through wich the diffraction occurs
2) An abandonment of quantum phenomena in favor of a geometric interpretation of scattering within a crystal, which leads to a loss of information in quantum correction effects

By employing a physics-based framework for the simulation, these major flaws are handeled within the calculations and no approximations are needed.

Currently, the programs is capable of reporducing experimentally backed results for simple mineral lattices (crystaline materials falling into one of the 6 families of Bravais latticies), as well as special user-defined "unique" crystals that require a bit more finese in modeling (such as B-form DNA) 

Some main goals of the project that are yet to be handeled fully, but are vital to propogating the usefulness of the simulation: 
1) Speed up the process of light propogation. This can be achieved in a few ways, the most promising is to employ a combination of quantum mechanic based scale factors (to account for physical corrections that will lead to the shifting/broadening of diffraction lines), while using FFT for the bulk of the transform
  a) This will need to be done in "layers" in planes parallel the projections screen to preserve 3 Dimensional information in the diffraction pattern
  b) Propogation of the light can be handeled via FFT, however overal intensity scales will be set by calculations gathered via QFT, briding the two methods
2) Investigation of more "exoctic" dynamic diffraction effects. The perfect periodicity of the simulated crstals is actually quite detrimental in real experiments, as it leads to relatively strong secondary effects. Usually, imperfect sections in a crystal average these effects towards the median of the diffraction intensity, however this simulation provides a perfect sandbox to study such effects
3) Addition of crystal impurities, to simulate more 
