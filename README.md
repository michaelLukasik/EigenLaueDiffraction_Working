# EigenLaueDiffraction_Working

A working directory for an xray diffraction simulation for periodic crystaline structures, ranging from simple minerals such as table salt to biological structures such as DNA. Originally created to produce [Laue Diffraction Patterns](https://en.wikipedia.org/wiki/Laue_equations), the framework and tools provided here aims to to simulate single crystal (or near single) diffraction for any crystaline structure. The project takes a quantum mechanical based approach to calculating the far (and less precicely, near) field diffraction patterns for vairous materials. Currently, the simulation is suited for crystals <100 units in size (or roughly <5000 atoms) at a resolution scale down to ~1nm at the imaging plane. 
<p align="center">
  <img width="460" height="300" src="https://github.com/michaelLukasik/EigenLaueDiffraction_Working/assets/138163589/df58a72c-e56b-48de-8039-1291df445a60">
</p>
<p align="center">
  A Laue Diffraction pattern captured by me at Brown University around 2019
</p>


## Objectives and Motivation behind the project
Typical approaches to opical-path tracing employ a fourier-optics approach to solve the wave equation via approximation to first order and using a fourier series to simulate "propogation" of the wave over some distance. 

This approach has three glaring flaws that are often overlooked in favour of simplicity of calculations: 
1) A complete loss of phase information at the point of direct measurement (plotting the intensity). This leads to the eponymously dubbed "phase problem" frequently dealt with in modern day cystolography.
2) The implicit "squashing" of the third spacial directions when approximating the radiating stucture as a 2-Dimensional Appeture through which the diffraction occurs. 
3) An abandonment of quantum phenomena in favor of a geometric interpretation of scattering within a crystal, which leads to a loss of information in quantum correction effects

By employing a physics-based framework for the simulation, these flaws are handeled within the calculations and the only approximation needed is handled by the user at run time (via terminating the partial wave analysis summation at an early index)

## Current Status and Immediate Goals

Currently, the program is capable of reporducing experimentally backed results for simple mineral lattices (crystaline materials falling into one of the 6 families of Bravais latticies), as well as special user-defined "unique" crystals that require a bit more finese in modeling (such as B-form DNA).

Some main goals of the project that are yet to be handeled fully, but are vital to expanding the usefulness of the simulation: 
1) **Speed up the process of light propogation** (This project was created with this roadbloack eternally in mind). This can be achieved in a few ways, the most promising is to employ a combination of quantum mechanic based scale factors (to account for physical corrections that will lead to the shifting/broadening of diffraction lines), while using FFT for the bulk of the transform from the diffraction to imaging planes.
  a) This will need to be done in "layers" in planes parallel the projection screen to preserve 3-Dimensional information in the diffraction pattern.
  b) Propogation of the light can be handeled via FFT, however overall intensity scales will be set by calculations gathered via QFT, bridging the two methods.
  c) Employing Bluestein's Method can allow for a "zoom" feature that will be beneficial for XRD applications. This will allow for special calculations for "zones of interest" for a closer study into the simulation's performance for certain crystaline structures. [For example, the "missing 4th layer" in Franklin's Photo 51 of B-Form DNA](https://scripts.iucr.org/cgi-bin/paper?S0365110X53001939) should be able to be studied in detail. [See this paper published in nature](https://www.nature.com/articles/s41377-020-00362-z) for more inspiration behind this thought.
3) Investigation of more "exoctic" dynamic diffraction effects, such as the [pendellösung effect](https://opg.optica.org/oe/fulltext.cfm?uri=oe-16-12-9097&id=163314). The perfect periodicity of the simulated crstals is actually quite detrimental in real experiments as it leads to relatively strong secondary effects. Usually, imperfect sections in a crystal average these effects towards the median of the diffraction intensity, however this simulation provides a perfect sandbox to study such effects
4) Addition of crystal impurities, to simulate more accurate real world senarios. These do not need to be large, just an extra addition of a Cystal::cellStructure tacked onto an existing structure. Real world extinction lengths limit these effects heavily in effective distance; running these tests on the current iteration of the simulation could be done before propogation optimization.


Some Preliminary plots are shown below. The orignial "goal" of the project was to accurately reproduce Laue Spots. These spots occur at precise angles and locations for a given crystal within a stereographic projection. By "growing" the crystal and using simulated band-pass filters and other neat tricks allowed by the convolution thereom we can isolate these peak determine their location. This gives us a way to calibrate the simulation to real world experiments. 

Shown below is an average output for two different crystals (olivine and salt) of the simulation and the effects proper filtering has on the pattern. The broad regions of intenisty represent a set of "solutions" to the wave equation rather than precise values determined by crystal geometry. These extra solutions are usually filtered out via large amounts of repetition in a real world crystal on the scale of ~10 microns (or thought of another way, higher sampling frequency in the cyrstal space leads to more precise measurements in diffraction space). However, our simulated crystals of ~10 units across are at a scale of about ~10-100 Angstroms, thus leading to a broadening of the peaks. 

![OlivinevsSodiumFilterTestLarge](https://github.com/michaelLukasik/EigenLaueDiffraction_Working/assets/138163589/005c3129-ea74-4888-84e7-805b9019f396)

After image processing and filtering the **heavily** broadened diffraction lines, pattern recognition is ran on a binary-thresholded representation of the diffraction pattern, and the pixel location of the simulated Laue Spots can be translated into real-space coordinates. These calculations can then be compared with experiment, as well as used as a nice self-check to make sure the geometry of the simulated crystal is accurately telegraphed by the Laue Pattern. Shown below is an example of a N<sub>x,y,z</sub> = 6 (total of 6<sup>3</sup> = 216 unit cells) crystal of Salt, desribed by its primative orientation, and the red circles are a result of OpenCV's blob detection reature tuned for a while circle on a dark background. 

![LaueDetectionKeypointsDark](https://github.com/michaelLukasik/EigenLaueDiffraction_Working/assets/138163589/773b8de9-aba9-4290-80dc-371da6213a45)

And overlaying the locations of the detected Laue spots (red circule) onto the direct observation at the projection screen within the simulation, we see nice agreement. The underyling structure revealing within the illumnated laue spots is merely an effect of aliasing due to improper sampling at run time. The sensitivity to thresholding provides a good analouge to the variable of time in real world experiments -- a single simulation run for a single orientation of a crystal can take ~10 seconds, real world experiments take on the order of days depending on the equipment. Increased false-positive rates increase with radial distance, this is adjustable via proper treshholding and is analogous to pre-maturely removing a photographic plate during the experiement. 

![LaueDetectionKeypoints](https://github.com/michaelLukasik/EigenLaueDiffraction_Working/assets/138163589/d4ce6cd2-858b-4a25-a47a-321c251c192c)

## More complicated crystaline structure

To be a crystal in a physical sense speaks only to the structure of an object defined by its highly predicably and periodic nature. In this sense, anything mathching this periodic condition can be treated (albeit with a bit more tact) as a crystal. This has been used heavily in the fields of chemsitry and biology, where XRay Diffraction and Crystalography data are cruitial to understanding chemical structure and function. Ear;y pioners of this thinking included Rosalin Franklin, who produced the now famous Photo 51, characterizing the double-helix shape of DNA using nothing but the principles explored earlier in this project. Thus, recreation of her results should be possible within the framework of the simulation. 

Realistically, there are some effects that appear within the simulation that make an exact 1:1 comparison with experimental results difficult, however proper data analysis can ease the discrepencies. The most glaring difficulty presented currently is the inability to scale the crystal to arbitrary sizes without sacrificing enormous emounts of time. This leads to difficulty in modeling the "congregate effect" of large-crystal systems; that is to say the physical edges of our simulated crystal impede on the "idealized" diffraction pattern that would be created in aggregate. To illustrate this point, we examine the simulated diffraction pattern of B-Form DNA, the most common form of DNA we all likely think of in our minds.





