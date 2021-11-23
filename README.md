# UNIFAC-Thermodinamics_II
A project of mine from the Thermodinamics II for Chemical Engineering course.

This project is from back in 2016, but I have never been more proud of a code of mine, and I will like to boast adding that I made created the whole algorithm and code in less than 6 hours.
Sadly the course was on spanish so some of the files will be on this language I´ll try to sumarize the problem to solve and translate the comments on the code.
The work was made in a group of four classmates, but it wasn't the only project that semester so I focused on this one and my friends made the other ones that didn't required code.
The program was coded on SciLab and there was no necesity for a graphic interface so none was done. First I'll sumarize the objectives of the project and then the work phases.

Summary (the whole instructions for the work can be found on the file "Work_Instructions.pdf" in spanish):

Imagine a separation system, on it you have variuos fluids with different boiling points (no azeotropes) and you wish to separe them via a flashing method. Obviously, in real life, you can't just extract each fliud on a 100% purity on a single step.
So the norm in the industry is to have a recycle line to return some of the destilled fluid (the one that goes out from above and is principally composed of the most volatile fluid) so the system can extract an even more pure substance when the stationary state is reached on the system (when things stop to vary that much).

The work has the following information: 

The system has a Feed line (F) with a ternary mixture with a composition z_i. My system, set by the profesor for our group, was 0.5 bencene, 0.2 toluene and 0.3 of 1-butanol. 
A set temperature for the separator of 100 °C (the fluid enters at that same temperature).
The Vapour exit (V) must have a composition y_i of 80% of the most volatile substance (in our case, bencene). And in some cases a Recycle line (R) from the upper exit to the feed with a proportion of R_R, defined as the flow of R/S, where S is the vapour exit line before the recicle line.
The R line should be considered as condensated before adding it to the Feed line.

Work phases:

a. Calculate V{y_i} using Raoult's law and the operating pressure as the mean of the bubble and dew temperatures. Without recycle line, R_R = 0.
b. Same as a. using R_R = 0.95.
c. Same as b. (with R_R = 0.95) but using a model for the calculation of the activity coeficients using the Wilson equation. The code should calculate Gamma 12 and Gamma 21.
d. Same as a. via a model for the calculation of the activity coeficients using the UNIFAC methodology (my favorite).

The files on this repository will include, at least, the work instructions PDF, the separate code for each work phase on a .sce file, a compiled code for the work on a single .sce file and the project written report in spanish.
