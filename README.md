## Description

These Matlab files are simulate the agent-based model as described in:
"Spatio-temporal modelling of human type 1 diabetes," Kyle C. A.
Wedgwood, Sarah J. Richardson, Noel G. Morgan and Krasimira
Tsaneva-Atanasova, <i> Frontiers in Physiology</i>, 2017.

## Platform
The code has been developed for Matlab R2014b and tested on Matlab
R2016a and comes with no guarantee that it will work on earlier
releases.

## Running the code
Simply run the Driver.m script to start a simulation. Various
model-based parameters can be set in this file. Note that spatial
parameters have been scaled by a factor of 4 to increase code speed.

## Files
<b> Driver.m </b> - Run the full code. Also contains model parameter settings.

<b> IBM.m </b>    - Main function file defining the full SDE model for cell
movement and interactions.

<b> modelForcesAllCellTypesFaster.m </b> - Function file defining forces on
immune cells arising from each other, the islet membrane and the beta
cells.

<b> modelForcesBetaCells.m </b> - Function file defining forces on beta cells
during initial placement (before full agent-based model starts).

<b> plot*.m </b> - Contains plotting routines to generate frames.

<b> plotResults.m </b> - Script to plot figures as in paper.

<b> plotIBMTimeCourse.m </b> - Script to produce immune cell invasion profiles as in
paper.

## Output
The code produces a subfolder named Anim which itself has a subdirectory
that can be named in the parameters structure. This subfolder will
contains frames of the simulation at specific time points (this can be
adjusted in the parameters structure also).

Frames can be put together into an animation using freely available
software (e.g. <i> ffmpeg </i>.

## Distribution
This code is distributed under the GNU General Licence V3.0. Users are
free to share, use and modify at will.
