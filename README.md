# electromagnetism
Simulation of electromagnetic waves written with Julia.

Many of the files here are reduntant, as it is a live reflection of my progress in testing and refining the program.
The most useful are **App/2d_sim.jl** and **App/3d_sim.jl**, for runs in 2D and 3D respectively. I will document my general method here.

## Method
I start with Coulomb's law to define the electric field **E** based on the positions and charges of point charges while accounting for time delay from the speed of light. From this I calculate the curl of the electric field **∇xE**, and combine it with the Maxwell-Faraday equation (and approximating *dt* as Δt). In this way I approximate the change in the magnetic field **dB** as

**dB** = -(**∇xE**)Δt

The vector field **dB** is then added to the magnetic field **B** from the previous time step (assuming that **B** = 0 everywhere when t = 0)

Every time step a new **E** and **B** are calculated and the display is updated.

## Current Work
As I am attempting to implement the Ampere-Maxwell equation, I am running into problems, of which there are many. Here are a few I am working on now, with their possible solutions below.

**1)** Currently point charges have a pre-defined position as a function of time, regardless of what other charges are doing.

> Give point charges a mass and approximate F = ma for each charge to find their velocity at each time step, then use that to iterate their positions

**2)** Implementation of the Ampere-Maxwell equation

> This will probably require a major revamping of how I find the electric field. WIP

## Other Considerations
Here are other items I am looking at but not currently focused on.

**1)** There are better ways to approximate than Euler's method

**2)** Multithreading

**3)** Utilizing the GPU (Metal.jl)
