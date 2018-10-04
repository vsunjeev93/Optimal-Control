# Sensitvity Analysis for Runaway
This repository consists of code written during my first doctoral project. 

The objective is to obtain crticial parameters for runaway at steady state and this is done by applying parametric sensitvity analysis (PSA) using the direct differential method. Since the equations are highly unstable a time relexation approach is used to solve the boundary value problem, i.e the corresponding unsteady state equations are solved for a very long time until steady state is reached. 

A finite difference discretization is done in the spatial domain. The resulting discretized ordinary differential equations are stiff and require stiff solvers. One such solver is ode15s in MATLAB 

About the model:-
The model consists of 3 semi-linear PDE. They consist of one dimensional mass, fluid phase temeprature and wall temperature conservation equations. 

