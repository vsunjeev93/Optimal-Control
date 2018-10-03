# Sensitvity Analysis for Runaway
This repository consists of code written during my first doctoral project. 

The objective is to obtain crticial parameters for runaway and this is doen using parameteric sensitvity analysis (PSA) using the direct differential method. Since the equations are highly unstable a time relexation approach is used to solve the boundary value problem, i.e the corresponding unsteady state equations are solved for a very long time until steady state is reached. The a finite difference discretization is done in the spatial domain. The resulting discretized ordinary differential equations can be solved used standard ODE techniques (eg. Runge-Kutta methods). 

About the model:-
The model consists of 3 semi-linear PDE. They consist of mass fluid phase temeprature and wall temperature conservation equations in one dimensional equations. 

