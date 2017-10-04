# ElmerMAN
Asymptotic Numerical Method (ANM) developpements as ELMER FEM Modules

Based on several works we first decided to test ANM algorithms in the framework of ELMER FEM multi-physic software.
You will find an algorithm that performs steady bifurcation analysis and branch switching in the case of simple bifurcation.

In this project, we first provide ELMER Modules for steady bifurcation analysis in an "as is" state.
It is a first validation step toward a complete modular toolbox.

You will need ELMER V7 compiled with MUMPS direct Solver in order to use those modules.
The preprint article describe the methods, the algorithm, the implementation and test cases.
The code is provided with no support at that time.

The present Modules are based on those works:

 + Guevel, Y., Boutyour, H., & Cadou, J. M. (2011). Automatic detection and branch switching methods for steady bifurcation in fluid mechanics. Journal of Computational Physics, 230(9), 3614 – 3629.

 + Cochelin, B. & Medale, M. (2013). Power series analysis as a major breakthrough to improve the efficiency of Asymptotic Numerical Method in the vicinity of bifurcations. Journal of Computational Physics, 236, 594 – 607.

 + Jawadi, A., Boutyour, H., & Cadou, J. M. (2013). Asymptotic Numerical Method for steady flow of power-law fluids. Journal of Non-Newtonian Fluid Mechanics, 202, 22 – 31.


For further details on the methods please read:

 + Cochelin, B. (1994). A path-following technique via an asymptotic-numerical method. Computers and Structures, 53(5), 1181 – 1192.

 + Cadou, J. M., Potier-Ferry, M., Cochelin, B., & Damil, N. (2001). ANM for stationary Navier–Stokes equations and with Petrov–Galerkin formulation. International Journal for Numerical Methods in Engineering, 50(4), 825–845.

 + Cadou, J. M., Potier-Ferry, M., & Cochelin, B. (2006). A numerical method for the computation of bifurcation points in fluid mechanics. European Journal of Mechanics - B/Fluids, 25(2), 234 – 254.

+ Guevel, Y., Girault, G., & Cadou, J. M. (2014). Parametric analysis of steady bifurcations in 2d incompressible viscous flow with high order algorithm. Computers & Fluids, 100, 185 – 195.

+ Medale, M. & Cochelin, B. (2015). High performance computations of steady-state bifurcations in 3D incompressible fluid flows by Asymptotic Numerical Method. Journal of Computational Physics, 299, 581 – 596.



Hopf part to come in ELMER :

+ Brezillon, A., Girault, G., & Cadou, J. M. (2010). A numerical algorithm coupling a bifurcating indicator and a direct method for the computation of Hopf bifurcation points in fluid mechanics. Computers & Fluids, 39(7), 1226 – 1240.

+ Girault, G., Guevel, Y., & Cadou, J. M. (2012). An algorithm for the computation of multiple Hopf bifurcation points based on Padé approximants. International Journal for Numerical Methods in Fluids, 68(9), 1189–1206.

+ Heyman, J., Girault, G., Guevel, Y., Allery, C., Hamdouni, A., & Cadou, J. M. (2013). Computation of Hopf bifurcations coupling reduced order models and the Asymptotic Numerical Method. Computers & Fluids, 76, 73 – 85.




A cheap non linear solver:

+ Cadou, J. M. & Potier-Ferry, M. (2010). A solver combining reduced basis and convergence acceleration with applications to non-linear elasticity. International Journal for Numerical Methods in Biomedical Engineering, 26(12), 1604–1617.
