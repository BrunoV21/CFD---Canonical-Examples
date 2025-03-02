3 files are presented in this directory each one targeted towards a specific regime, subsonic, transonic, and supersonic.

At the beginning of each file, the user can change the reference values for the physics of the problem.

A structured grid is used where the cell size is kept constant and equal between the x and y directions. 

Under parameters set up the user can change the dimensions of the volume, a 0.75 relation to x is recommended by the code for the y coordinate but this can be changed if the user wishes to. 

Cell size, the x coordinate where the airfoil starts, Mach number, alpha(relaxation coefficient), and the number of iterations can all be set in this section.

Once the code has been run, the time taken for setting up and solving the algebraic system of equations is displayed in the terminal. 

Plots of potential, u and v velocity, and density fields are then made available to the user, as well as a plot with the 'Residuals' which represent the evolution of the sum of the modulus of the difference of the value calculated in the current and previous iteration across all cells of the domain.


The directory No_Uniform_Grid_Spacing contains 2 files, the Mesh file was meant to be a script for structured non-uniform mesh generation. A multiblock approach was used where growth ratios, initial spacings, and maximum lengths. The version presented also includes a uniform grid with equal spacings that was being used as troubleshooting for the FiniteVolumes_Supersonic_V3. The coefficients for the interior points that were also included in the report appear to be properly set, however, the boundary coefficients were not forced properly and so the code does not converge. Unfortunately the tight time schedule did not help with solving the issue. The authors did however feel it was productive to include this code in the files sent for further learning.