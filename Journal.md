# Progress

## Goals

- Control position of topological defects in 2D
- Get Chiral LC director

## Questions (to Roland and Stephan)

- Spikes at the 3D Boundary
- Behavior of defect location? (see last point in Saturday, 11th October 2025 and Tuesday, 16th December 2025)
- Initial Stepsize determination
- Criterion to reject Stepsize because of intersecting mesh cells

## Things to discuss with Santiago

- In [Instabilty Hedgehog](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.59.563) and [scaling analysis](https://arxiv.org/pdf/1512.08164) they talk about the elastic parameter and how it affects the Ring Defect e.g., for our ellipse to circle example, this might be interesting. Also in the second paper, they divide by R^3 or the volume of the system. Should we do this too? In out examples, volume change was never a problem, so this should be fine

## Log book

### Thursday, 9th October 2025

- Read *LangtangenFEM2019* chapter 4,5,6 for better understanding of FEM


### Friday, 10th October 2025

- Organized the GitHub
- Analyzed the behavior of derivatives of DG Elements
    - Using DG1 and DG0 elements the DOFs give senseful non zero values, altough expected to be 0? (TestDG0.py)
    - project() of derivative of DG0 Function does not work? 
    - In the LC-Code, inner(grad(Q), grad(Q)) gives (unexpected) non-zero values
- Started Implementation of the target function of LC defects at arbitrary positions and top. charge
    - Use the Frank Oseen approx from *CorteseDefects2018* eq. (3) and (7) for local description of defects
    - First Approach: define small circles around the defects as subdomains, and impose target objective on those subdomains 

### Saturday, 11th October 2025

- Fixed a bug, where in normalizing the Eigenvectors to find the directorfield a division by 0 could happen
- Inital numerical experiment to find a target shape, varying width of rectangle to see influence of the position of the defects. 
- Analyzing the defects in the numerical solution for the Rectangle (*ResultsThesis/4.3/LCRectangle*), there seems to be two defects of +1/2 charge 
    - Changing the aspect ratio changes the location of the defects (see *Results/Test10to1*, *Results/Test5to1*, *Results/Test2to1*)
    - Making the mesh anisotropic changes the location of the defects, they locate at the edge of the fine mesh --> Indicate unphysical reason for the location of the defects? (see *Results/TestAnisotropicBig*, *Results/TestAnisotropicSmall*)
    - Different Meshsizes change the position of the defect --> Indicate unphysical reason for the location of the defects? (see *Results/Test5to1EvenFinerMesh*, *Results/Test5to1FinerMesh*,*Results/Test5to1*)

### Tuesday, 14th October 2025 / Wednesday, 22th October 2025

- Generalized saving of results to 3D
- Fixed a bug in compute_orientation() using high quadrature degree *form_compiler_parameters={'quadrature_degree': 2}*
- Created OneSidedWallBox.geo for a 3D Test: 3D Box with corner in (0,0,0) and (1,1,5), and one long surface (***corners (0,0,0) and (0,1,5)***) as a physical group. The goal is to get a twist of the box to achieve a pseudo chiral nematic
- Bug: center of mass calculation only 2D
- **TODO**: Rewrite all objective Functions 
- Added distinguished *ds* measures for boundary surfaces and movable surfaces

The first test does not yield good results. The shape derivative looks fine, but the resulting deformation is not as expected. --> New approach instead of only using one anchoring wall, also prescribe the opposite wall

### Monday, 27th October 2025

- Added the opposite wall as an anchoring surface
- Getting different results, with combinations of using *H1* or Elasticity inner product, tangential smoothing or normal projections
- Best result is found with the elasticity inner product, no tangential smoothing, with normal projection. (see Results/TestChiral)
- **Issue:** Spikes at the boundary of the mesh, which is bad, tangentual smoothing does not do much
- Next goal is to go back to look at the defects in 2d, try some cases to control the defect. Perspectivly also fix 2 boundaries and see what the other 2 boundaries do

### Tuesday, 28th October 2025

- **Problem**: How to prescribe the Q-Tensor for defects? We only have a local description of the defects.


### Meeting: Thursday, 30th October 2025

- Look at absolute error in sol of State equ
- Target funtion is not a background mesh, but attached to the mesh
    - Solution: Project after every mesh movement
    - **Declare the target function with spatial coordinates as an UFL expression**
        - X= SpatialCoordinate(mesh) *Declareing its dependant on the mesh*
        - UFL Expression
- Revisit the 2D case
- project cuts of the dependency
- First track a spatial constant field

- Create the geometry, which is a solution. -> 

### Friday, 31st October 2025

- Removed the projection() from the definition of the target functions
    - In 2D this results in a much better solution, see *Results/TestEllipseOptimal*, almost a perfect circle
- Created a twisted box to analyze the intuitive solution to the chiral problem see *Results/TestIrreg*
    - There are a few numerical issues with the forward solution, first some unregular directors are visible and also the solution does not seem to be consistent, as you can see in the debugging file from the 0th to first iteration.
- Tried an easier testcase of a homogenous director field, with an azimuthal angle, trying to achieve rotation of the box (prelim result *Results/TestRotateBox*), good beginning but many numerical issues.
- Done **TODO:** Check all instances of the different *ds* 
- **TODO:** Implment inital alpha from HerzogMeshQuality
- **TODO:** Run the programm to rotate with really small initialstepsize for 100 iterations.

### Monday, 3rd November 2025

- Implemented the initial stepsize calculation method proposed in *HerzogMeshQuality*
    - The stepsizes are too small
- Problem: There are many rigid movements it seems
    - Tuning the delta parameter higher somehow adds rigid movements

### Meeting: Thursday, 6th November 2025

**TODO:**
- Analyze the circle case in more detail: Try different initial shapes
- Prescribe only part of the geometry as target: start with half of the circle
- The issue of spikes at interfaces of anchoring surfaces: Try to smooth the transition 

### Friday, 7th November 2025

- Trying different Ellipses
- **Problem:** Meshing wide Ellipse does not work (e.g. radius 5 to 1)
- Trying to implement perturbations of the initial domain 

### Week of 10-14th November 2025

- Implemented representation of Shape derivative
- Tried to round the corners of the 3D Box --> No improvement
- Tried to smooth out the transition from BC to no BC --> No improvment
- Added Rectangle with only 2 Sides as BC, still same issue with spikes at the boundary of BC    to no BC

### Meeting, Wednesday 12th of November 2025

- For the 2D case, try to fix the non BC boundaries 

### Tuesday, 18th of November 2025

- First succesfull run of the 2D Two Sided BC Test Case (see Results/TestPrelim2WallsWithCut and Results/TestPrelim2WallsNoCut)
    - Modified the objective function by a normalization of the area/volume of the integration area
    - Added BCs on the shape gradient, so that the walls which fixes the non anchoring walls **TODO** To make this implementation and the implementation of the ds measures more elegent, make it necessary to put the non anchoring walls into a physical group at mesh generation
    - Added a regularize method, to smooth out spikes of the shape derivative at the interface of anchoring to non anchoring **TODO** Make the regularization function better (cleaner in a mathematical way, e.g. analyze the histogram of shape derivativevalues)
- Tried only integrating on a sub area: No good results but **TODO** explore that more, to possible analyze the different components of the shape derivative and where they come from
    - Looking at the shape derivative you can see the effect of trying to shrink the integration area at the interface of the integration to non-integration, showing that the shape derivative wants to use non physical information to decrease the objective functional. 

- **TODO** Standardize the physical group numbering

### Friday, 21st of November 2025

- Successfull run with a 3D Ellipsioid deforming into a sphere (see Results/TestEllipsiod3D)
    - Initial Stepsize need to be large enough to deform the sharp ends into round caps
    - Added the PlotGeometricalInformation method to plot radii etc.
    - The variance of the radius goes down significantly
    - The center of mass has some unwanted movement
    - The Beltrami Smoothing needed to be turned up in order to avoid spikes

### Thursday, 27th of November 2025

- Succesfull run with the rigid rotation test case
    - First tried with the standard LET Parameter (*E = 1, nu = 0.4*): Good Results, but the volume is getting compressed (see *Results/TestRigidRotationNormalElasticity*)
    - Then tried to subtract the compression part of the strain tensor, but that did not change much (see *Results/TestRigidRotationOnlyShear*)
    - Then increased the Youngs Modulus from *1* to *100* to make the volume less compressable, that yielded very good results (see *Results/TestRigidRotationGoodResult*) Note that in order to remove the rigid motions, one also has to increase the delta paremeter
    - **Question:** Does the subtraction of the trace do anything in LET Theory?
    - **Question:** Here we used the target Alignment as initial guess, and one can see, that the Energy cost ist not strong enough to align the LCs in the middle

- Revisited the pseudo chiral LC Case:
    - Used Same parameters as in the prior test case. Results look promising (see *Results/TestPseudoChiralPrelim*), further investigations needed, e.g.:
        - Fix the bottom boundary
        - Different LET Parameters and deltas
        - Center the box around the origin
        - Rewrite the target function without the Expression (Already implemented, need to run it as next step)
    

### Tuesday, 2nd December 2025

- Succesfull pseudo chiral LC Case, see *Results/TestChiralSuccess*:
    - With the Expression Method in the target function
    - Both deltas set to 0
    - Changed the forword solver, adding a relaxation step

### Thursday, 4th December 2025

- To improve the result from Tuesday, ran it with delta = 0.2, which makes the result better, since it removes rigid movements better see *Results/TestChiralSuccessWithDelta*
- Started to work on Saturn Ring defect
    - Found analytical expression in *AlamaSaturnDefect*
    - Tested forward solution: Looks alright
    - **TODO:** Implement objective

### Monday, 15th Decemember 2025

- New *Results/TestRigidRotationVeryGoodResult*, where I removed the Laplace Beltrami smoothing, which is has better results.
- Started with the Saturn Defect
    - In 2D everything is running, but there is a problem in the forward solution, since it introduces additional defects. Needs more analyzing

### Tuesday, 16th December 2025

- Got a good result for the saturn ring defect *Results/TestSaturnGood*
    - The forward solution looks good. The initialization of the forward solver is important! Otherwise there might be other defects introduces, see the 0-th iteration of *Results/TestSaturnMultDefects*
    - In the iteration 9 the stepsize is too large, such that the mesh intercepts itself, which makes the iterations after it invalid.
    - Mesh size is really important! Take care of the Choose of *L_c*, which influences the coherence length and thus the minimum mesh size 

### Friday, 26th December 2025

- Major Code Changes
- **Collecting Final Results in *PreliminaryResults/***
    - 2D *EllipseToCircle/*
    - 2D *Rectangle2Sides/*
- In 3D the S0 parameter still is unclear??? 

### Sunday, 28th December 2025

- Tried to reproduce the Rigid Rotation result with new code, a bug appeared where the rotation is exactly the wrong direction
- Found the Bug, it is due to the use of the normal projection, on the corners that is not well defined and causes this issue
- Also because of Mesh Size, we increased L_c (Choosen the value *2e-1* arbitraraly, important that it is high enough!)
- **Collected Final Results in *PreliminaryResults/***
    - 3D *RigidRotation/*
    - 3D *RigidRotationShowCaseLowYoungsModulus* Same paramters, just showing that role the Young's Modulus plays
    - 3D *RigidRotationShowCaseProjectionGoneWrong/* Showcase of the Bug mentioned above
    - 3D *RigidRotationShowCaseHighLc/* Showcase of why the higher L_c: Newtonsolver needs longer to converge and alignment not secured (but normally achieved through the initial guess)

### Monday, 29th December 2025

- **Collected Final Results in *PreliminaryResults/***
    - 3D *Chiral/* 
- Tried the 3D Ellipoid into circle for a little bit, but there is always an issue of meshquality and compression.
- Re-Done the Armijo Linesearch reset of alpha after each iteration (small mistake before) and the relative error calculation (also small mistake before)