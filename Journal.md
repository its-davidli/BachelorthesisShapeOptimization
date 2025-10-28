# Progress

## Goals

- Control position of topological defects in 2D
- Get Chiral LC director

## Questions (to Roland and Stephan)

- Behavior of derivatives of DG Elements, see 2nd point 10th october 
- Behavior of defect location? (see last point in Saturday, 11th October 2025)
- Spikes at the 3D Boundary
- How to shape optimize for only a local objective

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
