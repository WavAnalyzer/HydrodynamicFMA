================================================================================
FMA_LIBRARY USER MANUAL (V5.0)
PART 1: THE BASE CLASS (FMA)
================================================================================

1. DESCRIPTION AND INHERITANCE HIERARCHY
--------------------------------------------------------------------------------
The Fourier Modal Analysis (FMA) Library is a stateful Object-Oriented Programming (OOP) framework for MATLAB. It computes wave scattering, hydrodynamics, and wave energy conversion for periodic arrays of structures.

The "FMA" class is the master Base Class. It does not run simulations on its own. Instead, it acts as the mathematical and physical engine room. It stores the shared fluid parameters, numerical tuning settings, global wave states, and complex matrix algebra tools required by its Child classes.

Class Inheritance and Composition Diagram:

      [FMA] (Base Class - Core Math & Fluid Physics)
        |
        |---> [FMA_Breakwater] (Child - Bottom-mounted/submerged structures)
        |       |---> geometry (FMA_Breakwater_Geometry)
        |       |---> coeff    (FMA_Breakwater_Coeffs)
        |       |---> solver   (FMA_Breakwater_Solver)
        |
        |---> [FMA_WEC]        (Child - Surface-piercing arrays)
                |---> geometry (FMA_WEC_Geometry)
                |---> coeff    (FMA_WEC_Coeffs)
                |---> solver   (FMA_WEC_Solver)


2. PROPERTIES (CLASS ATTRIBUTES)
--------------------------------------------------------------------------------
When an FMA object (or its children) is created, it stores the following internal data. They can be accessed via dot notation (e.g., obj.lambda).

--- Wave and Environmental Properties ---
* lambda 
  Meaning: Incident wave length [meters].
  Data Type: Double (1x1)
  Range: > 0

* period 
  Meaning: Array periodicity or domain width (d_array) [meters].
  Data Type: Double (1x1)
  Range: > 0

* theta 
  Meaning: Incident wave angle [radians]. 
  Data Type: Double (1x1)
  Range: [0, pi/2). (Note: Use 0.0001 for exact normal incidence).

* h 
  Meaning: Global water depth [meters].
  Data Type: Double (1x1)
  Range: > 0

--- Numerical Tuning Parameters ---
* Num_ord 
  Meaning: Fourier truncation order.
  Data Type: Integer (1x1)
  Range: >= 0 (Total propagating modes evaluated = 2 * Num_ord + 1).

* Np 
  Meaning: Number of vertical evanescent modes to resolve.
  Data Type: Integer (1x1)
  Range: > 0 (Typically 2 to 5 for standard accuracy).

* damp 
  Meaning: Artificial numerical damping coefficient.
  Data Type: Double (1x1)
  Range: >= 0.0 (Defaults to 0.0).

* delta 
  Meaning: Tanh smoothing parameter to prevent Gibbs ringing at corners.
  Data Type: Double (1x1)
  Range: ~ 1e-6

* condno 
  Meaning: Condition number tolerance for robust SVD pseudo-inversion.
  Data Type: Double (1x1)
  Range: ~ 1e-8

--- Array and Structural Grid Parameters ---
* N_Slice 
  Meaning: Total number of discrete slices per structural body.
  Data Type: Integer (1x1)
  Range: >= 1

* N_grating 
  Meaning: Total number of array rows (gratings).
  Data Type: Integer (1x1)
  Range: >= 1

* Gap 
  Meaning: Center-to-center spacing distance between two structural rows [meters].
  Data Type: Double (1x1)
  Range: >= 0 (Use 0 for single-row arrays).

--- Global Wavenumber States (Dynamically Updated) ---
* k0: Free-space incident wavenumber (2 * pi / lambda). (Double, 1x1)
* kh: Damped vertical wavenumbers in the open water region. (Complex Array, 1 x (Np+1))
* kyiI: Transverse y-axis diffraction wavenumbers across diffraction orders. (Double Array, (2*Num_ord+1) x 1)


3. METHODS (FUNCTIONS) - CORE MATHEMATICS & MATRIX ALGEBRA
--------------------------------------------------------------------------------
These methods handle the heavy mathematical lifting. Since they belong to the base class, they are called using the instantiated object (e.g., obj.MethodName).

Note on Dimensions: 
- r = 2 * Num_ord represents the Fourier harmonic span.
- D_tot = Np * (2 * Num_ord + 1) represents the total dimension of the coupled evaluating modes.
- N_cases represents the number of simultaneous incident wave or radiation scenarios being evaluated.
- M and N represent arbitrary generic matrix dimensions for utility solvers.

[obj] = FMA(lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap)
> Function: Base class constructor. Initializes base physical constants, numerical tuning parameters, and structural grid parameters.
> Inputs:
  - lambda, period, theta, h: Wave and environmental parameters. (Double, 1x1)
  - Num_ord, Np, damp, delta, condno: Numerical tuning parameters.
  - N_Slice, N_grating, Gap: Optional structural grid parameters (nargin > 9).
> Outputs:
  - obj: Instantiated FMA base object.

[tG] = obj.FFTC(r, f, fshift, ef, e1_f)
> Function: Calculates the Fourier transform coefficients of a rectangular geometric profile using a 'tanh' smoothing function (N_fft = 65536) to prevent Gibbs ringing at sharp boundaries.
> Inputs:
  - r: Fourier truncation limit (2 * Num_ord). (Integer, 1x1)
  - f: Normalized width of the block (f2 - f1). (Double, 1x1)
  - fshift: Normalized center shift of the block ((f1 + f2) / 2). (Double, 1x1)
  - ef: Step height inside the block. (Complex/Double, 1x1)
  - e1_f: Step height outside the block. (Complex/Double, 1x1)
> Outputs:
  - tG: Complex Fourier coefficients. (Complex Array, 1 x (2r+1))

[T] = obj.TTMatrix(Q)
> Function: Constructs a block Toeplitz matrix from a 3D tensor of Fourier coefficients, used to map convolutions into matrix multiplications.
> Inputs:
  - Q: Tensor of Fourier coefficients. (Complex Array, Np x Np x (2r+1))
> Outputs:
  - T: 2D block Toeplitz matrix. (Complex Matrix, D_tot x D_tot)

[KXI] = obj.GetWavenumbers(alpha)
> Function: Calculates open water vertical/transverse wavenumbers and k-vectors, updates internal states (obj.kh, obj.k0, obj.kyiI), and returns the diagonal KXI matrix.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)
> Outputs:
  - KXI: Diagonal matrix of transverse x-axis wavenumbers. (Complex Diagonal Matrix, D_tot x D_tot)
  - Updates obj.kh, obj.k0, obj.kyiI internally.

[AiT] = obj.GetAiT(Ah, KXI)
> Function: Assembles the global orthogonal inner product matrix (AiT).
> Inputs:
  - Ah: Local open water orthogonal matrix. (Complex Matrix, Np x Np)
  - KXI: Diagonal transverse x-axis wavenumber matrix. (Complex Diagonal Matrix, D_tot x D_tot)
> Outputs:
  - AiT: Global matching matrix for velocity matching. (Complex Matrix, D_tot x D_tot)

[Q, V, W] = obj.EigenModes(kd, f1, f2)
> Function: Sets up and solves the general eigenvalue problem inside the periodic fluid domain using internal wave properties (k0, kh, kyiI).
> Inputs:
  - kd: Grating interior vertical wavenumbers. (Complex Array, 1 x (Np+1))
  - f1, f2: Slice filling factors. (Double, 1x1)
> Outputs:
  - Q: Eigenvalues / complex wavenumbers. (Complex Diagonal Matrix, D_tot x D_tot)
  - V: Velocity eigenvectors / amplitudes. (Complex Matrix, D_tot x D_tot)
  - W: Potential eigenvectors / amplitudes. (Complex Matrix, D_tot x D_tot)

[Qg_sorted, Wg_sorted, Vg_sorted] = obj.SortModes(Qg_raw, Wg_raw, Vg_raw, kd, total_dim)
> Function: Sorts the numerically derived eigenvalues/eigenvectors to strictly match the physical order of propagating and evanescent vertical modes.
> Inputs:
  - Qg_raw, Wg_raw, Vg_raw: Raw matrices from the eigenvalue solver. (Complex Matrices, total_dim x total_dim)
  - kd: Grating interior vertical wavenumbers. (Complex Array, 1 x (Np+1))
  - total_dim: Total dimension of the mode matrices (D_tot). (Integer, 1x1)
> Outputs:
  - Qg_sorted, Wg_sorted, Vg_sorted: Physically sorted eigenmodes. (Complex Matrices, total_dim x total_dim)

[X] = obj.RobustSolve(A, B)
> Function: A fail-safe linear system solver (X = A \ B). If the reciprocal condition number of A falls below the object's 'condno' tolerance, it safely falls back to an SVD-based pseudo-inverse.
> Inputs:
  - A: Left-hand side system matrix. (Complex Matrix, M x M)
  - B: Right-hand side matrix/vector. (Complex Matrix, M x N)
> Outputs:
  - X: Solved matrix/vector. (Complex Matrix, M x N)

[S] = obj.StarProduct(Sa, Sb)
> Function: Combines two global scattering matrices using the Redheffer Star Product to cascade adjacent geometric domains.
> Inputs:
  - Sa: Upstream scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)
  - Sb: Downstream scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)
> Outputs:
  - S: Combined global scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)

[R_in] = obj.ApplyReflection(Sg, R_load)
> Function: Applies a terminating boundary condition (reflection load) to a scattering matrix, collapsing it into a single reflection matrix.
> Inputs:
  - Sg: Local scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)
  - R_load: Downstream reflection matrix load. (Complex Matrix, D_tot x D_tot)
> Outputs:
  - R_in: Input reflection matrix. (Complex Matrix, D_tot x D_tot)

[R_in, Gen_out] = obj.CascadeLeftward(S, R_load, Gen_in, src_L, src_R)
> Function: Cascades a scattering matrix with a downstream (right) reflection load, while simultaneously propagating internally generated waves leftwards. Used primarily for calculating radiation potentials.
> Inputs:
  - S: Local scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)
  - R_load: Downstream (right) reflection matrix load. (Complex Matrix, D_tot x D_tot)
  - Gen_in: Incident generated wave amplitude from the right. (Complex Matrix, D_tot x N_cases)
  - src_L: Leftward source vector generated within the local domain. (Complex Matrix, D_tot x N_cases)
  - src_R: Rightward source vector generated within the local domain. (Complex Matrix, D_tot x N_cases)
> Outputs:
  - R_in: Combined input reflection matrix. (Complex Matrix, D_tot x D_tot)
  - Gen_out: Total leftward-propagating generated wave amplitude. (Complex Matrix, D_tot x N_cases)

[R_in, Gen_out] = obj.CascadeRightward(S, R_load, Gen_in, src_R, src_L)
> Function: Cascades a scattering matrix with an upstream (left) reflection load, while simultaneously propagating internally generated waves rightwards. Used primarily for calculating radiation potentials.
> Inputs:
  - S: Local scattering matrix. (Complex Matrix, 2*D_tot x 2*D_tot)
  - R_load: Upstream (left) reflection matrix load. (Complex Matrix, D_tot x D_tot)
  - Gen_in: Incident generated wave amplitude from the left. (Complex Matrix, D_tot x N_cases)
  - src_R: Rightward source vector generated within the local domain. (Complex Matrix, D_tot x N_cases)
  - src_L: Leftward source vector generated within the local domain. (Complex Matrix, D_tot x N_cases)
> Outputs:
  - R_in: Combined input reflection matrix. (Complex Matrix, D_tot x D_tot)
  - Gen_out: Total rightward-propagating generated wave amplitude. (Complex Matrix, D_tot x N_cases)

[mroots] = obj.DFS(alpha, N, h)
> Function: Calculates the roots of the standard water wave dispersion equation using a robust Homotopy tracking method.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)
  - N: Number of evanescent roots to find. (Integer, 1x1)
  - h: Water depth [m]. (Double, 1x1)
> Outputs:
  - mroots: Vertical wavenumbers (1 propagating, N evanescent). (Complex Array, 1 x (N+1))

--- Private Helper Methods (Access = private) ---
[out] = obj.oneroot(alpha, guess)
> Function: Newton-Raphson internal solver for converging a single wave dispersion root to a tolerance of 1e-9.

[mroot] = obj.homotopy(alpha, N)
> Function: Robust homotopy tracking method for calculating deep/intermediate water vertical roots from initial unit/imaginary guesses.

================================================================================
FMA_LIBRARY USER MANUAL 
PART 2: THE BREAKWATER CLASS (FMA_Breakwater)
================================================================================

1. DESCRIPTION AND COMPOSITION
--------------------------------------------------------------------------------
The "FMA_Breakwater" class inherits all the base physics and mathematical tools 
from the master "FMA" class. It is specifically tailored for evaluating wave 
scattering and vertical wave forces over bottom-mounted or submerged breakwater 
structures (where the structure's clearance is measured from the seabed up).

To keep the architecture clean and stateful, this class uses Composition. It 
acts as the orchestration hub and owns three "Grandchild" objects that handle 
the specific geometric and physical data:
  1. obj.geometry -> Handles the spatial discretization of the breakwater array.
  2. obj.coeff    -> Handles analytical hydrodynamic integrals and effective matrices.
  3. obj.solver   -> Handles the global scattering matrices and state outputs.


2. PROPERTIES (CLASS ATTRIBUTES)
--------------------------------------------------------------------------------
The stateful architecture stores variables internally rather than passing them 
continuously between functions. Below are the properties for the master class 
and its nested subclasses.

--- A. Master Class: FMA_Breakwater ---
* Inherits all base wave, numerical, and grid parameters from the FMA base class 
  (lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap).
* geometry: Instance pointer to the FMA_Breakwater_Geometry class.
* coeff: Instance pointer to the FMA_Breakwater_Coeffs class.
* solver: Instance pointer to the FMA_Breakwater_Solver class.
* kd: Wavenumbers of dispersion mode in the breakwater region. (Complex Array, 1 x (Np+1))

--- B. Geometry Subclass: FMA_Breakwater_Geometry ---
* N_Slice: Total number of discrete slices per structural body. (Integer, 1x1)
* f1: Vector of starting normalized filling factors for each slice. (Double Array, 1 x N_Slice)
* f2: Vector of ending normalized filling factors for each slice. (Double Array, 1 x N_Slice)
* d: Vector of structural clearances measured from the seabed to the structure [m]. (Double Array, 1 x N_Slice)
* L: Vector of discrete structural slice thicknesses in the x-direction [m]. (Double Array, 1 x N_Slice)

--- C. Coefficient Subclass: FMA_Breakwater_Coeffs ---
* Parent: Instance pointer back to the master FMA_Breakwater object.
* Ah: Local open water orthogonal matrix. (Complex Matrix, Np x Np)
* KXI: Diagonal matrix of transverse x-axis wavenumbers. (Complex Diagonal Matrix, D_tot x D_tot)
* AiT: Global matching matrix for velocity matching. (Complex Matrix, D_tot x D_tot)
* AdEff: Effective inner product tensor for the structural domain. (Complex Tensor, Np x Np x (2r+1))
* BdEff: Effective cross-coupling tensor matching open water to the structure. (Complex Tensor, Np x Np x (2r+1))
* BhEff: True transpose of the effective cross-coupling tensor. (Complex Tensor, Np x Np x (2r+1))

--- D. Solver Subclass: FMA_Breakwater_Solver ---
* Parent: Instance pointer back to the master FMA_Breakwater object.
* Sin_cyl: Input scattering block matrices per slice. (Cell Array, 1 x N_Slice)
* Sout_cyl: Output scattering block matrices per slice. (Cell Array, 1 x N_Slice)
* Q_cyl: Eigenvalues (wavenumbers) per slice. (Cell Array, 1 x N_Slice)
* W_cyl: Potential eigenvectors per slice. (Cell Array, 1 x N_Slice)
* X_cyl: Phase shift transmission matrices per slice. (Cell Array, 1 x N_Slice)
* kd_all: Local structural vertical wavenumbers per slice. (Cell Array, 1 x N_Slice)
* A_L: Global leftward-propagating internal wave amplitudes. (Cell Array, 1 x (N_blocks+1))
* A_R: Global rightward-propagating internal wave amplitudes. (Cell Array, 1 x (N_blocks+1))
* R: Global Reflection Matrix for the incident wave. (Complex Matrix, D_tot x 1)
* T: Global Transmission Matrix for the incident wave. (Complex Matrix, D_tot x 1)
* F_ex: Excitation force complex amplitudes for each row. (Complex Array, N_grating x 1)


3. CLASS CONSTRUCTOR & GEOMETRY METHODS
--------------------------------------------------------------------------------
To use this class, you must instantiate it with the specific physical and 
numerical environment parameters. 

[bw] = FMA_Breakwater(lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap)
> Function: Builds the Breakwater simulation environment and allocates internal states.
> Inputs: 
  - lambda: Incident wave length [m]. (Double, 1x1)
  - period: Array periodicity [m]. (Double, 1x1)
  - theta: Incident wave angle [radians]. (Double, 1x1)
  - h: Water depth [m]. (Double, 1x1)
  - Num_ord: Fourier truncation order. (Integer, 1x1)
  - Np: Number of vertical evanescent modes. (Integer, 1x1)
  - damp: Artificial numerical damping. (Double, 1x1)
  - delta: Tanh smoothing parameter. (Double, 1x1)
  - condno: SVD pseudo-inversion tolerance. (Double, 1x1)
  - N_Slice: Number of discrete slices per structural body. (Integer, 1x1)
  - N_grating: Number of structural rows in the array. (Integer, 1x1)
  - Gap: Center-to-center spacing between structural rows [m]. (Double, 1x1)
> Outputs:
  - bw: The instantiated Breakwater object.

bw.init_kd(alpha, d_draft)
> Function: Initiates the dispersion mode wavenumbers (bw.kd) in the breakwater region using obj.DFS(alpha, obj.Np, d_draft).
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)
  - d_draft (Optional): Structural clearance from the seabed [m]. Defaults to bw.geometry.d(1) if omitted. (Double, 1x1)

Once instantiated, define the structural geometry utilizing the nested geometry 
object. This automatically populates the structural vectors (f1, f2, d, L).

bw.geometry.SlicedCircle(P, r, dL, clearance_val)
> Function: Discretizes a circular cross-section into a series of rectangular staircase slices.
> Inputs:
  - P: Periodicity of the array [m]. (Double, 1x1)
  - r: Radius of the cylinder [m]. (Double, 1x1)
  - dL: Thickness of each discrete slice [m]. (Double, 1x1)
  - clearance_val: Clearance from the seabed to the structure [m]. (Double, 1x1)

bw.geometry.SlicedRectangle(P, width, total_length, N_slice, clearance_val)
> Function: Discretizes a rectangular cross-section into slices.
> Inputs:
  - P: Periodicity of the array [m]. (Double, 1x1)
  - width: Width of the rectangle in the y-direction (2a) [m]. (Double, 1x1)
  - total_length: Total length of the rectangle in the x-direction [m]. (Double, 1x1)
  - N_slice: Number of discrete slices. (Integer, 1x1)
  - clearance_val: Clearance from the seabed to the structure [m]. (Double, 1x1)


4. METHODS (FUNCTIONS) - HYDRODYNAMIC COEFFICIENTS
--------------------------------------------------------------------------------
These analytical integrals and matrix initializations resolve the depth matching 
conditions for submerged structures. They are handled internally by the solver 
but remain accessible via the coefficient handle (e.g., bw.coeff.Method).

obj.coeff.initGlobalMatrices()
> Function: Initializes the global open water orthogonal and wavenumber matching matrices (KXI, Ah, AiT) using the current parent wave state (kh, kyiI, h).

obj.coeff.initEffectiveMatrices(ind)
> Function: Assembles the 2D modal inner product and cross-coupling matrices (Ad, Bd, and true transpose Bh = Bd.') for slice 'ind', then applies FFTC across the slice filling factors to populate AdEff, BdEff, and BhEff.
> Inputs:
  - ind: Index of the discrete structural slice being evaluated. (Integer, 1x1)

[out] = obj.coeff.IntA(kh, h, d)
> Function: Evaluates the orthogonal inner product for the fluid domain located directly above the breakwater (clearance d) or full depth (h == d).
> Inputs:
  - kh: Vertical wavenumber in the open water region. (Complex, 1x1)
  - h: The water depth [m]. (Double, 1x1)
  - d: Clearance from the seabed [m]. (Double, 1x1)
> Outputs:
  - out: The integrated analytical coefficient. (Complex, 1x1)

[out] = obj.coeff.IntB(kh, kd, h, d)
> Function: Evaluates the cross-coupling integral matching the full depth (h) fluid domain to the reduced depth (d) fluid domain over the structure.
> Inputs:
  - kh: Vertical wavenumber in the open water region. (Complex, 1x1)
  - kd: Vertical wavenumber in the structural region. (Complex, 1x1)
  - h: The water depth [m]. (Double, 1x1)
  - d: Clearance from the seabed [m]. (Double, 1x1)
> Outputs:
  - out: The integrated analytical cross-coupling coefficient. (Complex, 1x1)


5. METHODS (FUNCTIONS) - SIMULATION SOLVERS
--------------------------------------------------------------------------------
These solvers evaluate the global array physics. The stateful OOP architecture 
means the solver requires minimal inputs, updating object properties dynamically.

Note on Dimensions: 
- D_tot = Np * (2 * Num_ord + 1) represents the total dimension of the coupled evaluating modes.
- N_blocks = 3 * N_Slice * N_grating + N_grating - 1 represents the total cascade blocks.
- N_x and N_y represent the number of points in the spatial evaluation grids.

obj.solver.clearStates()
> Function: Flushes high-memory cell arrays and matrices (Sin_cyl, Sout_cyl, X_cyl, Q_cyl, W_cyl, kd_all, A_L, A_R, R, T, F_ex) to prevent memory leaks during iterative frequency loops.

obj.solver.FMA_Stratify(alpha)
> Function: Main orchestrator for computing the complete scattering and hydrodynamic parameters for single or multi-row breakwater arrays. Calls GetWavenumbers, initGlobalMatrices, clearStates, SMatrices, FieldCoeffs, and Hydrodyms.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)
> Outputs:
  - Results are stored directly to obj.solver.R, obj.solver.T, obj.solver.F_ex, obj.solver.A_L, and obj.solver.A_R.

obj.solver.SMatrices(alpha)
> Function: Evaluates the local eigenmodes (Q_cyl, W_cyl), phase shift matrices (X_cyl), and input/output block scattering matrices (Sin_cyl, Sout_cyl) for each structural slice.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)

obj.solver.FieldCoeffs()
> Function: Solves the global matrix cascade block system (including inter-row gap blocks) via rightward and leftward sweeps to determine internal wave amplitudes (A_R, A_L) and global scattering matrices (R, T).

obj.solver.Hydrodyms()
> Function: Integrates internal wave fields against spatial slice geometries to evaluate the vertical excitation force complex amplitudes (obj.solver.F_ex) for each breakwater row.

[eta_amp] = obj.solver.WaveField(x_grid, y_grid)
> Function: Computes the normalized absolute free-surface elevation |eta| across a 2D spatial grid, accounting for scattering.
> Inputs:
  - x_grid: 1D array defining the spatial X mesh. (Double Array, 1 x N_x)
  - y_grid: 1D array defining the spatial Y mesh. (Double Array, 1 x N_y)
> Outputs:
  - eta_amp: Matrix of absolute free-surface elevations mapped to the input grid. (Double Matrix, N_y x N_x)

================================================================================
FMA_LIBRARY USER MANUAL 
PART 3: THE WEC CLASS (FMA_WEC)
================================================================================

1. DESCRIPTION AND COMPOSITION
--------------------------------------------------------------------------------
The "FMA_WEC" class inherits all the base physics and mathematical tools from 
the master "FMA" class. It is specifically tailored for evaluating wave 
scattering, hydrodynamics (added mass and radiation damping), and wave fields 
for surface-piercing Wave Energy Converter (WEC) arrays (where the structure's 
draft is measured from the free surface downwards).

To keep the architecture clean and stateful, this class uses Composition. It 
acts as the orchestration hub and owns three "Grandchild" objects that handle 
the specific geometric and physical data:
  1. obj.geometry -> Handles the spatial discretization of the WEC array.
  2. obj.coeff    -> Handles analytical draft and radiation integrals.
  3. obj.solver   -> Handles the coupled scattering and radiation matrices.


2. PROPERTIES (CLASS ATTRIBUTES)
--------------------------------------------------------------------------------
The stateful architecture stores variables internally rather than passing them 
continuously between functions. Below are the properties for the master class 
and its nested subclasses.

--- A. Master Class: FMA_WEC ---
* Inherits all base wave, numerical, and grid parameters from the FMA base class 
  (lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap).
* rho: Fluid density [kg/m^3]. (Double, 1x1)
* g: Gravitational acceleration [m/s^2]. (Double, 1x1)
* kd: Wavenumbers of dispersion mode in the WEC structural region. (Double Array, 1 x (Np+1))
* geometry: Instance pointer to the FMA_WEC_Geometry class.
* coeff: Instance pointer to the FMA_WEC_Coeffs class.
* solver: Instance pointer to the FMA_WEC_Solver class.

--- B. Geometry Subclass: FMA_WEC_Geometry ---
* N_Slice: Total number of discrete slices per structural body. (Integer, 1x1)
* f1: Vector of starting normalized filling factors for each slice. (Double Array, 1 x N_Slice)
* f2: Vector of ending normalized filling factors for each slice. (Double Array, 1 x N_Slice)
* d: Vector of structural drafts measured from the free surface down to the structure [m]. (Double Array, 1 x N_Slice)
* L: Vector of discrete structural slice thicknesses in the x-direction [m]. (Double Array, 1 x N_Slice)

--- C. Coefficient Subclass: FMA_WEC_Coeffs ---
* Parent: Instance pointer back to the master FMA_WEC object.
* Ah: Local open water orthogonal matrix. (Complex Matrix, Np x Np)
* KXI: Diagonal matrix of transverse x-axis wavenumbers. (Complex Diagonal Matrix, D_tot x D_tot)
* AiT: Global matching matrix for velocity matching. (Complex Matrix, D_tot x D_tot)
* AdEff: Effective inner product tensor for the structural domain. (Complex Tensor, Np x Np x (2r+1))
* BdEff: Effective cross-coupling tensor matching open water to the structure. (Complex Tensor, Np x Np x (2r+1))
* BhEff: True transpose of the effective cross-coupling tensor. (Complex Tensor, Np x Np x (2r+1))
* F1: Fourier-transformed psi-component radiation array. (Complex Array, D_tot x 1)
* F2: Fourier-transformed chi-component radiation array. (Complex Array, D_tot x 1)

--- D. Solver Subclass: FMA_WEC_Solver ---
* Parent: Instance pointer back to the master FMA_WEC object.
* Sin_cyl: Input scattering block matrices per slice. (Cell Array, 1 x N_Slice)
* Sout_cyl: Output scattering block matrices per slice. (Cell Array, 1 x N_Slice)
* Q_cyl: Eigenvalues (wavenumbers) per slice. (Cell Array, 1 x N_Slice)
* W_cyl: Potential eigenvectors per slice. (Cell Array, 1 x N_Slice)
* X_cyl: Phase shift transmission matrices per slice. (Cell Array, 1 x N_Slice)
* kd_all: Local structural vertical wavenumbers per slice. (Cell Array, 1 x N_Slice)
* src_in_cyl: Leftward radiation source potentials per slice. (Cell Array, 1 x N_Slice)
* src_out_cyl: Rightward radiation source potentials per slice. (Cell Array, 1 x N_Slice)
* A_L: Global leftward-propagating internal wave amplitudes. (Cell Array, 1 x (N_blocks+1))
* A_R: Global rightward-propagating internal wave amplitudes. (Cell Array, 1 x (N_blocks+1))
* R: Global Reflection Matrix (scattering + radiation cases). (Complex Matrix, D_tot x N_cases)
* T: Global Transmission Matrix (scattering + radiation cases). (Complex Matrix, D_tot x N_cases)
* F_ex: Excitation force complex amplitudes for each row. (Complex Array, N_grating x 1)
* Hydro_Matrix: Complex hydrodynamic matrix (Real = Added Mass, Imag = Damping). (Complex Matrix, N_grating x N_grating)


3. CLASS CONSTRUCTOR & GEOMETRY METHODS
--------------------------------------------------------------------------------
To use this class, you must instantiate it with the specific physical and 
numerical environment parameters. 

[wec] = FMA_WEC(lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap, rho, g)
> Function: Builds the WEC Array simulation environment and allocates internal states.
> Inputs: 
  - lambda, period, theta, h, Num_ord, Np, damp, delta, condno, N_Slice, N_grating, Gap: (Same as Breakwater Base Parameters).
  - rho: Fluid density [kg/m^3]. (Double, 1x1)
  - g: Gravitational acceleration [m/s^2]. (Double, 1x1)
> Outputs:
  - wec: The instantiated WEC Array object.

wec.init_kd(d_draft)
> Function: Initiates the dispersion mode wavenumbers (wec.kd = (0:Np) * pi / (h - d_draft)) in the fluid region beneath the WEC draft.
> Inputs:
  - d_draft (Optional): Structural draft from the free surface [m]. Defaults to wec.geometry.d(1) if omitted. (Double, 1x1)

Once instantiated, define the structural geometry utilizing the nested geometry 
object. This automatically populates the structural vectors (f1, f2, d, L).

wec.geometry.SlicedCircle(P, r, dL, draft_val)
> Function: Discretizes a circular cross-section into a series of rectangular staircase slices.
> Inputs:
  - P: Periodicity of the array [m]. (Double, 1x1)
  - r: Radius of the cylinder [m]. (Double, 1x1)
  - dL: Thickness of each discrete slice [m]. (Double, 1x1)
  - draft_val: Draft from the free surface to the structure [m]. (Double, 1x1)

wec.geometry.SlicedRectangle(P, width, total_length, N_slice, draft_val)
> Function: Discretizes a rectangular cross-section into slices.
> Inputs:
  - P: Periodicity of the array [m]. (Double, 1x1)
  - width: Width of the rectangle in the y-direction (2a) [m]. (Double, 1x1)
  - total_length: Total length of the rectangle in the x-direction [m]. (Double, 1x1)
  - N_slice: Number of discrete slices. (Integer, 1x1)
  - draft_val: Draft from the free surface to the structure [m]. (Double, 1x1)


4. METHODS (FUNCTIONS) - HYDRODYNAMIC COEFFICIENTS
--------------------------------------------------------------------------------
These analytical integrals and matrix initializations resolve the depth matching 
conditions and radiation potentials for surface-piercing structures. They are 
handled internally by the solver but remain accessible via the coefficient 
handle (e.g., wec.coeff.Method).

obj.coeff.initGlobalMatrices()
> Function: Initializes the global open water orthogonal and wavenumber matching matrices (KXI, Ah, AiT) using the current parent wave state.

obj.coeff.initEffectiveMatrices(ind)
> Function: Assembles the 3D effective modal tensors (AdEff, BdEff, BhEff) and the Fourier-transformed radiation forcing vectors (F1, F2) for slice 'ind'.
> Inputs:
  - ind: Index of the discrete structural slice being evaluated. (Integer, 1x1)

[out] = obj.coeff.IntAh(kh, h)
> Function: Evaluates the orthogonal inner product for the full water depth.
> Inputs: 
  - kh: Vertical wavenumber in the open water region. (Complex, 1x1)
  - h: The water depth [m]. (Double, 1x1)
> Outputs: 
  - out: The integrated analytical coefficient. (Complex, 1x1)

[out] = obj.coeff.IntAd(kd, h, d)
> Function: Evaluates the orthogonal inner product for the fluid domain located beneath the WEC (depth h-d).
> Inputs: 
  - kd: Vertical wavenumber in the structural grating region. (Complex, 1x1)
  - h: The water depth [m]. (Double, 1x1)
  - d: Draft of the structure from the free surface [m]. (Double, 1x1)
> Outputs: 
  - out: The integrated analytical coefficient. (Complex, 1x1)

[out] = obj.coeff.IntB(kh, kd, h, d)
> Function: Evaluates the cross-coupling integral matching the full depth (h) to the reduced depth (h-d) beneath the structure.
> Inputs: 
  - kh, kd: Open water and structural vertical wavenumbers. (Complex, 1x1)
  - h, d: Water depth and structural draft [m]. (Double, 1x1)
> Outputs: 
  - out: The integrated analytical cross-coupling coefficient. (Complex, 1x1)

[out] = obj.coeff.F11_psi(kd, h, d)
> Function: Evaluates the radiation forcing integral (translating body motion to fluid pressure).
> Inputs: 
  - kd, h, d: Structural wavenumber, water depth, structural draft.
> Outputs: 
  - out: The radiation forcing coefficient. (Complex, 1x1)

[out] = obj.coeff.F12_chi(kh, h, d, L)
[out] = obj.coeff.F12_psi(kd, h, d, L)
[out] = obj.coeff.F13_chi(kh, h, d, L)
> Functions: Evaluate the cross-radiation integral components across the slice half-length L.
> Inputs: 
  - kh / kd: Respective regional vertical wavenumbers. (Complex, 1x1)
  - h, d, L: Water depth, draft, slice half-thickness (dL/2) [m]. (Double, 1x1)
> Outputs: 
  - out: The integrated analytical coefficients. (Complex, 1x1)


5. METHODS (FUNCTIONS) - SIMULATION SOLVERS
--------------------------------------------------------------------------------
These solvers evaluate the global array physics, integrating scattering and 
radiation seamlessly.

Note on Dimensions: 
- D_tot = Np * (2 * Num_ord + 1) represents the total dimension of the coupled evaluating modes.
- N_cases = 1 + N_grating represents the number of simultaneous incident wave + radiation scenarios.
- N_x and N_y represent the number of points in the spatial evaluation grids.

obj.solver.clearStates()
> Function: Flushes high-memory cell arrays and matrices (Sin_cyl, Sout_cyl, X_cyl, Q_cyl, W_cyl, kd_all, src_in_cyl, src_out_cyl, A_L, A_R, R, T, F_ex, Hydro_Matrix) to prevent memory leaks during iterative frequency loops.

obj.solver.FMA_Stratify(alpha)
> Function: Main orchestrator that computes the complete scattering, radiation, and hydrodynamic parameters for single or multi-row WEC arrays. Calls GetWavenumbers, initGlobalMatrices, clearStates, SMatrices, FieldCoeffs, and Hydrodyms.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)
> Outputs (Stored Natively in Solver Properties):
  - obj.solver.R: Global Reflection Matrix. Column 1 is pure scattering; Columns 2:end are radiation potentials. (Complex Matrix, D_tot x N_cases)
  - obj.solver.T: Global Transmission Matrix. (Complex Matrix, D_tot x N_cases)
  - obj.solver.F_ex: Excitation force complex amplitudes for each row. (Complex Array, N_grating x 1)
  - obj.solver.Hydro_Matrix: Complex hydrodynamic matrix (Real = Added Mass, Imag = Damping). (Complex Matrix, N_grating x N_grating)
  - obj.solver.A_L, obj.solver.A_R: Internal wave amplitudes. (Cell Arrays of Complex Matrices, D_tot x N_cases)

obj.solver.SMatrices(alpha)
> Function: Evaluates the local eigenmodes, phase shift matrices, input/output scattering matrices (Sin_cyl, Sout_cyl), and radiation source vectors (src_in_cyl, src_out_cyl) for each structural slice.
> Inputs:
  - alpha: Dispersion parameter (k * tanh(kh)). (Double, 1x1)

obj.solver.FieldCoeffs()
> Function: Assembles the global multi-case cascade system (1 scattering case + N_grating radiation cases) and sweeps rightward and leftward to solve A_R, A_L, R, and T.

obj.solver.Hydrodyms()
> Function: Integrates the internal wave fields and analytical particular solutions across all slices to populate the excitation force vector (F_ex) and radiation hydrodynamic matrix (Hydro_Matrix).

[eta_amp] = obj.solver.WaveField(Xi, Rad_Multiplier, x_grid, y_grid)
> Function: Computes the normalized absolute free-surface elevation |eta| across a 2D spatial grid, accounting for combined scattering and radiation forces.
> Inputs:
  - Xi: Complex array of solved body motions (Response Amplitude Operators). (Complex Array, N_grating x 1)
  - Rad_Multiplier: Radiation multiplier scaling (omega^2 / g). (Double, 1x1)
  - x_grid: 1D array defining the spatial X mesh. (Double Array, 1 x N_x)
  - y_grid: 1D array defining the spatial Y mesh. (Double Array, 1 x N_y)
> Outputs:
  - eta_amp: Matrix of absolute free-surface elevations mapped to the input grid. (Double Matrix, N_y x N_x)
================================================================================
================================================================================