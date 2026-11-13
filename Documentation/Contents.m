% FMA Hydrodynamics Library
% Version 5.0 (Stateful OOP Architecture)
%
% A Fourier Modal Analysis OOP framework for wave scattering, 
% hydrodynamics, and wave energy conversion.
%
% Base Classes
%   FMA                     - Shared wave states, physical constants, and robust matrix algebra.
%
% Geometry Classes
%   FMA_Breakwater_Geometry - Spatial discretization for bottom-mounted/submerged breakwaters.
%   FMA_WEC_Geometry        - Spatial discretization for surface-piercing WEC arrays.
%
% Breakwater Modules
%   FMA_Breakwater          - Master orchestration hub for breakwater arrays.
%   FMA_Breakwater_Coeffs   - Analytical integrals and effective matrix assemblies.
%   FMA_Breakwater_Solver   - Global cascade, scattering matrices, and hydrodynamics.
%
% WEC Modules
%   FMA_WEC                 - Master orchestration hub for surface-piercing WEC arrays.
%   FMA_WEC_Coeffs          - Analytical integrals, effective matrices, and radiation vectors.
%   FMA_WEC_Solver          - Global cascade, radiation/scattering, and WEC hydrodynamics.