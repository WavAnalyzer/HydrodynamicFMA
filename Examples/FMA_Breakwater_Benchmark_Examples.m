% Dissipasive Fourier Modal Analysis Model Validation (OOP Version)
% Benchmarked in 27/05/2026 using refactored stateful FMA library[cite: 7]
% To use the examples, uncomment the block you want to benchmark
clc; clear all; close all;

disp('===================================================');
disp('   FMA_Library Model Validations (Breakwater)      ');
disp('===================================================');

%% Case 1: A G Abul-Azm Ocean Engineering 1994, 21(7), 683-706. Figure 7-10
% disp('Running Case 1.1 (Normal Incidence, L=0.5h)...');
% h = 1;
% kh_vec = 0.1:0.1:6;
% N = 5; 
% f1 = 0;
% f2 = 1;
% Num_ord = 0;
% L = 0.5 * h;
% S = 1;
% d = 0.5 * h;
% theta = 0;
% 
% delta = 1e-6; 
% condno = 1e-8; 
% damp = 0;
% Gap=0;
% N_slice=1;
% N_grating=1;
% 
% R1 = zeros(size(kh_vec));
% T1 = zeros(size(kh_vec));
% F1 = zeros(size(kh_vec));
% 
% for i = 1:length(kh_vec)
%     k = kh_vec(i) / h;
%     lambda = 2 * pi / k;
%     alpha = k * tanh(k * h);
%     period = 2 * lambda;
%     
%     % 1. Instantiate the FMA_Breakwater Object with grid/environment variables
%     bw = FMA_Breakwater(lambda, period, theta, h, Num_ord, N, damp, delta, condno, N_slice, N_grating, Gap);
%     
%     % 2. Assign structural geometry directly to the geometry instance
%     bw.geometry.N_Slice = N_slice;
%     bw.geometry.f1 = f1;
%     bw.geometry.f2 = f2;
%     bw.geometry.d = d;
%     bw.geometry.L = L;
%     
%     % 3. Run the internal solver (no outputs needed, state is stored internally)
%     bw.solver.FMA_Stratify(alpha);
%     
%     % 4. Extract results directly from the solver's properties
%     R1(i) = bw.solver.R(Num_ord+1);
%     T1(i) = bw.solver.T(Num_ord+1);
%     F1(i) = -1i * bw.solver.F_ex / (h * period);
%     fprintf('----------kh = %.3f ---------\n', kh_vec(i));
% end
% 
% figure('Name', 'Case 1.1', 'Position', [100, 100, 900, 400]);
% subplot(1,2,1);
% plot(kh_vec, abs(R1), 'r', 'LineWidth', 1.5); hold on;
% plot(kh_vec, abs(T1), 'k', 'LineWidth', 1.5);
% title('Case 1.1: Magnitudes'); xlabel('kh'); ylabel('|R|, |T|'); legend('|R|', '|T|');
% xlim([0 6]); ylim([0 1.2]); grid on;
% 
% subplot(1,2,2);
% plot(kh_vec, phase(R1) + 2*pi, 'r', 'LineWidth', 1.5); hold on;
% plot(kh_vec, phase(T1) - 2*pi, 'k', 'LineWidth', 1.5);
% title('Case 1.1: Phases'); xlabel('kh'); ylabel('Phase'); legend('Phase(R)', 'Phase(T)');
% xlim([0 6]); ylim([-6.5 9.5]); grid on;
% figure;
% subplot(1,2,1)
% plot(kh_vec,abs(F1),'b');
% xlim([0 6]);
% subplot(1,2,2)
% plot(kh_vec,phase(F1),'k');

%% Case 2: Bragg test (Il Nuovo Cimento C 13, 823-833, 1990)
% disp('Running Case 2 (Bragg Test - ResoSmatrix_Bragg)...');
% h = 1;
% D = 0.16 * h;
% Lt = 32 * h;
% ff = 1 / 8.0;
% theta = 0;
% Ng = 10;
% Np = 5;
% f1=0;f2=1;
% 
% Lp = Lt / Ng;
% Lg = Lp * ff;
% d = 1 - pi * D / 2;
% Sg = 1;
% dL = Lg / Sg;
% L=ones(1,Sg)*dL;
% Lh = Lp * (1 - ff);
% f1=ones(1,Sg)*f1;
% f2=ones(1,Sg)*f2;
% Num_ord=0;
% 
% kappa = 0.5 : 0.02 : 2.52;
% delta = 1e-8;
% condno = 1e-8;
% damp = 0;
% Gap = Lp; % Array spacing parameter
% 
% R4 = zeros(size(kappa));
% A4 = zeros(size(kappa));
% 
% for i = 1:length(kappa)
%     k0 = kappa(i) * Ng * pi / Lt;
%     lambda = 2 * pi / k0;
%     period = 2 * lambda;
%     alpha = k0 * tanh(k0 * h);
%     
%     % Initialize object and geometry vectors
%     bw = FMA_Breakwater(lambda, period, theta, h, Num_ord, Np, damp, delta, condno, Sg, Ng, Gap);
%     bw.geometry.N_Slice = Sg;
%     bw.geometry.f1 = f1;
%     bw.geometry.f2 = f2;
%     bw.geometry.L = L;
%     bw.geometry.d = ones(1,Sg)*d; % vectorized clearance
%     
%     bw.solver.FMA_Stratify(alpha);
%     
%     R_val = bw.solver.R(Num_ord+1);
%     R4(i) = abs(R_val);
%     A4(i) = angle(R_val);
%     
%     fprintf('---------- kappa = %.2f ---------\n', kappa(i));
% end
% 
% figure('Name', 'Case 2', 'Position', [250, 250, 900, 400]);
% subplot(1,2,1);
% plot(kappa, R4, 'k', 'LineWidth', 1.5);
% title('Case 2: Bragg Test Magnitude'); xlabel('\kappa'); ylabel('|R|');
% xlim([0.5 2.5]); ylim([0 1.1]); grid on;
% 
% subplot(1,2,2);
% plot(kappa, A4, 'r', 'LineWidth', 1.5);
% title('Case 2: Bragg Test Phase'); xlabel('\kappa'); ylabel('Phase(R)');
% xlim([0.5 2.5]); grid on;

%% Case 3: Hydrodynamic BIC- Single-row breakwater grating
clc; clear; close all;
fprintf('Step 1: High-Resolution Sweep to Pinpoint the Fano Resonance...\n');

% =========================================================================
% 1. Global Parameters
% =========================================================================
period = 1.0;            
L = 0.2;                 
L_total = L;             
f1_vec = -0.38;            
f2_vec = 0.38;  
N_slice = 1;             
N_grating = 1;           
Gap = 0;                 

h = 1.0;                 
d_clearance = 0.1;       

Np = 2; Num_ord = 15; delta = 1e-4; damp = 0.0; condno = 1e-10;          
G = 2 * pi / period;

ky_norm = 0.05; 
KY = ky_norm * G; 

% =========================================================================
% 2. Ultra-High Resolution Sweep (Optimized for OOP)
% =========================================================================
N_k = 1000;
k_norm_vec = linspace(0.8995, 0.9005, N_k);
R_mag = zeros(1, N_k);

% PRE-ALLOCATE: Instantiate the FMA object ONCE outside the loop to prevent overhead[cite: 7]
k_inc_init = k_norm_vec(1) * G;
lambda_init = 2 * pi / k_inc_init;
theta_init = asin(KY / k_inc_init);

bw = FMA_Breakwater(lambda_init, period, theta_init, h, Num_ord, Np, damp, delta, condno, N_slice, N_grating, Gap);
bw.geometry.N_Slice = N_slice;
bw.geometry.f1 = f1_vec;
bw.geometry.f2 = f2_vec;
bw.geometry.L = L;
bw.geometry.d = d_clearance;

for i = 1:N_k
    k_norm = k_norm_vec(i);
    k_inc = k_norm * G;
    
    % Dynamically update the state properties rather than rebuilding the object
    bw.lambda = 2 * pi / k_inc;
    bw.theta = asin(KY / k_inc);
    alpha = k_inc * tanh(k_inc * h);
    
    bw.solver.FMA_Stratify(alpha);       
    R_mag(i) = abs(bw.solver.R(Num_ord + 1));
    
    if mod(i, 200) == 0
        fprintf('Find the quasi-BIC %d/%d\n', i, N_k);
    end
end

[~, peak_idx] = max(R_mag);
k_peak_norm = k_norm_vec(peak_idx);
fprintf('-> Fano Peak found exactly at k/G = %.6f\n', k_peak_norm);

% =========================================================================
% 3. Calculate 2D Near-Field Map at the Exact Peak
% =========================================================================
fprintf('\nStep 2: Calculating 2D Near-Field Map at the Exact Peak...\n');

k_inc = k_peak_norm * G;
bw.lambda = 2 * pi / k_inc;
bw.theta = asin(KY / k_inc);
alpha = k_inc * tanh(k_inc * h);

% Evaluate fields internally
bw.solver.FMA_Stratify(alpha);     

x_grid = linspace(-1.0, 1.0, 400);
y_grid = linspace(-1 * period, 1 * period, 400);

% Direct call leveraging internal state arrays 
eta_amp = bw.solver.WaveField(x_grid, y_grid);

% Fix for original script's missing variable error: calculate max_val, max_x, max_y[cite: 7]
max_val = max(eta_amp(:));
[r_idx, c_idx] = find(eta_amp == max_val, 1);
max_x = x_grid(c_idx);
max_y = y_grid(r_idx);

% =========================================================================
% 4. Visualization & Hotspot Detection
% =========================================================================
fprintf('Rendering Plots...\n');
figure('Color', 'w', 'Position', [100, 150, 1000, 450], 'Name', 'Quasi-BIC Hotspot');

ax1 = subplot(1, 2, 1);
plot(k_norm_vec, R_mag, 'b-', 'LineWidth', 2); hold on;
plot(k_peak_norm, R_mag(peak_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('High-Resolution Fano Profile', 'FontSize', 14);
xlabel('Frequency (k / G)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|R| (Reflection)', 'FontSize', 12, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 11);

ax2 = subplot(1, 2, 2);
[X, Y] = meshgrid(x_grid, y_grid);
pcolor(X, Y, eta_amp); shading interp;
colormap(jet); 
c = colorbar; 
ylabel(c, 'Enhancement |u|/|u_{inc}|', 'FontSize', 12, 'FontWeight', 'bold');

hold on;
patch([-L_total/2, L_total/2, L_total/2, -L_total/2], [f1_vec, f1_vec, f2_vec, f2_vec], ...
      [0.6 0.6 0.6], 'EdgeColor', 'w', 'LineWidth', 1.5);
xlabel('x / Period', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y / Period', 'FontSize', 12, 'FontWeight', 'bold');
axis equal; 
xlim([min(x_grid) max(x_grid)]); 
ylim([min(y_grid) max(y_grid)]);
set(gca, 'FontSize', 11, 'Layer', 'top');

fprintf('\n=== HOTSPOT DETECTED ===\n');
fprintf('Maximum Field Enhancement: %.2f x\n', max_val);
fprintf('Exact Coordinate: x = %.3f, y = %.3f\n', max_x, max_y);

disp('===================================================');
disp('           All Validations Complete!               ');
disp('===================================================');