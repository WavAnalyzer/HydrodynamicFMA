clc; clear; close all;
%% Case 1: Benchmark Ning et al. (2017) 2D case
% % Effect of Structure Drafts (d1/h = d2/h = 0.08, 0.125, 0.15)
% % Validated Sign Rules using stateful OOP FMA library[cite: 11]
% 
% h = 10.0;            
% Bg = 0.2 * h;        % Fixed clear spacing D/h = 0.2
% a_width = 0.6 * h;   % Fixed breadth a1/h = a2/h = 0.6 (This is L)
% 
% theta = 0.0001; Np = 2; Num_ord = 0; f1 = 0; f2 = 1; 
% delta = 1e-6; condno = 1e-8; damp = 0.0;         
% 
% rho = 1000; g = 9.81; 
% kh_list = 0.1 : 0.1 : 10.0; 
% draft_ratios = [0.08, 0.125, 0.15]; 
% N_slice = 1;
% 
% K_R_data = zeros(length(draft_ratios), length(kh_list));
% K_T_data = zeros(length(draft_ratios), length(kh_list));
% CWR_data = zeros(length(draft_ratios), length(kh_list));
% 
% fprintf('Computing Fig. 4 Benchmark (Draft Variations)...\n');
% 
% % =========================================================================
% % PRE-ALLOCATE OOP OBJECTS OUTSIDE THE LOOPS TO PREVENT INSTANTIATION OVERHEAD
% % =========================================================================
% % Use dummy lambda/period/draft for initial setup
% dummy_lambda = 100; dummy_period = 500; dummy_draft = 0.1 * h;
% 
% % 1. Object for Coupled Array (N_grating = 2)
% Gap_coupled = Bg + a_width; % Center-to-center spacing[cite: 11]
% wec_coupled = FMA_WEC(dummy_lambda, dummy_period, theta, h, Num_ord, Np, damp, delta, condno, N_slice, 2, Gap_coupled, rho, g);
% wec_coupled.geometry.N_Slice = N_slice;
% wec_coupled.geometry.f1 = f1;
% wec_coupled.geometry.f2 = f2;
% wec_coupled.geometry.L = a_width;
% wec_coupled.geometry.d = dummy_draft;
% 
% % 2. Object for Isolated Array (N_grating = 1) for PTO calculations
% wec_iso = FMA_WEC(dummy_lambda, dummy_period, theta, h, Num_ord, Np, damp, delta, condno, N_slice, 1, 0, rho, g);
% wec_iso.geometry.N_Slice = N_slice;
% wec_iso.geometry.f1 = f1;
% wec_iso.geometry.f2 = f2;
% wec_iso.geometry.L = a_width;
% wec_iso.geometry.d = dummy_draft;
% % =========================================================================
% 
% for j = 1:length(draft_ratios)
%     d_draft = draft_ratios(j) * h;
%     
%     % Update geometry drafts for current outer loop
%     wec_coupled.geometry.d = d_draft;
%     wec_iso.geometry.d = d_draft;
%     
%     fprintf('--> Processing Draft Ratio: d/h = %.3f\n', draft_ratios(j));
%     
%     for i = 1:length(kh_list)
%         kh_val = kh_list(i); k = kh_val / h; lambda = 2 * pi / k; 
%         period = 5 * lambda; alpha = k * tanh(k * h); omega = sqrt(g * k * tanh(k * h));
%         slice_Ly = (f2 - f1) * period; 
%         
%         %% STEP 1 & 2: SOLVE COUPLED ARRAY HYDRODYNAMICS
%         % Update wave state for coupled WEC
%         wec_coupled.lambda = lambda;
%         wec_coupled.period = period;
%         wec_coupled.solver.FMA_Stratify(alpha);
%         
%         % Extract values directly from internal states[cite: 11]
%         F_ex = wec_coupled.solver.F_ex;
%         Hydro_Matrix = wec_coupled.solver.Hydro_Matrix;
%         R_all = wec_coupled.solver.R;
%         T_all = wec_coupled.solver.T;
%             
%         F_ex_2D_conj = (rho * g) * F_ex / slice_Ly;
%         Ref_Scat = R_all(Num_ord+1,1);
%         Trans_Scat = T_all(Num_ord+1,1);
%         Rad_R_A = T_all(Num_ord+1,2);
%         Rad_R_B = T_all(Num_ord+1,3);
%         Rad_L_A = R_all(Num_ord+1,2);
%         Rad_L_B = R_all(Num_ord+1,3);
%         Hydro_2D = Hydro_Matrix / slice_Ly;
%                     
%         % Validated rules: Exact real/imag parts, NO abs() to preserve negative cross-coupling[cite: 11]
%         mu_mat = rho * real(Hydro_2D);
%         lambda_mat = rho * omega * imag(Hydro_2D); 
%         
%         %% STEP 3: ISOLATED HYDRODYNAMICS FOR PTO
%         % Update wave state for isolated WEC
%         wec_iso.lambda = lambda;
%         wec_iso.period = period;
%         wec_iso.solver.FMA_Stratify(alpha);
%             
%         % The 4th output in original script (Fez_iso_A) was the Hydro_Matrix[cite: 11]
%         Hydro_iso = wec_iso.solver.Hydro_Matrix / slice_Ly;
%         mu_iso = rho * real(Hydro_iso);
%         lambda_iso = rho * omega * imag(Hydro_iso);
%                               
%         %% STEP 4: EQUATION OF MOTION (Unchanged)[cite: 11]
%         M_2D = [rho * a_width * d_draft, 0; 
%                 0, rho * a_width * d_draft];
%              
%         K_2D = [rho * g * a_width, 0; 
%                 0, rho * g * a_width];
%                 
%         lambda_PTO_single = sqrt((K_2D(1,1)/omega - omega*(M_2D(1,1) + mu_iso))^2 + lambda_iso^2);
%         PTO_Matrix = [lambda_PTO_single, 0; 0, lambda_PTO_single];
%              
%         % Validated Conjugated EOM
%         System_Matrix_Conj = -omega^2 * (M_2D + mu_mat) - 1i * omega * (lambda_mat + PTO_Matrix) + K_2D;
%         Xi_conj = System_Matrix_Conj \ F_ex_2D_conj; 
%         
%         %% STEP 5: POWER EXTRACTION (CWR) (Unchanged)[cite: 11]
%         P_1 = 0.5 * omega^2 * lambda_PTO_single * abs(Xi_conj(1))^2;
%         P_2 = 0.5 * omega^2 * lambda_PTO_single * abs(Xi_conj(2))^2;
%         
%         C_group = 0.5 * (omega / k) * (1 + 2*k*h / sinh(2*k*h));
%         P_incident = 0.5 * rho * g * C_group;
%         
%         CWR_data(j, i) = (P_1 + P_2) / P_incident;
%         
%         %% STEP 6: SUPERPOSITION (Unchanged)[cite: 11]
%         Rad_Multiplier = omega^2 / g; 
%         
%         T_total = Trans_Scat + Rad_Multiplier * (Rad_R_A * Xi_conj(1) + Rad_R_B * Xi_conj(2)); 
%         R_total = Ref_Scat   + Rad_Multiplier * (Rad_L_A * Xi_conj(1) + Rad_L_B * Xi_conj(2)); 
%         
%         K_T_data(j, i) = abs(T_total);
%         K_R_data(j, i) = abs(R_total);
%         if mod(i, 20) == 0
%             fprintf('----------kh=%.2f/%.2f--------------\n',kh_val,max(kh_list));
%         end
%     end
% end
% 
% fprintf('Computation Complete. Plotting Figure 4...\n');
% 
% %% Plotting Fig 4 (Replicating Ning 2017) (Unchanged)[cite: 11]
% figure('Position', [100, 100, 1200, 400], 'Name', 'Ning 2017 - Figure 4');
% 
% % (a) Reflection Coefficient
% subplot(1,3,1);
% plot(kh_list, K_R_data(1,:), 'k-', 'LineWidth', 2); hold on;
% plot(kh_list, K_R_data(2,:), 'r--', 'LineWidth', 2);
% plot(kh_list, K_R_data(3,:), 'b-.', 'LineWidth', 2);
% title('(a)'); xlabel('kh'); ylabel('K_R'); 
% xlim([0 10]); ylim([0 1.0]); 
% legend('d_1/h = d_2/h = 0.08', 'd_1/h = d_2/h = 0.125', 'd_1/h = d_2/h = 0.15', 'Location', 'SouthEast');
% grid on;
% 
% % (b) Transmission Coefficient
% subplot(1,3,2);
% plot(kh_list, K_T_data(1,:), 'k-', 'LineWidth', 2); hold on;
% plot(kh_list, K_T_data(2,:), 'r--', 'LineWidth', 2);
% plot(kh_list, K_T_data(3,:), 'b-.', 'LineWidth', 2);
% title('(b)'); xlabel('kh'); ylabel('K_T'); 
% xlim([0 10]); ylim([0 1.0]); 
% legend('d_1/h = d_2/h = 0.08', 'd_1/h = d_2/h = 0.125', 'd_1/h = d_2/h = 0.15');
% grid on;
% 
% % (c) Capture Width Ratio (CWR)
% subplot(1,3,3);
% plot(kh_list, CWR_data(1,:), 'k-', 'LineWidth', 2); hold on;
% plot(kh_list, CWR_data(2,:), 'r--', 'LineWidth', 2);
% plot(kh_list, CWR_data(3,:), 'b-.', 'LineWidth', 2);
% title('(c)'); xlabel('kh'); ylabel('\eta_{total}'); 
% xlim([0 10]); ylim([0 1.0]); 
% legend('d_1/h = d_2/h = 0.08', 'd_1/h = d_2/h = 0.125', 'd_1/h = d_2/h = 0.15');
% grid on;
%% Case 2:Benchmark Figure 8(a,b) of Yue 2019. Equivalent for Rectangular WECs using FMA_WEC solver
% h = 10.0;
% a = 0.3 * h;         % Half-width a/h = 0.3
% H_draft = 0.2 * h;   % Draft H/h = 0.2
% 
% rho = 1000; g = 9.81;
% Area = 4 * a^2;      % Rectangular area (2a x 2a)
% M_3D = rho * Area * H_draft;
% K_3D = rho * g * Area;
% 
% Np = 2;
% Num_ord = 10; % Increased slightly for better high-frequency stability
% delta = 1e-6; condno = 1e-8; damp = 0.0; theta_I = 0.0001;
% 
% N_slice = 1; % 1 slice for a perfect rectangle!
% N_grating = 1;
% dL = 2 * a;
% Gap = 0; % No gap for a single-row grating
% 
% % % 2. Calculate optimal PTO from isolated response
% kr_a = 0.88;
% k_res = kr_a / a;
% alpha_res = k_res * tanh(k_res * h);
% omega_res = sqrt(g * alpha_res);
% X3_bar_peak = 0.35;
% X3_dim = X3_bar_peak * (rho * g * Area);
% C_g = (omega_res / (2 * k_res)) * (1 + (2 * k_res * h) / sinh(2 * k_res * h));
% b_PTO = (k_res / (4 * rho * g * C_g)) * (X3_dim^2);
% 
% % % 3. Pre-calculate Isolated WEC Power Curve vs k*a
% fprintf('Pre-calculating Isolated WEC response...\n');
% ka_list = linspace(0.05, 5.5, 150);
% P_iso_list = zeros(size(ka_list));
% 
% d_iso = 30 * (2*a);
% fy_iso = (2*a / 2) / d_iso;
% f1_iso = 0.5 - fy_iso;
% f2_iso = 0.5 + fy_iso;
% 
% % =========================================================================
% % PRE-ALLOCATE OOP OBJECTS OUTSIDE THE LOOPS
% % =========================================================================
% dummy_lambda = 100;
% % Instantiate Isolated WEC Object
% wec_iso = FMA_WEC(dummy_lambda, d_iso, theta_I, h, Num_ord, Np, damp, delta, condno, N_slice, N_grating, Gap, rho, g);
% wec_iso.geometry.N_Slice = N_slice;
% wec_iso.geometry.f1 = f1_iso;
% wec_iso.geometry.f2 = f2_iso;
% wec_iso.geometry.L = dL;
% wec_iso.geometry.d = H_draft;
% 
% % Instantiate Array WEC Object (will dynamically update period and filling factors)
% wec_arr = FMA_WEC(dummy_lambda, d_iso, theta_I, h, Num_ord, Np, damp, delta, condno, N_slice, N_grating, Gap, rho, g);
% wec_arr.geometry.N_Slice = N_slice;
% wec_arr.geometry.L = dL;
% wec_arr.geometry.d = H_draft;
% % =========================================================================
% 
% for i = 1:length(ka_list)
%     ka = ka_list(i);
%     k = ka / a;
%     lambda = 2 * pi / k;
%     alpha = k * tanh(k * h);
%     omega = sqrt(g * alpha);
% 
%     % Update wave property and run solver
%     wec_iso.lambda = lambda;
%     wec_iso.solver.FMA_Stratify(alpha);
%     
%     % Extract internal states
%     F_ex_solver = wec_iso.solver.F_ex;
%     Hydro_Matrix_solver = wec_iso.solver.Hydro_Matrix;
% 
%     F_ex = (rho * g) * F_ex_solver;
%     mu_mat = rho * real(Hydro_Matrix_solver);
%     lambda_mat = rho * omega * imag(Hydro_Matrix_solver);
% 
%     Sys_Conj = -omega^2 * (M_3D + mu_mat) - 1i * omega * (lambda_mat + b_PTO) + K_3D;
%     Xi = Sys_Conj \ F_ex;
%     P_iso_list(i) = 0.5 * b_PTO * omega^2 * abs(Xi)^2;
% end
% 
% % Smoothing
% windowSize = 5;
% b_filt = (1/windowSize)*ones(1,windowSize);
% P_iso_list = filter(b_filt, 1, P_iso_list);
% 
% % % 4. Grid Setup for Figure 8(a)
% d_2a_list = linspace(1.05, 7.0, 100);
% kd_pi_list = linspace(0.1, 8.0, 200);
% 
% [D_2A, KD_PI] = meshgrid(d_2a_list, kd_pi_list);
% Q_gain = zeros(size(D_2A));
% KA_MAT = zeros(size(D_2A));
% 
% fprintf('Computing Array Gain Contour (Full Grid)...\n');
% total_pts = numel(D_2A);
% count = 0;
% 
% for col = 1:length(d_2a_list)
%     d_ratio = d_2a_list(col);
%     d_array = d_ratio * (2 * a);
% 
%     % Dynamically update the array geometry properties
%     wec_arr.period = d_array;
%     fy = (2*a / 2) / d_array;
%     wec_arr.geometry.f1 = 0.5 - fy;
%     wec_arr.geometry.f2 = 0.5 + fy;
% 
%     for row = 1:length(kd_pi_list)
%         count = count + 1;
%         kd_pi = kd_pi_list(row);
% 
%         kd = kd_pi * pi;
%         k = kd / d_array;
%         ka = k * a;
%         KA_MAT(row, col) = ka; 
% 
%         lambda = 2 * pi / k;
%         alpha = k * tanh(k * h);
%         omega = sqrt(g * alpha);
% 
%         % Update wavelength and run solver
%         wec_arr.lambda = lambda;
%         wec_arr.solver.FMA_Stratify(alpha);
%         
%         % Extract internal states
%         F_ex_solver = wec_arr.solver.F_ex;
%         Hydro_Matrix_solver = wec_arr.solver.Hydro_Matrix;
%         
%         F_ex = (rho * g) * F_ex_solver;
%         mu_mat = rho * real(Hydro_Matrix_solver);
%         lambda_mat = rho * omega * imag(Hydro_Matrix_solver);
% 
%         Sys_Conj = -omega^2 * (M_3D + mu_mat) - 1i * omega * (lambda_mat + b_PTO) + K_3D;
%         Xi = Sys_Conj \ F_ex;
% 
%         P_array = 0.5 * b_PTO * omega^2 * abs(Xi)^2;
%         P_iso_interp = interp1(ka_list, P_iso_list, ka, 'linear', 'extrap');
% 
%         Q_gain(row, col) = P_array / P_iso_interp;
% 
%         if mod(count, 10) == 0
%             fprintf('Progress: %.1f%%\n', 100*count/total_pts);
%         end
%     end
% end
% 
% % % 5. Plotting Figure 8(a) Equivalent
% fprintf('Computation Complete. Rendering Clean Plot...\n');
% figure('Position', [100, 100, 600, 500], 'Name', 'Figure 8(a) - Rectangular WECs');
% 
% Q_plot = Q_gain;
% Q_plot(Q_plot > 8) = 8.0;
% Q_plot(Q_plot < 0.0) = 0.0;
% Q_plot(KA_MAT > 2) = NaN;
% 
% contourf(D_2A, KD_PI, Q_plot, 60, 'LineStyle', 'none');
% colormap(jet);
% cb = colorbar;
% ylabel(cb, 'Array Gain, q');
% caxis([0 8]);
% 
% hold on;
% 
% mask_d = [min(d_2a_list), max(d_2a_list), min(d_2a_list)];
% mask_kd = [2.0 * 2 * min(d_2a_list) / pi, 2.0 * 2 * max(d_2a_list) / pi, max(kd_pi_list)*2];
% fill(mask_d, mask_kd, 'w', 'EdgeColor', 'none'); 
% 
% kd_pi_max = 2.0 * 2 * d_2a_list / pi;
% plot(d_2a_list, kd_pi_max, 'r-', 'LineWidth', 2, 'DisplayName', 'ka = 3.0 (Plot Limit)');
% 
% xlabel('Periodicity, d / 2a');
% ylabel('Wavenumber, kd / \pi');
% title('Figure 8(a) Eqv: Rectangular Arrays, \theta_I = 0^\circ');
% legend('Location', 'NorthWest');
% axis([1.05 7 0.8 8]);
% grid on;

%% Case3: Benchmark for Tokic & Yue (2019) - Figure 15 (a-e)
% % Adapted for stateful FMA_WEC architecture
% fprintf('Initializing Fig. 15(a-e) Single Row Benchmark...\n');
% 
% % Shared Physical Parameters[cite: 14]
% h = 10.0; a = 0.3 * h; H_draft = 0.2 * h;
% rho = 1000; g = 9.81;
% M_3D = rho * pi * a^2 * H_draft;
% K_3D = rho * g * pi * a^2;
% 
% Np = 5; Num_ord = 12; delta = 1e-6; condno = 1e-8; damp = 0.0;
% N_slice = 20; N_grating = 1; % Single Row
% Gap = 0; % No gap for single row
% 
% % Fig 15 Array Configuration[cite: 14]
% d_over_2a = 3.0;
% d_array = d_over_2a * (2 * a);
% dL = 2 * a / N_slice;
% 
% % Calculate the constant optimal b_PTO from isolated body at k_r*a = 0.88[cite: 14]
% kr_a = 0.88; k_res = kr_a / a; alpha_res = k_res * tanh(k_res * h); omega_res = sqrt(g * alpha_res);
% X3_bar_peak = 0.35; X3_dim = X3_bar_peak * (rho * g * pi * a^2);
% C_g = (omega_res / (2 * k_res)) * (1 + (2 * k_res * h) / sinh(2 * k_res * h));
% b_PTO_opt = (k_res / (4 * rho * g * C_g)) * (X3_dim^2);
% 
% % Define the 5 cases from Figure 15[cite: 14]
% cases = struct('kd_pi', {1.67, 1.67, 2.01, 3.72, 2.01}, ...
%                'theta_I_deg', {0, 0, 0, 0, 30}, ...
%                'PTO_on', {false, true, true, true, true}, ...
%                'label', {'a', 'b', 'c', 'd', 'e'});
% 
% x_grid = linspace(-5*h, 5*h, 600);
% y_grid = linspace(-d_array, d_array, 400);
% [X, Y] = meshgrid(x_grid, y_grid);
% 
% % =========================================================================
% % PRE-ALLOCATE OOP OBJECT OUTSIDE THE LOOP
% % =========================================================================
% dummy_lambda = 100;
% dummy_theta = 0.0001;
% 
% wec = FMA_WEC(dummy_lambda, d_array, dummy_theta, h, Num_ord, Np, damp, delta, condno, N_slice, N_grating, Gap, rho, g);
% wec.geometry.SlicedCircle(d_array, a, dL, H_draft);
% 
% figure('Position', [100, 100, 1200, 800], 'Name', 'Fig 15(a-e)');
% 
% for c_idx = 1:length(cases)
%     kd_pi = cases(c_idx).kd_pi;
%     theta_I = max(cases(c_idx).theta_I_deg * pi / 180, 0.0001); % Avoid exact 0 singularity[cite: 14]
%     b_PTO = cases(c_idx).PTO_on * b_PTO_opt;
% 
%     k = (kd_pi * pi) / d_array;
%     lambda = 2 * pi / k;
%     alpha = k * tanh(k * h); omega = sqrt(g * alpha);
% 
%     fprintf('Solving 15(%c): kd/pi=%.2f, theta=%.0f, PTO=%d\n', cases(c_idx).label, kd_pi, cases(c_idx).theta_I_deg, cases(c_idx).PTO_on);
% 
%     % Update dynamic wave state
%     wec.lambda = lambda;
%     wec.theta = theta_I;
%     
%     % Execute internal solver cascade
%     wec.solver.FMA_Stratify(alpha);
%     
%     % Extract internal physics states
%     F_ex = (rho * g) * wec.solver.F_ex;
%     mu_mat = rho * real(wec.solver.Hydro_Matrix);
%     lambda_mat = rho * omega * imag(wec.solver.Hydro_Matrix);
% 
%     % Solve Equation of Motion[cite: 14]
%     Sys_Conj = -omega^2 * (M_3D + mu_mat) - 1i * omega * (lambda_mat + b_PTO) + K_3D;
%     Xi = Sys_Conj \ F_ex;
%     Rad_Multiplier = omega^2 / g;
% 
%     % Generate Wave Field leveraging internal wave states
%     eta_amp = wec.solver.WaveField(Xi, Rad_Multiplier, x_grid, y_grid);
% 
%     subplot(5, 1, c_idx);
%     contourf(X/h, Y/h, eta_amp, 80, 'LineStyle', 'none');
%     colormap(jet); colorbar; caxis([0, 2.0]); hold on;
% 
%     theta_cyl = linspace(0, 2*pi, 50);
%     fill((a*cos(theta_cyl))/h, (a*sin(theta_cyl))/h, [0.2 0.2 0.2], 'EdgeColor', 'k');
% 
%     title(sprintf('(%c) kd=%.2f\\pi, PTO Active: %d', cases(c_idx).label, kd_pi, cases(c_idx).PTO_on));
%     ylabel('y/h'); axis equal; xlim([-5, 5]); ylim([-1.8, 1.8]); set(gca, 'Layer', 'top'); box on;
% end


%% Case 4: 2D Field Maps of the Symmetric and Antisymmetric Trapped Modes


fprintf('Calculating 2D Field Maps for Double-Row WEC...\n');

% =========================================================================
% 1. Shared Geometric & Physical Parameters
% =========================================================================
h_water = 10.0;
L = 0.3 * h_water;         
H_draft = 0.2 * h_water;   
d_array = 8.0 * L;         

Np = 2; Num_ord = 12; delta = 1e-6; condno = 1e-8; damp = 0.0;
rho = 1000; g = 9.81;

N_grating = 2;       
N_slice = 20;        
dL = 2 * L / N_slice;

f1 = -L/d_array; f2 = L/d_array;

L_vec = repmat(dL, 1, N_slice);
d_vec = repmat(H_draft, 1, N_slice);
f1_vec = repmat(f1, 1, N_slice);
f2_vec = repmat(f2, 1, N_slice);

M_3D = rho * H_draft * (2*L) * (2*L);
K_3D = rho * g * (2*L) * (2*L);
M_mat = eye(N_grating) * M_3D;  
K_mat = eye(N_grating) * K_3D;
b_PTO = 0; 

G = 2 * pi / d_array;
kx_norm = 0.35;

% Define the plotting grid, and a SHIFTED evaluation grid
y_grid = linspace(-1.5 * d_array, 1.5 * d_array, 600);    
y_eval = y_grid + d_array / 2; % Shifts the WEC center to y=0

% =========================================================================
% PRE-ALLOCATE OOP OBJECT OUTSIDE THE SECTIONS
% =========================================================================
dummy_lambda = 100;
dummy_theta = 0.0001;
dummy_Gap = d_array;

wec = FMA_WEC(dummy_lambda, d_array, dummy_theta, h_water, Num_ord, Np, damp, delta, condno, N_slice, N_grating, dummy_Gap, rho, g);
wec.geometry.N_Slice = N_slice;
wec.geometry.L = L_vec;
wec.geometry.d = d_vec;
wec.geometry.f1 = f1_vec;
wec.geometry.f2 = f2_vec;

% =========================================================================
% 2. State A: Antisymmetric Mode 
% =========================================================================
h_paper_A = 2.7700;       
D_c2c_A = h_paper_A * d_array; 
x_grid_A = linspace(-2 * d_array, 6*d_array, 600);
d_x_A = D_c2c_A;                  

fprintf('\n--- Auto-tuning Antisymmetric Mode (D_c2c = %.4f d_array) ---\n', h_paper_A);
k_fine_vec = 0.64151; 
k_test = k_fine_vec;

k_inc = k_test * G;
lam = 2 * pi / k_inc;
alp = k_inc * tanh(k_inc * h_water);
om = sqrt(g * alp);
th = asin(kx_norm / k_test);
    
% Dynamically update states
wec.lambda = lam;
wec.theta = th;
wec.Gap = D_c2c_A;

% Execute Solver
wec.solver.FMA_Stratify(alp);
    
% Extract results
F_ex = (rho * g) * wec.solver.F_ex;
mu_mat = rho * real(wec.solver.Hydro_Matrix);
lambda_mat = rho * om * imag(wec.solver.Hydro_Matrix);
Sys_Conj = -om^2 * (M_mat + mu_mat) - 1i * om * (lambda_mat + b_PTO) + K_mat;
best_Xi_A = Sys_Conj \ F_ex;

best_k_A = k_test;
best_Rad_A = om^2 / g;

k_inc_A = best_k_A * G;
phase_diff_A = abs(rad2deg(angle(best_Xi_A(1))) - rad2deg(angle(best_Xi_A(2))));

fprintf('-> Peak Locked: k/G = %.5f\n', best_k_A);
fprintf('-> WEC Phase Difference: %.1f degrees\n', phase_diff_A);

% Generate wave field natively
eta_amp_A = wec.solver.WaveField(best_Xi_A, best_Rad_A, x_grid_A, y_eval);

% =========================================================================
% 3. State B: Symmetric Mode 
% =========================================================================
h_paper_S = 3.6400;       
D_c2c_S = h_paper_S * d_array; 
x_grid_S = linspace(-2 * d_array, 6*d_array, 600);
d_x_S = D_c2c_S; 

fprintf('\n--- Auto-tuning Symmetric Mode (D_c2c = %.4f d_array) ---\n', h_paper_S);
k_fine_vec_S = 0.6457; 
k_test = k_fine_vec_S;

k_inc = k_test * G;
lam = 2 * pi / k_inc;
alp = k_inc * tanh(k_inc * h_water);
om = sqrt(g * alp);
th = asin(kx_norm / k_test);
    
% Dynamically update states
wec.lambda = lam;
wec.theta = th;
wec.Gap = D_c2c_S;

% Execute Solver
wec.solver.FMA_Stratify(alp);
    
% Extract results
F_ex = (rho * g) * wec.solver.F_ex;
mu_mat = rho * real(wec.solver.Hydro_Matrix);
lambda_mat = rho * om * imag(wec.solver.Hydro_Matrix);
Sys_Conj = -om^2 * (M_mat + mu_mat) - 1i * om * (lambda_mat + b_PTO) + K_mat;
best_Xi_S = Sys_Conj \ F_ex;

best_k_S = k_test;
best_Rad_S = om^2 / g;

k_inc_S = best_k_S * G;
phase_diff_S = abs(rad2deg(angle(best_Xi_S(1))) - rad2deg(angle(best_Xi_S(2))));

fprintf('-> Peak Locked: k/G = %.5f\n', best_k_S);
fprintf('-> WEC Phase Difference: %.1f degrees\n', phase_diff_S);

% Generate wave field natively
eta_amp_S = wec.solver.WaveField(best_Xi_S, best_Rad_S, x_grid_S, y_eval);

% =========================================================================
% 4. Visualization (Normalized to d_array)
% =========================================================================
try turbo_map = turbo(); catch; turbo_map = 'jet'; end
figure('Color','w','Position',[100, 100, 1000, 900]);

w_y_norm = abs(f2_vec(1)-f1_vec(1));     
L_norm = L / d_array;                    

% The WECs are now centered at integers!
y_centers_norm = [0, 1, -1]; 

% --- Subplot 1: Antisymmetric Mode ---
subplot(2,1,1);
imagesc(x_grid_A / d_array, y_grid / d_array, eta_amp_A);
set(gca, 'YDir', 'normal'); 
colormap(turbo_map); 
c1 = colorbar;
ylabel(c1, 'Hydro |\eta|'); % Fixed for pre-R2014b MATLAB compatibility
try caxis([0 5]); catch; clim([0 5]); end 
hold on;

for yc = y_centers_norm
    rectangle('Position', [-L_norm, yc - w_y_norm/2, 2*L_norm, w_y_norm], 'EdgeColor', 'w', 'LineWidth', 1.5);
    rectangle('Position', [h_paper_A - L_norm, yc - w_y_norm/2, 2*L_norm, w_y_norm], 'EdgeColor', 'w', 'LineWidth', 1.5);
end
hold off;

ylim([-1.5, 1.5]);
xlim([-2, 6]); 
axis equal;
xlabel('x / d_{array}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y / d_{array}', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Antisymmetric Mode (D_{c2c} = %.4f d_{array}, k/G = %.5f)', h_paper_A, best_k_A), 'FontSize', 13);

% --- Subplot 2: Symmetric Mode ---
subplot(2,1,2);
imagesc(x_grid_S / d_array, y_grid / d_array, eta_amp_S);
set(gca, 'YDir', 'normal'); 
colormap(turbo_map);
c2 = colorbar;
ylabel(c2, 'Hydro |\eta|'); % Fixed for pre-R2014b MATLAB compatibility
try caxis([0 5]); catch; clim([0 5]); end 
hold on;

for yc = y_centers_norm
    rectangle('Position', [-L_norm, yc - w_y_norm/2, 2*L_norm, w_y_norm], 'EdgeColor', 'w', 'LineWidth', 1.5);
    rectangle('Position', [h_paper_S - L_norm, yc - w_y_norm/2, 2*L_norm, w_y_norm], 'EdgeColor', 'w', 'LineWidth', 1.5);
end
hold off;

ylim([-1.5, 1.5]);
xlim([-2, 6]);
axis equal;
xlabel('x / d_{array}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y / d_{array}', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Symmetric Mode (D_{c2c} = %.4f d_{array}, k/G = %.5f)', h_paper_S, best_k_S), 'FontSize', 13);

fprintf('\nField Mapping Complete!\n');



