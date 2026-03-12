%% =======================================================================
%  NEMOH TIME-DOMAIN SIMULATION - Cummins Equation with Prony Method
%  =======================================================================
%  Solves the 6-DOF Cummins equation using Nemoh hydrodynamic outputs:
%
%    [M + A_inf] x_ddot(t) + int_0^t K(t-tau) x_dot(tau) dtau
%        + (C_h + C_moor) x(t) = F_exc(t) + F_ext(t)
%
%  The radiation convolution integral is replaced by a Prony state-space
%  approximation:  K_ij(t) ~ sum_k alpha(i,j,k) * exp(beta(i,j,k) * t)
%
%  Wave excitation forces are computed in the frequency domain from
%  Nemoh transfer functions, then reconstructed in time domain.
%
%  References:
%    Cummins, W.E. (1962) "The Impulse Response Function and Ship Motions"
%    Perez & Fossen (2009) "A Matlab Toolbox for Parametric Identification
%       of Radiation-Force Models of Ships and Offshore Structures"
%  =======================================================================

clear; clc; close all;

%% ===== USER CONFIGURATION =============================================

% --- Nemoh Case Directory -----------------------------------------------
%  Must contain: results/, Mechanics/, Nemoh.cal
nemoh_dir = fullfile('C:','Users','Fonta','OneDrive','Desktop', ...
    'Florid_Scale_Model_OpenFAST','Nemoh','SimplePontoon');

% --- Simulation Time ----------------------------------------------------
t_end   = 360;       % Simulation duration [s]
dt_out  = 0.01;     % Output time step [s] (ODE solver adapts internally)

% --- Wave Type ----------------------------------------------------------
wave_type = 'irregular';   % 'regular' or 'irregular'
wave_dir  = 0;           % Wave direction [deg] (must match Nemoh run)

% Regular wave parameters
reg_H     = 0.02;       % Wave height [m]
reg_T     = 1.0;        % Wave period [s]
reg_phase = 0;           % Initial phase [rad]

% Irregular wave parameters (JONSWAP / Pierson-Moskowitz)
irr_Hs       = 0.014;    % Significant wave height [m]
irr_Tp       = 3.2;     % Peak period [s]
irr_gamma    = 3.3;     % Peak enhancement factor (1.0 = Pierson-Moskowitz)
irr_ncomp    = 128;     % Number of wave components
irr_wmin     = 1;     % Min frequency [rad/s]
irr_wmax     = 15.0;    % Max frequency [rad/s]
irr_seed     = 42;      % Random seed for phase generation
% Wave components are uniformly spaced in wavenumber (k) space.
% The dispersion relation w^2 = g*k*tanh(k*d) maps each k_n to a
% non-uniformly spaced w_n. Incommensurate frequencies → signal does
% not repeat, giving good statistics with few components (128 typical).

% --- Environment -------------------------------------------------------
water_depth  = 1.0;     % Water depth [m] (for dispersion relation)
g_acc        = 9.81;    % Gravitational acceleration [m/s^2]

% --- Retardation Function & Prony Fit -----------------------------------
n_prony    = 12;          % Number of exponential terms per DOF pair
t_irf_max  = 3;         % Max time for retardation function [s]
dt_irf     = 0.01;       % Time step for retardation function [s]
prony_tol  = 1e-3;       % Skip DOF pairs with max|K|/max|K_all| below this

% --- Mooring Stiffness (Additional Linear 6x6) -------------------------
%  Added to hydrostatic stiffness: C_total = C_hydrostatic + C_mooring
%  DOF order: [2urge, sway, heave, roll, pitch, yaw]
K_mooring = zeros(6);
% Examples:
K_mooring(1,1) = 50;    % Surge restoring [N/m]
% K_mooring(2,2) = 50;    % Sway restoring  [N/m]
% K_mooring(3,3) = 0;     % Heave (usually already in hydrostatics)
% K_mooring(4,4) = 5;     % Roll  [N.m/rad]
% K_mooring(5,5) = 10;    % Pitch [N.m/rad]
% K_mooring(6,6) = 20;    % Yaw   [N.m/rad]
% Off-diagonal coupling terms can also be specified.
% ---- Option A: 3-line CATENARY spread at 120 deg ----
% depth = 1.0;  R_anch = 4.0;  R_fair = 0.15;
% for i = 1:3
%     ang = (i-1)*120 * pi/180;
%     lines(i).type     = 'catenary';
%     lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
%     lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); 0];
%     lines(i).L0 = 2.5;   lines(i).w = 0.2;   lines(i).EA = 5000;
% end

% ---- Option B: 3-line TAUT WIRE spread at 120 deg ----
%  Taut lines run at an angle from seabed anchors to near-waterline
%  fairleads. Pretension is set by making L0 shorter than the actual
%  anchor-to-fairlead distance. No weight-per-length needed.
%
%  Geometry:  anchor radius = 1.5 m, fairlead radius = 0.15 m
%             depth = 1.0 m, fairlead z = -0.01 m (just below waterline)
%             chord length ≈ sqrt(1.35^2 + 0.99^2) ≈ 1.67 m
%             L0 = 1.58 m  → ~5% pretension
%             EA = 2500 N   (elastic cord / nylon scale-model line)
depth = water_depth;  R_anch = 1.5;  R_fair = 0.15;  z_fair = -0.01;
for i = 1:3
    ang = (i-1)*120 * pi/180;
    lines(i).type     = 'taut';
    lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
    lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); z_fair];
    lines(i).L0 = 1.6;     % unstretched length [m]  (< chord → pretension)
    lines(i).EA = 1000;      % axial stiffness [N]
end

% ---- Option C: 4-line TAUT WIRE spread at 90 deg ----
%  Four orthogonal lines for symmetric restoring in surge & sway.
% depth = 1.0;  R_anch = 1.5;  R_fair = 0.12;  z_fair = -0.01;
% for i = 1:4
%     ang = (i-1)*90 * pi/180;
%     lines(i).type     = 'taut';
%     lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
%     lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); z_fair];
%     lines(i).L0 = 1.55;     lines(i).EA = 2500;
% end

[K_mooring, F0] = compute_mooring_stiffness(lines, zeros(6,1));
% --- External Force with Application Point --------------------------------
%  Specify a 3D force vector and where it acts. Moments about CoG are
%  computed automatically:  M = cross(r, F),  r = application_point - CoG.
%
%  ext_force: applied force [Fx; Fy; Fz] in [N]
%             - constant 3x1 vector, OR
%             - function handle ext_force(t) returning 3x1 [N]
%  ext_point: [x; y; z] force application point in body-fixed frame [m]
%  CoG:       [x; y; z] centre of gravity in body-fixed frame [m]
%
%  The resulting 6-DOF vector passed to the solver is:
%    F_external(t) = [Fx; Fy; Fz; Mx; My; Mz]
%
CoG       = [0; 0; -0.01];
ext_force = [0.36; 0; 0];
ext_point = [0; 0; -0.2];
% Examples:
%  Horizontal surge force at z = -0.35 m (creates pitch moment):
%    ext_force = [0.3; 0; 0];
%    ext_point = [0; 0; -0.35];
%    -> moment = cross([0;0;-0.35]-CoG, [0.3;0;0]) = [0; -0.105; 0] N.m
%
%  Time-varying force at a point:
%    ext_force = @(t) [0.5*sin(2*pi*t/5); 0; 0];
%    ext_point = [0; 0; -0.35];
%
%  Constant heave preload (at CoG, so no moment):
%    ext_force = [0; 0; -5.0];
%    ext_point = CoG;

% Build the 6-DOF external force function (forces + moments about CoG)
F_external = build_ext_force(ext_force, ext_point, CoG);

% --- Initial Conditions -------------------------------------------------
q0 = zeros(6,1);   % Displacement [m, m, m, rad, rad, rad]
v0 = zeros(6,1);   % Velocity     [m/s, m/s, m/s, rad/s, rad/s, rad/s]

% --- ODE Solver Options -------------------------------------------------
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.02);

% --- Active Degrees of Freedom ------------------------------------------
%  Specify which DOFs to simulate. Disabled DOFs are locked at zero.
%  1=Surge, 2=Sway, 3=Heave, 4=Roll, 5=Pitch, 6=Yaw
active_dofs = [1, 3, 5];   % Surge, Heave, Pitch

% --- Labels for Plotting ------------------------------------------------
dof_labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
dof_units  = {'m','m','m','deg','deg','deg'};
% Conversion factor: ODE works in radians internally, display in degrees
dof_scale  = [1, 1, 1, 180/pi, 180/pi, 180/pi];
vel_units  = {'m/s','m/s','m/s','deg/s','deg/s','deg/s'};

%% ===== LOAD NEMOH DATA ================================================
fprintf('========================================\n');
fprintf(' Nemoh Time-Domain Simulation\n');
fprintf('========================================\n\n');
fprintf('Loading Nemoh data from:\n  %s\n\n', nemoh_dir);

% --- Mass Matrix (6x6) ---
M = load_matrix(fullfile(nemoh_dir, 'Mechanics', 'Inertia.dat'));
fprintf('Mass matrix loaded. Diagonal: [');
fprintf(' %.4e', diag(M));
fprintf(' ]\n');

% --- Hydrostatic Stiffness (6x6) ---
K_h = load_matrix(fullfile(nemoh_dir, 'Mechanics', 'Kh.dat'));
fprintf('Hydrostatic stiffness loaded. K_33=%.2f, K_44=%.2f, K_55=%.2f\n', ...
    K_h(3,3), K_h(4,4), K_h(5,5));

% --- Radiation Coefficients: A(omega), B(omega) ---
[omega_rad, A_rad, B_rad] = load_radiation( ...
    fullfile(nemoh_dir, 'results', 'RadiationCoefficients.tec'));
n_freq = length(omega_rad);
fprintf('Radiation coefficients: %d frequencies [%.2f - %.2f] rad/s\n', ...
    n_freq, omega_rad(1), omega_rad(end));

% --- Excitation Force Transfer Function ---
[omega_exc, F_exc_mag, F_exc_phase] = load_excitation( ...
    fullfile(nemoh_dir, 'results', 'ExcitationForce.tec'));
fprintf('Excitation force loaded: %d frequencies, 6 DOFs\n', length(omega_exc));

%% ===== APPLY DOF CONSTRAINTS ===========================================
dof_mask = false(6,1);
dof_mask(active_dofs) = true;
disabled_dofs = find(~dof_mask);
n_active = length(active_dofs);

fprintf('\nActive DOFs (%d/6):', n_active);
for ii = 1:n_active
    fprintf(' %s', dof_labels{active_dofs(ii)});
end
fprintf('\n');
if ~isempty(disabled_dofs)
    fprintf('Disabled DOFs:');
    for ii = 1:length(disabled_dofs)
        fprintf(' %s', dof_labels{disabled_dofs(ii)});
    end
    fprintf(' (locked at zero)\n');
end

% Zero out disabled DOF rows/cols in mass and stiffness
M(disabled_dofs, :) = 0;  M(:, disabled_dofs) = 0;
K_h(disabled_dofs, :) = 0;  K_h(:, disabled_dofs) = 0;
K_mooring(disabled_dofs, :) = 0;  K_mooring(:, disabled_dofs) = 0;

% Set diagonal to identity for disabled DOFs so M+A_inf stays invertible
for dd = disabled_dofs'
    M(dd, dd) = 1;
end

% Zero out disabled DOFs in radiation coefficients
A_rad(disabled_dofs, :, :) = 0;  A_rad(:, disabled_dofs, :) = 0;
B_rad(disabled_dofs, :, :) = 0;  B_rad(:, disabled_dofs, :) = 0;

% Zero out disabled DOFs in excitation force
F_exc_mag(disabled_dofs, :) = 0;
F_exc_phase(disabled_dofs, :) = 0;

% Zero out initial conditions for disabled DOFs
q0(disabled_dofs) = 0;
v0(disabled_dofs) = 0;

%% ===== COMPUTE RETARDATION FUNCTION K(t) FROM B(omega) ================
fprintf('\nComputing retardation functions K(t) from B(omega)...\n');

t_irf = (0 : dt_irf : t_irf_max)';
n_t_irf = length(t_irf);

% K_ij(t) = (2/pi) * integral_0^inf B_ij(w) * cos(w*t) dw
K_ret = zeros(6, 6, n_t_irf);

cos_mat = cos(omega_rad(:) * t_irf(:)');  % [n_freq x n_t_irf]

for i = 1:6
    for j = 1:6
        B_ij = squeeze(B_rad(i, j, :));
        integrand = B_ij(:) .* cos_mat;  % [n_freq x n_t_irf]
        K_ret(i, j, :) = (2/pi) * trapz(omega_rad, integrand);
    end
end

% Infinite-frequency added mass (highest frequency approximation)
A_inf = A_rad(:, :, end);
fprintf('A_inf diagonal: [');
fprintf(' %.4e', diag(A_inf));
fprintf(' ]\n');
fprintf('Retardation functions computed on [0, %.1f] s with dt=%.3f s (%d points)\n', ...
    t_irf_max, dt_irf, n_t_irf);

%% ===== FIT PRONY COEFFICIENTS ==========================================
fprintf('\nFitting Prony coefficients (%d terms per DOF pair)...\n', n_prony);

alpha_prony = zeros(6, 6, n_prony);
beta_prony  = zeros(6, 6, n_prony);

K_max_global = max(abs(K_ret(:)));
n_fitted = 0;

for i = 1:6
    for j = 1:6
        K_ij = squeeze(K_ret(i, j, :));

        % Skip negligible retardation functions
        if max(abs(K_ij)) < prony_tol * K_max_global
            continue;
        end

        [a, b] = fit_prony(t_irf, K_ij, n_prony);
        alpha_prony(i, j, :) = a;
        beta_prony(i, j, :)  = b;
        n_fitted = n_fitted + 1;
    end
end

fprintf('  Fitted %d / 36 DOF pairs (rest below threshold)\n', n_fitted);

% Report fit quality
fprintf('  Prony fit quality:\n');
worst_nrmse = 0;
for i = 1:6
    for j = 1:6
        K_ij = squeeze(K_ret(i, j, :));
        if max(abs(K_ij)) < prony_tol * K_max_global; continue; end

        K_fit = zeros(n_t_irf, 1);
        for k = 1:n_prony
            K_fit = K_fit + alpha_prony(i,j,k) * exp(beta_prony(i,j,k) * t_irf);
        end
        nrmse = sqrt(mean((K_ij(:) - K_fit).^2)) / max(abs(K_ij));
        worst_nrmse = max(worst_nrmse, nrmse);
        if nrmse > 0.10
            fprintf('    K(%d,%d): NRMSE = %.1f%% ** consider increasing n_prony **\n', ...
                i, j, nrmse*100);
        end
    end
end
fprintf('  Worst NRMSE: %.1f%%\n', worst_nrmse*100);

%% ===== GENERATE WAVE INPUT =============================================
fprintf('\nWave input: %s\n', wave_type);

switch lower(wave_type)
    case 'regular'
        wave_omega = 2*pi / reg_T;
        wave_amp   = reg_H / 2;
        wave_eps   = reg_phase;
        n_waves    = 1;
        S_wave     = [];
        fprintf('  H = %.4f m, T = %.2f s, omega = %.3f rad/s\n', ...
            reg_H, reg_T, wave_omega);

    case 'irregular'
        rng(irr_seed);

        % k-space discretisation:
        %   Uniform spacing in wavenumber k, mapped to non-uniform omega
        %   via the dispersion relation: w^2 = g*k*tanh(k*d).
        %   Incommensurate frequencies → signal does not repeat.
        k_min = solve_dispersion(irr_wmin, water_depth, g_acc);
        k_max = solve_dispersion(irr_wmax, water_depth, g_acc);
        dk_irr = (k_max - k_min) / irr_ncomp;
        k_vec  = k_min + ((1:irr_ncomp)' - 0.5) * dk_irr;  % midpoints

        % Dispersion: omega_n and group velocity Cg_n
        wave_omega = sqrt(g_acc * k_vec .* tanh(k_vec * water_depth));
        Cg = (wave_omega ./ (2*k_vec)) .* ...
             (1 + 2*k_vec*water_depth ./ sinh(2*k_vec*water_depth));
        dw_vec = Cg * dk_irr;             % non-uniform dw per component
        n_waves = irr_ncomp;
        wave_eps = 2*pi * rand(n_waves, 1);

        % JONSWAP spectrum (gamma=1 gives Pierson-Moskowitz)
        S_wave = jonswap_spectrum(wave_omega, irr_Hs, irr_Tp, irr_gamma);
        wave_amp = sqrt(2 * S_wave .* dw_vec);

        Hs_check = 4 * sqrt(sum(S_wave .* dw_vec));
        fprintf('  Hs = %.4f m, Tp = %.2f s, gamma = %.1f\n', irr_Hs, irr_Tp, irr_gamma);
        fprintf('  k-space: %d components, dk = %.4f rad/m\n', n_waves, dk_irr);
        fprintf('  k in [%.4f, %.4f] rad/m  ->  omega in [%.3f, %.3f] rad/s\n', ...
            k_vec(1), k_vec(end), wave_omega(1), wave_omega(end));
        fprintf('  Hs (reconstructed) = %.4f m\n', Hs_check);

    otherwise
        error('Unknown wave_type: "%s". Use ''regular'' or ''irregular''.', wave_type);
end

%% ===== INTERPOLATE EXCITATION FORCE TO WAVE FREQUENCIES ================
F_mag_w   = zeros(6, n_waves);
F_phase_w = zeros(6, n_waves);

for i = 1:6
    F_mag_w(i,:)   = interp1(omega_exc, F_exc_mag(i,:), wave_omega(:), 'linear', 0);
    F_phase_w(i,:) = interp1(omega_exc, F_exc_phase(i,:), wave_omega(:), 'linear', 0);
end

% Pre-compute products for fast ODE evaluation:
%   F_exc_i(t) = sum_n wave_amp(n) |F_i(wn)| cos(wn*t + eps_n + phase_i(wn))
%              = F_cos(i,:) * cos(base_phase) - F_sin(i,:) * sin(base_phase)
F_amp_wave = F_mag_w .* wave_amp(:)';       % [6 x n_waves]
F_cos_pre  = F_amp_wave .* cos(F_phase_w);  % [6 x n_waves]
F_sin_pre  = F_amp_wave .* sin(F_phase_w);  % [6 x n_waves]

%% ===== SET UP AND SOLVE STATE-SPACE ODE ================================
fprintf('\nSolving Cummins equation...\n');

n_dof    = 6;
n_prony_states = n_dof * n_dof * n_prony;
n_states = 2*n_dof + n_prony_states;
fprintf('  State vector: %d states (6 pos + 6 vel + %d Prony)\n', ...
    n_states, n_prony_states);

% Pre-invert effective mass matrix
M_total = M + A_inf;
M_inv   = M_total \ eye(6);

% Total stiffness
K_total = K_h + K_mooring;

% Pack parameters for ODE function
params = struct();
params.M_inv       = M_inv;
params.K_total     = K_total;
params.n_prony     = n_prony;
params.alpha_prony = alpha_prony;   % [6, 6, n_prony]
params.beta_prony  = beta_prony;    % [6, 6, n_prony]
params.n_waves     = n_waves;
params.wave_omega  = wave_omega(:);
params.wave_eps    = wave_eps(:);
params.F_cos_pre   = F_cos_pre;     % [6 x n_waves]
params.F_sin_pre   = F_sin_pre;     % [6 x n_waves]
params.F_external  = F_external;
params.disabled_dofs = disabled_dofs;

% Initial state vector
x0 = zeros(n_states, 1);
x0(1:6)  = q0;
x0(7:12) = v0;

% Time span
tspan = 0 : dt_out : t_end;

% Solve
tic;
[t_sol, x_sol] = ode45(@(t, x) cummins_rhs(t, x, params), tspan, x0, ode_opts);
elapsed = toc;
fprintf('  Solved in %.2f s (%d output steps)\n', elapsed, length(t_sol));

%% ===== EXTRACT AND POST-PROCESS RESULTS ================================
fprintf('\nPost-processing results...\n');

q_sol = x_sol(:, 1:6);      % Displacements [m or rad]
v_sol = x_sol(:, 7:12);     % Velocities    [m/s or rad/s]

% Reconstruct wave elevation
eta = zeros(length(t_sol), 1);
for n = 1:n_waves
    eta = eta + wave_amp(n) * cos(wave_omega(n) * t_sol + wave_eps(n));
end

% Reconstruct excitation forces
F_exc_t = zeros(length(t_sol), 6);
for it = 1:length(t_sol)
    base_phase = wave_omega(:) * t_sol(it) + wave_eps(:);
    cv = cos(base_phase);
    sv = sin(base_phase);
    F_exc_t(it, :) = (F_cos_pre * cv - F_sin_pre * sv)';
end

% Compute radiation memory force (for output only)
F_rad_t = zeros(length(t_sol), 6);
for it = 1:length(t_sol)
    I_3d = reshape(x_sol(it, 13:end)', [6, 6, n_prony]);
    F_rad_t(it, :) = sum(sum(I_3d, 3), 2)';
end

% Reconstruct external force time series (6-DOF, incl. moments)
F_ext_t = zeros(length(t_sol), 6);
for it = 1:length(t_sol)
    F_ext_t(it, :) = F_external(t_sol(it))';
end

% Compute net force and accelerations
%   F_net = F_exc + F_ext - F_rad - K_total*q  (= [M+A_inf]*a)
F_net_t = zeros(length(t_sol), 6);
a_sol   = zeros(length(t_sol), 6);
for it = 1:length(t_sol)
    rhs = F_exc_t(it,:)' + F_ext_t(it,:)' - K_total * q_sol(it,:)' - F_rad_t(it,:)';
    F_net_t(it,:) = rhs';
    a_sol(it,:)   = (M_inv * rhs)';
end

% Compute net hydrodynamic force (excl. mooring restoring)
%   F_hydro_net = F_exc + F_ext - F_rad - K_h*q  (what mooring must resist)
F_hydro_net_t = F_exc_t + F_ext_t - F_rad_t - (K_h * q_sol')';

% Compute total hydrodynamic + external force on body for each DOF:
%   F_total = F_exc - F_rad + F_ext
%   (restoring force is the body's reactive response, not an applied load)
F_total_t = F_exc_t - F_rad_t + F_ext_t;

% Save results to struct
results = struct();
results.t           = t_sol;
results.q           = q_sol;
results.v           = v_sol;
results.a           = a_sol;
results.eta         = eta;
results.F_exc       = F_exc_t;
results.F_rad       = F_rad_t;
results.F_ext       = F_ext_t;
results.F_total     = F_total_t;
results.F_net       = F_net_t;
results.F_hydro_net = F_hydro_net_t;
results.wave_type   = wave_type;
results.dof_labels  = {dof_labels};
results.params      = params;

% Print external force info
F_ext_6dof = F_external(0);
fprintf('External force (6-DOF): [');
fprintf(' %.4e', F_ext_6dof);
fprintf(' ]\n');
fprintf('  Force  = [%.4f, %.4f, %.4f] N  at point [%.3f, %.3f, %.3f] m\n', ...
    ext_force(1), ext_force(2), ext_force(3), ext_point(1), ext_point(2), ext_point(3));
fprintf('  Moment = [%.4e, %.4e, %.4e] N.m  (about CoG [%.3f, %.3f, %.3f])\n', ...
    F_ext_6dof(4), F_ext_6dof(5), F_ext_6dof(6), CoG(1), CoG(2), CoG(3));

% Save to .mat file
save(fullfile(nemoh_dir, 'timeseries_results.mat'), 'results', ...
    'K_ret', 't_irf', 'alpha_prony', 'beta_prony', 'A_inf', 'M', 'K_h', ...
    'K_mooring', 'F_total_t', 'F_hydro_net_t', 'F_ext_t');
fprintf('Results saved to: timeseries_results.mat\n');

%% ===== PLOT RESULTS =====================================================
fprintf('Plotting results...\n');

% --- Figure 1: Wave Elevation ---
figure('Name', 'Wave Elevation', 'Position', [50 600 900 250]);
plot(t_sol, eta*100, 'b', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('\eta [cm]');
title(sprintf('Wave Elevation (%s)', wave_type)); grid on;
xlim([0 t_end]);

% --- Figure 2: Body Displacements (active DOFs only) ---
n_plot_rows = ceil(n_active / 2);
figure('Name', 'Body Displacements', 'Position', [50 50 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);
    plot(t_sol, q_sol(:,di) * dof_scale(di), 'b', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel(sprintf('%s [%s]', dof_labels{di}, dof_units{di}));
    title(dof_labels{di}); grid on;
    xlim([0 t_end]);
end
sgtitle('Body Displacements', 'FontWeight', 'bold');

% --- Figure 3: Body Velocities (active DOFs only) ---
figure('Name', 'Body Velocities', 'Position', [100 50 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);
    plot(t_sol, v_sol(:,di) * dof_scale(di), 'r', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel(sprintf('[%s]', vel_units{di}));
    title([dof_labels{di}, ' Velocity']); grid on;
    xlim([0 t_end]);
end
sgtitle('Body Velocities', 'FontWeight', 'bold');

% --- Figure 4: Excitation, Radiation & External Forces (active DOFs only) ---
figure('Name', 'Hydrodynamic Forces', 'Position', [150 50 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);
    plot(t_sol, F_exc_t(:,di), 'Color', [0 0.5 0], 'LineWidth', 1); hold on;
    plot(t_sol, -F_rad_t(:,di), 'Color', [0.8 0 0], 'LineWidth', 0.8);
    plot(t_sol, F_ext_t(:,di), 'Color', [0.6 0 0.8], 'LineWidth', 1.2, ...
        'LineStyle', '--');
    xlabel('Time [s]');
    if di <= 3; ylabel('[N]'); else; ylabel('[N.m]'); end
    title(dof_labels{di}); grid on;
    legend('F_{exc}', '-F_{rad}', 'F_{ext}', 'Location', 'best');
    xlim([0 t_end]);
end
sgtitle('Excitation, Radiation and External Forces', 'FontWeight', 'bold');

% --- Figure 4b: Net Hydro Force (excl. mooring, active DOFs only) ---
%   F_hydro_net = F_exc + F_ext - F_rad - K_h*q  (what mooring must resist)
figure('Name', 'Net Hydro Force (excl. mooring)', ...
    'Position', [175 50 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);
    plot(t_sol, F_hydro_net_t(:,di), 'Color', [0.85 0.33 0], 'LineWidth', 1);
    hold on;
    plot(t_sol, F_ext_t(:,di), 'Color', [0.6 0 0.8], 'LineWidth', 1.2, ...
        'LineStyle', '--');
    xlabel('Time [s]');
    if di <= 3; ylabel('[N]'); else; ylabel('[N.m]'); end
    title(sprintf('%s  |  max = %.3e', dof_labels{di}, max(abs(F_hydro_net_t(:,di)))));
    grid on;
    legend('F_{hydro,net}', 'F_{ext}', 'Location', 'best');
    xlim([0 t_end]);
end
sgtitle('Net Hydro Force  (F_{exc} + F_{ext} - F_{rad} - K_h q)  [excl. mooring]', ...
    'FontWeight', 'bold');

% --- Figure 4c: Net Force on Body (incl. mooring, active DOFs only) ---
%   F_net = F_exc + F_ext - F_rad - K_total*q  (= [M+A_inf]*a)
figure('Name', 'Net Force (incl. mooring)', ...
    'Position', [200 50 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);
    plot(t_sol, F_net_t(:,di), 'Color', [0.1 0.1 0.6], 'LineWidth', 1);
    hold on;
    plot(t_sol, F_ext_t(:,di), 'Color', [0.6 0 0.8], 'LineWidth', 1.2, ...
        'LineStyle', '--');
    xlabel('Time [s]');
    if di <= 3; ylabel('[N]'); else; ylabel('[N.m]'); end
    title(sprintf('%s  |  max = %.3e', dof_labels{di}, max(abs(F_net_t(:,di)))));
    grid on;
    legend('F_{net}', 'F_{ext}', 'Location', 'best');
    xlim([0 t_end]);
end
sgtitle('Net Force  (F_{exc} + F_{ext} - F_{rad} - K_{total} q)  [incl. mooring]', ...
    'FontWeight', 'bold');

% --- Figure 5: Prony Fit Verification (active DOFs only) ---
figure('Name', 'Prony Fit Verification', 'Position', [200 100 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    K_ii = squeeze(K_ret(di, di, :));

    subplot(n_plot_rows, min(n_active,2), pi);
    K_fit = zeros(n_t_irf, 1);
    for k = 1:n_prony
        K_fit = K_fit + alpha_prony(di,di,k) * exp(beta_prony(di,di,k) * t_irf);
    end
    plot(t_irf, K_ii, 'b-', 'LineWidth', 1.5); hold on;
    plot(t_irf, K_fit, 'r--', 'LineWidth', 1.2);
    xlabel('Time [s]'); ylabel(sprintf('K_{%d%d}(t)', di, di));
    title(dof_labels{di}); grid on;
    legend('Nemoh', 'Prony', 'Location', 'best');
end
sgtitle('Retardation Function: Prony Fit vs Cosine Transform', 'FontWeight', 'bold');

% --- Figure 6: Radiation Damping Verification (active DOFs only) ---
figure('Name', 'B(w) Prony Verification', 'Position', [250 100 1100 250*n_plot_rows]);
omega_plot = linspace(omega_rad(1), omega_rad(end), 500)';
for pi = 1:n_active
    di = active_dofs(pi);
    B_ii = squeeze(B_rad(di, di, :));

    B_fit = zeros(size(omega_plot));
    for k = 1:n_prony
        ak = alpha_prony(di, di, k);
        bk = beta_prony(di, di, k);
        B_fit = B_fit + ak * (-bk) ./ (bk^2 + omega_plot.^2);
    end

    subplot(n_plot_rows, min(n_active,2), pi);
    plot(omega_rad, B_ii, 'b-', 'LineWidth', 1.5); hold on;
    plot(omega_plot, B_fit, 'r--', 'LineWidth', 1.2);
    xlabel('\omega [rad/s]'); ylabel(sprintf('B_{%d%d}(\\omega)', di, di));
    title(dof_labels{di}); grid on;
    legend('Nemoh', 'Prony', 'Location', 'best');
end
sgtitle('Radiation Damping: Prony Reconstruction vs Nemoh', 'FontWeight', 'bold');

% --- Figure 7: Spectrum (irregular waves only) ---
if strcmpi(wave_type, 'irregular')
    figure('Name', 'Wave Spectrum', 'Position', [300 200 900 400]);

    subplot(1, 2, 1);
    plot(wave_omega, S_wave, 'b', 'LineWidth', 1.5);
    xlabel('\omega [rad/s]'); ylabel('S(\omega) [m^2 s/rad]');
    if irr_gamma > 1
        title(sprintf('JONSWAP Spectrum (\\gamma = %.1f)', irr_gamma));
    else
        title('Pierson-Moskowitz Spectrum');
    end
    grid on;

    subplot(1, 2, 2);
    dt_avg = mean(diff(t_sol));
    fs = 1 / dt_avg;
    nfft = min(2^nextpow2(length(t_sol)/4), length(t_sol));
    [pxx, f_psd] = pwelch(q_sol(:,3), hanning(nfft), nfft/2, nfft, fs);
    plot(2*pi*f_psd, pxx/(2*pi), 'b', 'LineWidth', 1);
    xlabel('\omega [rad/s]'); ylabel('S_x(\omega) [m^2 s/rad]');
    title('Heave Response Spectrum'); grid on;
end

%% ===== CLIP TRANSIENT & COMPUTE STATISTICS =============================
%  Discard the first 20% of the time series (transient build-up) before
%  computing PDFs, exceedance probabilities, and summary statistics.

transient_frac = 0.20;   % fraction of time series to discard
i_start = find(t_sol >= transient_frac * t_end, 1, 'first');
if isempty(i_start); i_start = 1; end

t_ss           = t_sol(i_start:end);
q_ss           = q_sol(i_start:end, :);
v_ss           = v_sol(i_start:end, :);
F_total_ss     = F_total_t(i_start:end, :);
F_net_ss       = F_net_t(i_start:end, :);
F_hydro_net_ss = F_hydro_net_t(i_start:end, :);
F_ext_ss       = F_ext_t(i_start:end, :);

fprintf('\nStatistics computed on steady-state window:\n');
fprintf('  Discarded first %.0f%% (t < %.1f s), using %d of %d samples\n', ...
    transient_frac*100, t_sol(i_start), length(t_ss), length(t_sol));

force_units = {'N','N','N','N.m','N.m','N.m'};

%% ===== PDF OF DISPLACEMENTS AND TOTAL FORCES ===========================
fprintf('Computing probability density functions...\n');

% --- Figure 8: Displacement PDFs ---
fig_pdf_disp = figure('Name', 'Displacement PDF', ...
    'Position', [350 100 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);

    q_i = q_ss(:, di) * dof_scale(di);
    [f_kde, xi_kde] = ksdensity(q_i);
    max_q = max(abs(q_i));
    [~, idx_max] = max(abs(q_i));
    max_q_signed = q_i(idx_max);

    histogram(q_i, 60, 'Normalization', 'pdf', ...
        'FaceColor', [0.7 0.85 1], 'EdgeColor', [0.4 0.6 0.8]); hold on;
    plot(xi_kde, f_kde, 'b-', 'LineWidth', 2);
    yl = ylim;
    plot([max_q_signed max_q_signed], yl, 'r--', 'LineWidth', 1.5);
    plot([-max_q_signed -max_q_signed], yl, 'r--', 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    xlabel(sprintf('%s [%s]', dof_labels{di}, dof_units{di}));
    ylabel('Probability Density');
    title(sprintf('%s  |  max = %.4e %s', dof_labels{di}, max_q, dof_units{di}));
    legend('Histogram', 'KDE', sprintf('Max = %.3e', max_q_signed), ...
        'Location', 'best');
    grid on;
end
sgtitle(sprintf('PDF of Displacements  (t > %.1f s)', t_sol(i_start)), ...
    'FontWeight', 'bold');

% --- Figure 9: Total Force PDFs ---
fig_pdf_force = figure('Name', 'Total Force PDF', ...
    'Position', [400 100 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);

    F_i = F_total_ss(:, di);
    [f_kde, xi_kde] = ksdensity(F_i);
    max_F = max(abs(F_i));
    [~, idx_max] = max(abs(F_i));
    max_F_signed = F_i(idx_max);

    histogram(F_i, 60, 'Normalization', 'pdf', ...
        'FaceColor', [1 0.85 0.7], 'EdgeColor', [0.8 0.5 0.3]); hold on;
    plot(xi_kde, f_kde, 'Color', [0.8 0.2 0], 'LineWidth', 2);
    yl = ylim;
    plot([max_F_signed max_F_signed], yl, 'k--', 'LineWidth', 1.5);

    xlabel(sprintf('%s [%s]', dof_labels{di}, force_units{di}));
    ylabel('Probability Density');
    title(sprintf('%s  |  max = %.4e %s', dof_labels{di}, max_F, force_units{di}));
    legend('Histogram', 'KDE', sprintf('Max = %.3e', max_F_signed), ...
        'Location', 'best');
    grid on;
end
sgtitle(sprintf('PDF of Total Force  (t > %.1f s)', t_sol(i_start)), ...
    'FontWeight', 'bold');

%% ===== EXCEEDANCE PROBABILITY PLOTS ====================================
fprintf('Computing exceedance probabilities...\n');

% --- Figure 10: Displacement Exceedance ---
fig_exc_disp = figure('Name', 'Displacement Exceedance', ...
    'Position', [450 100 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);

    q_i = abs(q_ss(:, di)) * dof_scale(di);
    q_sorted = sort(q_i, 'ascend');
    n_pts = length(q_sorted);
    P_exc = 1 - (1:n_pts)' / (n_pts + 1);  % Weibull plotting position

    semilogy(q_sorted, P_exc, 'b-', 'LineWidth', 1.5); hold on;

    % Mark key exceedance levels
    levels = [0.01, 0.05, 0.10, 0.50];
    colors_exc = {[0.8 0 0], [1 0.4 0], [0.9 0.6 0], [0.5 0.5 0.5]};
    for li = 1:length(levels)
        idx = find(P_exc <= levels(li), 1, 'first');
        if ~isempty(idx)
            plot(q_sorted(idx), levels(li), 'o', ...
                'MarkerSize', 8, 'MarkerFaceColor', colors_exc{li}, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
            text(q_sorted(idx)*1.05, levels(li), ...
                sprintf(' %.0f%%: %.3e', levels(li)*100, q_sorted(idx)), ...
                'FontSize', 8, 'Color', colors_exc{li});
        end
    end

    xlabel(sprintf('|%s| [%s]', dof_labels{di}, dof_units{di}));
    ylabel('Exceedance Probability P(X > x)');
    title(dof_labels{di}); grid on;
    ylim([1/(2*n_pts) 1]);
end
sgtitle(sprintf('Displacement Exceedance Probability  (t > %.1f s)', ...
    t_sol(i_start)), 'FontWeight', 'bold');

% --- Figure 11: Total Force Exceedance ---
fig_exc_force = figure('Name', 'Force Exceedance', ...
    'Position', [500 100 1100 250*n_plot_rows]);
for pi = 1:n_active
    di = active_dofs(pi);
    subplot(n_plot_rows, min(n_active,2), pi);

    F_i = abs(F_total_ss(:, di));
    F_sorted = sort(F_i, 'ascend');
    n_pts = length(F_sorted);
    P_exc = 1 - (1:n_pts)' / (n_pts + 1);

    semilogy(F_sorted, P_exc, 'Color', [0.8 0.2 0], 'LineWidth', 1.5); hold on;

    for li = 1:length(levels)
        idx = find(P_exc <= levels(li), 1, 'first');
        if ~isempty(idx)
            plot(F_sorted(idx), levels(li), 'o', ...
                'MarkerSize', 8, 'MarkerFaceColor', colors_exc{li}, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
            text(F_sorted(idx)*1.05, levels(li), ...
                sprintf(' %.0f%%: %.3e', levels(li)*100, F_sorted(idx)), ...
                'FontSize', 8, 'Color', colors_exc{li});
        end
    end

    xlabel(sprintf('|%s| [%s]', dof_labels{di}, force_units{di}));
    ylabel('Exceedance Probability P(X > x)');
    title(dof_labels{di}); grid on;
    ylim([1/(2*n_pts) 1]);
end
sgtitle(sprintf('Total Force Exceedance Probability  (t > %.1f s)', ...
    t_sol(i_start)), 'FontWeight', 'bold');

% --- Print summary statistics table (steady-state only) ---
fprintf('\n==================== RESPONSE STATISTICS (steady-state) ====================\n');
fprintf('  Transient clip: first %.0f%% discarded (t < %.1f s)\n', ...
    transient_frac*100, t_sol(i_start));
fprintf('%-8s  %12s  %12s  %12s  %12s\n', ...
    'DOF', 'Max |Disp|', 'Std Disp', 'Max |Force|', 'Std Force');
fprintf('----------------------------------------------------------------------------\n');
for pi = 1:n_active
    di = active_dofs(pi);
    sc = dof_scale(di);

    % Exceedance values at 1% and 5%
    q_i_abs = sort(abs(q_ss(:,di)) * sc, 'descend');
    F_i_abs = sort(abs(F_total_ss(:,di)), 'descend');
    n_ss = length(q_i_abs);
    q_1pct = q_i_abs(max(1, round(0.01*n_ss)));
    q_5pct = q_i_abs(max(1, round(0.05*n_ss)));
    F_1pct = F_i_abs(max(1, round(0.01*n_ss)));
    F_5pct = F_i_abs(max(1, round(0.05*n_ss)));

    fprintf('%-8s  %10.4e %s  %10.4e  %10.4e %s  %10.4e\n', ...
        dof_labels{di}, ...
        max(abs(q_ss(:,di))) * sc, dof_units{di}, ...
        std(q_ss(:,di)) * sc, ...
        max(abs(F_total_ss(:,di))), force_units{di}, ...
        std(F_total_ss(:,di)));
    fprintf('          1%%exc: %.3e      5%%exc: %.3e   1%%exc: %.3e      5%%exc: %.3e\n', ...
        q_1pct, q_5pct, F_1pct, F_5pct);
end
fprintf('============================================================================\n');

% --- Print mean values (to verify external force effect) ---
fprintf('\n==================== MEAN VALUES (steady-state) ============================\n');
fprintf('%-8s  %12s  %12s  %12s  %12s  %12s\n', ...
    'DOF', 'Mean Disp', 'Mean F_ext', 'Mean F_hnet', 'Mean F_net', 'F_ext(t=0)');
fprintf('----------------------------------------------------------------------------\n');
F_ext_check = F_external(0);
for pi = 1:n_active
    di = active_dofs(pi);
    fprintf('%-8s  %10.4e %s  %10.4e  %10.4e  %10.4e  %10.4e\n', ...
        dof_labels{di}, ...
        mean(q_ss(:,di)) * dof_scale(di), dof_units{di}, ...
        mean(F_ext_ss(:,di)), ...
        mean(F_hydro_net_ss(:,di)), ...
        mean(F_net_ss(:,di)), ...
        F_ext_check(di));
end
fprintf('  Expected: Mean F_net ~ 0 (= M*a).  Mean displacement ~ F_ext / K_total.\n');
fprintf('  F_hydro_net mean ~ F_ext * K_moor/K_total  (surge), or ~ 0 if K_h dominates (pitch).\n');
fprintf('============================================================================\n');

% --- Export PDF figures to file ---
pdf_dir = fullfile(nemoh_dir, 'pdf_outputs');
if ~exist(pdf_dir, 'dir'); mkdir(pdf_dir); end

exportgraphics(fig_pdf_disp, fullfile(pdf_dir, 'displacement_pdf.pdf'), ...
    'ContentType', 'vector', 'Resolution', 300);
exportgraphics(fig_pdf_force, fullfile(pdf_dir, 'total_force_pdf.pdf'), ...
    'ContentType', 'vector', 'Resolution', 300);
exportgraphics(fig_exc_disp, fullfile(pdf_dir, 'displacement_exceedance.pdf'), ...
    'ContentType', 'vector', 'Resolution', 300);
exportgraphics(fig_exc_force, fullfile(pdf_dir, 'force_exceedance.pdf'), ...
    'ContentType', 'vector', 'Resolution', 300);
fprintf('PDF files saved to:\n  %s\n', pdf_dir);

fprintf('\n========================================\n');
fprintf(' Simulation complete.\n');
fprintf('========================================\n');

%% =======================================================================
%                        LOCAL FUNCTIONS
%  =======================================================================

function M = load_matrix(filepath)
%LOAD_MATRIX  Load a 6x6 matrix from Nemoh Mechanics file
%   Reads Inertia.dat, Kh.dat, or similar whitespace-delimited 6x6 files.
    M = zeros(6);
    fid = fopen(filepath, 'r');
    if fid == -1
        warning('load_matrix: Cannot open %s. Returning zeros.', filepath);
        return;
    end
    for row = 1:6
        line = fgetl(fid);
        if ~ischar(line); break; end
        vals = sscanf(line, '%e');
        n = min(length(vals), 6);
        M(row, 1:n) = vals(1:n)';
    end
    fclose(fid);
end


function [omega, A, B] = load_radiation(filepath)
%LOAD_RADIATION  Load radiation coefficients from Nemoh RadiationCoefficients.tec
%   Returns:
%     omega : [n_freq x 1]    frequencies [rad/s]
%     A     : [6, 6, n_freq]  added mass A_ij(omega)
%     B     : [6, 6, n_freq]  radiation damping B_ij(omega)
%
%   File has 6 zones (one per excitation DOF). Each zone contains 13 columns:
%     omega, A(resp1), B(resp1), A(resp2), B(resp2), ..., A(resp6), B(resp6)

    n_dof = 6;
    fid = fopen(filepath, 'r');
    if fid == -1; error('Cannot open %s', filepath); end

    raw_data = {};
    zone_count = 0;
    current_data = [];

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line); break; end
        line = strtrim(line);
        if isempty(line); continue; end

        % Skip VARIABLES header lines
        if strncmpi(line, 'VARIABLES', 9) || line(1) == '"'
            continue;
        end

        % Detect Zone header
        if strncmpi(line, 'Zone', 4)
            if zone_count > 0
                raw_data{zone_count} = current_data; %#ok<AGROW>
            end
            zone_count = zone_count + 1;
            current_data = [];
            continue;
        end

        % Parse data
        vals = sscanf(line, '%e');
        if ~isempty(vals)
            current_data = [current_data; vals']; %#ok<AGROW>
        end
    end
    % Store last zone
    if zone_count > 0
        raw_data{zone_count} = current_data;
    end
    fclose(fid);

    % Extract omega from first zone
    omega = raw_data{1}(:, 1);
    n_f = length(omega);
    A = zeros(n_dof, n_dof, n_f);
    B = zeros(n_dof, n_dof, n_f);

    for j = 1:min(zone_count, n_dof)  % j = excitation DOF (zone)
        data = raw_data{j};
        for k = 1:n_dof  % k = response DOF
            col_A = 2*(k-1) + 2;
            col_B = 2*(k-1) + 3;
            if col_B <= size(data, 2)
                A(k, j, :) = data(:, col_A);
                B(k, j, :) = data(:, col_B);
            end
        end
    end
end


function [omega, F_mag, F_phase] = load_excitation(filepath)
%LOAD_EXCITATION  Load excitation force from Nemoh ExcitationForce.tec
%   Returns:
%     omega   : [n_freq x 1]   frequencies [rad/s]
%     F_mag   : [6, n_freq]    excitation force magnitudes
%     F_phase : [6, n_freq]    excitation force phases [rad]
%
%   File contains columns:
%     omega, |F1|, angle(F1), |F2|, angle(F2), ..., |F6|, angle(F6)

    n_dof = 6;
    fid = fopen(filepath, 'r');
    if fid == -1; error('Cannot open %s', filepath); end

    data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line); break; end
        line = strtrim(line);
        if isempty(line); continue; end
        if strncmpi(line, 'VARIABLES', 9) || line(1) == '"'; continue; end
        if strncmpi(line, 'Zone', 4); continue; end

        vals = sscanf(line, '%e');
        if ~isempty(vals)
            data = [data; vals']; %#ok<AGROW>
        end
    end
    fclose(fid);

    omega = data(:, 1);
    n_f = length(omega);
    F_mag   = zeros(n_dof, n_f);
    F_phase = zeros(n_dof, n_f);

    for i = 1:n_dof
        col_mag   = 2*(i-1) + 2;
        col_phase = 2*(i-1) + 3;
        if col_phase <= size(data, 2)
            F_mag(i, :)   = data(:, col_mag);
            F_phase(i, :) = data(:, col_phase);  % radians
        end
    end
end


function [alpha, beta] = fit_prony(t, K_vals, n_terms)
%FIT_PRONY  Fit retardation function with real decaying exponentials
%   K(t) ~ sum_k alpha(k) * exp(beta(k) * t)
%
%   Uses fixed logarithmically-spaced poles (time constants) with linear
%   least-squares for the coefficients. All beta values are negative real,
%   ensuring stability.
%
%   Inputs:
%     t       : [n x 1] time vector (starting from 0 or near 0)
%     K_vals  : [n x 1] retardation function values
%     n_terms : number of exponential terms
%
%   Outputs:
%     alpha : [n_terms x 1] coefficients
%     beta  : [n_terms x 1] exponents (all negative real)

    t = t(:);
    K_vals = K_vals(:);

    % Logarithmically spaced time constants from dt to t_max
    t_max = t(end);
    t_min = max(t(2) - t(1), t_max / 1000);
    tau = logspace(log10(t_min), log10(t_max), n_terms);
    beta = -1 ./ tau(:);

    % Build basis matrix: Phi(i,k) = exp(beta(k) * t(i))
    Phi = exp(t * beta');  % [n_t x n_terms]

    % Solve via least squares (may use non-negative or unconstrained)
    alpha = Phi \ K_vals;

    alpha = alpha(:);
    beta  = beta(:);
end


function S = jonswap_spectrum(omega, Hs, Tp, gamma)
%JONSWAP_SPECTRUM  JONSWAP (or Pierson-Moskowitz) wave spectrum
%   S(omega) in [m^2 s/rad]
%   Set gamma = 1.0 for Pierson-Moskowitz spectrum.
%
%   S(w) = alpha * (5/16) * Hs^2 * wp^4 / w^5 * exp(-5/4*(wp/w)^4) * gamma^r
%   where r = exp(-(w-wp)^2 / (2*sigma^2*wp^2))
%         sigma = 0.07 (w <= wp), 0.09 (w > wp)
%         alpha = 1 - 0.287*ln(gamma)

    omega = omega(:);
    wp = 2*pi / Tp;

    % Spectral width parameter
    sigma = 0.09 * ones(size(omega));
    sigma(omega <= wp) = 0.07;

    % Peak enhancement exponent
    r = exp(-(omega - wp).^2 ./ (2 * sigma.^2 * wp^2));

    % Normalizing coefficient for Hs
    alpha_norm = 1 - 0.287 * log(gamma);

    % Pierson-Moskowitz base spectrum
    S_pm = (5/16) * Hs^2 * wp^4 ./ omega.^5 .* exp(-5/4 * (wp./omega).^4);

    % JONSWAP = PM * peak enhancement
    S = alpha_norm * S_pm .* gamma.^r;

    % Avoid numerical issues at very low frequencies
    S(omega < wp/10) = 0;
    S(~isfinite(S)) = 0;
end


function dxdt = cummins_rhs(t, x, p)
%CUMMINS_RHS  Right-hand side of Cummins equation with Prony states
%
%   State vector:  x = [q(6); v(6); I_prony(6*6*n_p)]
%
%   Equations:
%     dq/dt = v
%     [M+A_inf] dv/dt = F_exc(t) + F_ext(t) - K_total*q - mu
%     dI(i,j,k)/dt = beta(i,j,k)*I(i,j,k) + alpha(i,j,k)*v(j)
%     mu(i) = sum_j sum_k I(i,j,k)

    n_dof = 6;
    n_p   = p.n_prony;

    q = x(1:6);
    v = x(7:12);

    % --- Wave excitation force (frequency-domain superposition) ---
    % F_exc_i(t) = sum_n A_n |F_i(wn)| cos(wn*t + eps_n + phi_i(wn))
    %            = F_cos(i,:) * cos(base_phase) - F_sin(i,:) * sin(base_phase)
    base_phase = p.wave_omega * t + p.wave_eps;  % [n_waves x 1]
    cv = cos(base_phase);
    sv = sin(base_phase);
    F_exc = p.F_cos_pre * cv - p.F_sin_pre * sv;  % [6 x 1]

    % --- External force ---
    F_ext = p.F_external(t);

    % --- Radiation memory from Prony states ---
    I_3d = reshape(x(13:end), [n_dof, n_dof, n_p]);  % I_3d(i,j,k)
    mu = sum(sum(I_3d, 3), 2);  % sum over k then j -> [6 x 1]

    % --- Velocity equation ---
    dv = p.M_inv * (F_exc + F_ext - p.K_total * q - mu);

    % --- Enforce disabled DOF constraints ---
    dv(p.disabled_dofs) = 0;
    v(p.disabled_dofs) = 0;

    % --- Prony state derivatives ---
    % dI(i,j,k)/dt = beta(i,j,k)*I(i,j,k) + alpha(i,j,k)*v(j)
    v_exp = repmat(reshape(v, [1, n_dof, 1]), [n_dof, 1, n_p]);
    dI_3d = p.beta_prony .* I_3d + p.alpha_prony .* v_exp;

    dxdt = [v; dv; dI_3d(:)];
end


function k = solve_dispersion(omega, d, g)
%SOLVE_DISPERSION  Solve the linear dispersion relation for wavenumber k
%   omega^2 = g * k * tanh(k * d)
%
%   Uses Newton-Raphson iteration with deep-water initial guess.
%   Handles scalar or vector omega input.

    omega = omega(:);
    k = omega.^2 / g;   % deep-water initial guess

    for iter = 1:50
        kd  = k * d;
        thk = tanh(kd);
        f   = omega.^2 - g * k .* thk;
        df  = -g * (thk + k * d .* (1 - thk.^2));   % sech^2 = 1 - tanh^2
        dk  = -f ./ df;
        k   = k + dk;
        if max(abs(dk ./ max(k, 1e-12))) < 1e-12
            break;
        end
    end
end


function F_ext_fn = build_ext_force(force_input, point, cog)
%BUILD_EXT_FORCE  Build a 6-DOF external force function from force + application point
%   Given a 3D force and its application point, returns a function handle
%   that produces the full [Fx;Fy;Fz;Mx;My;Mz] vector including moments
%   about the centre of gravity:  M = cross(r, F),  r = point - cog.
%
%   force_input : 3x1 constant vector [N], OR function handle f(t)->3x1
%   point       : 3x1 application point in body-fixed frame [m]
%   cog         : 3x1 centre of gravity in body-fixed frame [m]

    r = point(:) - cog(:);   % moment arm from CoG to application point
    if isa(force_input, 'function_handle')
        F_ext_fn = @(t) force_and_moment(force_input(t), r);
    else
        F6 = force_and_moment(force_input(:), r);
        F_ext_fn = @(t) F6;   % constant — avoid repeated cross product
    end
end


function F6 = force_and_moment(F3, r)
%FORCE_AND_MOMENT  Combine 3D force with moment about CoG
%   F6 = [Fx; Fy; Fz; Mx; My; Mz]  where  M = cross(r, F3)
    F3 = F3(:);
    M  = cross(r, F3);
    F6 = [F3; M];
end
