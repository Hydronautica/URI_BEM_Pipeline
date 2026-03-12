%% =======================================================================
%  PARAMETRIC Hs-Tp SWEEP: 1% Exceedance Force Contour Maps
%  =======================================================================
%  Loops over a grid of (Hs, Tp) sea-state combinations, runs a Cummins
%  equation time-domain simulation for each, and plots 2D contour maps of
%  the 1% exceedance probability for two force quantities:
%
%  1) F_hydro_net = F_exc + F_ext - F_rad - K_h * q
%     Net hydrodynamic loading (excl. mooring restoring).
%     This is what the mooring system must resist.
%
%  2) F_net = F_exc + F_ext - F_rad - K_total * q   (K_total = K_h + K_moor)
%           = [M + A_inf] * a
%     True net (unbalanced) force on the body, including mooring.
%
%  The 1% exceedance value is the threshold exceeded by only 1% of the
%  steady-state |F| samples (after discarding the initial transient).
%  =======================================================================

clear; clc; close all;

%% ===== SWEEP CONFIGURATION =============================================

% --- Sea-state grid ---------------------------------------------------
Hs_vec = linspace(0.005, 0.050, 5);   % Significant wave height [m]
Tp_vec = linspace(0.6,   3.0,   5);   % Peak period [s]

% --- JONSWAP parameters (fixed across sweep) --------------------------
irr_gamma = 1.0;       % Peak enhancement factor (1 = PM)
irr_ncomp = 128;       % Number of wave components
irr_wmin  = 0.8;       % Min frequency [rad/s]
irr_wmax  = 15.0;      % Max frequency [rad/s]
irr_seed  = 42;        % Random seed for phase generation
% Wave components are uniformly spaced in wavenumber (k) space.
% The dispersion relation w^2 = g*k*tanh(k*d) maps each k_n to a
% non-uniformly spaced w_n. Because these frequencies are incommensurate,
% the wave signal does not repeat, giving good statistics with few
% components (128 is typically sufficient).

% --- Simulation -------------------------------------------------------
t_end   = 1800;          % Simulation duration [s]
dt_out  = 0.01;         % Output time step [s]
transient_frac = 0.20;  % Fraction of time series to discard (20%)

% --- Prony fit --------------------------------------------------------
n_prony   = 12;
t_irf_max = 3;
dt_irf    = 0.01;
prony_tol = 1e-3;

% --- Active DOFs ------------------------------------------------------
active_dofs = [1, 3, 5];   % 1=Surge, 3=Heave, 5=Pitch

% --- ODE Solver Options -----------------------------------------------
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.02);

% --- Nemoh Case Directory ---------------------------------------------
nemoh_dir = fullfile('C:','Users','Fonta','OneDrive','Desktop', ...
    'Florid_Scale_Model_OpenFAST','Nemoh','SimplePontoon');

% --- Environment ------------------------------------------------------
water_depth = 1.0;      % Water depth [m] (used for dispersion & mooring)
g_acc       = 9.81;     % Gravitational acceleration [m/s^2]

% --- Mooring Configuration --------------------------------------------
K_mooring = zeros(6);
K_mooring(1,1) = 50;    % Surge restoring [N/m]

depth = water_depth;  R_anch = 1.5;  R_fair = 0.15;  z_fair = -0.01;
for i = 1:3
    ang = (i-1)*120 * pi/180;
    lines(i).type     = 'taut';
    lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
    lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); z_fair];
    lines(i).L0 = 1.6;
    lines(i).EA = 1000;
end
[K_mooring, ~] = compute_mooring_stiffness(lines, zeros(6,1));

% --- External Force with Application Point ----------------------------
%  Specify a 3D force and where it acts. Moments about CoG are computed
%  automatically:  M = cross(r, F),  r = application_point - CoG.
%
%  ext_force: 3x1 constant [N], or function handle f(t)->3x1
%  ext_point: [x; y; z] application point [m]
%  CoG:       [x; y; z] centre of gravity [m]
CoG       = [0; 0; -0.01];
ext_force = [0.36; 0; 0];
ext_point = [0; 0; -0.2];
% Example: ext_force = [0.3; 0; 0]; ext_point = [0; 0; -0.35];
F_external = build_ext_force(ext_force, ext_point, CoG);

% --- Labels -----------------------------------------------------------
dof_labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
dof_units  = {'m','m','m','deg','deg','deg'};
dof_scale  = [1, 1, 1, 180/pi, 180/pi, 180/pi];
vel_units  = {'m/s','m/s','m/s','deg/s','deg/s','deg/s'};
force_units = {'N','N','N','N.m','N.m','N.m'};

%% ===== ONE-TIME SETUP (Nemoh data, Prony fit) ==========================
fprintf('========================================\n');
fprintf(' Hs-Tp Parametric Sweep\n');
fprintf('========================================\n\n');
fprintf('Grid: %d Hs x %d Tp = %d simulations\n', ...
    length(Hs_vec), length(Tp_vec), length(Hs_vec)*length(Tp_vec));

% Print external force info
F_ext_6dof = F_external(0);
fprintf('External force (6-DOF): [');
fprintf(' %.4e', F_ext_6dof);
fprintf(' ]\n');
fprintf('  Force  = [%.4f, %.4f, %.4f] N  at point [%.3f, %.3f, %.3f] m\n', ...
    ext_force(1), ext_force(2), ext_force(3), ext_point(1), ext_point(2), ext_point(3));
fprintf('  Moment = [%.4e, %.4e, %.4e] N.m  (about CoG [%.3f, %.3f, %.3f])\n', ...
    F_ext_6dof(4), F_ext_6dof(5), F_ext_6dof(6), CoG(1), CoG(2), CoG(3));

fprintf('Loading Nemoh data from:\n  %s\n\n', nemoh_dir);

% --- Load Nemoh outputs ---
M = load_matrix(fullfile(nemoh_dir, 'Mechanics', 'Inertia.dat'));
K_h = load_matrix(fullfile(nemoh_dir, 'Mechanics', 'Kh.dat'));
[omega_rad, A_rad, B_rad] = load_radiation( ...
    fullfile(nemoh_dir, 'results', 'RadiationCoefficients.tec'));
[omega_exc, F_exc_mag, F_exc_phase] = load_excitation( ...
    fullfile(nemoh_dir, 'results', 'ExcitationForce.tec'));

fprintf('  Frequencies: %d [%.2f - %.2f] rad/s\n', ...
    length(omega_rad), omega_rad(1), omega_rad(end));

% --- Apply DOF constraints ---
dof_mask = false(6,1);
dof_mask(active_dofs) = true;
disabled_dofs = find(~dof_mask);
n_active = length(active_dofs);

M(disabled_dofs, :) = 0;  M(:, disabled_dofs) = 0;
K_h(disabled_dofs, :) = 0;  K_h(:, disabled_dofs) = 0;
K_mooring(disabled_dofs, :) = 0;  K_mooring(:, disabled_dofs) = 0;
for dd = disabled_dofs'
    M(dd, dd) = 1;
end
A_rad(disabled_dofs, :, :) = 0;  A_rad(:, disabled_dofs, :) = 0;
B_rad(disabled_dofs, :, :) = 0;  B_rad(:, disabled_dofs, :) = 0;
F_exc_mag(disabled_dofs, :) = 0;
F_exc_phase(disabled_dofs, :) = 0;

% --- Retardation function K(t) from B(omega) ---
fprintf('Computing retardation functions...\n');
t_irf = (0 : dt_irf : t_irf_max)';
n_t_irf = length(t_irf);
K_ret = zeros(6, 6, n_t_irf);
cos_mat = cos(omega_rad(:) * t_irf(:)');
for i = 1:6
    for j = 1:6
        B_ij = squeeze(B_rad(i, j, :));
        integrand = B_ij(:) .* cos_mat;
        K_ret(i, j, :) = (2/pi) * trapz(omega_rad, integrand);
    end
end

% --- Infinite-frequency added mass ---
A_inf = A_rad(:, :, end);

% --- Prony fit ---
fprintf('Fitting Prony coefficients (%d terms)...\n', n_prony);
alpha_prony = zeros(6, 6, n_prony);
beta_prony  = zeros(6, 6, n_prony);
K_max_global = max(abs(K_ret(:)));
for i = 1:6
    for j = 1:6
        K_ij = squeeze(K_ret(i, j, :));
        if max(abs(K_ij)) < prony_tol * K_max_global; continue; end
        [a, b] = fit_prony(t_irf, K_ij, n_prony);
        alpha_prony(i, j, :) = a;
        beta_prony(i, j, :)  = b;
    end
end

% --- Pre-compute shared matrices ---
M_total = M + A_inf;
M_inv   = M_total \ eye(6);
K_total = K_h + K_mooring;

n_dof = 6;
n_prony_states = n_dof * n_dof * n_prony;
n_states = 2*n_dof + n_prony_states;

% --- Wave discretisation info (k-space) ---
k_min_info = solve_dispersion(irr_wmin, water_depth, g_acc);
k_max_info = solve_dispersion(irr_wmax, water_depth, g_acc);
dk_info    = (k_max_info - k_min_info) / irr_ncomp;
fprintf('Setup complete. State vector: %d states\n', n_states);
fprintf('Wave discretisation: %d components in k-space [%.4f, %.4f] rad/m  (dk = %.4f)\n', ...
    irr_ncomp, k_min_info, k_max_info, dk_info);
fprintf('  -> non-uniform omega in [%.3f, %.3f] rad/s  (depth = %.2f m)\n\n', ...
    irr_wmin, irr_wmax, water_depth);

%% ===== PARAMETRIC SWEEP ================================================
n_Hs = length(Hs_vec);
n_Tp = length(Tp_vec);

% Result matrices: 1% exceedance for each active DOF
%   F_hydro_net = F_exc + F_ext - F_rad - K_h*q   (what mooring must resist)
%   F_net       = F_exc + F_ext - F_rad - K_total*q (true net force = [M+A_inf]*a)
exc_1pct_hydro = zeros(n_Hs, n_Tp, 6);   % 1% exceedance of |F_hydro_net|
exc_1pct_net   = zeros(n_Hs, n_Tp, 6);   % 1% exceedance of |F_net|

% Also store max and std for reference
max_hydro = zeros(n_Hs, n_Tp, 6);
std_hydro = zeros(n_Hs, n_Tp, 6);
max_net   = zeros(n_Hs, n_Tp, 6);
std_net   = zeros(n_Hs, n_Tp, 6);

total_runs = n_Hs * n_Tp;
run_count = 0;
t_sweep_start = tic;

fprintf('Starting sweep: %d runs\n', total_runs);
fprintf('%-5s  %-8s  %-8s', 'Run', 'Hs [m]', 'Tp [s]');
for pi = 1:n_active
    fprintf('  %21s', [dof_labels{active_dofs(pi)} ' hydro/net']);
end
fprintf('  %8s\n', 'Time [s]');
fprintf('%s\n', repmat('-', 1, 35 + 23*n_active));

tspan = 0 : dt_out : t_end;
i_start_clip = find(tspan >= transient_frac * t_end, 1, 'first');
if isempty(i_start_clip); i_start_clip = 1; end

for iH = 1:n_Hs
    for iT = 1:n_Tp
        run_count = run_count + 1;
        t_run = tic;

        Hs = Hs_vec(iH);
        Tp = Tp_vec(iT);

        % --- Generate JONSWAP spectrum (k-space discretisation) ---
        %  Uniform spacing in wavenumber k, mapped to non-uniform omega
        %  via the dispersion relation: w^2 = g*k*tanh(k*d).
        %  Incommensurate frequencies → signal does not repeat.
        rng(irr_seed);
        k_min = solve_dispersion(irr_wmin, water_depth, g_acc);
        k_max = solve_dispersion(irr_wmax, water_depth, g_acc);
        dk    = (k_max - k_min) / irr_ncomp;
        k_vec = k_min + ((1:irr_ncomp)' - 0.5) * dk;   % midpoints

        % Dispersion: omega_n and group velocity Cg_n
        wave_omega = sqrt(g_acc * k_vec .* tanh(k_vec * water_depth));
        Cg = (wave_omega ./ (2*k_vec)) .* ...
             (1 + 2*k_vec*water_depth ./ sinh(2*k_vec*water_depth));
        dw_vec = Cg * dk;              % non-uniform dw per component
        n_waves = irr_ncomp;
        wave_eps = 2*pi * rand(n_waves, 1);

        S_wave   = jonswap_spectrum(wave_omega, Hs, Tp, irr_gamma);
        wave_amp = sqrt(2 * S_wave .* dw_vec);

        % --- Interpolate excitation force to wave frequencies ---
        F_mag_w   = zeros(6, n_waves);
        F_phase_w = zeros(6, n_waves);
        for i = 1:6
            F_mag_w(i,:)   = interp1(omega_exc, F_exc_mag(i,:), wave_omega(:), 'linear', 0);
            F_phase_w(i,:) = interp1(omega_exc, F_exc_phase(i,:), wave_omega(:), 'linear', 0);
        end

        F_amp_wave = F_mag_w .* wave_amp(:)';
        F_cos_pre  = F_amp_wave .* cos(F_phase_w);
        F_sin_pre  = F_amp_wave .* sin(F_phase_w);

        % --- Pack params and solve ODE ---
        params = struct();
        params.M_inv       = M_inv;
        params.K_total     = K_total;
        params.n_prony     = n_prony;
        params.alpha_prony = alpha_prony;
        params.beta_prony  = beta_prony;
        params.n_waves     = n_waves;
        params.wave_omega  = wave_omega(:);
        params.wave_eps    = wave_eps(:);
        params.F_cos_pre   = F_cos_pre;
        params.F_sin_pre   = F_sin_pre;
        params.F_external  = F_external;
        params.disabled_dofs = disabled_dofs;

        x0 = zeros(n_states, 1);

        [t_sol, x_sol] = ode45(@(t, x) cummins_rhs(t, x, params), ...
            tspan, x0, ode_opts);

        % --- Extract solution ---
        q_sol = x_sol(:, 1:6);

        % --- Reconstruct forces (steady-state portion only) ---
        idx_ss = i_start_clip : length(t_sol);
        n_ss   = length(idx_ss);

        F_exc_ss  = zeros(n_ss, 6);
        F_rad_ss  = zeros(n_ss, 6);

        for it_local = 1:n_ss
            it = idx_ss(it_local);
            % Excitation
            base_phase = wave_omega(:) * t_sol(it) + wave_eps(:);
            cv = cos(base_phase);
            sv = sin(base_phase);
            F_exc_ss(it_local, :) = (F_cos_pre * cv - F_sin_pre * sv)';
            % Radiation memory
            I_3d = reshape(x_sol(it, 13:end)', [6, 6, n_prony]);
            F_rad_ss(it_local, :) = sum(sum(I_3d, 3), 2)';
        end

        q_ss = q_sol(idx_ss, :);

        % --- Compute force quantities ---
        %  F_hydro_net = F_exc + F_ext - F_rad - K_h*q
        %               (what mooring must resist, excludes mooring restoring)
        %  F_net       = F_exc + F_ext - F_rad - K_total*q
        %               (true net force on body = [M+A_inf]*a)
        F_hydro_net = zeros(n_ss, 6);
        F_net_ss    = zeros(n_ss, 6);
        for it_local = 1:n_ss
            F_ext = F_external(t_sol(idx_ss(it_local)));
            rhs_hydro = F_exc_ss(it_local,:)' + F_ext ...
                - F_rad_ss(it_local,:)' - K_h * q_ss(it_local,:)';
            F_hydro_net(it_local, :) = rhs_hydro';
            F_net_ss(it_local, :)    = (rhs_hydro - K_mooring * q_ss(it_local,:)')';
        end

        % --- Compute 1% exceedance for each DOF ---
        for di = 1:6
            % F_hydro_net (force minus mooring)
            Fh_abs = sort(abs(F_hydro_net(:, di)), 'descend');
            idx_1pct = max(1, round(0.01 * n_ss));
            exc_1pct_hydro(iH, iT, di) = Fh_abs(idx_1pct);
            max_hydro(iH, iT, di) = Fh_abs(1);
            std_hydro(iH, iT, di) = std(F_hydro_net(:, di));

            % F_net (true net force on body)
            Fn_abs = sort(abs(F_net_ss(:, di)), 'descend');
            exc_1pct_net(iH, iT, di) = Fn_abs(idx_1pct);
            max_net(iH, iT, di) = Fn_abs(1);
            std_net(iH, iT, di) = std(F_net_ss(:, di));
        end

        % --- Progress ---
        elapsed_run = toc(t_run);
        fprintf('%-5d  %-8.4f  %-8.3f', run_count, Hs, Tp);
        for pi = 1:n_active
            di = active_dofs(pi);
            fprintf('  %10.3e/%10.3e', ...
                exc_1pct_hydro(iH, iT, di), exc_1pct_net(iH, iT, di));
        end
        fprintf('  %8.2f\n', elapsed_run);
    end
end

elapsed_sweep = toc(t_sweep_start);
fprintf('\nSweep complete: %d runs in %.1f s (avg %.2f s/run)\n', ...
    total_runs, elapsed_sweep, elapsed_sweep/total_runs);

%% ===== SAVE SWEEP RESULTS ==============================================
F_ext_6dof_save = F_external(0);
save(fullfile(nemoh_dir, 'sweep_results.mat'), ...
    'Hs_vec', 'Tp_vec', 'exc_1pct_hydro', 'exc_1pct_net', ...
    'max_hydro', 'std_hydro', 'max_net', 'std_net', ...
    'active_dofs', 'dof_labels', 'force_units', ...
    'K_mooring', 'K_h', 'K_total', 'irr_gamma', 'transient_frac', ...
    'ext_force', 'ext_point', 'CoG', 'F_ext_6dof_save');
fprintf('Results saved to: sweep_results.mat\n');

%% ===== CONTOUR PLOTS ===================================================
fprintf('Plotting contour maps...\n');

[Tp_grid, Hs_grid] = meshgrid(Tp_vec, Hs_vec);

pdf_dir = fullfile(nemoh_dir, 'pdf_outputs');
if ~exist(pdf_dir, 'dir'); mkdir(pdf_dir); end

% --- Figures A: Net Hydrodynamic Force (F_hydro_net, mooring excluded) ---
fig_hydro = gobjects(n_active, 1);
for pi = 1:n_active
    di = active_dofs(pi);
    Z = exc_1pct_hydro(:, :, di);

    fig_hydro(pi) = figure('Name', ...
        sprintf('1%% Exc Hydro - %s', dof_labels{di}), ...
        'Position', [50 + 50*pi, 150, 700, 550]);

    [C_data, h_contour] = contourf(Tp_grid, Hs_grid, Z, 12);
    clabel(C_data, h_contour, 'FontSize', 8, 'Color', 'k');
    colormap(jet);
    cb = colorbar;
    cb.Label.String = sprintf('1%% Exceedance |F_{hydro,net}| [%s]', force_units{di});
    cb.Label.FontSize = 11;

    xlabel('T_p [s]', 'FontSize', 12);
    ylabel('H_s [m]', 'FontSize', 12);
    title(sprintf('%s: 1%% Exc. — Net Hydro Force (excl. mooring)', ...
        dof_labels{di}), 'FontSize', 12);
    grid on;
    set(gca, 'Layer', 'top');

    fname = sprintf('sweep_1pct_hydro_%s.pdf', lower(dof_labels{di}));
    exportgraphics(fig_hydro(pi), fullfile(pdf_dir, fname), ...
        'ContentType', 'vector', 'Resolution', 300);
end

% --- Figures B: True Net Force on Body (F_net, mooring included) ---
fig_net = gobjects(n_active, 1);
for pi = 1:n_active
    di = active_dofs(pi);
    Z = exc_1pct_net(:, :, di);

    fig_net(pi) = figure('Name', ...
        sprintf('1%% Exc Net - %s', dof_labels{di}), ...
        'Position', [80 + 50*pi, 100, 700, 550]);

    [C_data, h_contour] = contourf(Tp_grid, Hs_grid, Z, 12);
    clabel(C_data, h_contour, 'FontSize', 8, 'Color', 'k');
    colormap(parula);
    cb = colorbar;
    cb.Label.String = sprintf('1%% Exceedance |F_{net}| [%s]', force_units{di});
    cb.Label.FontSize = 11;

    xlabel('T_p [s]', 'FontSize', 12);
    ylabel('H_s [m]', 'FontSize', 12);
    title(sprintf('%s: 1%% Exc. — Net Force on Body (incl. mooring)', ...
        dof_labels{di}), 'FontSize', 12);
    grid on;
    set(gca, 'Layer', 'top');

    fname = sprintf('sweep_1pct_net_%s.pdf', lower(dof_labels{di}));
    exportgraphics(fig_net(pi), fullfile(pdf_dir, fname), ...
        'ContentType', 'vector', 'Resolution', 300);
end

fprintf('Contour PDFs saved to:\n  %s\n', pdf_dir);

% --- Print summary tables ---
fprintf('\n===== 1%% EXCEEDANCE: NET HYDRO FORCE (excl. mooring) =====\n');
fprintf('%-8s  %-8s', 'Hs [m]', 'Tp [s]');
for pi = 1:n_active
    fprintf('  %12s', dof_labels{active_dofs(pi)});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 18 + 14*n_active));
for iH = 1:n_Hs
    for iT = 1:n_Tp
        fprintf('%-8.4f  %-8.3f', Hs_vec(iH), Tp_vec(iT));
        for pi = 1:n_active
            fprintf('  %12.4e', exc_1pct_hydro(iH, iT, active_dofs(pi)));
        end
        fprintf('\n');
    end
end

fprintf('\n===== 1%% EXCEEDANCE: NET FORCE ON BODY (incl. mooring) ====\n');
fprintf('%-8s  %-8s', 'Hs [m]', 'Tp [s]');
for pi = 1:n_active
    fprintf('  %12s', dof_labels{active_dofs(pi)});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 18 + 14*n_active));
for iH = 1:n_Hs
    for iT = 1:n_Tp
        fprintf('%-8.4f  %-8.3f', Hs_vec(iH), Tp_vec(iT));
        for pi = 1:n_active
            fprintf('  %12.4e', exc_1pct_net(iH, iT, active_dofs(pi)));
        end
        fprintf('\n');
    end
end
fprintf('==========================================================\n');

fprintf('\n========================================\n');
fprintf(' Parametric sweep complete.\n');
fprintf('========================================\n');

%% =======================================================================
%                        LOCAL FUNCTIONS
%  =======================================================================

function M = load_matrix(filepath)
%LOAD_MATRIX  Load a 6x6 matrix from Nemoh Mechanics file
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
        if strncmpi(line, 'VARIABLES', 9) || line(1) == '"'; continue; end
        if strncmpi(line, 'Zone', 4)
            if zone_count > 0
                raw_data{zone_count} = current_data; %#ok<AGROW>
            end
            zone_count = zone_count + 1;
            current_data = [];
            continue;
        end
        vals = sscanf(line, '%e');
        if ~isempty(vals)
            current_data = [current_data; vals']; %#ok<AGROW>
        end
    end
    if zone_count > 0
        raw_data{zone_count} = current_data;
    end
    fclose(fid);

    omega = raw_data{1}(:, 1);
    n_f = length(omega);
    A = zeros(n_dof, n_dof, n_f);
    B = zeros(n_dof, n_dof, n_f);

    for j = 1:min(zone_count, n_dof)
        data = raw_data{j};
        for k = 1:n_dof
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
            F_phase(i, :) = data(:, col_phase);
        end
    end
end


function [alpha, beta] = fit_prony(t, K_vals, n_terms)
%FIT_PRONY  Fit retardation function with real decaying exponentials
    t = t(:);
    K_vals = K_vals(:);
    t_max = t(end);
    t_min = max(t(2) - t(1), t_max / 1000);
    tau = logspace(log10(t_min), log10(t_max), n_terms);
    beta = -1 ./ tau(:);
    Phi = exp(t * beta');
    alpha = Phi \ K_vals;
    alpha = alpha(:);
    beta  = beta(:);
end


function S = jonswap_spectrum(omega, Hs, Tp, gamma)
%JONSWAP_SPECTRUM  JONSWAP (or Pierson-Moskowitz) wave spectrum
    omega = omega(:);
    wp = 2*pi / Tp;
    sigma = 0.09 * ones(size(omega));
    sigma(omega <= wp) = 0.07;
    r = exp(-(omega - wp).^2 ./ (2 * sigma.^2 * wp^2));
    alpha_norm = 1 - 0.287 * log(gamma);
    S_pm = (5/16) * Hs^2 * wp^4 ./ omega.^5 .* exp(-5/4 * (wp./omega).^4);
    S = alpha_norm * S_pm .* gamma.^r;
    S(omega < wp/10) = 0;
    S(~isfinite(S)) = 0;
end


function dxdt = cummins_rhs(t, x, p)
%CUMMINS_RHS  Right-hand side of Cummins equation with Prony states
    n_dof = 6;
    n_p   = p.n_prony;

    q = x(1:6);
    v = x(7:12);

    base_phase = p.wave_omega * t + p.wave_eps;
    cv = cos(base_phase);
    sv = sin(base_phase);
    F_exc = p.F_cos_pre * cv - p.F_sin_pre * sv;

    F_ext = p.F_external(t);

    I_3d = reshape(x(13:end), [n_dof, n_dof, n_p]);
    mu = sum(sum(I_3d, 3), 2);

    dv = p.M_inv * (F_exc + F_ext - p.K_total * q - mu);

    dv(p.disabled_dofs) = 0;
    v(p.disabled_dofs) = 0;

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
    r = point(:) - cog(:);
    if isa(force_input, 'function_handle')
        F_ext_fn = @(t) force_and_moment(force_input(t), r);
    else
        F6 = force_and_moment(force_input(:), r);
        F_ext_fn = @(t) F6;
    end
end


function F6 = force_and_moment(F3, r)
%FORCE_AND_MOMENT  Combine 3D force with moment about CoG
    F3 = F3(:);
    M  = cross(r, F3);
    F6 = [F3; M];
end
