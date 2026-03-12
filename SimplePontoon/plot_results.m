%% plot_results.m
%  Visualise Nemoh BEM results for the SimplePontoon case.
%  Reads all output files from results/ and creates publication-quality plots.
%
%  Outputs plotted:
%    1) Added mass & radiation damping (diagonal terms)
%    2) Excitation force magnitude & phase
%    3) Froude-Krylov vs Diffraction vs Total excitation
%    4) Impulse response functions (radiation)
%    5) Excitation force IRFs
%    6) Full 6x6 added mass & damping matrices (from CM.dat / CA.dat)

clear; close all; clc;

%% ---- User settings --------------------------------------------------------
resDir  = fullfile(fileparts(mfilename('fullpath')), 'results');
dofNames = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
dofUnits_force  = {'N/m','N/m','N/m','N\cdotm/m','N\cdotm/m','N\cdotm/m'};
dofUnits_A = {'kg','kg','kg','kg\cdotm','kg\cdotm','kg\cdotm^2'};
dofUnits_B = {'kg/s','kg/s','kg/s','kg\cdotm/s','kg\cdotm/s','kg\cdotm^2/s'};
nDof = 6;

%% ========================================================================
%  0) MESH VISUALISATION
%  ========================================================================
caseDir = fileparts(mfilename('fullpath'));
meshFile = fullfile(caseDir, 'SimplePontoon.dat');
[meshNodes, meshPanels] = read_nemoh_mesh(meshFile);

% Classify panels: lid (all z==0) vs wetted (at least one z<0)
isLid = false(size(meshPanels,1),1);
for ip = 1:size(meshPanels,1)
    zv = meshNodes(meshPanels(ip,:), 3);
    if all(abs(zv) < 1e-8)
        isLid(ip) = true;
    end
end
nLid = sum(isLid);
nWet = sum(~isLid);

fprintf('Mesh: %d nodes, %d panels (%d wetted, %d lid)\n', ...
    size(meshNodes,1), size(meshPanels,1), nWet, nLid);

% --- Figure: 3D mesh with colour-coded panels ---
figure('Name','Mesh Visualisation','Position',[50 50 1400 600]);

% (a) Full 3D view
subplot(1,2,1);
% Wetted panels in blue
if nWet > 0
    patch('Vertices', meshNodes, ...
          'Faces', meshPanels(~isLid,:), ...
          'FaceColor', [0.3 0.6 1.0], ...
          'EdgeColor', 'k', ...
          'FaceAlpha', 0.7, ...
          'LineWidth', 0.3);
end
hold on;
% Lid panels in red
if nLid > 0
    patch('Vertices', meshNodes, ...
          'Faces', meshPanels(isLid,:), ...
          'FaceColor', [1.0 0.3 0.3], ...
          'EdgeColor', 'k', ...
          'FaceAlpha', 0.5, ...
          'LineWidth', 0.3);
end
% Waterline
plot3(meshNodes(:,1), meshNodes(:,2), zeros(size(meshNodes,1),1), ...
      '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 1);
axis equal; grid on; view(35, 25);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('3D Mesh (%d wetted + %d lid panels)', nWet, nLid));
legend({'Wetted surface','Lid (IRR)'},'Location','best');
colormap(gca, 'parula');

% (b) Cross-section view (y-z plane at x~0)
subplot(1,2,2);
% Find panels near x=0
xc_panels = zeros(size(meshPanels,1),1);
for ip = 1:size(meshPanels,1)
    xc_panels(ip) = mean(meshNodes(meshPanels(ip,:), 1));
end
% Plot all panels coloured by z-coordinate of centroid
zc_panels = zeros(size(meshPanels,1),1);
for ip = 1:size(meshPanels,1)
    zc_panels(ip) = mean(meshNodes(meshPanels(ip,:), 3));
end
patch('Vertices', meshNodes, ...
      'Faces', meshPanels(~isLid,:), ...
      'FaceVertexCData', zc_panels(~isLid), ...
      'FaceColor', 'flat', ...
      'EdgeColor', [0.3 0.3 0.3], ...
      'FaceAlpha', 0.85, ...
      'LineWidth', 0.3);
hold on;
if nLid > 0
    patch('Vertices', meshNodes, ...
          'Faces', meshPanels(isLid,:), ...
          'FaceColor', [1.0 0.3 0.3], ...
          'EdgeColor', [0.3 0.3 0.3], ...
          'FaceAlpha', 0.5, ...
          'LineWidth', 0.3);
end
axis equal; grid on; view(0, 0);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Side View (y-z)');
cb = colorbar; ylabel(cb, 'Panel centroid z [m]');

sgtitle('SimplePontoon BEM Mesh');

% --- Figure: Multiple views ---
figure('Name','Mesh Views','Position',[100 100 1400 900]);

viewAngles = {[0 90], [0 0], [90 0], [35 25]};
viewNames  = {'Top (x-y)', 'Front (x-z)', 'Side (y-z)', 'Isometric'};

for iv = 1:4
    subplot(2,2,iv);
    if nWet > 0
        patch('Vertices', meshNodes, ...
              'Faces', meshPanels(~isLid,:), ...
              'FaceColor', [0.3 0.6 1.0], ...
              'EdgeColor', 'k', ...
              'FaceAlpha', 0.7, ...
              'LineWidth', 0.3);
    end
    hold on;
    if nLid > 0
        patch('Vertices', meshNodes, ...
              'Faces', meshPanels(isLid,:), ...
              'FaceColor', [1.0 0.3 0.3], ...
              'EdgeColor', 'k', ...
              'FaceAlpha', 0.5, ...
              'LineWidth', 0.3);
    end
    axis equal; grid on;
    view(viewAngles{iv});
    xlabel('x'); ylabel('y'); zlabel('z');
    title(viewNames{iv});
end
sgtitle('SimplePontoon Mesh - Multiple Views');

%% ========================================================================
%  1) RADIATION COEFFICIENTS  (RadiationCoefficients.tec)
%  ========================================================================
fprintf('Reading RadiationCoefficients.tec ...\n');
fid = fopen(fullfile(resDir,'RadiationCoefficients.tec'),'r');
% Skip 7 header lines (VARIABLES + column names)
for k = 1:7, fgetl(fid); end

A = zeros(0, nDof, nDof);  % (nFreq, i, j)
B = zeros(0, nDof, nDof);
omega_rad = [];

for iDof = 1:nDof
    zoneLine = fgetl(fid);  % Zone header
    % Parse I=  value for number of frequencies
    tok = regexp(zoneLine, 'I=\s*(\d+)', 'tokens');
    nFreq = str2double(tok{1}{1});
    block = zeros(nFreq, 1 + 2*nDof);
    for row = 1:nFreq
        block(row,:) = sscanf(fgetl(fid), '%e')';
    end
    if iDof == 1
        omega_rad = block(:,1);
        A = zeros(nFreq, nDof, nDof);
        B = zeros(nFreq, nDof, nDof);
    end
    for jDof = 1:nDof
        A(:, iDof, jDof) = block(:, 2*(jDof-1)+2);
        B(:, iDof, jDof) = block(:, 2*(jDof-1)+3);
    end
end
fclose(fid);
nFreq = length(omega_rad);
fprintf('  %d frequencies, %d DoFs\n', nFreq, nDof);

%% ---- Plot diagonal added mass & damping -----------------------------------
figure('Name','Added Mass & Radiation Damping (diagonal)','Position',[50 50 1200 800]);
for i = 1:nDof
    subplot(nDof, 2, 2*(i-1)+1);
    plot(omega_rad, A(:,i,i), 'b-', 'LineWidth', 1.2);
    ylabel(sprintf('A_{%d%d} [%s]', i, i, dofUnits_A{i}));
    if i == 1, title('Added Mass'); end
    if i == nDof, xlabel('\omega [rad/s]'); end
    grid on;
    text(0.02, 0.92, dofNames{i}, 'Units','normalized', ...
         'FontWeight','bold', 'FontSize',9);

    subplot(nDof, 2, 2*(i-1)+2);
    plot(omega_rad, B(:,i,i), 'r-', 'LineWidth', 1.2);
    ylabel(sprintf('B_{%d%d} [%s]', i, i, dofUnits_B{i}));
    if i == 1, title('Radiation Damping'); end
    if i == nDof, xlabel('\omega [rad/s]'); end
    grid on;
end
sgtitle('Hydrodynamic Radiation Coefficients (diagonal)');

%% ========================================================================
%  2) EXCITATION FORCE  (ExcitationForce.tec)
%  ========================================================================
[omega_exc, Fexc_abs, Fexc_ang] = read_force_tec(...
    fullfile(resDir,'ExcitationForce.tec'), nDof);
fprintf('Reading ExcitationForce.tec ... %d freqs\n', length(omega_exc));

figure('Name','Excitation Force','Position',[100 100 1200 800]);
for i = 1:nDof
    subplot(nDof, 2, 2*(i-1)+1);
    plot(omega_exc, Fexc_abs(:,i), 'b-', 'LineWidth', 1.2);
    ylabel(sprintf('|F_{exc,%d}| [%s]', i, dofUnits_force{i}));
    if i == 1, title('Magnitude'); end
    if i == nDof, xlabel('\omega [rad/s]'); end
    grid on;
    text(0.02, 0.92, dofNames{i}, 'Units','normalized', ...
         'FontWeight','bold', 'FontSize',9);

    subplot(nDof, 2, 2*(i-1)+2);
    plot(omega_exc, rad2deg(Fexc_ang(:,i)), 'r-', 'LineWidth', 1.2);
    ylabel(sprintf('\\angle F_{exc,%d} [deg]', i));
    if i == 1, title('Phase'); end
    if i == nDof, xlabel('\omega [rad/s]'); end
    grid on; ylim([-180 180]);
end
sgtitle('Wave Excitation Force / Moment');

%% ========================================================================
%  3) FK vs DIFFRACTION vs TOTAL
%  ========================================================================
[omega_fk,  Ffk_abs,  Ffk_ang]  = read_force_tec(...
    fullfile(resDir,'FKForce.tec'), nDof);
[omega_dif, Fdif_abs, Fdif_ang] = read_force_tec(...
    fullfile(resDir,'DiffractionForce.tec'), nDof);

figure('Name','FK vs Diffraction vs Total','Position',[150 150 1200 800]);
plotDofs = [1, 3, 5];  % Surge, Heave, Pitch (most relevant for head seas)
for idx = 1:length(plotDofs)
    i = plotDofs(idx);

    subplot(length(plotDofs), 2, 2*(idx-1)+1);
    plot(omega_fk,  Ffk_abs(:,i),  'g-',  'LineWidth', 1.2); hold on;
    plot(omega_dif, Fdif_abs(:,i), 'm--', 'LineWidth', 1.2);
    plot(omega_exc, Fexc_abs(:,i), 'b-',  'LineWidth', 1.5);
    ylabel(sprintf('|F_%d| [%s]', i, dofUnits_force{i}));
    legend('Froude-Krylov','Diffraction','Total','Location','best');
    if idx == 1, title('Force Magnitude'); end
    if idx == length(plotDofs), xlabel('\omega [rad/s]'); end
    grid on;
    text(0.02, 0.92, dofNames{i}, 'Units','normalized', ...
         'FontWeight','bold', 'FontSize',9);

    subplot(length(plotDofs), 2, 2*(idx-1)+2);
    plot(omega_fk,  rad2deg(Ffk_ang(:,i)),  'g-',  'LineWidth', 1.2); hold on;
    plot(omega_dif, rad2deg(Fdif_ang(:,i)), 'm--', 'LineWidth', 1.2);
    plot(omega_exc, rad2deg(Fexc_ang(:,i)), 'b-',  'LineWidth', 1.5);
    ylabel(sprintf('\\angle F_%d [deg]', i));
    legend('FK','Diffr','Total','Location','best');
    if idx == 1, title('Phase'); end
    if idx == length(plotDofs), xlabel('\omega [rad/s]'); end
    grid on; ylim([-180 180]);
end
sgtitle('Froude-Krylov / Diffraction / Total Excitation (head seas)');

%% ========================================================================
%  4) RADIATION IRF  (IRF.tec)
%  ========================================================================
fprintf('Reading IRF.tec ...\n');
fid = fopen(fullfile(resDir,'IRF.tec'),'r');
for k = 1:7, fgetl(fid); end  % skip headers

Ainf = zeros(0, nDof, nDof);
Kirf = zeros(0, nDof, nDof);
t_irf = [];

for iDof = 1:nDof
    zoneLine = fgetl(fid);
    tok = regexp(zoneLine, 'I=\s*(\d+)', 'tokens');
    nT = str2double(tok{1}{1});
    block = zeros(nT, 1 + 2*nDof);
    for row = 1:nT
        block(row,:) = sscanf(fgetl(fid), '%e')';
    end
    if iDof == 1
        t_irf = block(:,1);
        Ainf = zeros(nT, nDof, nDof);
        Kirf = zeros(nT, nDof, nDof);
    end
    for jDof = 1:nDof
        Ainf(:, iDof, jDof) = block(:, 2*(jDof-1)+2);
        Kirf(:, iDof, jDof) = block(:, 2*(jDof-1)+3);
    end
end
fclose(fid);
fprintf('  %d time steps\n', length(t_irf));

figure('Name','Radiation IRF (diagonal)','Position',[200 200 1200 700]);
for i = 1:nDof
    subplot(2, 3, i);
    yyaxis left
    plot(t_irf, Kirf(:,i,i), 'b-', 'LineWidth', 1.2);
    ylabel(sprintf('K_{%d%d}(t)', i, i));
    yyaxis right
    plot(t_irf, Ainf(:,i,i), 'r--', 'LineWidth', 1.0);
    ylabel(sprintf('A_{%d%d}^\\infty', i, i));
    xlabel('t [s]'); title(dofNames{i}); grid on;
end
sgtitle('Radiation Impulse Response Functions');

%% ========================================================================
%  5) EXCITATION FORCE IRF  (IRF_excForce.tec)
%  ========================================================================
fprintf('Reading IRF_excForce.tec ...\n');
fid = fopen(fullfile(resDir,'IRF_excForce.tec'),'r');
for k = 1:7, fgetl(fid); end  % skip headers

zoneLine = fgetl(fid);
tok = regexp(zoneLine, 'I=\s*(\d+)', 'tokens');
nT_exc = str2double(tok{1}{1});
block = zeros(nT_exc, 1 + nDof);
for row = 1:nT_exc
    block(row,:) = sscanf(fgetl(fid), '%e')';
end
fclose(fid);

t_exc_irf = block(:,1);
Kexc_irf  = block(:, 2:end);

figure('Name','Excitation Force IRF','Position',[250 250 1200 500]);
for i = 1:nDof
    subplot(2, 3, i);
    plot(t_exc_irf, Kexc_irf(:,i), 'b-', 'LineWidth', 1.0);
    xlabel('t [s]'); ylabel(sprintf('K_{exc,%d}(t)', i));
    title(dofNames{i}); grid on;
end
sgtitle('Excitation Force Impulse Response Functions');

%% ========================================================================
%  6) FULL 6x6 ADDED MASS & DAMPING  (CM.dat, CA.dat)
%  ========================================================================
[omega_cm, CM_full] = read_matrix_dat(fullfile(resDir,'CM.dat'), nDof);
[omega_ca, CA_full] = read_matrix_dat(fullfile(resDir,'CA.dat'), nDof);
fprintf('Reading CM.dat & CA.dat ... %d freqs\n', length(omega_cm));

% Plot full 6x6 added mass matrix
figure('Name','Full Added Mass Matrix A(w)','Position',[300 50 1400 900]);
for i = 1:nDof
    for j = 1:nDof
        subplot(nDof, nDof, (i-1)*nDof + j);
        plot(omega_cm, squeeze(CM_full(:,i,j)), 'b-', 'LineWidth', 0.8);
        if i == 1, title(dofNames{j}); end
        if j == 1, ylabel(dofNames{i}); end
        if i == nDof, xlabel('\omega'); end
        set(gca, 'FontSize', 7); grid on;
    end
end
sgtitle('Added Mass A_{ij}(\omega)');

% Plot full 6x6 damping matrix
figure('Name','Full Damping Matrix B(w)','Position',[350 50 1400 900]);
for i = 1:nDof
    for j = 1:nDof
        subplot(nDof, nDof, (i-1)*nDof + j);
        plot(omega_ca, squeeze(CA_full(:,i,j)), 'r-', 'LineWidth', 0.8);
        if i == 1, title(dofNames{j}); end
        if j == 1, ylabel(dofNames{i}); end
        if i == nDof, xlabel('\omega'); end
        set(gca, 'FontSize', 7); grid on;
    end
end
sgtitle('Radiation Damping B_{ij}(\omega)');

%% ========================================================================
%  7) SUMMARY FIGURE - KEY RESULTS
%  ========================================================================
figure('Name','Summary - Key Hydrodynamic Coefficients','Position',[50 50 1400 900]);

% Heave added mass
subplot(3,3,1);
plot(omega_rad, A(:,3,3), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('A_{33} [kg]');
title('Heave Added Mass'); grid on;

% Heave damping
subplot(3,3,2);
plot(omega_rad, B(:,3,3), 'r-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('B_{33} [kg/s]');
title('Heave Damping'); grid on;

% Heave excitation
subplot(3,3,3);
plot(omega_exc, Fexc_abs(:,3), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('|F_{exc,3}| [N/m]');
title('Heave Excitation'); grid on;

% Roll added mass
subplot(3,3,4);
plot(omega_rad, A(:,4,4), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('A_{44} [kg\cdotm^2]');
title('Roll Added Mass'); grid on;

% Roll damping
subplot(3,3,5);
plot(omega_rad, B(:,4,4), 'r-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('B_{44} [kg\cdotm^2/s]');
title('Roll Damping'); grid on;

% Surge excitation
subplot(3,3,6);
plot(omega_exc, Fexc_abs(:,1), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('|F_{exc,1}| [N/m]');
title('Surge Excitation'); grid on;

% Pitch added mass
subplot(3,3,7);
plot(omega_rad, A(:,5,5), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('A_{55} [kg\cdotm^2]');
title('Pitch Added Mass'); grid on;

% Pitch damping
subplot(3,3,8);
plot(omega_rad, B(:,5,5), 'r-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('B_{55} [kg\cdotm^2/s]');
title('Pitch Damping'); grid on;

% Pitch excitation
subplot(3,3,9);
plot(omega_exc, Fexc_abs(:,5), 'b-', 'LineWidth', 1.5);
xlabel('\omega [rad/s]'); ylabel('|F_{exc,5}| [N\cdotm/m]');
title('Pitch Excitation'); grid on;

sgtitle('SimplePontoon - Key Hydrodynamic Results');

fprintf('\nAll plots generated.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [omega, F_abs, F_ang] = read_force_tec(filepath, nDof)
%READ_FORCE_TEC  Read Nemoh force .tec file (Excitation/Diffraction/FK).
%  Returns omega (rad/s), magnitude array, phase array (radians).
    fid = fopen(filepath, 'r');
    for k = 1:7, fgetl(fid); end   % skip header lines
    zoneLine = fgetl(fid);          % zone header
    tok = regexp(zoneLine, 'I=\s*(\d+)', 'tokens');
    nF = str2double(tok{1}{1});
    block = zeros(nF, 1 + 2*nDof);
    for row = 1:nF
        block(row,:) = sscanf(fgetl(fid), '%e')';
    end
    fclose(fid);
    omega = block(:,1);
    F_abs = zeros(nF, nDof);
    F_ang = zeros(nF, nDof);
    for i = 1:nDof
        F_abs(:,i) = block(:, 2*(i-1)+2);
        F_ang(:,i) = block(:, 2*(i-1)+3);
    end
end

function [omega, M] = read_matrix_dat(filepath, nDof)
%READ_MATRIX_DAT  Read Nemoh CM.dat or CA.dat (frequency-dependent 6x6 matrix).
    fid = fopen(filepath, 'r');
    hdr = fgetl(fid);
    tok = regexp(hdr, '(\d+)', 'tokens');
    nF = str2double(tok{1}{1});
    omega = zeros(nF, 1);
    M = zeros(nF, nDof, nDof);
    for k = 1:nF
        omega(k) = sscanf(fgetl(fid), '%e');
        for i = 1:nDof
            vals = sscanf(fgetl(fid), '%e');
            M(k, i, :) = vals(:)';
        end
    end
    fclose(fid);
end

function [nodes, panels] = read_nemoh_mesh(filepath)
%READ_NEMOH_MESH  Read Nemoh .dat mesh file.
%  Returns:
%    nodes:  Nx3 array of (x,y,z) coordinates
%    panels: Mx4 array of node indices (1-based)
    fid = fopen(filepath, 'r');
    % Header line: magic number and symmetry flag
    fgetl(fid);

    % Read nodes until terminator (index == 0)
    nodeList = [];
    while true
        line = strtrim(fgetl(fid));
        vals = sscanf(line, '%f');
        if isempty(vals), continue; end
        idx = vals(1);
        if idx == 0, break; end
        nodeList(end+1, :) = vals(2:4)'; %#ok<AGROW>
    end

    % Read panels until terminator (all zeros)
    panelList = [];
    while true
        line = strtrim(fgetl(fid));
        if feof(fid), break; end
        vals = sscanf(line, '%d');
        if isempty(vals), continue; end
        if all(vals == 0), break; end
        panelList(end+1, :) = vals(1:4)'; %#ok<AGROW>
    end
    fclose(fid);

    nodes  = nodeList;
    panels = panelList;
end
