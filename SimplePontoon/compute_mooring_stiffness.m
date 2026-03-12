function [K_moor, F_equil, info] = compute_mooring_stiffness(lines, body_pos, options)
%COMPUTE_MOORING_STIFFNESS  Linearized 6x6 mooring stiffness matrix
%
%  Computes the linearized stiffness matrix of a mooring system by
%  numerically differentiating the quasi-static mooring force about a
%  given body position.  Supports catenary and taut-wire line models.
%
%  Usage:
%    [K_moor, F_equil, info] = compute_mooring_stiffness(lines, body_pos)
%    [K_moor, F_equil, info] = compute_mooring_stiffness(lines, body_pos, opts)
%
%  Inputs:
%    lines : struct array (one element per mooring line) with fields:
%        .type     : 'catenary' or 'taut'
%        .anchor   : [3x1] anchor position in global frame [m]
%        .fairlead : [3x1] fairlead position on body, relative to CG [m]
%        .L0       : unstretched line length [m]
%        .EA       : axial stiffness [N]
%        .w        : submerged weight per unit length [N/m] (catenary only)
%                    w = (rho_line - rho_water) * g * A_cross
%                    (positive = line sinks)
%
%    body_pos : [6x1] body position at linearization point
%               [surge; sway; heave; roll; pitch; yaw] in [m; rad]
%
%    options  : (optional) struct with fields:
%        .delta_t : perturbation for translations [m]  (default 1e-5)
%        .delta_r : perturbation for rotations [rad]   (default 1e-5)
%        .plot    : true/false, plot mooring layout     (default true)
%
%  Outputs:
%    K_moor  : [6x6] linearized stiffness matrix
%              Plug directly into run_nemoh_timeseries.m as K_mooring
%    F_equil : [6x1] mooring force/moment at body_pos [N; N.m]
%    info    : struct with per-line results (tensions, angles, etc.)
%
%  Examples:
%    %% --- 3-line catenary spread at 120 deg ---
%    depth = 1.0;  R_anch = 2.0;  R_fair = 0.15;
%    for i = 1:3
%        ang = (i-1)*120 * pi/180;
%        lines(i).type     = 'catenary';
%        lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
%        lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); 0];
%        lines(i).L0 = 2.5;   lines(i).w = 0.2;   lines(i).EA = 5000;
%    end
%    [K, F0] = compute_mooring_stiffness(lines, zeros(6,1));
%
%    %% --- 4-line taut spread at 90 deg ---
%    for i = 1:4
%        ang = (i-1)*90 * pi/180;
%        lines(i).type     = 'taut';
%        lines(i).anchor   = [1.5*cos(ang); 1.5*sin(ang); -1.0];
%        lines(i).fairlead = [0.1*cos(ang); 0.1*sin(ang); -0.01];
%        lines(i).L0 = 1.3;   lines(i).EA = 2000;
%    end
%    [K, F0] = compute_mooring_stiffness(lines, zeros(6,1));
%
%  Reference:
%    Jonkman (2009) "Definition of the Floating System for Phase IV of OC3"
%    Faltinsen (1990) "Sea Loads on Ships and Offshore Structures"

    % --- Defaults ---
    if nargin < 2 || isempty(body_pos); body_pos = zeros(6,1); end
    if nargin < 3; options = struct(); end
    if ~isfield(options, 'delta_t'); options.delta_t = 1e-5; end
    if ~isfield(options, 'delta_r'); options.delta_r = 1e-5; end
    if ~isfield(options, 'plot');    options.plot = true; end

    body_pos = body_pos(:);
    n_lines  = length(lines);

    % ---- Force at equilibrium ----
    [F_equil, line_info] = total_mooring_force(lines, body_pos);

    % ---- Numerical linearization (central differences) ----
    delta = [options.delta_t * ones(3,1); options.delta_r * ones(3,1)];
    K_moor = zeros(6);

    for j = 1:6
        bp = body_pos;  bp(j) = bp(j) + delta(j);
        Fp = total_mooring_force(lines, bp);

        bm = body_pos;  bm(j) = bm(j) - delta(j);
        Fm = total_mooring_force(lines, bm);

        % K = -dF/dq  (restoring stiffness convention)
        K_moor(:,j) = -(Fp - Fm) / (2*delta(j));
    end

    % Symmetrize to remove numerical noise
    K_moor = 0.5 * (K_moor + K_moor');

    % ---- Collect output info ----
    info.lines   = line_info;
    info.F_equil = F_equil;

    % ---- Print summary ----
    dof_labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
    fprintf('\n============ MOORING STIFFNESS ============\n');
    fprintf('Lines: %d    Linearized at: [', n_lines);
    fprintf(' %.4f', body_pos); fprintf(' ]\n\n');
    for i = 1:n_lines
        fprintf('  Line %d (%s): T = %.4f N, angle = %.1f deg\n', ...
            i, lines(i).type, line_info(i).tension, line_info(i).angle_deg);
    end
    fprintf('\nEquilibrium force [N, N.m]:\n  [');
    fprintf(' %+.4e', F_equil); fprintf(' ]\n');
    fprintf('\nK_moor [N/m  N/rad  N.m/m  N.m/rad]:\n');
    fprintf('          ');
    for j = 1:6; fprintf('%10s', dof_labels{j}); end; fprintf('\n');
    for i = 1:6
        fprintf('  %-6s', dof_labels{i});
        for j = 1:6
            fprintf(' %+9.3f', K_moor(i,j));
        end
        fprintf('\n');
    end
    fprintf('===========================================\n');

    % ---- Plot ----
    if options.plot
        plot_mooring(lines, body_pos, line_info);
    end
end


%% ===================================================================
%                       INTERNAL FUNCTIONS
%  ===================================================================

function [F, li] = total_mooring_force(lines, bp)
%TOTAL_MOORING_FORCE  Sum forces/moments from all mooring lines

    F = zeros(6,1);
    n = length(lines);
    pos  = bp(1:3);
    R    = rot_zyx(bp(4), bp(5), bp(6));

    for i = 1:n
        fl_glob = pos + R * lines(i).fairlead(:);
        anch    = lines(i).anchor(:);

        switch lower(lines(i).type)
            case 'catenary'
                [Fl, T, ai] = catenary_force(anch, fl_glob, ...
                    lines(i).L0, lines(i).w, lines(i).EA);
            case 'taut'
                [Fl, T, ai] = taut_force(anch, fl_glob, ...
                    lines(i).L0, lines(i).EA);
            otherwise
                error('Unknown line type: %s', lines(i).type);
        end

        r_arm = R * lines(i).fairlead(:);
        M     = cross(r_arm, Fl);

        F(1:3) = F(1:3) + Fl;
        F(4:6) = F(4:6) + M;

        if nargout > 1
            li(i).F      = Fl;
            li(i).tension     = T;
            li(i).fairlead_g  = fl_glob;
            li(i).angle_deg   = ai.angle_deg;
        end
    end
end


function [Fl, T, ai] = taut_force(anch, fl, L0, EA)
%TAUT_FORCE  Force from a taut (pre-tensioned) mooring line

    d  = fl - anch;
    Lc = norm(d);

    if Lc < 1e-12
        Fl = zeros(3,1); T = 0;
        ai.angle_deg = 0; ai.stretch = 0;
        return
    end

    u = d / Lc;           % unit vector anchor → fairlead
    stretch = Lc - L0;

    if stretch > 0
        T = EA * stretch / L0;
    else
        T = 0;            % line slack
    end

    Fl = -T * u;          % pulls fairlead toward anchor

    ai.stretch   = stretch;
    ai.angle_deg = atan2d(abs(d(3)), sqrt(d(1)^2 + d(2)^2));
end


function [Fl, T, ai] = catenary_force(anch, fl, L0, w, EA)
%CATENARY_FORCE  Force from a catenary mooring line
%   2D catenary solved in the vertical plane containing anchor & fairlead.

    dx = fl(1) - anch(1);
    dy = fl(2) - anch(2);
    dz = fl(3) - anch(3);

    D = sqrt(dx^2 + dy^2);      % horizontal span
    H = dz;                      % vertical span (positive up)

    % Solve 2D catenary
    [FH, FV] = solve_catenary(D, H, L0, w, EA);

    T = sqrt(FH^2 + FV^2);

    % Horizontal unit vector: fairlead → anchor (pull direction)
    if D > 1e-12
        ux = -dx / D;
        uy = -dy / D;
    else
        ux = 0;  uy = 0;
    end

    % Force on fairlead: horizontal toward anchor, vertical downward
    Fl = [FH*ux;  FH*uy;  -FV];

    ai.FH = FH;  ai.FV = FV;
    ai.angle_deg = atan2d(FV, FH);
    ai.D = D;    ai.H = H;
end


function [FH, FV] = solve_catenary(D, H, L0, w, EA)
%SOLVE_CATENARY  2D extensible catenary equation solver
%
%  Given horizontal span D, vertical span H, line properties (L0, w, EA),
%  find horizontal and vertical tension at the fairlead (FH, FV).
%
%  Handles both seabed-contact and fully-suspended configurations.
%  Uses Newton-Raphson with analytic Jacobian.

    % --- Near-weightless → taut wire fallback ---
    if w < 1e-12
        Lg = sqrt(D^2 + H^2);
        T  = EA * max(0, Lg - L0) / max(L0, 1e-12);
        FH = T * D / max(Lg, 1e-12);
        FV = T * H / max(Lg, 1e-12);
        return
    end

    % --- Initial guess (taut-wire based) ---
    Lg = sqrt(D^2 + H^2);
    T0 = max(EA * max(0, Lg - L0) / L0,  w * L0 * 0.1);
    FH = max(T0 * D / max(Lg, 1e-10),  0.1 * w * L0);
    FV = max(w * L0 * 0.5,  T0 * abs(H) / max(Lg, 1e-10));

    % --- Newton-Raphson ---
    for iter = 1:200
        [r, J] = cat_residual(FH, FV, D, H, L0, w, EA);

        if norm(r) < 1e-12;  return;  end

        det_J = J(1,1)*J(2,2) - J(1,2)*J(2,1);
        if abs(det_J) < 1e-30;  break;  end          % singular

        dx = -[J(2,2)*r(1) - J(1,2)*r(2);
               -J(2,1)*r(1) + J(1,1)*r(2)] / det_J;

        % Back-tracking line search
        alpha = 1.0;
        for ls = 1:30
            FH_try = FH + alpha*dx(1);
            FV_try = FV + alpha*dx(2);
            if FH_try > 0 && FV_try > 0
                r_try = cat_residual(FH_try, FV_try, D, H, L0, w, EA);
                if norm(r_try) < norm(r)
                    FH = FH_try;
                    FV = FV_try;
                    break
                end
            end
            alpha = alpha * 0.5;
        end

        if alpha < 1e-10;  break;  end
    end

    if norm(r) > 1e-6
        warning('Catenary solver did not converge (|r|=%.2e). D=%.4f H=%.4f L0=%.4f', ...
            norm(r), D, H, L0);
    end
end


function [r, J] = cat_residual(FH, FV, D, H, L0, w, EA)
%CAT_RESIDUAL  Residual and analytic Jacobian for catenary equations

    ratio = FV / FH;
    s1    = sqrt(1 + ratio^2);           % sqrt(1 + (FV/FH)^2)
    a1    = asinh(ratio);                % asinh(FV/FH)

    if FV < w * L0
        % ---- Seabed contact ----
        Ls = FV / w;

        D_calc = (L0 - Ls) + FH/w * a1 + FH*L0/EA;
        H_calc = FH/w * (s1 - 1) + FV^2 / (2*w*EA);

        r = [D_calc - D;  H_calc - H];

        % Jacobian
        dDdFH = a1/w  - FV/(w*FH*s1)  + L0/EA;
        dDdFV = -1/w  + 1/(w*s1);
        dHdFH = (s1 - 1)/w  - FV^2/(w*FH^2*s1);
        dHdFV = FV/(w*FH*s1) + FV/(w*EA);
    else
        % ---- Fully suspended ----
        ratio2 = (FV - w*L0) / FH;
        s2     = sqrt(1 + ratio2^2);
        a2     = asinh(ratio2);

        D_calc = FH/w * (a1 - a2) + FH*L0/EA;
        H_calc = FH/w * (s1 - s2) + (FV*L0 - w*L0^2/2)/EA;

        r = [D_calc - D;  H_calc - H];

        dDdFH = (a1-a2)/w - FV/(w*FH*s1) + (FV-w*L0)/(w*FH*s2) + L0/EA;
        dDdFV = 1/(w*s1) - 1/(w*s2);
        dHdFH = (s1-s2)/w - FV^2/(w*FH^2*s1) + (FV-w*L0)^2/(w*FH^2*s2);
        dHdFV = FV/(w*FH*s1) - (FV-w*L0)/(w*FH*s2) + L0/EA;
    end

    J = [dDdFH, dDdFV;  dHdFH, dHdFV];
end


function R = rot_zyx(roll, pitch, yaw)
%ROT_ZYX  ZYX Euler rotation matrix (body to global)

    cr = cos(roll);   sr = sin(roll);
    cp = cos(pitch);  sp = sin(pitch);
    cy = cos(yaw);    sy = sin(yaw);

    R = [ cy*cp,   cy*sp*sr - sy*cr,   cy*sp*cr + sy*sr;
          sy*cp,   sy*sp*sr + cy*cr,   sy*sp*cr - cy*sr;
          -sp,     cp*sr,              cp*cr            ];
end


function plot_mooring(lines, bp, li)
%PLOT_MOORING  Visualize mooring layout (plan view + side view)

    n = length(lines);
    pos = bp(1:3);
    R   = rot_zyx(bp(4), bp(5), bp(6));

    figure('Name','Mooring Layout','Position',[100 200 1100 450]);
    colors = lines2colors(n);

    % --- Plan view ---
    subplot(1,2,1); hold on; axis equal; grid on;
    title('Plan View (top)');
    xlabel('X [m]'); ylabel('Y [m]');

    % Body marker
    plot(pos(1), pos(2), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.8]);

    for i = 1:n
        fl = pos + R * lines(i).fairlead(:);
        an = lines(i).anchor(:);
        plot([an(1) fl(1)], [an(2) fl(2)], '-', 'Color', colors(i,:), 'LineWidth', 2);
        plot(an(1), an(2), 'v', 'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
        plot(fl(1), fl(2), 'o', 'Color', colors(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        text(an(1), an(2), sprintf('  A%d',i), 'Color', colors(i,:), 'FontWeight','bold');
    end
    legend_entries = arrayfun(@(i) sprintf('Line %d (T=%.2fN)', i, li(i).tension), 1:n, 'UniformOutput', false);
    % Don't call legend with too many handles; just annotate

    % --- Side view ---
    subplot(1,2,2); hold on; grid on;
    title('Side View (elevation)');
    xlabel('Horizontal distance [m]'); ylabel('Z [m]');

    plot(0, pos(3), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.8]);

    for i = 1:n
        fl = pos + R * lines(i).fairlead(:);
        an = lines(i).anchor(:);
        D_i = sqrt((fl(1)-an(1))^2 + (fl(2)-an(2))^2);

        if strcmpi(lines(i).type, 'catenary') && isfield(li(i), 'F')
            % Draw catenary shape
            draw_catenary_shape(D_i, fl(3)-an(3), an(3), lines(i), li(i), colors(i,:));
        else
            % Straight line
            plot([0 D_i], [an(3) fl(3)], '-', 'Color', colors(i,:), 'LineWidth', 2);
        end
        plot(0, an(3), 'v', 'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
        plot(D_i, fl(3), 'o', 'Color', colors(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        text(D_i, fl(3), sprintf('  L%d: %.2fN',i,li(i).tension), 'Color', colors(i,:));
    end

    % Draw seabed
    xl = xlim;
    plot(xl, [min(arrayfun(@(l) l.anchor(3), lines)) min(arrayfun(@(l) l.anchor(3), lines))], ...
        'k--', 'LineWidth', 1);

    sgtitle('Mooring Layout', 'FontWeight', 'bold');
end


function draw_catenary_shape(D, H, z_anchor, line, li, col)
%DRAW_CATENARY_SHAPE  Plot the catenary curve in the side view

    FH = li.F(1)^2 + li.F(2)^2;
    FH = sqrt(FH);  % horizontal tension magnitude
    FV_fair = abs(li.F(3));
    w  = line.w;
    L0 = line.L0;
    EA = line.EA;

    if w < 1e-12
        plot([0 D], [z_anchor, z_anchor + H], '-', 'Color', col, 'LineWidth', 2);
        return
    end

    % Determine if seabed contact
    if FV_fair < w * L0
        Ls = FV_fair / w;
        Lb = L0 - Ls;
    else
        Ls = L0;
        Lb = 0;
    end

    n_pts = 100;

    if Lb > 0
        % Seabed portion
        x_bed = linspace(0, Lb, 20);
        z_bed = z_anchor * ones(size(x_bed));

        % Suspended portion (s from 0 = touchdown to Ls = fairlead)
        s = linspace(0, Ls, n_pts);
        x_cat = Lb + FH/w * asinh(w*s / FH);
        z_cat = z_anchor + FH/w * (sqrt(1 + (w*s/FH).^2) - 1);

        % Add elastic stretch (approximate shift)
        x_stretch = FH * s / EA;
        x_cat = x_cat + x_stretch;

        plot([x_bed, x_cat], [z_bed, z_cat], '-', 'Color', col, 'LineWidth', 2);
    else
        % Fully suspended
        s = linspace(0, L0, n_pts);
        FV_anch = FV_fair - w * L0;   % vertical tension at anchor

        x_cat = FH/w * (asinh((FV_anch + w*s)/FH) - asinh(FV_anch/FH));
        z_cat = z_anchor + FH/w * (sqrt(1+((FV_anch+w*s)/FH).^2) - sqrt(1+(FV_anch/FH)^2));

        plot(x_cat, z_cat, '-', 'Color', col, 'LineWidth', 2);
    end
end


function c = lines2colors(n)
%LINES2COLORS  Generate distinct colors for n mooring lines

    base = [0.8 0.2 0.2;   % red
            0.2 0.6 0.2;   % green
            0.2 0.2 0.8;   % blue
            0.8 0.6 0.0;   % orange
            0.6 0.2 0.8;   % purple
            0.0 0.7 0.7;   % cyan
            0.8 0.4 0.6;   % pink
            0.4 0.4 0.4];  % gray

    c = base(mod((1:n)-1, size(base,1)) + 1, :);
end
