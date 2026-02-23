% Modify the file names below: 1st: gri mesh file, 2nd: U state
viz('mesh_refined_2394.gri','coarse_mesh_unsteady_firstorder_RK3_roe_t_171.000000.txt')

% viz.m - Main script equivalent to viz.py
% Usage: run from command line with arguments or call functions directly
% Command line equivalent: viz(meshFile, uFile, p, field, fname)

function viz(varargin)
    if nargin < 1
        disp('Pass at least one argument: Mesh filename');
        return;
    end
    meshFile = varargin{1};
    Mesh     = readgri(meshFile);
    U     = [];
    p     = 0;
    field = 'Mach';
    fname = '';

    if nargin >= 2, U     = readU(varargin{2});    end
    if nargin >= 3, p     = str2double(varargin{3}); end
    if nargin >= 4, field = varargin{4};             end
    if nargin >= 5, fname = varargin{5};             end

    if ~isempty(U)
        frange = [];
        if strcmpi(field, 'mach')
            frange = [0.0, 0.9];
        end
        plotstate(Mesh, U, p, field, frange, fname);
        plotCoefficients(Mesh, U);
    else
        plotmesh(Mesh, 'mesh.png');
    end
end


%-----------------------------------------------------------
function U = readU(fname)
    U = load(fname);
end


%-----------------------------------------------------------
function [cpVal, cx, cy] = cp(P,x1,y1,x2,y2,normalX,normalY)
    % Calculate the pressure coefficient and x/y force coefficients
    g    = 1.4;
    gmi  = g - 1;
    P0   = 1/g;
    Pout = P0 * 0.7;
    c=18.804;
    l = sqrt((x2-x1).^2+(y2-y1).^2);
    Mout2 = (2/gmi) * ((P0/Pout)^(gmi/g) - 1);
    qout  = 0.5 * g * Pout * Mout2;
    cpVal = (P - Pout) / qout;
    cx = -P.*l.*normalX./(qout.*c);
    cy = -P.*l.*normalY./(qout.*c);
end

function [nx, ny] = edgeNormals(x1, y1, x2, y2, direction)
    % Step 1: Compute edge tangent vectors
    dx = x2 - x1;
    dy = y2 - y1;
    % Step 2: Rotate tangent 90 degrees to get normal
    % Two choices: (-dy, dx) or (dy, -dx)
    nx = -dy;
    ny =  dx;
    % Step 3: Normalize to unit length
    len = sqrt(nx.^2 + ny.^2);
    nx  = nx ./ len;
    ny  = ny ./ len;
    % Step 4: Ensure normals point "above" the curve (positive y direction)
    % If the normal points downward (ny < 0), flip it
    if direction == 0 % want normals pointing above curve
        flip    = ny < 0;
        % colorPlot = 'red';
    elseif direction == 1 % normals pointing below curve
        flip = ny > 0;
        % colorPlot = 'blue';
    end
    nx(flip) = -nx(flip);
    ny(flip) = -ny(flip);
    % figure(45); hold on;
    % plot(nx,ny,'Color',colorPlot)
end

%-----------------------------------------------------------
function [Pup, Plow, xup, xlow] = getEdgeP(Mesh, U)
    V     = Mesh.V;
    B_up  = Mesh.B_up;
    B_low = Mesh.B_low;

    U_mod = getField(U, 'mach');

    Pup  = U_mod(B_up(:,3),  6);
    Plow = U_mod(B_low(:,3), 6);

    % Average x coordinate of the two edge nodes
    xup  = (V(B_up(:,1),  1) + V(B_up(:,2),  1)) / 2;
    xlow = (V(B_low(:,1), 1) + V(B_low(:,2), 1)) / 2;
end


%-----------------------------------------------------------
function COEFF = plotCoefficients(Mesh, U)
    V     = Mesh.V;
    B_up  = Mesh.B_up;
    B_low = Mesh.B_low;

    U_mod = getField(U, 'mach');

    Pup  = U_mod(B_up(:,3),  6);
    Plow = U_mod(B_low(:,3), 6);

    % Average x coordinate of the two edge nodes
    xup  = (V(B_up(:,1),  1) + V(B_up(:,2),  1)) / 2;
    xlow = (V(B_low(:,1), 1) + V(B_low(:,2), 1)) / 2;
    % [Pup, Plow, xup, xlow] = getEdgeP(Mesh, U);
    [nx_up, ny_up] = edgeNormals(V(B_up(:,1),  1),V(B_up(:,1),  2),V(B_up(:,2),  1),V(B_up(:,2),  2), 0);
    [nx_low, ny_low] = edgeNormals(V(B_low(:,1),  1),V(B_low(:,1),  2),V(B_low(:,2),  1),V(B_low(:,2),  2), 1);

    [cp_up, cx_up, cy_up]  = cp(Pup,V(B_up(:,1),  1),V(B_up(:,1),  2),V(B_up(:,2),  1),V(B_up(:,2),  2),nx_up, ny_up);
    [cp_low, cx_low, cy_low] = cp(Plow,V(B_low(:,1),  1),V(B_low(:,1),  2),V(B_low(:,2),  1),V(B_low(:,2),  2),nx_low, ny_low);
    CX_low = sum(cx_low);
    CY_low = sum(cy_low);
    CX_up = sum(cx_up);
    CY_up = sum(cy_up);
    CX = CX_up+CX_low % outputs the force coefficient in x
    CY = CY_up+CY_low % outputs the force coefficient in y
    COEFF = [CX_low, CY_low; CX_up, CY_up];

    fig = figure(1);%'Position', [100, 100, 800, 600]);
    hold on;
    plot(xup,  cp_up,  'm-', 'LineWidth', 1.5, 'DisplayName', 'Upper Surface');
    plot(xlow, cp_low, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Lower Surface');
    hold off;

    title('Pressure Coefficient Distribution', 'FontSize', 14);
    xlabel('$x$ [mm]',  'FontSize', 12);
    ylabel('$C_p$',  'FontSize', 12);
    legend('FontSize', 12);
    set(gca, 'YDir', 'reverse');  % invert y-axis (aerodynamic convention)
    grid on;

    saveas(fig, 'pressureCoeff.png');
    set(fig, 'visible', 'on');
end


%-----------------------------------------------------------
function U_mod = getField(U, field)
    r  = U(:,1); ru = U(:,2); rv = U(:,3); rE = U(:,4);
    g  = 1.4;
    V  = sqrt(ru.^2 + rv.^2) ./ r;
    p  = (g-1.) .* (rE - 0.5.*r.*V.^2);
    c  = sqrt(g.*p./r);
    M  = V ./ c;
    S  = p ./ r.^g;

    U_mod = [U, V, p, c, M, S];
    % Columns: r(1), ru(2), rv(3), rE(4), V(5), p(6), c(7), M(8), S(9)

    if ~strcmpi(field, 'mach')
        U_mod = [];
    end
end


%-----------------------------------------------------------
function UN = getnodestate(Mesh, U)
    V  = Mesh.V; E = Mesh.E;
    Nv = size(V, 1);
    Ne = size(E, 1);
    UN    = zeros(Nv, 1);
    count = zeros(Nv, 1);
    for e = 1:Ne
        for i = 1:3
            n = E(e, i);
            UN(n)    = UN(n)    + U(e);
            count(n) = count(n) + 1;
        end
    end
    UN = UN ./ count;
end


%-----------------------------------------------------------
function plotmesh(Mesh, fname)
    V  = Mesh.V; E = Mesh.E; BE = Mesh.BE;
    f  = figure('Position', [100, 100, 1200, 1200]);
    triplot(triangulation(E, V(:,1), V(:,2)), 'k-');
    hold on;
    for i = 1:size(BE, 1)
        plot(V(BE(i,1:2), 1), V(BE(i,1:2), 2), '-', 'LineWidth', 1, 'Color', 'black');
    end
    hold off;
    axis equal; axis off;
    if ~isempty(fname)
        saveas(f, fname);
    else
        set(f, 'visible', 'on');
    end
    close(f);
end


%-----------------------------------------------------------
function plotstate(Mesh, U, p, field, frange, fname)

    V  = Mesh.V; E = Mesh.E; BE = Mesh.BE;
    f  = figure(5);%'Position', [100, 100, 1200, 1200]);

    U_mod = getField(U, field);
    F     = U_mod(:, 9);  % Mach number (column 8)

    if p == 0
        trisurf(E, V(:,1), V(:,2), zeros(size(V,1),1), ...
            'FaceColor', 'flat', 'FaceVertexCData', F, ...
            'EdgeColor', 'none');
        view(2);
    else
        if ~isempty(frange)
            vc = linspace(frange(1), frange(2), 21);
        else
            vc = 20;
        end
        UN = getnodestate(Mesh, F);
        tricontourf_matlab(V, E, UN, vc);  % see note below
    end

    hold on;
    for i = 1:size(BE, 1)
        plot(V(BE(i,1:2), 1), V(BE(i,1:2), 2), '-', 'LineWidth', 2, 'Color', 'black');
    end
    hold off;

    axis equal; axis off;
    colormap('jet');
    % if ~isempty(frange)
    %     caxis(frange);
    % end
    cb = colorbar('Location', 'southoutside');
    clim([0.71,0.75]);
    cb.FontSize = 16;

    if ~isempty(fname)
        saveas(f, fname);
    else
        set(f, 'visible', 'on');
    end

end


%-----------------------------------------------------------
% Helper: contour fill on unstructured triangular mesh (equivalent to tricontourf)
function tricontourf_matlab(V, E, UN, vc)
    if isscalar(vc)
        tricontourf(V(:,1), V(:,2), E-1, UN, vc);  % requires Mapping Toolbox or R2020b+
    else
        tricontourf(V(:,1), V(:,2), E-1, UN, vc);
    end
end


function Mesh = readgri(fname)
% READGRI  Read a .gri mesh file and return mesh structure.
%   Mesh = readgri(fname) reads the mesh file and returns a struct with:
%     V      - vertices (Nn x 2)
%     E      - elements (Ne x 3), 1-indexed
%     IE     - interior edges (ni x 4): n1, n2, elem1, elem2
%     BE     - boundary edges (nb x 4): n1, n2, elem, group
%     Bname  - cell array of boundary names
%     B_up   - boundary edges for Curve1 (upper surface)
%     B_low  - boundary edges for Curve5 (lower surface)
    fname
    f = fopen(fname, 'r');

    % Read header: Nn, Ne, dim
    header = sscanf(fgetl(f), '%d %d %d');
    Nn = header(1); Ne = header(2);

    % Read vertices
    V = zeros(Nn, 2);
    for n = 1:Nn
        V(n,:) = sscanf(fgetl(f), '%f %f')';
    end

    % Read boundaries
    NB    = sscanf(fgetl(f), '%d');
    B     = {};
    Bname = {};
    B_upper = {};
    B_lower = {};

    for i = 1:NB
        s     = strsplit(strtrim(fgetl(f)));
        Nb    = str2double(s{1});
        bname = s{3};
        Bname{end+1} = bname;

        Bi = zeros(Nb, 2);
        for n = 1:Nb
            row    = sscanf(fgetl(f), '%d %d')';
            Bi(n,:) = row;  % keep 1-indexed (MATLAB convention)
        end
        B{end+1} = Bi;

        if strcmp(bname, 'Curve1')
            B_upper{end+1} = Bi;
        elseif strcmp(bname, 'Curve5')
            B_lower{end+1} = Bi;
        end
    end

    % Read elements
    E   = [];
    Ne0 = 0;
    while Ne0 < Ne
        s  = sscanf(fgetl(f), '%d');
        ne = s(1);
        Ei = zeros(ne, 3);
        for n = 1:ne
            Ei(n,:) = sscanf(fgetl(f), '%d %d %d')';
        end
        E   = [E; Ei];  % keep 1-indexed
        Ne0 = Ne0 + ne;
    end
    fclose(f);

    % Convert boundaries to 0-indexed equivalent for edgehash,
    % then convert back — here we keep everything 1-indexed in MATLAB
    [IE, BE]      = edgehash(E, B);
    [~,  B_up]   = edgehash(E, B_upper);
    [~,  B_low]  = edgehash(E, B_lower);

    Mesh.V     = V;
    Mesh.E     = E;
    Mesh.IE    = IE;
    Mesh.BE    = BE;
    Mesh.Bname = Bname;
    Mesh.B_up  = B_up;
    Mesh.B_low = B_low;
end


function [IE, BE] = edgehash(E, B)
% EDGEHASH  Identify interior and boundary edges.
%   IE contains [n1, n2, elem1, elem2] for each interior edge (1-indexed)
%   BE contains [n1, n2, elem,  group] for each boundary edge (1-indexed)

    Ne = size(E, 1);
    Nn = max(E(:));

    % Sparse hash matrix to track which element owns each edge
    H  = sparse(Nn, Nn);

    IE = zeros(ceil(Ne * 1.5), 4);
    ni = 0;

    for e = 1:Ne
        for i = 1:3
            n1 = E(e, i);
            n2 = E(e, mod(i, 3) + 1);  % wraps: 1->2, 2->3, 3->1

            if H(n2, n1) == 0
                H(n1, n2) = e;
            else
                eR = H(n2, n1);
                ni = ni + 1;
                IE(ni, :) = [n1, n2, e, eR];
                H(n2, n1) = 0;
            end
        end
    end
    IE = IE(1:ni, :);

    % Count total boundary edges
    nb0 = 0;
    for g = 1:length(B)
        nb0 = nb0 + size(B{g}, 1);
    end

    BE = zeros(nb0, 4);
    nb = 0;

    for g = 1:length(B)
        Bi = B{g};
        for b = 1:size(Bi, 1)
            n1 = Bi(b, 1);
            n2 = Bi(b, 2);
            if H(n1, n2) == 0
                tmp = n1; n1 = n2; n2 = tmp;  % swap if reversed
            end
            nb = nb + 1;
            BE(nb, :) = [n1, n2, H(n1, n2), g];
        end
    end
end