function [Floe, Nb] = initial_concentration_kobuk(c2_boundary,target_concentration,height, NumFloes, min_floe_size)
    % This function is used to generate the initial floe field
    % target_concentration is a scalar used to determine what percentage of the region should be covered
    [Ny, Nx] = size(target_concentration);
    c = flipud(target_concentration);
    x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/1:max(c2_boundary(1,:));
    y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/1:max(c2_boundary(2,:));
    c2_boundary_poly = polyshape(c2_boundary');
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    Lx = max(x)*2; Ly= max(y)*2;
    
    %Create floes that act as boundaries and wont move
    x1 = [-Lx/2 0 0 -Lx/2]; y1 = [Ly/2 Ly/2 5e3 5e3]; B1 = polyshape(x1,y1);
    x2 = [Lx/2 0 0 Lx/2]; y2 = [Ly/2 Ly/2 5e3 5e3]; B2 = polyshape(x2,y2);
    x3 = [-Lx/2 0 0 -Lx/2]; y3 = [-Ly/2 -Ly/2 5e3 5e3]; B3 = polyshape(x3,y3);
    x4 = [Lx/2 0 0 Lx/2]; y4 = [-Ly/2 -Ly/2 5e3 5e3]; B4 = polyshape(x4,y4);
    load('kobuk_poly.mat','upperlake')
    B1 = subtract(B1,upperlake); B2 = subtract(B2,upperlake); B3 = subtract(B3,upperlake); B4 = subtract(B4,upperlake);
    Floe1 = initialize_floe_values(B1,height, 0);
    Floe2 = initialize_floe_values(B2,height, 0);
    Floe3 = initialize_floe_values(B3,height, 0);
    Floe4 = initialize_floe_values(B4,height, 0);
    bound = subtract(c2_boundary_poly, upperlake);
    Floe = [Floe1 Floe2 Floe3 Floe4];
    
    Nb = 4;

    % Assume polyIn is your input polyshape
    polyIn = ...  % your full polygon

    % Get the bounding box of the polyshape
    [xmin, xmax] = bounds(polyIn.Vertices(:,1));
    [ymin, ymax] = bounds(polyIn.Vertices(:,2));
    ymid = (ymin + ymax)/2;

    % Define a rectangle that spans the bottom half
    bottomRect = polyshape([xmin xmax xmax xmin], [ymin ymin ymid ymid]);

    % Take intersection
    bottomHalf = intersect(polyIn, bottomRect);

    polynya = polyshape()
    mask_polygon = subtract(upperlake,);  % Manually define this (e.g., from shapefile or polygon vertices)

    % Initialize one Floe with all subfloes
    Floe = initialize_floe_values(subfloe_polys, height, 300);  % ← pass an array of polyshapes

end