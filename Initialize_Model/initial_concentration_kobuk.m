function [Floe, Nb] = initial_concentration_kobuk(c2_boundary,target_concentration,height, NumFloes, min_floe_size)
%% This function is used to generate the initial floe field

%Identify the grids to align with the concentrations specified by the input
[Ny, Nx] = size(target_concentration);
c = flipud(target_concentration);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
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
%Loop through all the regions of the domain to create new floes
for jj = 1:Ny
    for ii = 1:Nx
        if c(jj,ii)>0
            boundary = polyshape([x(ii) x(ii) x(ii+1) x(ii+1)], [y(jj) y(jj+1) y(jj+1) y(jj)]);
            polyB = intersect(c2_boundary_poly,boundary);
            N = 4*ceil(NumFloes*area(polyB)/area(c2_boundary_poly)/c(jj,ii));
            X = 0.95*dx/2*(2*rand(N,1)-1)+(x(ii)+x(ii+1))/2;
            Y = 0.95*dy/2*(2*rand(N,1)-1)+(y(jj)+y(jj+1))/2;
            in = inpolygon(X,Y,polyB.Vertices(:,1),polyB.Vertices(:,2));
            [~, b,~,~,~] = polybnd_voronoi([X(in) Y(in)],boundary.Vertices);

            clear poly
            for kk = 1:length(b)
                poly(kk) = polyshape(b{kk});
            end
            polyout=subtract([poly],bound);
            clear polyfloe; polyfloe = [];
            for kk = 1:length(polyout)
                R = regions(polyout(kk));
                polyfloe = [polyfloe; R];
            end
            if isempty(polyfloe)
                continue;  % skip this grid cell (ii, jj) if no valid floe shapes were created
            end
            Nf = randperm(length(polyfloe));
            Atot = 0;
            count = 1;

            while Atot/area(polyB)<= c(jj,ii)
                    floenew = initialize_floe_values(polyfloe(Nf(count)),height, 0);
                    Floe = [Floe floenew];
                    count = count+1;
                    Atot = Atot+floenew.area;
                    if count > length(Nf)
                        Atot = area(polyB)+1;
                    end
            end
        end
    end
end
areas = cat(1,Floe.area);
Floe(areas<min_floe_size)=[];

disp('initial_concentration_kobuk was called')
end

