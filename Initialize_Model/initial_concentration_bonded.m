function [Floe,bonds,Nbound,Nbond] = initial_concentration(c2_boundary,target_concentration,height, NumFloes, NumSubfloes, min_floe_size)
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
load('kobuk_poly.mat','wholelake', 'lake_with_polynya');
R(1) = subtract(B1,wholelake); R(2) = subtract(B2,wholelake); R(3) = subtract(B3,wholelake); R(4) = subtract(B4,wholelake);

nx=40; ny=4;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;
[xx,yy] = meshgrid(Xc,Yc); X = xx(:); Y = yy(:);%[yc(1) yc(2) yc(2) yc(1)];

Nbound = 0;
Floe = [];
FloeBnds = [];
FloeNum = 1;

for ii = 1:length(R)
    FloeNEW = initialize_floe_values(R(ii),height,NumSubfloes);
    Floe = [Floe FloeNEW];
end
Nbound = length(Floe);
for ii = 1:Nbound
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

%Loop through all the regions of the domain to create new floes
for jj = 1:Ny
    for ii = 1:Nx
        if c(jj,ii)>0
            boundary = polyshape([x(ii) x(ii) x(ii+1) x(ii+1)], [y(jj) y(jj+1) y(jj+1) y(jj)]);
            boundary = intersect(boundary,c2_boundary_poly);
            frac = area(boundary)/(area(boundary)-sum(area(R)));
            N = ceil(NumFloes*frac/c(jj,ii));
            if N == 1
                Floe = initialize_floe_values(c2_boundary_poly,height,NumSubfloes);
            else
                X = 0.975*dx/2*(2*rand(N,1)-1)+(x(ii)+x(ii+1))/2;
                Y = 0.975*dy/2*(2*rand(N,1)-1)+(y(jj)+y(jj+1))/2;
                in = inpolygon(X,Y,boundary.Vertices(:,1),boundary.Vertices(:,2));
                X = X(in); Y = Y(in);
                [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices);

                count = 1; clear polyfloe
                for m = 1:length(b)
                    pnew = polyshape(b{m});
                    pnow = intersect(pnew,lake_with_polynya);
                    if area(pnow)>0
                        polyfloe(count) = pnow;
                        count = count + 1;
                    end
                end
                
                poly = polyfloe;
                
                for iii = 1:length(polyfloe)
                    bonds(iii).bond = [];
                end
                r = sqrt(min_floe_size);

                Nf = 1:length(poly);%randperm(length(b));
                for iii = 1:length(poly)-1
                    for jjj = iii+1:length(poly)
                        [k,dist] = dsearchn(poly(iii).Vertices,poly(jjj).Vertices);
                        if length(k(dist<0.1)) > 1
                            verts = poly(iii).Vertices; verts = verts(k(dist<0.1),:);
                            Xi = (verts(1,1)+verts(2,1))/2; Yi = (verts(1,2)+verts(2,2))/2;
                            Xi1 = (verts(1,1)+Xi)/2; Yi1 = (verts(1,2)+Yi)/2;
                            Xi2 = (verts(2,1)+Xi)/2; Yi2 = (verts(2,2)+Yi)/2;
                            d = sqrt((verts(1,1)-verts(2,1))^2+(verts(1,2)-verts(2,2))^2);
                            dx=(verts(1,1)-verts(2,1))/d; dy = (verts(1,2)-verts(2,2))/d;
                            bonds_i.Num = jjj; bonds_j.Num = iii;
                            vertsy = [Yi1, Yi1+r/2*(dy*sqrt(3)-dx),Yi1+r/2*(-dy*sqrt(3)-dx)];
                            vertsx = [Xi1, Xi1+r/2*(dx*sqrt(3)+dy),Xi1+r/2*(-dx*sqrt(3)+dy)];
                            vertsy2 = [Yi1, Yi1+r/2*(dy*sqrt(3)+dx),Yi1+r/2*(-dy*sqrt(3)+dx)];
                            vertsx2 = [Xi1, Xi1+r/2*(dx*sqrt(3)-dy),Xi1+r/2*(-dx*sqrt(3)-dy)];
                            poly_bnd1 = polyshape(vertsx,vertsy);
                            poly_bnd1 = translate(poly_bnd1,-r/10*[dy -dx]);
                            poly_bnd2 = polyshape(vertsx2,vertsy2);
                            poly_bnd2 = translate(poly_bnd2,r/10*[dy -dx]);
                            vertsy = [Yi2, Yi2+r/2*(dy*sqrt(3)-dx),Yi2+r/2*(-dy*sqrt(3)-dx)];
                            vertsx = [Xi2, Xi2+r/2*(dx*sqrt(3)+dy),Xi2+r/2*(-dx*sqrt(3)+dy)];
                            vertsy2 = [Yi2, Yi2+r/2*(dy*sqrt(3)+dx),Yi2+r/2*(-dy*sqrt(3)+dx)];
                            vertsx2 = [Xi2, Xi2+r/2*(dx*sqrt(3)-dy),Xi2+r/2*(-dx*sqrt(3)-dy)];
                            poly_bnd3 = polyshape(vertsx,vertsy);
                            poly_bnd3 = translate(poly_bnd3,-r/10*[dy -dx]);
                            poly_bnd4 = polyshape(vertsx2,vertsy2);
                            poly_bnd4 = translate(poly_bnd4,r/10*[dy -dx]);
                            poly_bnd = union(poly_bnd1,poly_bnd2);
                            if area(intersect(poly_bnd,poly(iii)))>area(intersect(poly_bnd,poly(jjj)))
                                bonds_i.c(1) = poly_bnd1; bonds_i.c(2) = poly_bnd3;
                                bonds_j.c(1) = poly_bnd2; bonds_j.c(2) = poly_bnd4;
                                bonds_i.r_bnd = sqrt(area(poly_bnd1)); bonds_j.r_bnd = sqrt(area(poly_bnd2));
                            else
                                bonds_i.c(1) = poly_bnd2; bonds_i.c(2) = poly_bnd4;
                                bonds_j.c(1) = poly_bnd1; bonds_j.c(2) = poly_bnd3;
                                bonds_i.r_bnd = sqrt(area(poly_bnd2)); bonds_j.r_bnd = sqrt(area(poly_bnd1));
                            end

                            bonds_i.L = d; bonds_j.L = d;
                            bonds_i.Xb = []; bonds_j.Xb = [];
                            bonds_i.Yb = []; bonds_j.Yb = [];
                            bonds_i.Xc(1) = Xi1; bonds_j.Xc(1) = Xi1;
                            bonds_i.Xc(2) = Xi2; bonds_j.Xc(2) = Xi2;
                            bonds_i.Yc(1) = Yi1; bonds_j.Yc(1) = Yi1;
                            bonds_i.Yc(2) = Yi2; bonds_j.Yc(2) = Yi2;
                            bonds_i.Fx_p = 0; bonds_j.Fx_p = 0;
                            bonds_i.Fy_p = 0; bonds_j.Fy_p = 0;
                            bonds_i.Stress = [0;0]; bonds_j.Stress = [0;0];
                            bonds(iii).bond = [bonds(iii).bond bonds_i];
                            bonds(jjj).bond = [bonds(jjj).bond bonds_j];

                        end
                    end
                end

                Atot = 0;
                count = 1;
                polyFloes = [];

                polyFloes = poly;
                for kk= 1:length(polyFloes)
                    polyFloesNew = polyFloes(kk);
                    bond_tmp = bonds(kk).bond;
                    new_bnds = [];

                    if polyFloesNew.NumRegions > 1
                        polyout = sortregions(polyFloesNew,'area','descend');
                        R = regions(polyout);
                        polyFloesNew = R(1);
                    end
                    floenew = initialize_floe_values(polyFloesNew,height,NumSubfloes);
                    floenew.num = Nf(kk);
                    if ~isempty(bond_tmp)
                        for iii = 1:length(bond_tmp)
                            for jjj = 1:length(bond_tmp(iii).c)
                                bond_tmp2 = bond_tmp(iii);
                                [Xb,Yb] = centroid(bond_tmp(iii).c(jjj));  bond_tmp2.Xb = Xb-floenew.Xi; bond_tmp2.Yb = Yb-floenew.Yi;
                                bond_tmp2.Xc = bond_tmp(iii).Xc(jjj)-floenew.Xi; bond_tmp2.Yc = bond_tmp(iii).Yc(jjj)-floenew.Yi;
                                new_bnds = [new_bnds; bond_tmp2];
                            end

                        end
                        new_bnds=rmfield(new_bnds,{'c'});

                        floenew.bonds = new_bnds;
                    else
                        xx = 1; xx(1) =[1 2];
                    end
                    if abs(Nf(kk)-kk)>0
                        xx = 1; xx(1) =[1 2];
                    end
                    Floe = [Floe floenew];
                end
            end
        end
    end
end
Nums = cat(1,Floe.num);
for ii = Nbound+1:length(Floe)
    floe = Floe(ii);
    bnds = unique(cat(1,Floe(ii).bonds.Num));
    bonds1 = cat(1,Floe(ii).bonds.Num);
    for jj = 1:length(bnds)
        Lia = ismember(Nums,bnds(jj)); floe2 = Floe(Lia); floeNum = floe2.num;
        Num1 = sum(ismember(bonds1,floeNum)); BondNum2 = cat(1,floe2.bonds.Num); Num2 = sum(ismember(BondNum2,Floe(ii).num));
        if ~Num1
            xx = 1; xx(1) =[1 2];
        elseif ~Num2
            xx = 1; xx(1) =[1 2];
        elseif sum(Num1)~=sum(Num2)
            xx = 1; xx(1) =[1 2];
        end
    end
end

Nbond = 0;

end

