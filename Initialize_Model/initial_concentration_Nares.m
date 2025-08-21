function [Floe, Nb] = initial_concentration_Nares(c2_boundary,target_concentration,height, NumFloes, min_floe_size,ISLANDS)
%% This function is used to generate the initial floe field

%Identify the grids to align with the concentrations specified by the input
[Ny, Nx] = size(target_concentration);
c = flipud(target_concentration);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
dx = x(2)-x(1);
dy = y(2)-y(1);
Lx = abs(min(x));Ly = abs(min(y));
Bx = [-Lx -Lx Lx    Lx];
By = [-Ly  -3e5 -3e5  -Ly];
B = polyshape(Bx',By');

%Create floes that act as boundaries and wont move
load( './Initialize_Model/Kennedy','R')
%load( './Nares/Nares_Strait_segments','R')
% polyout = sortregions(Nares,'area','ascend');
% R = regions(polyout);

%Remove islands
if ~ISLANDS
 %R(1:4) = [];
 R(1:3) = [];
end

Floe = []; bound = c2_boundary_poly;
for ii = 1:length(R)
    FloeNEW = initialize_floe_values(R(ii),height,1);
    bound = subtract(bound, FloeNEW.poly);
    Floe = [Floe FloeNEW];
end
N = length(Floe);
for ii = 1:N
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

% Floes1=fracture_floe(Floe(N-2),10,Floe);
% Floes2=fracture_floe(Floe(N-1),10,Floe);
% Floes3=fracture_floe(Floe(N),10,Floe);
% Floes4=fracture_floe(Floe(N-3),3,Floe);
% Floes5=fracture_floe(Floe(N-4),3,Floe);
% FloeNEW = [Floes1 Floes2 Floes3 Floes4 Floes5];
% for ii = 1:length(FloeNEW)
%     FloeNEW(ii).poly = polyshape(FloeNEW(ii).c_alpha'+[FloeNEW(ii).Xi FloeNEW(ii).Yi]);
% end
% Floe = [Floe(1:N-5) FloeNEW];
Bu = union([Floe.poly]);
Bu = union(Bu,B);

Nb = length(Floe);
%Loop through all the regions of the domain to create new floes
for jj = 1:Ny
    for ii = 1:Nx
        if c(jj,ii)>0
            boundary = polyshape([x(ii) x(ii) x(ii+1) x(ii+1)], [y(jj) y(jj+1) y(jj+1) y(jj)]);
            polyB = intersect(bound,boundary); %Use these when having
%             boundaries
            N = 4*ceil(NumFloes*area(polyB)/area(bound)/c(jj,ii)); %Use these when having
            
%             boundaries
            X = 0.95*dx/2*(2*rand(N,1)-1)+(x(ii)+x(ii+1))/2;
            Y = 0.95*dy/2*(2*rand(N,1)-1)+(y(jj)+y(jj+1))/2;
            in = inpolygon(X,Y,polyB.Vertices(:,1),polyB.Vertices(:,2));
            [~, b,~,~,~] = polybnd_voronoi([X(in) Y(in)],boundary.Vertices);
%             [~, b,~,~,~] = polybnd_voronoi([X(in) Y(in)],polyB.Vertices); %%Use these when having
%             boundaries

            clear nanb
             for kk = 1:length(b)
                 nanb(kk) = isnan(max(max(b{kk})));
                 if length(b{kk})<3
                     nanb(kk) = true;
                 end
             end
             b(nanb == 1) = [];
             for kk = 1:length(b)
                 poly(kk) = polyshape(b{kk});
             end

            %load('FloeInitWide.mat','poly');
            poly = intersect([poly],polyB);
            polyout=subtract([poly],Bu);
            clear polyfloe; polyfloe = [];
            for kk = 1:length(polyout)
                R = regions(polyout(kk));
                polyfloe = [polyfloe; R];
            end
            plot(polyfloe)
            error()
            polyIce = [];
            for count = 1:length(polyfloe)
                poly = polyfloe(count);
                R = regions(poly);
                for kk = 1:length(R)
                    if R(kk).NumHoles > 0
                        p = holes(R(kk));
                        if length(p)>1
                            xx = 1; xx(1) =[1 2];
                        end
                        for iii = 1:length(p)
                            pfull = rmholes(R(kk));
                            [Xi2,Yi2] = centroid(p);
                            L = [ Xi2 Yi2; Xi2+1 Yi2];
                            PcT = cutpolygon([pfull.Vertices; pfull.Vertices(1,:)], L, 1);
                            pnew = polyshape(PcT);
                            pnew = subtract(pnew,union(p));
                            PcB = cutpolygon([pfull.Vertices; pfull.Vertices(1,:)], L, 2);
                            pnew2 = polyshape(PcB);
                            pnew2 = subtract(pnew2,union(p));
                            poly = [regions(pnew); regions(pnew2)];
                        end
                    else
                        poly = R(kk);
                    end
                    polyIce = [polyIce; poly];
                end
            end
            Atot = 0;
            count = 1;

            Nf = 1:length(polyIce);%randperm(length(b));
            while Atot/area(polyB)<=c(jj,ii)
                    floenew = initialize_floe_values(polyIce(Nf(count)),height);
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
areas(1:Nb) = min_floe_size;
Floe(areas<min_floe_size)=[];
end

% bonds.Num = [];
% bonds.r_bnd = [];
% bonds.L = [];
% bonds.Xb = [];
% bonds.Yb = [];
% bonds.Xc = [];
% bonds.Yc = [];
% bonds.Fx_p = 0;
% bonds.Fy_p = 0;
% bonds.Stress = [0;0];
% for ii = 1:Nb
%     Floe(ii).bonds = bonds;
% end

% polyout = sortregions(kennedy,'area','ascend');
% R = regions(polyout);
% 
% new_poly = [];
% for i = 4:length(R)
%     N = 5;
%     poly = R(i)
%     [Xi,Yi] = centroid(poly);
%     c_alpha = [(poly.Vertices-[Xi Yi])' [poly.Vertices(1,1)-Xi; poly.Vertices(1,2)-Yi]];
%     rmax = sqrt(max(sum((poly.Vertices' - [Xi;Yi]).^2,1)));
%     X = rmax*(2*rand(N,1) - 1);
%     Y = rmax*(2*rand(N,1) - 1);
%     boundary = [-1 -1 1 1; -1 1 1 -1]*rmax; boundary = polyshape(boundary');
%     [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices);
%     for kk = 1:length(b)
%         poly(kk) = polyshape(b{kk});
%     end
%     polynew = intersect([poly],polyshape(c_alpha'));
%     A = area(polynew);
%     polynew = polynew(A>0);
%     polynew = translate(polynew,[Xi,Yi]);
%     new_poly = [new_poly polynew];
% end

% polysnew = polyfloe;
% load( 'Kennedy','R')
% Bu = union(R);
% polyfloe = [];
% for ii=1:length(polysnew)
%     polynew = subtract(polysnew(ii),Bu);
%     polyout = sortregions(polynew,'perimeter','descend');
%     R = regions(polyout);
%     poly = [];
%     for kk = 1:length(R)
%         if R(kk).NumHoles > 0
%             p = holes(R(kk));
%             for iii = 1:length(p)
%                 pfull = rmholes(R(kk));
%                 [Xi2,Yi2] = centroid(p);
%                 L = [ Xi2 Yi2; Xi2+1 Yi2];
%                 PcT = cutpolygon([pfull.Vertices; pfull.Vertices(1,:)], L, 1);
%                 pnew = polyshape(PcT);
%                 pnew = subtract(pnew,union(p));
%                 PcB = cutpolygon([pfull.Vertices; pfull.Vertices(1,:)], L, 2);
%                 pnew2 = polyshape(PcB);
%                 pnew2 = subtract(pnew2,union(p));
%                 poly = [regions(pnew); regions(pnew2)];
%             end
%         else
%             poly = R(kk);
%         end
%         polyfloe = [polyfloe; poly];
%     end
% end
