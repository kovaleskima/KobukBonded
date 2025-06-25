function [floes,Floe0] = FloeSimplify(floe,PACK,Floe0,polyboundary)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% Remap the main polygon to a shape with fewer vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
floes = [];
rho_ice = 920;
save('FloeS.mat','floe','Floe0')


% if isfield(floe,'poly')
%     if abs(centroid(floe.poly))<10
%         xx=1; xx(1) = [1 2];
%     end
% end

x=cat(1,Floe0.Xi);
y=cat(1,Floe0.Yi);
rmax=cat(1,Floe0.rmax);
alive=cat(1,Floe0.alive);
k=1;
floe.potentialInteractions = [];
for j=1:length(Floe0)
   % Floe0(j).poly = polyshape(Floe0(j).c_alpha'+[Floe0(j).Xi Floe0(j).Yi]);
    %FloeOld(j).poly = polyshape(Floe0(j).c_alpha'+[Floe0(j).Xi Floe0(j).Yi]);
    if alive(j) && sqrt((floe.Xi-x(j))^2 + (floe.Yi-y(j))^2) > 1 && sqrt((floe.Xi-x(j))^2 + (floe.Yi-y(j))^2)<(floe.rmax+rmax(j)) % if floes are potentially overlapping
        floe.potentialInteractions(k).floeNum=j;
        floe.potentialInteractions(k).c=[Floe0(j).c_alpha(1,:)+x(j); Floe0(j).c_alpha(2,:)+y(j)];
        floe.potentialInteractions(k).Ui=Floe0(j).Ui;
        floe.potentialInteractions(k).Vi=Floe0(j).Vi;
        floe.potentialInteractions(k).h=Floe0(j).h;
        floe.potentialInteractions(k).area=Floe0(j).area;
        floe.potentialInteractions(k).Xi=x(j);
        floe.potentialInteractions(k).Yi=y(j);
        floe.potentialInteractions(k).ksi_ice = Floe0(j).ksi_ice;
        k=k+1;
    end
end


%Set function to define the tolerance for Douglas-Peuker algorithm
SimpMin = @(A) log10(A)^3.5;

%Create new simplified polyshape
floenew = floe;
floenew.poly = polyshape(floenew.c_alpha'+[floenew.Xi floenew.Yi]);
%vertnew = DouglasPeucker(floe.poly.Vertices,SimpMin(floe.area));
[vertx,verty] = reducem(floe.c0(1,:)',floe.c0(2,:)',1000);
pnew = polyshape(vertx, verty);
if ~isempty(polyboundary)
    pnew = subtract(pnew,polyboundary); 
end
Atot = sum(area(pnew));
if isinf(Atot)
  save('floefail.mat','floe','pnew');
%    1
elseif isnan(Atot)
  save('floefail.mat','floe','pnew');
%   2
elseif Atot ==0
  save('floefail.mat','floe','pnew');
%   3
  R = [];
else
  polynew = scale(pnew,sqrt(floe.area/Atot));

  %Align center of old polygon with the enw one
  [x1,y1] = centroid(polynew);
  dx = floe.Xi-x1;
  dy = floe.Yi-y1;
  % if area(intersect(floe.poly,polynew))/floe.area > 0.95
  %     polynew = translate(polynew,[dx, dy]);
  % end

  %Check if simplifcation led to polygon having multiple regions
  polyout = sortregions(polynew,'area','descend');
  R = regions(polyout);
  R = R(area(R)>1e4);
  Atot = sum(area(R));
end

for jj = 1:length(R)
    for ii = 1:length(floe.potentialInteractions)
        if ~isinf(floe.potentialInteractions(ii).floeNum) & ~isempty(R)
            polyq = rmholes(R(jj));
            [Xi, Yi] = centroid(polyq);
            Xt = Xi+floe.Xi; Yt = Yi+floe.Yi;
            polyq = translate(polyq,[Xt,Yt]);
            poly = polyshape(floe.potentialInteractions(ii).c');
            polyI = intersect(poly,polyq);
            if area(polyI)/area(poly)>0.4
                %if isempty(Floe0(floe.potentialInteractions(ii).floeNum).poly)
                %    FloeO(floe.potentialInteractions(ii).floeNum).poly = polyshape(Floe0(floe.potentialInteractions(ii).floeNum).c_alpha'+[Floe0(floe.potentialInteractions(ii).floeNum).Xi Floe0(floe.potentialInteractions(ii).floeNum).Yi]);
                %save('FloeS.mat','floe','Floe0')
                %floeold = floenew;
                %end
                floenew = FuseFloes(floenew,Floe0(floe.potentialInteractions(ii).floeNum));
                Floe0(floe.potentialInteractions(ii).floeNum).alive = 0;
                if length(floenew)>1
                    m = cat(1,floenew.mass); u = cat(1,floenew.Ui); v = cat(1,floenew.Vi);
                    dksi = cat(1,floenew.dksi_ice_p); inertia_moment = cat(1,floenew.inertia_moment);
                    [~,I] = max(m);
                    height.mean = sum(m)/(area(floenew(I).poly)*rho_ice);
                    height.delta = 0;
                    tmp = initialize_floe_values(floenew(I).poly, height);
                    tmp.mass = sum(m); tmp.Ui = sum(m.*u)/tmp.mass; tmp.Vi = sum(m.*v)/tmp.mass;
                    tmp.dksi_ice_p = sum(dksi.*inertia_moment)/tmp.inertia_moment;
                    tmp.FxOA = [];%0;
                    tmp.FyOA = [];
                    floenew = tmp;
%                save('FloeS.mat','floe','Floe0')
%                    xx = 1; xx(1) =[1 2];
                end
            end
        end
    end
end


%% Calculate the new properties associated with this floe since it has a new shape
if length(R) == 1
    floes = floenew;
%     if R.NumHoles > 0
%         xx = 1;
%         xx(1)=[1 2];
%     end
    R = rmholes(R);
    floes.area = area(R);
    [Xi, Yi] = centroid(R);
    floes.Xi = Xi+floenew.Xi; floes.Yi = Yi+floenew.Yi;
    floes.h = floes.mass/(rho_ice*floes.area);
    floes.c0 = [R.Vertices(:,1)-Xi,R.Vertices(:,2)-Yi; R.Vertices(1,1)-Xi,R.Vertices(1,2)-Yi]';
    floes.angles = polyangles(R.Vertices(:,1),R.Vertices(:,2));
    A_rot=[cos(floes.alpha_i) -sin(floes.alpha_i); sin(floes.alpha_i) cos(floes.alpha_i)]; %rotation matrix
    floes.c_alpha=A_rot*floes.c0;
    floes.inertia_moment = PolygonMoments(floes.c_alpha',floes.h);
    floes.rmax = max(sqrt(floes.c_alpha(1,:).^2+floes.c_alpha(2,:).^2));
    floes.X = floes.rmax*(2*rand(1000,1) - 1);
    floes.Y = floes.rmax*(2*rand(1000,1) - 1);
    floes.A = inpolygon(floes.X,floes.Y,floes.c_alpha(1,:),floes.c_alpha(2,:));
    poly = polyshape(floes.c_alpha');
    if abs(centroid(poly))>100
        xx=1; xx(1) = [1 2];
    end
elseif length(R) > 1
    for ii = 1:length(R)
        floenew = floe;
        poly1new = R(ii);
        polya = rmholes(poly1new);
%         if poly1new.NumHoles > 0
%             xx = 1;
%             xx(1)=[1 2];
%         end
        if ~PACK
            poly1new = polya;
%             xx = 1; xx(1) =[1 2];
        end
        [Xi,Yi] = centroid(poly1new);
        floenew.area = area(poly1new);
        poly1new = translate(poly1new,[floe.Xi,floe.Yi]);
        floenew.poly = poly1new;
        floenew.mass = floe.mass*area(R(ii))/Atot;
        floenew.h = floenew.mass/(rho_ice*floenew.area);
        floenew.c_alpha = [(polya.Vertices-[Xi Yi])' [polya.Vertices(1,1)-Xi; polya.Vertices(1,2)-Yi]];
        floenew.angles = polyangles(polya.Vertices(:,1),polya.Vertices(:,2));
        floenew.c0 = floenew.c_alpha;
        floenew.inertia_moment = PolygonMoments(floenew.c0',floenew.h);
        floenew.rmax = max(sqrt(floenew.c_alpha(1,:).^2+floenew.c_alpha(2,:).^2));
        floenew.X = floenew.rmax*(2*rand(1000,1) - 1);
        floenew.Y = floenew.rmax*(2*rand(1000,1) - 1);
        floenew.A = inpolygon(floenew.X,floenew.Y,floenew.c_alpha(1,:),floenew.c_alpha(2,:));
%         floenew.Xg = floe.Xg;
%         floenew.Yg = floe.Yg;
%         floenew.X = floe.X;
%         floenew. Y = floe.Y;
%         
%         [in] = inpolygon(floenew.X(:)+Xi, floenew.Y(:)+Yi,floenew.poly.Vertices(:,1),floenew.poly.Vertices(:,2));
%         floenew.A=reshape(in,length(floenew.X),length(floenew.X));
        
        floenew.Xi = floe.Xi+Xi; floenew.Yi = floe.Yi+Yi; floenew.alive = 1;
        floenew.alpha_i = 0; floenew.Ui = floe.Ui; floenew.Vi = floe.Vi;
        floenew.dXi_p = floe.dXi_p; floenew.dYi_p = floe.dYi_p;
        floenew.dUi_p = floe.dUi_p; floenew.dVi_p = floe.dVi_p;
        floenew.dalpha_i_p = 0; floenew.ksi_ice = floenew.area/floe.area*floe.ksi_ice;
        floenew.dksi_ice_p = floe.dksi_ice_p;
        floenew.interactions = [];
        floenew.potentialInteractions = [];
        floenew.collision_force = 0;
        floenew.collision_torque = 0;
        floenew.OverlapArea = 0;
        
         poly = polyshape(floenew.c_alpha');
%         if abs(centroid(poly))>100
%             xx=1; xx(1) = [1 2];
%         end
%         if isfield(floe,'poly')
%             if area(intersect(floenew.poly,floe.poly)) < 10 && floenew.area > 67600000
%                 xx=1; xx(1) = [1 2];
%             end
%         end
        floes = [floes floenew];
    end
else
    floes = floe;
    floes.alive = 0;
end
if isfield(floes,'potentialInteractions')
    floes=rmfield(floes,{'potentialInteractions'});
end
%if floe.mass/sum(cat(1,floes.mass))-1 > 1e-3
%    xx = 1;
%    xx(1) = [1 2];
%end
%if sum(cat(1,floes.mass))/floe.mass-1 > 1e-3
%    xx = 1;
%    xx(1) = [1 2];
%end
% h = cat(1,floes.h);
% if max(h)/floe.h-1 > 0.01
%     xx = 1;
%     xx(1) = [1 2];
% end
% for ii = 1:length(floes)
%     if isempty(floes(ii).SubFloes.inertia)
%         xx=1;
%         xx(1) = [1 2];
%     end
% end
% if abs(floe.area/sum(cat(1,floes.area))-1)>0.05
%     xx = 1;
%     xx(1) =[1 2];
% end
% for ii = 1:length(floes)
%     if abs(floes(ii).area/area(polyshape(floes(ii).c_alpha'))-1)>1e-3
%         xx = 1;
%         xx(1) =[1 2];
%     end
% end


warning('on',id)
warning('on',id2)
warning('on',id3)
end

