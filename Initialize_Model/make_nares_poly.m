load latlonNares2

polyGL = polyshape(GL.lon,GL.lat) 
polyCA = polyshape(CA.lon,CA.lat) 

close all
plot(polyCA)
hold on
plot(polyGL)                                        

crozier = [80.500278, -67.287500;...
80.503333, -67.278333;...
80.505000, -67.262222;...
80.506111, -67.255833;...
80.506944, -67.234722;...
80.505000, -67.229722;...
80.505556, -67.222778;...
80.504167, -67.204167;...
80.501944, -67.193611;...
80.496944, -67.190000;...
80.491667, -67.170278;...
80.483056, -67.162500;...
80.481667, -67.182222;...
80.483056, -67.210556;...
80.485278, -67.217222;...
80.487222, -67.237778;...
80.493889, -67.256389;...
80.500556, -67.287222];


Lat_cr = crozier(:,1);
Lon_cr = crozier(:,2);
polyCR = polyshape(Lon_cr,Lat_cr);

Hans = [80.831667, -66.448333;...
80.828333, -66.438333;...
80.827500, -66.419722;...
80.823333, -66.416111;...
80.820833, -66.440278;...
80.822778, -66.479444;...
80.824167, -66.488611;...
80.826944, -66.490833;...
80.830278, -66.476389;...
80.831667, -66.457500;...
80.831389, -66.452500];

Lat_hans = Hans(:,1);
Lon_hans = Hans(:,2);
polyHans = polyshape(Lon_hans,Lat_hans);


close all
plot(polyCA)
hold on
plot(polyGL) 
plot(polyCR)
plot(polyHans);
fig = figure(1);
exportgraphics(fig,['polys.jpg'] ,'resolution',300); 

%% Convert from lat/lon

% Default parameters for NSIDC/SCAR polar stereographic projection (south pole)
earthradius = 6378137.0; % WGS84 radius
eccentricity = 0.08181919; % Earth's misshapenness
lat_true = -70; % Latitude of true scale (standard parallel)
lon_posy = -45; % Longitude along positive Y axis

R = regions(polyCA)
for ii = 1:length(R)
    clear Xca; clear Yca
    for jj = 1:length(R(ii).Vertices(:,2))
        [Xca(jj),Yca(jj)] = polarstereo_fwd(R(ii).Vertices(jj,2),R(ii).Vertices(jj,1), earthradius, eccentricity, lat_true, lon_posy);
    end
    polyCAps(ii) = polyshape(Xca,Yca);
end

R = regions(polyGL)
for ii = 1:length(R)
    clear Xgl; clear Ygl
    for jj = 1:length(R(ii).Vertices(:,2))
        [Xgl(jj),Ygl(jj)] = polarstereo_fwd(R(ii).Vertices(jj,2),R(ii).Vertices(jj,1), earthradius, eccentricity, lat_true, lon_posy);
    end
    polyGLps(ii) = polyshape(Xgl,Ygl);
end

for jj = 1:length(polyCR.Vertices(:,2))
    [Xcr(jj),Ycr(jj)] = polarstereo_fwd(polyCR.Vertices(jj,2),polyCR.Vertices(jj,1), earthradius, eccentricity, lat_true, lon_posy);
end
polyCRps = polyshape(Xcr,Ycr);

for jj = 1:length(polyHans.Vertices(:,2))
    [Xhans(jj),Yhans(jj)] = polarstereo_fwd(polyHans.Vertices(jj,2),polyHans.Vertices(jj,1), earthradius, eccentricity, lat_true, lon_posy);
end
polyHANSps = polyshape(Xhans,Yhans);

close all
plot(polyGLps)
hold on
plot(polyCAps)
plot(polyCRps)
plot(polyHANSps)


fig = figure(1);
exportgraphics(fig,['poly.jpg'] ,'resolution',300); 

Nares = union([polyGLps]);
Nares = union(Nares,union([polyCAps]));
Nares = union(Nares,polyCRps);
Nares = union(Nares, polyHANSps);

[Xc,Yc] = centroid(Nares);
Nares = translate(Nares,-[Xc,Yc]);

polynew = rotate(Nares,12.5);
polynew = scale(polynew,1/145);
polynew = translate(polynew,[-0.61e5,-0.325e5])
close all
%plot(polynew)
%hold on
%plot([Floe(1:Nb).poly])

xbox = [-30000 -30000 30000 30000];
ybox = [-50000 50000 50000 -50000];
pbox = polyshape(xbox, ybox);
%plot(pbox)

poly = intersect(pbox,polynew)
R = regions(poly);
R(area(R)<500)=[];
area(R)

ratio = 2; %Ly/Lx;%
fig=figure('Position',[10 10 200 200*ratio],'visible','on');  
set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
figure(fig)
clf(fig);

plot(R)

fig = figure(1);
exportgraphics(fig,['polynew.jpg'] ,'resolution',300); 
