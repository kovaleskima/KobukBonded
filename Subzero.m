%% Set Flags

RIDGING=false; 

FRACTURES=true;

PERIODIC=false;

PACKING = false;

WELDING = false;

CORNERS = false;

COLLISION = true;

AVERAGE = true;

RAFTING = false;

KEEP_MIN = true;

SIMPLIFY = false;

ifPlot = true; %Plot floe figures or not?

ifPlotStress = false;


%% Initialize model vars

%add paths
paths

dt=1; %Time step in sec
height.mean = 0.1;
height.delta = 0;

%Define ocean currents
nDTpack = 5500;
rho_ice=920;
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack);

%Define 10m winds
winds.x = ocean.Xocn; winds.y = ocean.Yocn;
U0 = 10; V0 = -10;
winds.u=U0*ones(size(ocean.Xocn));
winds.v=V0*ones(size(ocean.Xocn));

%Define boundaries and floes that act as boundaries
side = 1;
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[c2_boundary(1,:); c2_boundary(2,:)]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height,1);
uright = 0; uleft = 0; %Define speeds that boundaries might be moving with
min_floe_size = 2*Lx*Ly/10000;% Define the minimum floe size you want in initial configuration

target_concentration = 1;  % initialize all to open ocean

% Redefine initial floe state how I actually want them
NumFloes = 1200;
[Floe, bonds, Nb, Nbond] = initial_concentration_again(c2_boundary,target_concentration,height, NumFloes, 0, min_floe_size);
Floe0 = Floe;
Nums = cat(1,Floe.num);
for ii = 1+Nb:length(Floe)
    floe = Floe(ii);
    bnds = unique(cat(1,Floe(ii).bonds.Num));
    bonds1 = cat(1,Floe(ii).bonds.Num);
end


%Define Modulus for floe interactions
global Modulus r_mean L_mean
Modulus = 7.5e5*(mean(sqrt(cat(1,Floe.area)))+min(sqrt(cat(1,Floe.area))));
r_mean = mean(sqrt(cat(1,Floe.area)));
L = [];
for ii = 1:length(Floe)
    L_tmp=cat(1,Floe(ii).bonds.L);
    L = [L; L_tmp];
end
L_mean = median(L);
save('Modulus.mat','Modulus','r_mean','L_mean');

%%

dhdt = 1; %Set to 1 for ice to grow in thickness over time

nDTOut=10; %Output frequency (in number of time steps)

nSnapshots=100; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20;

tStart = tic; 

%set flags to be used for updating ice ocean interactions
doInt.flag = false; 
doInt.step = 10;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=10; Ny=10;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initiailize Eulearian Data
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);
SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
DSigX = 0; DSigY= 0; DSig1= 0; DSig2= 0;
DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
U = zeros(Ny, Nx); V = zeros(Ny, Nx);
dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
Fx = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
Sig = zeros(Ny, Nx); mass = zeros(Ny,Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);
% error()

%% Initialize time and other stuff to zero
if isempty(dir('Floes_bnds')); disp('Creating folder: Floes_bnds'); mkdir('Floes_bnds'); end
if isempty(dir('./Floes_bnds/figs')); disp('Creating folder: figs'); mkdir('./Floes_bnds/figs'); end
if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end

Time = i_step*dt;

%% Solving for floe trajectories
tic; 
% while side < 2.5
while im_num<nSnapshots

    if mod(i_step,10)==0        
        disp(' ');
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)]);
        disp(' ');
        tic
        doInt.flag=true;
    else
        doInt.flag=false;
    end

    %Simplify the floe shapes if there are too many vertices
    if mod(i_step,nSimp)==0 && SIMPLIFY
        FloeOld = Floe;
        parfor j=1:length(FloeOld)
            FloeOld(j).poly = polyshape(FloeOld(j).c_alpha'+[FloeOld(j).Xi FloeOld(j).Yi]);
        end
        floenew = [];
        for ii = 1:length(Floe)
            ParFloes(ii).floenew = [];
            ParFloes(ii).kill = [];
            ParFloes(ii).verts = 0;
        end
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        parfor ii = 1+Nb:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                [floe2,kill] = FloeSimplify(Floe(ii),0,FloeOld,polyAU);
                if isfield(floe2,'poly')
                    floe2=rmfield(floe2,{'poly'});
                end
                if isempty(kill)
                    ParFloes(ii).kill = [ParFloes(ii).kill kill];
                end
                for jj = 1:length(floe2)
                    if jj == 1
                        Floe(ii) = floe2(jj);
                        ParFloes(ii).verts=length(floe2(jj).c0(1,:));
                    else
                        ParFloes(ii).floenew = [ParFloes(ii).floenew floe2(jj)];
                    end
                end
            end
        end
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        floenew =[]; kill = [];
        clear FloeOld
        for ii = 1:length(ParFloes)
            floenew = [floenew ParFloes(ii).floenew];
            kill = [kill ParFloes(ii).kill];
        end
        kill = unique(kill(kill>Nb));
        live = cat(1,Floe.alive);
        live(kill) = 0;
        Floe(live==0)=[];
        Floe =[Floe floenew];
    end

    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        

        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
        if ifPlot
            [fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);
            exportgraphics(fig,['./Floes_bnds/figs/' num2str(im_num,'%03.f') '.jpg']);

        end
        

        if AVERAGE
            SigXXa = SigXX/fix(nDTOut); SigYXa = SigYX/fix(nDTOut);
            SigXYa = SigXY/fix(nDTOut); SigYYa = SigYY/fix(nDTOut);
            DivSigXa = DivSigX/fix(nDTOut); DivSig1a = DivSig1/fix(nDTOut);
            DivSigYa = DivSigY/fix(nDTOut); DivSig2a = DivSig2/fix(nDTOut);
            Eux = Eux/fix(nDTOut); Evx = Evx/fix(nDTOut);
            Euy = Euy/fix(nDTOut); Evy = Evy/fix(nDTOut);
            U = U/fix(nDTOut); V = V/fix(nDTOut);
            dU = dU/fix(nDTOut); dV = dV/fix(nDTOut);
            Fx = Fx/fix(nDTOut); Fy = Fy/fix(nDTOut);
            Sig = Sig/fix(nDTOut); 
            mass = mass/fix(nDTOut);
            if max(abs(U))>0.5
                xx = 1; xx(1) =[1 2];
            end
        else
            SigXXa = squeeze(eularian_data.stressxx); SigYXa = squeeze(eularian_data.stressyx);
            SigXYa = squeeze(eularian_data.stressxy); SigYYa = squeeze(eularian_data.stressyy);
            DivSigXa = DSigX; DivSig1a = DSig1;
            DivSigYa = DSigY; DivSig2a = DSig2;
            Eux = squeeze(eularian_data.strainux); Evx = squeeze(eularian_data.strainvx);
            Euy = squeeze(eularian_data.strainuy); Evy = squeeze(eularian_data.strainvy);
            U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
            dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv);
            Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_x);
            Sig = Sig+squeeze(eularian_data.stress);
        end
        
                
        if ifPlotStress
            [fig] =plot_basic_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./Stresses/figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
            figure(2)
            plot(Xc,SigXXa,'kx','linewidth',2); title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
            drawnow
        end
        
    end
    
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice_new(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean, height, min_floe_size, PERIODIC,3,3,Nb);
        end
    end
    
    
    if mod(i_step,nDTOut)==0
        save(['./output/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','Nb','Nbond','eularian_data','SigXXa','SigXYa', 'SigYYa','U','V','dU','Fx','mass','c2_boundary');
        SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
        SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
        DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
        DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
        Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
        Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
        U = zeros(Ny, Nx); V = zeros(Ny, Nx);
        dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
        Fx = U; Fy = U;
        Sig = zeros(Ny, Nx); mass = zeros(Ny, Nx);
        
        M = cat(1,Floe.mass);
        Mtot(im_num) = sum(M)+sum(Vdnew(:));
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);
    
    if AVERAGE
        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
        SigXX = SigXX+squeeze(eularian_data.stressxx); SigYX = SigYX+squeeze(eularian_data.stressyx);
        SigXY = SigXY+squeeze(eularian_data.stressxy); SigYY = SigYY+squeeze(eularian_data.stressyy);
        Eux = Eux+squeeze(eularian_data.strainux); Evx = Evx+squeeze(eularian_data.strainvx);
        Euy = Euy+squeeze(eularian_data.strainuy); Evy = Evy+squeeze(eularian_data.strainvy);
        U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
        dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv); 
        Sig = Sig+squeeze(eularian_data.stress);
        mass = mass+squeeze(eularian_data.Mtot);

    end
    
    %Perform welding if selected to at designated rate
    if WELDING && mod(i_step,25)==0
        weldrate = 150;%Set rate at which floes will meld
        A=cat(1,Floe.area);
        if max(A) > Amax
           Amax = max(A);
        end
        %Do welding over different sized grid boxes
        if mod(i_step,5000)==0
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/2,1,1);
        elseif mod(i_step,500)==0
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,2,2);
        else
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,3,3);
        end
    end


    if FRACTURES && mod(i_step,50)==0 %&& im_num > 40
        compactness = sum(cat(1,Floe.area))/area(c2_boundary_poly);
        [Floe,Princ] = FracMohr(Floe,Nb,min_floe_size,compactness);
    end
 
    if CORNERS && mod(i_step,10)==0
        keep = rand(length(Floe),1)>0.7;
        keep(1:Nb) = ones(Nb,1);
        if KEEP_MIN
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
         overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
         keep(overlapArea>0.15) = 0;
        fracturedFloes=corners(Floe(~keep),Nb);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end

    end   

    %Advect the dissolved mass
    Vdnew = Vd(:,:,1)+dissolvedNEW;
    dissolvedNEW=zeros(Ny,Nx);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
    if ~KEEP_MIN
        Area=cat(1,Floe.area);
        dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
        if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
        Floe=Floe(Area> min_floe_size);
    end

    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
    Time=Time+dt; i_step=i_step+1; %update time index
    
end

tEnd = toc(tStart)