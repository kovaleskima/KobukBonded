function [Floe,Princ] = FracMohr(Floe,Nb,min_floe_size,concentration)
%Use Mohr's cone to determine which floes are fractured
%   If Principal stresses are outside the cone then the floes are fractured
rho_ice=920;

        Pstar = 2.25e5; C = 20;
        h = mean(cat(1,Floe.h));
        P = Pstar*h*exp(-C*(1-concentration));
        t = linspace(0,2*pi) ;
        a = P*sqrt(2)/2 ; b = a/2 ;
        x = [0 -1.5e5 -5.5e5];%a*cos(t) ;
        y = [0 -5.5e5 -1.5e5];%b*sin(t) ;
        Mohr = polyshape(x,y);
       Mohr = rotate(Mohr,45);
       Mohr = translate(Mohr,[-P/2, -P/2]);
        A = cat(1,Floe.area);

q = 5.2; SigC = 250e3;
Sig1 = (1/q+1)*SigC/(1/q-q);
Sig2 = q*Sig1+SigC;
Sig11 = -3e4;
Sig22 = q*Sig11+SigC;
MohrX = [Sig1; Sig11; Sig22];
MohrY = [Sig2; Sig22; Sig11];
Mohr = polyshape(-MohrX,-MohrY);
        p1 = -59520; p2 = 4e8;
        for ii = 1:length(Floe) %for every floe, get principal stresses
            Stress = eig(Floe(ii).Stress);
            Princ1(ii) = max(Stress); Princ2(ii) = min(Stress); 
            Princ(ii,1) = max(Stress);
            Princ(ii,2) = min(Stress);
            Princ(ii,3) = Floe(ii).area;
            Stresses(ii) = vecnorm([Princ1 Princ2]);
        end
        Princ1 = Princ(:,1); Princ2 = Princ(:,2); 
        [in,out] = inpolygon(Princ1,Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2)); % determine if the principal stresses are in the Mohr Cone

count = 1; Sig1 = []; Sig2 = [];
FloeNums = cat(1,Floe.num);
for ii = 1+Nb:length(Floe) %from 1+# of bonds to the number of floes (??)

    if ~isempty(Floe(ii).bonds) %if there are no bonds on a floe
        bnds = unique(cat(1,Floe(ii).bonds.Num));

        for jj = length(bnds):-1:1
            Lia = ismember(cat(1,Floe(ii).bonds.Num),bnds(jj));
            bnd_tmp = Floe(ii).bonds(Lia); %add a temporary bond to a floe?
            Sig1(count) = 0; Sig2(count) = 0;
            for kk = 1:length(bnd_tmp) %for every temporary bond added
                Sig1(count) = Sig1(count)+bnd_tmp(kk).Stress(1)/50; Sig2(count) = Sig2(count)+bnd_tmp(kk).Stress(2)/50; %add temporary bond stress to principal stresses
                bnd_tmp(kk).Stress = [0; 0]; %clear temporary bond stress
            end
            Sig1(count) = Sig1(count)/kk; Sig2(count) = Sig2(count)/kk; %something physics related, check sigma matrix
            if abs(Sig1(count))>0.75e4

                Floe(ii).bonds(Lia) = [];
                Lia = ismember(FloeNums,bnds(jj)); 
                BndNums = cat(1,Floe(Lia).bonds.Num); 
                Lia2 = ismember(BndNums,Floe(ii).num);
                Floe(Lia).bonds(Lia2) = []; %if the floe is bonded and the bonds exist, clear it (break it)
            elseif abs(Sig2(count))>0.75e4 
                Floe(ii).bonds(Lia) = [];
                Lia = ismember(FloeNums,bnds(jj));
                BndNums = cat(1,Floe(Lia).bonds.Num);
                Lia2 = ismember(BndNums,Floe(ii).num);
                Floe(Lia).bonds(Lia2) = []; %do the same if the other principal stress threshold is exceeded.
            else
                Floe(ii).bonds(Lia) = bnd_tmp; %otherwise keep the temporary bond value
            end
            count = count+1;
        end
    end
end
        
        p = rand(length(Floe),1);
        StressNorm = Stresses/max(Stresses);%rmoutliers(Stresses));

        %dead code here?
        keep = zeros(length(Floe),1);
        keep(in) = 1; %keep floes in mohr cone
        keep(A<min_floe_size)=1;
        keep(1:Nb) = ones(Nb,1);
        
        
        %just to only fracture bonds
        keep = ones(length(Floe),1); %(Nbond+1:length(Floe)) = ones(length(Floe)-Nbond,1);
        
        
        keep = logical(keep);
        for ii = 1:length(Floe)
            FracFloes(ii).floenew = [];
        end

        parfor ii = 1:length(keep)
            if ~keep(ii)
                FracFloes(ii).floenew=fracture_floe(Floe(ii),3,Floe,1);
            end
        end
        fracturedFloes =[];
        for ii = 1:length(FracFloes)
            fracturedFloes = [fracturedFloes FracFloes(ii).floenew];
        end
        if isfield(fracturedFloes,'potentialInteractions')
            fracturedFloes=rmfield(fracturedFloes,{'potentialInteractions'});
        end
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
end

