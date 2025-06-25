close all
FloeNumbers = 1:length(Floe);
Floes = Floe(Nbound+1:end);
Nums = cat(1,Floe.num);
clear f1; clear f2
count = 1;
for ii = 1:length(Floes)
    BondNum = unique(cat(1,Floes(ii).bonds.Num));
    for jj = 1:length(BondNum)
        f1(count) = Floes(ii).num;
        f2(count) = BondNum(jj);
        count = count+1;
    end
end
f2(f1 ==0) = [];f1(f1 ==0) = [];
G = digraph(f1,f2);
plot(G,'Layout','layered')
[bins,binsizes] = conncomp(G,'Type','weak');
Lia = ismember(FloeNumbers,Nums);Lia(1:Nbond) = 0;
clear poly
for iii = 1:length(Floe)
    poly(iii) = polyshape(Floe(iii).c_alpha'+[Floe(iii).Xi Floe(iii).Yi]);
end
count = 1;
clear p
for ii = 1+Nbond:length(binsizes)
    if Lia(ii)
        Lia1 = ismember(bins,ii);
        Lia2 = ismember(Nums,FloeNumbers(Lia1));
        ptmp = union([poly(Lia2)]);
        p(count) = ptmp;%rmholes(ptmp);
        count = count+1;  
    end
end

close all
figure
hold on
colors = distinguishable_colors(length(p));
for ii = 1:length(p)
    plot(p(ii),'FaceColor',colors(ii,:) )
end

Abox = area(c2_boundary_poly);
a = sqrt(area(p));
hFSD = histogram(a/1e3);
edges = hFSD.BinEdges;
val = hFSD.Values;
clear FSD
count3 = 1;
%Loop through to find floe sizes per km^2
for ii = 1:length(val)
    FSD(count3) = sum(val(ii:end))/Abox*1e6;
    count3 = count3+1;
end
%Plot FSD
bins = (edges(2:end)+edges(1:end-1))/2;
bins(FSD<1e-4) = []; FSD(FSD<1e-4) = [];
figure;
fig(1) = loglog(bins,FSD,'linewidth',2);
min_size = 1;
binsUpper = bins(bins>min_size); slopes1 = -2;
hold on
fig(2) = loglog([binsUpper(1) 50], 5*1e-1*[binsUpper(1) 50].^(slopes1),'k--','linewidth',2);
set(gca, 'YScale', 'log')
set(gca, 'xScale', 'log')
xlabel('Floe Size (km)','fontsize',20,'interpreter','latex')
ylabel({'FSD (floes per km$^2$)'},'fontsize',20,'interpreter','latex')
fig = figure(2);

%% 
for kk = 1:427
    load(['./Floes_bnds1/Floe' num2str(kk,'%07.f') '.mat'])
    c2_boundary_poly = polyshape(c2_boundary');
    [fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);
    exportgraphics(fig,['./figs/' num2str(kk,'%03.f') '.jpg'] ,'resolution',300);
end
%% 
fig = 0; num_chunks = zeros(1,427);Nbound = 0;edges = 0:2:102;
for kk = 138:423
    load(['./Floes_bnds2/Floe' num2str(kk,'%07.f') '.mat'])
    c2_boundary_poly = polyshape(c2_boundary'); fig = 0;
    %[fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);
    close all
    FloeNumbers = 1:length(Floe);
    Floes = Floe(Nbound+1:end);
    Nums = cat(1,Floe.num);
    clear f1; clear f2
    Xc = []; Xb = [];
    Yc = []; Yb = [];
    Xi = cat(1,Floe.Xi); Yi = cat(1,Floe.Yi);
    count = 1;
    for ii = 1:length(Floes)
        BondNum = [FloeNumbers(ii); unique(cat(1,Floes(ii).bonds.Num))];
        Xb = [Xb; cat(1,Floe(ii).bonds.Xb)+Xi(ii)];
        Xc = [Xc; cat(1,Floe(ii).bonds.Xc)+Xi(ii)];
        Yb = [Yb; cat(1,Floe(ii).bonds.Yb)+Yi(ii)];
        Yc = [Yc; cat(1,Floe(ii).bonds.Yc)+Yi(ii)];
        for jj = 1:length(BondNum)
            f1(count) = Floes(ii).num;
            f2(count) = BondNum(jj);
            count = count+1;
        end
    end
    f2(f1 ==0) = [];f1(f1 ==0) = [];
    G = digraph(f1,f2);
    plot(G,'Layout','layered')
    [bins,binsizes] = conncomp(G,'Type','weak');
    Lia = ismember(FloeNumbers,Nums);Lia(1:Nbond) = 0;
    clear poly
    for iii = 1:length(Floe)
        poly(iii) = polyshape(Floe(iii).c_alpha'+[Floe(iii).Xi Floe(iii).Yi]);
    end
    count = 1;
    clear p
    for ii = 1+Nbond:length(binsizes)
        if Lia(ii)
            Lia1 = ismember(bins,ii);
            Lia2 = ismember(Nums,FloeNumbers(Lia1));
            ptmp = union([poly(Lia2)]);
            p(count) = ptmp;%rmholes(ptmp);
            count = count+1;
        end
    end
    num_chunks(kk) = length(p);

    close all
    ratio=2.1;
    fig=figure('Position',[100 100 1000*ratio 1000],'visible','on');
    set(fig,'PaperSize',12*[ratio 1],'PaperPosition',12*[0 0 ratio 1]);
    figure(fig)
    clf(fig);
    subplot(1,2,1)
    plot(poly(1+Nbound:length(Floe)),'FaceColor','none','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
    hold on
    plot(Xc,Yc,'r.','linewidth',3);
    colors = distinguishable_colors(length(p));
    for ii = 1:length(p)
        plot(p(ii),'FaceColor',colors(ii,:) )
    end

    set(0,'CurrentFigure',fig);
    xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
    yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
    plot(xb,yb, 'k-','linewidth',2);

    colormap('gray'); caxis([0 1]);
    axis equal
    axis([-65000 65000 -65000 65000])
    set(gca,'Ydir','normal');
    box on

    subplot(1,2,2)
    plot(p,'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
    hold on
    plot(xb,yb, 'k-','linewidth',2);
    colormap('gray'); caxis([0 1]);
    axis equal
    axis([-65000 65000 -65000 65000])
    set(gca,'Ydir','normal');
    box on
    fig = figure(1);
    exportgraphics(fig,['./figs2/' num2str(kk,'%03.f') '.jpg'] ,'resolution',300);

    close all
    Abox = area(c2_boundary_poly);
    a = sqrt(area(p));
    hFSD = histogram(a/1e3,edges);
%     edges = hFSD.BinEdges;
    val = hFSD.Values;
    clear FSD
    count3 = 1;
    %Loop through to find floe sizes per km^2
    for ii = 1:length(val)
        FSD(count3) = sum(val(ii:end))/Abox*1e6;
        count3 = count3+1;
    end
    %Plot FSD
    bins = (edges(2:end)+edges(1:end-1))/2;
    bins(FSD<1e-6) = []; FSD(FSD<1e-6) = [];
    figure;
    fig(1) = loglog(bins,FSD,'linewidth',2);
    min_size = 1;
    binsUpper = bins(bins>min_size); slopes1 = -2;
    hold on
    fig(2) = loglog([binsUpper(1) 50], 5*1e-1*[binsUpper(1) 50].^(slopes1),'k--','linewidth',2);
    set(gca, 'YScale', 'log')
    set(gca, 'xScale', 'log')
    xlabel('Floe Size (km)','fontsize',20,'interpreter','latex')
    ylabel({'FSD (floes per km$^2$)'},'fontsize',20,'interpreter','latex')
    fig = figure(2);
    saveas(fig,['./figs2/FSD' num2str(kk,'%03.f') '.jpg'],'jpg');
end
plot(num_chunks)

%% 
network=zeros(length(Floe));
nums = cat(1,Floe.num);
bincolors = zeros(length(Floe),3);
count = 1;
for ii = 1:length(Floe)
    bnds = cat(1,Floe(ii).bonds.Num);
    Lia = ismember(nums,bnds);
    network(Lia,ii) = bins(Lia);
end
count = count + 1;
count2 = 1;
for kk = 1:length(binsizes)
    for ii = 1:binsizes(kk)
        bincolors(count2,:) = colors(kk,:);
        count2 = count2+1;
    end
end
%circularGraph(network,'Colormap',bincolors)
