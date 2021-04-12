%%% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 1 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/fschneider/Documents/GitHub/FigureGround_Ephys_Analysis

% Initialise figure
f                   = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);
ax0                 = axes('Position',[0 0 1 1],'Visible','off');
col                 = [0 .9 0; .9 0 0; 0 .9 0; .9 0 0];
row                 = linspace(.05,.84,5);
clm                 = linspace(.1,.76,4);
dim                 = [.15 .2];

%%% Plot SFG example stimulus %%%
load('PATH/DataStruct_2019-05-20.mat')

freqMat             = data.stimSpecs.freq_mat;
axA                 = axes('Position',[row(1) clm(4) dim]); hold on
no                  = 4;
idx                 = data.behaviour.stimNrPool == data.behaviour.stimID(no);
figOn               = unique(data.behaviour.figOn(idx));
figElem             = data.stimSpecs.fig{no};
plotSFGstim(freqMat, no, figOn, figElem);
axA.XAxis.Visible   = 'off';
tx1                 = text(31,10000, 'Figure', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
tx2                 = text(5,10000, 'Ground', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
axA.Position(1)     = row(1);
axA.Position(3:4)   = [.15 .2];
axA.XAxis.Visible   = 'off';

axA                 = axes('Position',[row(2) clm(4) dim]); hold on
idx                 = data.behaviour.stimNrPool == data.behaviour.stimID(16);
figOn               = unique(data.behaviour.figOn(idx));
plotSFGstim(freqMat, 16, figOn, figElem);
axA.Position(1)     = row(2);
axA.Position(3:4)   = [.15 .2];
axA.XAxis.Visible   = 'off';
axA.YAxis.Visible   = 'off';

%%% Plot paradigm %%%
axP                 = axes('Position',[row(1) clm(3)+.125 dim]); hold on
axP.FontSize        = 14;
axP.XLim            = [250 3750];
axP.XTick           = [500 2000 3500];
axP.XTickLabel      = {'0','1500', '3000'};
axP.XLabel.String   = 'Time [ms]';
axP.YLim            = [0 6];
axP.YTick           = [0 2 4];
axP.YTickLabel      = {'Sound','Figure', 'Response'};
axP.YAxis.FontSize  = 12;
p1                  = plot([zeros(1,500),ones(1,3000),zeros(1,500)], 'LineWidth', 2, 'Color','k');
p2                  = plot([zeros(1,1900)+2,zeros(1,1000)+3,zeros(1,1200)+2], 'LineWidth', 2, 'Color','k');
p3                  = plot([zeros(1,2000)+4,zeros(1,900)+5,zeros(1,1200)+4], 'LineWidth', 2, 'Color','k');
tx1                 = text(300,5.5, 'Test trial', 'FontSize', 12, 'Color', 'k');
tx2                 = text(2325,4.5, 'HI', 'FontSize', 10, 'Color', [.3 .3 .3]);
tx3                 = text(1300,4.5, 'MI', 'FontSize', 10, 'Color', [.3 .3 .3]);
tx4                 = text(3100,4.5, 'MI', 'FontSize', 10, 'Color', [.3 .3 .3]);
axP.Position(1)     = row(1);
axP.Position(3:4)   = [.15 .1];

axP                 = axes('Position',[row(2) clm(3)+.125 dim]); hold on
axP.FontSize        = 14;
axP.XLim            = [250 3750];
axP.XTick           = [500 2000 3500];
axP.XTickLabel      = {'0','1500', '3000'};
axP.XLabel.String   = 'Time [ms]';
axP.YLim            = [0 6];
axP.YAxis.Visible   = 'off';
p1                  = plot([zeros(1,500),ones(1,3000),zeros(1,500)], 'LineWidth', 2, 'Color','k');
p2                  = plot([zeros(1,4000)+2], 'LineWidth', 2, 'Color','k');
p3                  = plot([zeros(1,3500)+4,zeros(1,500)+5], 'LineWidth', 2, 'Color','k');
tx1                 = text(300,5.5, 'Control trial', 'FontSize', 12, 'Color', 'k');
tx2                 = text(3600,4.5, 'CR', 'FontSize', 10, 'Color', [.3 .3 .3]);
tx3                 = text(2100,4.5, 'FA', 'FontSize', 10, 'Color', [.3 .3 .3]);
axP.Position(1)     = row(2);
axP.Position(3:4)   = [.15 .1];


%%% Plot behaviour %%%
bdir                = 'PATH/';
load([bdir 'M1/dprime.mat'])
dpE                 = dp;
load([bdir 'M2/dprime.mat'])
dpD                 = dp;
load([bdir 'M1/mRT.mat'])
mRTE                = mRT;
load([bdir 'M2/mRT.mat'])
mRTD                = mRT;
load([bdir 'M1/RTsd.mat'])
sdRTE               = RTsd;
load([bdir 'M2/RTsd.mat'])
sdRTD               = RTsd;

of                  = 0.009;
bwdt                = 30;
lw                  = 2;
arrE                = [ones(size(dpE,1),1); ones(size(dpE,1),1)+1;];
arrD                = [ones(size(dpD,1)-10,1)+2; ones(size(dpD,1)-10,1)+3];

mat                 = vertcat(dpE(:,2), dpE(:,3), dpD(1:end-10,2), dpD(1:end-10,3));
arr                 = vertcat(arrE, arrD);

axB                 = axes('Position',[row(3)+of clm(3)+.125 dim]); hold on
ln                  = line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k');
bx                  = boxplot(mat, arr, 'Colors', 'k');
set(bx(end,:),'Visible','off')
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw})
axB.YLabel.String   = 'd prime';
axB.YLim            = [-.2 5.2];
axB.XTick           = [1 2 3 4];
axB.YTick           = [0 1 2 3 4 5];
axB.XTickLabel      = {'Coh8' 'Coh12' 'Coh8' 'Coh12'};
axB.XColor          = [0 0 0];
axB.YColor          = [0 0 0];
axB.FontSize        = 14;
axB.XTickLabelRotation = 30;
tx1                 = text(1.3,5.2, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
tx2                 = text(3.3,5.2, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
box off

x                   = [.8 1.8 2.8 3.8];
for i = 1:4
    xx              = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc              = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    
    if i == 1 || i == 3
        sc.MarkerFaceColor  = [0 1 0];
        sc.MarkerEdgeColor  = [0 1 0];
    else
        sc.MarkerFaceColor  = [1 0 0];
        sc.MarkerEdgeColor  = [1 0 0];
    end 
    sc.SizeData     = 20;
end
axB.Position(1)     = row(3)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1          = dpE(:,2); d2 = dpE(:,3);
        var         = 1;
    else
        d1          = dpD(1:end-10,2); d2 = dpD(1:end-10,3);
        var      	= 3;
    end
    
    pDP(ii)     	= signrank(d1,d2);
    if signrank(d1,d2) < .001
        star        = plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star        = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star        = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%% RT %%%%
arrE = []; arrD = []; mat = []; arr = [];
arrE                = [ones(size(mRTE,1),1); ones(size(mRTE,1),1)+1;];
arrD                = [ones(size(mRTD,1)-10,1)+2; ones(size(mRTD,1)-10,1)+3];
mat                 = vertcat(mRTE(:,2), mRTE(:,3), mRTD(1:end-10,2), mRTD(1:end-10,3));
arr                 = vertcat(arrE,[], arrD);

axB                 = axes('Position',[row(4)+of clm(3)+.125 dim]); hold on
ln                  = line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k');
bx                  = boxplot(mat, arr, 'Colors', 'k');
set(bx(end,:),'Visible','off')
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw})
axB.YLabel.String   = 'Reaction time [ms]';
axB.XTick           = [1 2 3 4];
axB.XTickLabel      = {'Coh8' 'Coh12' 'Coh8' 'Coh12'};
axB.YLim        	= [400 850];
axB.YTick           = [400 600 800];
axB.XColor          = [0 0 0];
axB.YColor          = [0 0 0];
axB.FontSize        = 14;
axB.XTickLabelRotation = 30;
tx1                 = text(1.3,850, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
tx2                 = text(3.3,850, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
box off

x                   = [.8 1.8 2.8 3.8];
for i = 1:4
    xx              = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc              = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    if i == 1 || i == 3
        sc.MarkerFaceColor  = [0 1 0];
        sc.MarkerEdgeColor  = [0 1 0];
    else
        sc.MarkerFaceColor  = [1 0 0];
        sc.MarkerEdgeColor  = [1 0 0];
    end
    sc.SizeData     = 20;
end

axB.Position(1)     = row(4)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1          = mRTE(:,2); d2 = mRTE(:,3);
        var       	= 1;
    else
        d1          = mRTD(1:end-10,2); d2 = mRTD(1:end-10,3);
        var         = 3;
    end
    
    pRT(ii)         = signrank(d1,d2);
    if signrank(d1,d2) < .001
        star        = plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star        = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star        = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%% CV %%%
CVE                 = sdRTE./mRTE;
CVD                 = sdRTD./mRTD;

arrE = []; arrD = []; mat = []; arr = [];
arrE                = [ones(size(CVE,1),1); ones(size(CVE,1),1)+1;];
arrD                = [ones(size(CVD,1)-10,1)+2; ones(size(CVD,1)-10,1)+3];

mat                 = vertcat(CVE(:,2), CVE(:,3),CVD(1:end-10,2), CVD(1:end-10,3));
arr                 = vertcat(arrE, arrD);
    
axB                 = axes('Position',[row(5)+of clm(3)+.125 dim]); hold on
ln                  = line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k');
bx                  = boxplot(mat, arr, 'Colors', 'k');
set(bx(end,:),'Visible','off')
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw})
axB.YLabel.String   = 'Coefficient of variation';
axB.YLim            = [0.1 .35];
axB.YTick           = [.1 .2 .3];
axB.XTick           = [1 2 3 4];
axB.XTickLabel      = {'Coh8' 'Coh12' 'Coh8' 'Coh12'};
axB.XColor          = [0 0 0];
axB.YColor          = [0 0 0];
axB.FontSize        = 14;
axB.XTickLabelRotation = 30;
tx1                 = text(1.3,.35, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
tx2                 = text(3.3,.35, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
box off

x                   = [.8 1.8 2.8 3.8];
for i = 1:4
    xx              = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc              = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    if i == 1 || i == 3
        sc.MarkerFaceColor  = [0 1 0];
        sc.MarkerEdgeColor  = [0 1 0];
    else
        sc.MarkerFaceColor  = [1 0 0];
        sc.MarkerEdgeColor  = [1 0 0];
    end
    sc.SizeData     = 20;
end

axB.Position(1)     = row(5)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1          = CVE(:,2); d2 = CVE(:,3);
        var         = 1;
    else
        d1          = CVD(1:end-10,2); d2 = CVD(1:end-10,3);
        var         = 3;
    end
    
    pCV(ii)         = signrank(d1,d2);
    if signrank(d1,d2) < .001
        star     	= plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star        = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star        = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BF + Latency + Phase locking maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ylim                = 17;
cmap                = flipud(jet(14));
lcmap               = [linspace(0,1,50)', linspace(0,1,50)', ones(50,1)];
alp                 = .9;
back                = [0 0 0];
typ                 = 'muae';

freqStart               = 180;                                              % Tuning low freq [Hz]
steps                   = 14;                                               % No of desired tones
PT(1)                   = freqStart;                                        % Starting frequency [Hz]
for i = 2:steps                                                           	% For no of tones...
    PT(i)               =  PT(i-1)*2^(1/2);                              	% 1/2 octave steps
end
frex                    = round(PT);

for iAn = 1:2
    
    if iAn == 1
        animalID    = 'M1';
    else
        animalID    = 'M2';
    end

    load(['PATH/tMap_' animalID '_' typ '.mat']);
    load(['PATH/lMap_' animalID '_' typ  '.mat']);
    load(['PATH/ccoord_' animalID '_' typ  '.mat']);
    
    AP            	= find(logical(sum(~isnan(mfr_mat),2)));              	% Boundaries of recording field
    ML              = find(logical(sum(~isnan(mfr_mat))));
    [xx,yy]         = coreBoundary(mfr_mat,AP,ML,false,animalID);         	% Get X & Y coordinates for field boundary
        
    x               = [];
    y               = [];
    
    for r = 1:length(xx)
        x           = horzcat(x, xx(r)-.5, xx(r)+.5);
    end
    
    for r = 1:length(yy)
        y           = horzcat(y, yy(r), yy(r));
    end
    
    switch animalID
        case 'M1'
            xax  	= [5 10 15];
            yax   	= [-17 -7];
            ytick   = {'+5' '+15'};
        case 'M2'
            xax    	= [5 10 15];
            yax   	= [-15 -6];
            ytick   = {'+4', '+14'};
    end
    
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iAn == 1
        
        %%% Tonotopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axT1        = axes('Position',[row(3)+of clm(2) dim]); hold on; axis equal
        im          = imagesc(1:size(mfr_mat,1),-18:-1, flipud(mfr_mat));
        p1          = plot(x,-y, 'Color', [0 0 0],'LineWidth', 4);
        axT1.YLim   = [-18 -5];
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        %%% Latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axL1        = axes('Position',[row(4)+of clm(2) dim]); hold on; axis equal
        im          = imagesc(1:size(mlat_mat,1),-18:-1, flipud(mlat_mat));
        hold on
        p1          = plot(x,-y, 'Color', [0 0 0],'LineWidth', 4);
        axL1.YLim   = [-18 -5];
        axL1.YAxis.Visible = 'off';
        caxis([15 80])
        
    else
        
        %%% Tonotopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axT1        = axes('Position',[row(3)+of clm(1) dim]); hold on; axis equal
        im          = imagesc([1:size(mfr_mat,1)], [-18:-1], flipud(mfr_mat));
        hold on
        p1          = plot(x,-y, 'Color', [0 0 0],'LineWidth', 4);
        axT1.YLim   = [-17 -4];
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        %%% Latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axL1        = axes('Position',[row(4)+of clm(1) dim]); hold on; axis equal
        im          = imagesc(1:size(mlat_mat,1),-18:-1, flipud(mlat_mat));
        hold on
        p1          = plot(x,-y, 'Color', [0 0 0],'LineWidth', 4);
        axL1.YLim   = [-17 -4];
        axL1.YAxis.Visible = 'off';
        caxis([15 80])
        
        %%% Phase-lock %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Load MUA data [M2]        
        %         co = [];
        %         for ii = 1:size(muaeD,2)
        %             if ~isempty(muaeD{ii}.phLock)
        %                 pL(ii,:) = muaeD{ii}.phLock;
        %                 co(ii,:) = muaeD{ii}.coord;
        %             end
        %         end
        %         save([dest_dir 'raw\PL_id_' animalID '_' typ  '.mat'], 'pL');
        %         save([dest_dir 'raw\PL_co_' animalID '_' typ  '.mat'], 'co');

        load('PATH/PL_co_M2_muae.mat')
        load('PATH/PL_id_M2_muae.mat')
        
        
        c                   = 1;
        for i = 1:size(co,1)
            MLr             = co(i,2);                              % ML penetration site on grid [mm]
            dep             = co(i,3);                              % Depth of recording from GT tip [mm]
            offset          = (dep/sind(90)) * sind(15);          	% Calculate offset [mm]
            adj             = MLr - offset;                      	% Adjusted ML value [mm]
            AP(c)           = co(i,1);
            ML(c)           = adj;
            c               = c+1;
        end
        
        axPL                = axes('Position',[row(5)+of clm(1) dim]); hold on; axis equal
        im                  = imagesc([1:size(mfr_mat,1)], [-18:-1], flipud(mfr_mat));
        p1                  = plot(x,-y, 'Color', [0 0 0],'LineWidth', 4);
        axPL.YLim           = [-17 -4];
        
        scatter(ML(sum(pL,2) == 4),-AP(sum(pL,2) == 4), '^r', 'filled', 'filled', 'MarkerFaceAlpha', 1)
        scatter(ML(sum(pL,2) == 3),-AP(sum(pL,2) == 3), '^r', 'filled', 'filled', 'MarkerFaceAlpha', .6)
        scatter(ML(sum(pL,2) == 2),-AP(sum(pL,2) == 2), '^r', 'filled', 'filled', 'MarkerFaceAlpha', .4)
        % scatter(ML(sum(pL,2) == 1),AP(sum(pL,2) == 1), '^r', 'filled', 'filled', 'MarkerFaceAlpha', .1)
        
        fill([12.5 13 13.5], [-8.5 -7.5 -8.5], 'r', 'FaceAlpha', .4, 'EdgeColor', 'none')
        fill([12.5 13 13.5], [-7.5 -6.5 -7.5], 'r', 'FaceAlpha', .6, 'EdgeColor', 'none')
        fill([12.5 13 13.5], [-6.5 -5.5 -6.5], 'r', 'FaceAlpha', 1, 'EdgeColor', 'none')
        
        text(14,-8, '2/4','FontSize', 12, 'Color', 'r', 'FontWeight', 'bold')
        text(14,-7, '3/4','FontSize', 12, 'Color', 'r', 'FontWeight', 'bold')
        text(14,-6, '4/4','FontSize', 12, 'Color', 'r', 'FontWeight', 'bold')
    end
    
    axT1.YTick              = yax;
    axT1.YTickLabel         = ytick;
    axT1.FontSize           = 14;
    cm                      = [back; flipud(jet(256))];
    colormap(axT1, cm)
    
    if iAn == 1
        axT1.YLabel.String  = {'M1';'Distance IAL [mm]'};
        axT1.XAxis.Visible  = 'off';
        axL1.XAxis.Visible  = 'off';
    else
        axT1.YLabel.String  = {'M2';'Distance IAL [mm]'};
        axT1.XTick          = xax;
        axL1.XTick          = xax;
    end
    
    axL1.YTick              = yax;
    axL1.YTickLabel         = ytick;
    axL1.FontSize           = 14;
    axL1.YAxis.Visible      = 'off';
    cm                      = [back; lcmap];
    colormap(axL1, cm)

    if iAn == 2
        axPL.YTick          = yax;
        axPL.YTickLabel     = ytick;
        axPL.YTick          = yax;
        axPL.XTick          = xax;
        axPL.FontSize       = 14;
        axPL.XLabel.String  = 'Grid position ML [mm]';
        axPL.YAxis.Visible  = 'off';
        cm                  = [back; gray(256)];
        colormap(axPL, cm)
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        axT1.XLabel.String  = 'Grid position ML [mm]';
        axL1.XLabel.String  = 'Grid position ML [mm]';
        
    end
end

cm                  = [back; flipud(jet(256))];
colormap(axT1, cm)
axCB                = axes('Position',[row(5)+.025 clm(2)+.05 .001 .1]);
axCB.Visible        = 'off';
colormap(axCB, cm)
cb                  = colorbar(axCB);
cb.Color            = [0 0 0];
cb.Position(3)    	= .01;
cb.Label.String   	= 'Best frequency [Hz]';
cb.FontSize         = 12;
caxis([floor(log(frex(1))) ceil(log(frex(11)))])
cb.Ticks            = [log(frex(1)),log(frex(7)),log(frex(11))];
cb.TickLabels       = {num2str(frex(1)), num2str(frex(7)), ['>' num2str(frex(11))]};

cm                  = [back; gray(256)];
axCB                = axes('Position',[row(5)+.01 clm(2)+.05 .001 .1]);
axCB.Visible        = 'off';
colormap(axCB, cm)
cb                  = colorbar(axCB);
cb.Color            = [0 0 0];
cb.Position(3)      = .01;
cb.Ticks            = [];

axCB                = axes('Position',[row(5)+.1 clm(2)+.05 .001 .1]);
axCB.Visible        = 'off';
colormap(axCB, lcmap)
cb                  = colorbar(axCB);
cb.Color            = [0 0 0];
cb.Position(3)      = .01;
cb.Label.String     = 'Latency [ms]';
cb.FontSize         = 12;
caxis([15 80])
cb.Ticks            = [15 48 80];
cb.TickLabels       = {'15' '48' '>80'};


offset              = 0.03;
text(0,.98, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(2)-offset,.98, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(3)-offset,.98, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(4)-offset,.98, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(5)-offset,.98, 'e', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0,.54, 'f', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(3)-offset,.54, 'g', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(4)-offset,.54, 'h', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(row(5)-offset,.32, 'i', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

text(row(3)+.06,.54, 'Tonotopy', 'Parent', ax0, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
text(row(4)+.05,.54, 'Peak Latency', 'Parent', ax0, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
text(row(5)+.06,.32, 'Phase Locking', 'Parent', ax0, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')

% addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
% dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
% export_fig([dest_dir 'FIG1'], '-r400',f);
% 
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f, [dest_dir 'FIG1'], '-dpdf', '-r400'); 