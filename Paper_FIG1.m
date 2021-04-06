%%% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 1 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Volumes/Felix_ExtDrive/Felix/Scripts/Stuff
clearvars -except muaeD muaeE data lfpD lfpE
f       = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);
ax0  	= axes('Position',[0 0 1 1],'Visible','off');
col     = [0 .9 0; .9 0 0; 0 .9 0; .9 0 0];

row     = linspace(.05,.84,5);
clm     = linspace(.1,.76,4);
dim     = [.15 .2];

%%% SFG stim %%%
% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Dollar\2019-05-20_12-01-40\Data\DataStruct_2019-05-20.mat')
load('/Volumes/Felix_ExtDrive/Rec/Dollar/2019-05-20_12-01-40/Data/DataStruct_2019-05-20.mat')

freqMat             = data.stimSpecs.freq_mat;
axA                 = axes('Position',[row(1) clm(4) dim]); hold on
no                  = 4;
idx                 = data.behaviour.stimNrPool == data.behaviour.stimID(no);
figOn               = unique(data.behaviour.figOn(idx));
figElem             = data.stimSpecs.fig{no};
plotSFGstim(freqMat, no, figOn, figElem)
axA.XAxis.Visible   = 'off';
text(31,10000, 'Figure', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold')
text(5,10000, 'Ground', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
axA.Position(1)     = row(1);
axA.Position(3:4)   = [.15 .2];
axA.XAxis.Visible   = 'off';

axA                 = axes('Position',[row(2) clm(4) dim]); hold on
idx                 = data.behaviour.stimNrPool == data.behaviour.stimID(16);
figOn               = unique(data.behaviour.figOn(idx));
plotSFGstim(freqMat, 16, figOn, figElem)
axA.Position(1)     = row(2);
axA.Position(3:4)   = [.15 .2];
axA.XAxis.Visible   = 'off';
axA.YAxis.Visible   = 'off';

%%% Paradigm %%%
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
plot([zeros(1,500),ones(1,3000),zeros(1,500)], 'LineWidth', 2, 'Color','k')
plot([zeros(1,1900)+2,zeros(1,1000)+3,zeros(1,1200)+2], 'LineWidth', 2, 'Color','k')
plot([zeros(1,2000)+4,zeros(1,900)+5,zeros(1,1200)+4], 'LineWidth', 2, 'Color','k')
text(300,5.5, 'Test trial', 'FontSize', 12, 'Color', 'k')
text(2325,4.5, 'HI', 'FontSize', 10, 'Color', [.3 .3 .3])
text(1300,4.5, 'MI', 'FontSize', 10, 'Color', [.3 .3 .3])
text(3100,4.5, 'MI', 'FontSize', 10, 'Color', [.3 .3 .3])
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
plot([zeros(1,500),ones(1,3000),zeros(1,500)], 'LineWidth', 2, 'Color','k')
plot([zeros(1,4000)+2], 'LineWidth', 2, 'Color','k')
plot([zeros(1,3500)+4,zeros(1,500)+5], 'LineWidth', 2, 'Color','k')
text(300,5.5, 'Control trial', 'FontSize', 12, 'Color', 'k')
text(3600,4.5, 'CR', 'FontSize', 10, 'Color', [.3 .3 .3])
text(2100,4.5, 'FA', 'FontSize', 10, 'Color', [.3 .3 .3])
axP.Position(1)     = row(2);
axP.Position(3:4)   = [.15 .1];

% bdir = 'Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Eric\Summary\';
bdir = '/Volumes/Felix_ExtDrive/Rec/';

%%% Behaviour %%%
load([bdir 'Eric/Summary/dprime.mat'])
dpE     = dp;
load([bdir 'Dollar/Summary/dprime.mat'])
dpD     = dp;
load([bdir 'Eric/Summary/mRT.mat'])
mRTE    = mRT;
load([bdir 'Dollar/Summary/mRT.mat'])
mRTD    = mRT;
load([bdir 'Eric/Summary/RTsd.mat'])
sdRTE   = RTsd;
load([bdir 'Dollar/Summary/RTsd.mat'])
sdRTD   = RTsd;

of      = 0.009;
bwdt    = 30;
lw      = 2;
arrE    = [ones(size(dpE,1),1); ones(size(dpE,1),1)+1;];
arrD    = [ones(size(dpD,1)-10,1)+2; ones(size(dpD,1)-10,1)+3];

% disp(mean([dpE(:,2), dpE(:,3)]))
% disp(nanmean([dpD(1:end-10,2), dpD(1:end-10,3)]))

mat     = vertcat(dpE(:,2), dpE(:,3), dpD(1:end-10,2), dpD(1:end-10,3));
arr     = vertcat(arrE, arrD);

axB     = axes('Position',[row(3)+of clm(3)+.125 dim]); hold on
line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k')
bx      = boxplot(mat, arr, 'Colors', 'k');
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
box off

text(1.3,5.2, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
text(3.3,5.2, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')

x = [.8 1.8 2.8 3.8];
for i = 1:4
    xx = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    if i == 1 || i == 3
        sc.MarkerFaceColor = [0 1 0];
        sc.MarkerEdgeColor = [0 1 0];
    else
        sc.MarkerFaceColor = [1 0 0];
        sc.MarkerEdgeColor = [1 0 0];
    end
    sc.SizeData = 20;
end
axB.Position(1)     = row(3)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1  = dpE(:,2); d2 = dpE(:,3);
        var = 1;
    else
        d1  = dpD(1:end-10,2); d2 = dpD(1:end-10,3);
        var = 3;
    end
    
    pDP(ii) = signrank(d1,d2);
    if signrank(d1,d2) < .001
        star = plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%% RT %%%%
arrE = []; arrD = []; mat = []; arr = [];
arrE    = [ones(size(mRTE,1),1); ones(size(mRTE,1),1)+1;];
arrD    = [ones(size(mRTD,1)-10,1)+2; ones(size(mRTD,1)-10,1)+3];
mat     = vertcat(mRTE(:,2), mRTE(:,3), mRTD(1:end-10,2), mRTD(1:end-10,3));
arr     = vertcat(arrE,[], arrD);

% disp(mean([mRTE(:,2), mRTE(:,3)]))
% disp(nanmean([mRTD(1:end-10,2), mRTD(1:end-10,3)]))

axB         = axes('Position',[row(4)+of clm(3)+.125 dim]); hold on
line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k')
bx          = boxplot(mat, arr, 'Colors', 'k');
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
box off

text(1.3,850, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
text(3.3,850, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')

x = [.8 1.8 2.8 3.8];
for i = 1:4
    xx = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    if i == 1 || i == 3
        sc.MarkerFaceColor = [0 1 0];
        sc.MarkerEdgeColor = [0 1 0];
    else
        sc.MarkerFaceColor = [1 0 0];
        sc.MarkerEdgeColor = [1 0 0];
    end
    sc.SizeData = 20;
end

axB.Position(1)     = row(4)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1  = mRTE(:,2); d2 = mRTE(:,3);
        var = 1;
    else
        d1  = mRTD(1:end-10,2); d2 = mRTD(1:end-10,3);
        var = 3;
    end
    
    pRT(ii) = signrank(d1,d2);
    if signrank(d1,d2) < .001
        star = plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%% CV %%%
CVE     = sdRTE./mRTE;
CVD     = sdRTD./mRTD;

% disp(mean(CVE))
% disp(nanmean(CVD(1:end-10,:)))

arrE = []; arrD = []; mat = []; arr = [];
arrE    = [ones(size(CVE,1),1); ones(size(CVE,1),1)+1;];
arrD    = [ones(size(CVD,1)-10,1)+2; ones(size(CVD,1)-10,1)+3];

mat     = vertcat(CVE(:,2), CVE(:,3),CVD(1:end-10,2), CVD(1:end-10,3));
arr     = vertcat(arrE, arrD);

axB     = axes('Position',[row(5)+of clm(3)+.125 dim]); hold on
line([2.5 2.5],[-1 1000], 'LineStyle', ':', 'LineWidth',2, 'Color','k')
bx      = boxplot(mat, arr, 'Colors', 'k');
set(bx(end,:),'Visible','off')
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw})
axB.YLabel.String = 'Coefficient of variation';
axB.YLim            = [0.1 .35];
axB.YTick           = [.1 .2 .3];
axB.XTick           = [1 2 3 4];
axB.XTickLabel      = {'Coh8' 'Coh12' 'Coh8' 'Coh12'};
axB.XColor          = [0 0 0];
axB.YColor          = [0 0 0];
axB.FontSize        = 14;
axB.XTickLabelRotation = 30;
box off

text(1.3,.35, 'M1', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')
text(3.3,.35, 'M2', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold')

x = [.8 1.8 2.8 3.8];
for i = 1:4
    xx = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(arr == i));
    sc = scatter(xx, mat(arr == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
    if i == 1 || i == 3
        sc.MarkerFaceColor = [0 1 0];
        sc.MarkerEdgeColor = [0 1 0];
    else
        sc.MarkerFaceColor = [1 0 0];
        sc.MarkerEdgeColor = [1 0 0];
    end
    sc.SizeData = 20;
end

axB.Position(1)     = row(5)+of;
axB.Position(2)     = clm(3)+.125;
axB.Position(3:4)   = [1/7 (clm(4)+.2)-(clm(3)+.125)];

y                   = diff(axB.YLim)*.06;
yy                  = axB.YLim(2) - y;
for ii = 1:2
    d1 = []; d2 = [];
    if ii == 1
        d1  = CVE(:,2); d2 = CVE(:,3);
        var = 1;
    else
        d1  = CVD(1:end-10,2); d2 = CVD(1:end-10,3);
        var = 3;
    end
    
    pCV(ii) = signrank(d1,d2);
    if signrank(d1,d2) < .001
        star = plot([.3 .5 .7]+var,[yy yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .01
        star = plot([.4 .6]+var,[yy yy], 'k*', 'LineWidth', 1.1);
    elseif signrank(d1,d2) < .05
        star = plot([.5]+var,[yy], 'k*', 'LineWidth', 1.1);
    end
    
    if signrank(d1,d2) < .05
        line([.2 .8]+var, [yy-y/2 yy-y/2], 'LineStyle', '-', 'LineWidth',2, 'Color','k')
        star.MarkerSize = 8;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BF + Latency + Phase locking maps
% ADD ARROW: arrow('Start',[0 0 0],'Stop',[1 0 0],'BaseAngle',90,'TipAngle',20,'Width',20,'Length',50)


% dest_dir = 'X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\';
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';

% Frequency pool
freq_pool       = 440 * 2 .^((-31:97)/24);                      % SFG frequency pool
freqStart       = 180;                                          % Tuning low freq [Hz]
steps           = 14;                                           % No of desired tones
PT(1)           = freqStart;                                    % Starting frequency [Hz]

for i = 2:steps                                                 % For no of tones...
    PT(i)       =  PT(i-1)*2^(1/2);                             % 1/2 octave steps
end
frex            = round(PT);
import          = 1;
flag            = 0;
typ             = 'muae';

if import == 0
    for iAn = 1:2
        par = []; mfr_mat = []; mlat_mat = []; cc = 0;
        
        if iAn == 1
            animalID = 'Eric';
        else
            animalID = 'Dollar';
        end
        
        % Chamber angle
        if strcmp(animalID, 'Eric')
            angle   = 10;
        elseif strcmp(animalID, 'Dollar')
            angle   = 15;
        end
        
        %         path    = ['Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\' animalID '\'];
        path    = ['/Volumes/Felix_ExtDrive/Rec/' animalID '/'];
        fr_mat  = nan(18);
        lat_mat = nan(18);
        
        all     = 0;
        incl    = 0;
        idx    	= 1;
        sz   	= 14;
        
        % Get recording dates and initialise output matrix
        par = getParam(animalID);
        for se = 1:length(fieldnames(par))-2
            rc                          = [];
            rc                          = fieldnames(par.(['sess' num2str(se)]));   	% Recording name, string
            for r = 1:length(fieldnames(par.(['sess' num2str(se)])))
                rec                     = par.(['sess' num2str(se)]).(rc{r});
                
                for i = 1:size(rec.nChan,2)
                    chan                = ['ch' num2str(rec.nChan(i))];
                    
                    if isfield(rec, chan)
                        coord               = rec.(chan).coord;
                    else
                        continue
                    end
                    
                    %                     source              = [path rec.fname '\Data\var\'];
                    source              = [path rec.fname '/Data/var/'];
                    
                    if flag == 1
                        typ             = 'spks';
                        fname           = [source 'Tun_Spks_' chan '_clus1.mat'];
                    else
                        typ             = 'muae';
                        fname        	= [source 'Tun_MUAe_' chan '.mat'];
                    end
                    
                    if exist(fname) == 2
                        load(fname) % load spiking data
                        
                        % Get BF from PSTH
                        [mFR, ~]    = getRF(mmat);
                        mmFR        = nanmean(mFR);
                        
                        if isnan(sum(mmFR)) || sum(sum(mmFR)) == 0
                            disp(['BF ' num2str(se) '-' num2str(r) '-' chan])
                            continue
                        end
                        
                        f           = fit((1:length(mmFR))',mmFR','smoothingspline','SmoothingParam',0.95);
                        ff          = f(1:sz);
                        bf          = find(ff == max(ff));
                        BF        	= round(frex(bf));
                        
                        % Load spike latency data
                        for iLat = 1:sz
                            for iVol = 1:3
                                [lat(iVol,iLat), latIdx(iVol,iLat)] = getLatency(mmat{idx, iLat, iVol}, 200, 200, 2, 3);
                                pktime = find(mmat{idx, iLat, iVol}(201:400) == max(mmat{idx, iLat, iVol}(201:400)));
                                if isempty(pktime) || numel(pktime) > 1 || max(mmat{idx, iLat, iVol}(201:400)) == 0
                                    pk(iVol,iLat) = nan;
                                else
                                    pk(iVol,iLat) = pktime;
                                end
                            end
                        end
                        
                        if sum(sum(isnan(pk))) == sz*3
                            disp(['CONTINUE ' num2str(se) '-' num2str(r) '-' chan])
                            continue
                        end
                        
                        %                         % Find min latency
                        %                         ons = lat(:,bf);
                        %
                        %                         if sum(ons < 10) > 0
                        %                             id          = ons < 10;
                        %                             diff        = 10 - ons;
                        %                             ons(id)     = 10;
                        %                             diff (~id)  = 0;
                        %                         else diff       = [0;0;0];
                        %                         end
                        %
                        %                         onLat    	= nanmean(ons);
                        pkLat    	= nanmean(pk(:,bf));
                        
                        if isnan(pkLat)
                            disp(['nanLAT ' num2str(se) '-' num2str(r) '-' chan])
                        elseif isnan(BF)
                            disp(['nanBF ' num2str(se) '-' num2str(r) '-' chan])
                        end
                        
                        matFR       = nan(18);
                        matLat      = nan(18);
                        
                        xxx           = round(coord(1));
                        yyy           = round(coord(2));
                        
%                         % Get corrected coordinates
%                         pen         = rec.(chan).coord(2);                              % ML penetration site on grid [mm]
%                         dep         = rec.(chan).coord(3);                              % depth of recording from GT tip [mm]
%                         offset      = (dep/sind(90)) * sind(angle);                     % calculate offset [mm]
%                         adj         = pen - offset;                                     % adjusted ML value [mm]
%                         xx          = rec.(chan).coord(1);
%                         yy          = adj;
%                         
%                         xxx           = round(xx);
%                         yyy           = round(yy);
                        
                        matFR(xxx, yyy)	= log(BF);
                        matLat(xxx, yyy)= pkLat;
                        % matLat(x, y)= onLat;

                        % Update matrix
                        fr_mat    	= cat(3, fr_mat, matFR);
                        lat_mat  	= cat(3, lat_mat, matLat);
  
                        cc          = cc+1;
                        bfmat(cc)   = find(ff == max(ff));
                        latmat(cc)  = pkLat;
                        ccoord(cc,:)= [xxx,yyy];
                        
                        incl     	= incl + 1;
                        all      	= all + 1;
                    end
                end
            end
        end
        
        % Average across cells
        nCells      = sum(~isnan(fr_mat),3);
        mfr_mat     = nanmean(fr_mat(:,:,2:end),3);
        mlat_mat    = nanmean(lat_mat(:,:,2:end),3);
        
        filtWidth   = 2;
        filtSigma   = 2;
        imageFilter = fspecial('gaussian',filtWidth,filtSigma);
        
        mfr_mat     = nanconv(mfr_mat,imageFilter, 'nanout');
        mlat_mat    = nanconv(mlat_mat,imageFilter, 'nanout');
        
        %         save([dest_dir 'raw\tMap_' animalID '_' typ '.mat'], 'mfr_mat');
        %         save([dest_dir 'raw\lMap_' animalID '_' typ  '.mat'], 'mlat_mat');
        %         save([dest_dir 'raw\ccoord_' animalID '_' typ  '.mat'], 'ccoord');
        %         save([dest_dir 'raw\bf_' animalID '_' typ  '.mat'], 'bfmat');
        %         save([dest_dir 'raw\lat_' animalID '_' typ  '.mat'], 'latmat');
        
        save(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/tMap_' animalID '_' typ '.mat'], 'mfr_mat');
        save(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/lMap_' animalID '_' typ  '.mat'], 'mlat_mat');
        save(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/ccoord_' animalID '_' typ  '.mat'], 'ccoord');
        save(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/bf_' animalID '_' typ  '.mat'], 'bfmat');
        save(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/lat_' animalID '_' typ  '.mat'], 'latmat');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ylim      = 17;
cmap      = flipud(jet(14));
lcmap     = [linspace(0,1,50)', linspace(0,1,50)', ones(50,1)];
alp       = .9;
back      = [0 0 0];

for iAn = 1:2
    
    if iAn == 1
        animalID = 'Eric';
    else
        animalID = 'Dollar';
    end
    
    %     load([dest_dir 'raw\tMap_' animalID '_' typ  '.mat']);
    %     load([dest_dir 'raw\lMap_' animalID '_' typ  '.mat']);
    %     load([dest_dir 'raw\ccoord_' animalID '_' typ  '.mat']);
    %     load([dest_dir 'raw\bf_' animalID '_' typ  '.mat']);
    %     load([dest_dir 'raw\lat_' animalID '_' typ  '.mat']);
    
    load(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/tMap_' animalID '_' typ '.mat']);
    load(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/lMap_' animalID '_' typ  '.mat']);
    load(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/ccoord_' animalID '_' typ  '.mat']);
    load(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/bf_' animalID '_' typ  '.mat']);
    load(['/Users/fschneider/ownCloud/NCL_revision/Figures/raw/lat_' animalID '_' typ  '.mat']);
    
    AP      = find(logical(sum(~isnan(mfr_mat),2)));                                           % Boundaries of recording field
    ML      = find(logical(sum(~isnan(mfr_mat))));
    [xx,yy] = coreBoundary(mfr_mat,AP,ML,false,animalID);                                    	% Get X & Y coordinates for field boundary
    
    x = [];
    for r = 1:length(xx)
        x = horzcat(x, xx(r)-.5, xx(r)+.5);
    end
    
    y = [];
    for r = 1:length(yy)
        y = horzcat(y, yy(r), yy(r));
    end
    
    switch animalID
        case 'Eric'
            xax       = [5 10 15];
            yax       = [-17 -7];
            ytick = {'+5' '+15'};
        case 'Dollar'
            %             xax       = [3 8 13];
            xax       = [5 10 15];
            yax       = [-15 -6];
            ytick = {'+4', '+14'};
    end
    
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iAn == 1
        
        %%% Tonotopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axT1 = axes('Position',[row(3)+of clm(2) dim]); hold on; axis equal
        imagesc(1:size(mfr_mat,1),-18:-1, flipud(mfr_mat));
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        axT1.YLim = [-18 -5];
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        %%% Latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axL1 = axes('Position',[row(4)+of clm(2) dim]); hold on; axis equal
        imagesc(1:size(mlat_mat,1),-18:-1, flipud(mlat_mat));
        hold on
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        axL1.YLim = [-18 -5];
        axPL.YAxis.Visible = 'off';
        caxis([15 80])
        
    else
        
        %%% Tonotopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axT1 = axes('Position',[row(3)+of clm(1) dim]); hold on; axis equal
        imagesc([1:size(mfr_mat,1)], [-18:-1], flipud(mfr_mat));
        hold on
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        axT1.YLim = [-17 -4];
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        %%% Latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axL1 = axes('Position',[row(4)+of clm(1) dim]); hold on; axis equal
        imagesc(1:size(mlat_mat,1),-18:-1, flipud(mlat_mat));
        hold on
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        axL1.YLim = [-17 -4];
        axPL.YAxis.Visible = 'off';
        caxis([15 80])
        
        %%% Phase-lock %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Dollar\Summary\muae.mat')
        %         muaeD = muae; clear muae
        %         co = [];
        %         for ii = 1:size(muaeD,2)
        %             if ~isempty(muaeD{ii}.phLock)
        %                 pL(ii,:) = muaeD{ii}.phLock;
        %                 co(ii,:) = muaeD{ii}.coord;
        %             end
        %         end
        %         save([dest_dir 'raw\PL_id_' animalID '_' typ  '.mat'], 'pL');
        %         save([dest_dir 'raw\PL_co_' animalID '_' typ  '.mat'], 'co');
        
        %         load('X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\raw\PL_co_Dollar_muae.mat')
        %         load('X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\raw\PL_id_Dollar_muae.mat')
        
        load('/Volumes/Felix_ExtDrive/Felix/Documents/Publications/FigGnd_Ephys/Figures/raw/PL_co_Dollar_muae.mat')
        load('/Volumes/Felix_ExtDrive/Felix/Documents/Publications/FigGnd_Ephys/Figures/raw/PL_id_Dollar_muae.mat')
        
        
        c = 1;
        for i = 1:size(co,1)
            MLr      	= co(i,2);                              % ML penetration site on grid [mm]
            dep         = co(i,3);                              % depth of recording from GT tip [mm]
            offset      = (dep/sind(90)) * sind(15);          	% calculate offset [mm]
            adj         = MLr - offset;                      	% adjusted ML value [mm]
            AP(c)      	= co(i,1);
            ML(c)       = adj;
            c = c+1;
        end
        
        axPL = axes('Position',[row(5)+of clm(1) dim]); hold on; axis equal
        imagesc([1:size(mfr_mat,1)], [-18:-1], flipud(mfr_mat));
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        axPL.YLim = [-17 -4];
        
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
    
    axT1.YTick = yax;
    axT1.YTickLabel = ytick;
    axT1.FontSize = 14;
    cm = [back; flipud(jet(256))];
    colormap(axT1, cm)
    
    if iAn == 1
        axT1.YLabel.String = {'M1';'Distance IAL [mm]'};
        axT1.XAxis.Visible = 'off';
        axL1.XAxis.Visible = 'off';
    else
        axT1.YLabel.String = {'M2';'Distance IAL [mm]'};
        axT1.XTick = xax;
        axL1.XTick = xax;
    end
    
    axL1.YTick      = yax;
    axL1.YTickLabel = ytick;
    axL1.FontSize   = 14;
    axL1.YAxis.Visible = 'off';
    cm              = [back; lcmap];
    colormap(axL1, cm)

    if iAn == 2
        axPL.YTick      = yax;
        axPL.YTickLabel = ytick;
        axPL.YTick      = yax;
        axPL.XTick      = xax;
        axPL.FontSize   = 14;
        axPL.XLabel.String = 'Grid position ML [mm]';
        axPL.YAxis.Visible = 'off';
        cm = [back; gray(256)];
        colormap(axPL, cm)
        caxis([floor(log(frex(1))) ceil(log(frex(11)))])
        
        axT1.XLabel.String = 'Grid position ML [mm]';
        axL1.XLabel.String = 'Grid position ML [mm]';
        
    end
end

cm = [back; flipud(jet(256))];
colormap(axT1, cm)
axCB = axes('Position',[row(5)+.025 clm(2)+.05 .001 .1]);
axCB.Visible = 'off';
colormap(axCB, cm)
cb = colorbar(axCB);
cb.Color = [0 0 0];
cb.Position(3) = .01;
cb.Label.String = 'Best frequency [Hz]';
cb.FontSize = 12;
caxis([floor(log(frex(1))) ceil(log(frex(11)))])
cb.Ticks = [log(frex(1)),log(frex(7)),log(frex(11))];
cb.TickLabels = {num2str(frex(1)), num2str(frex(7)), ['>' num2str(frex(11))]};

cm              = [back; gray(256)];
axCB            = axes('Position',[row(5)+.01 clm(2)+.05 .001 .1]);
axCB.Visible    = 'off';
colormap(axCB, cm)
cb              = colorbar(axCB);
cb.Color        = [0 0 0];
cb.Position(3)  = .01;
cb.Ticks        = [];

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


offset = 0.03;
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

% addpath X:\Felix\Scripts\Stuff\export_fig-master
% addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-8b0ba13
% export_fig([dest_dir 'FIG1'], '-r400',f);

%% SFIG Raw BF %%%%

f       = figInit('fig'); axis off; set(gcf,'color',[1 1 1])
ax0  	= axes('Position',[0 0 1 1],'Visible','off');
dim     = [.4 .4];

for iAn = 1:2
    
    if iAn == 1
        animalID = 'Eric';
    else
        animalID = 'Dollar';
    end
    
    switch animalID
        case 'Eric'
            xax       = [5 10 15];
            yax       = [-17 -7];
            ytick = {'+5' '+15'};
        case 'Dollar'
            %             xax       = [3 8 13];
            xax       = [5 10 15];
            yax       = [-15 -6];
            ytick = {'+4', '+14'};
    end
    
    load([dest_dir 'raw\ccoord_' animalID '_' typ  '.mat']);
    load([dest_dir 'raw\bf_' animalID '_' typ  '.mat']);
    load([dest_dir 'raw\lat_' animalID '_' typ  '.mat']);
    
    if iAn == 1
        axT2 = axes('Position',[.1 .55 dim]); hold on; axis equal
    else
        axT2 = axes('Position',[.1 .1 dim]); hold on; axis equal
    end
    
    fi      = fill([0 0 20 20],[-20 0 0 -20], back, 'LineStyle','none');
    C       = unique(ccoord(:,1:2),'rows');
    R       = .25 ;                                                           	% Radius of circle
    th      = linspace(0,2*pi) ;
    axT2.YLim = [-18 -4];
    
    for i = 1:size(C,1)
        rowfind     = (ccoord(:,1:2) - C(i,:)) == 0;
        idx         = sum(rowfind,2) == 2;
        BF          = mean(bfmat(idx));                                         % Show mean BF for overlapping coordinates
        col         = cmap(round(BF),:);                                        % Define colour code
        xxx         = R*cos(th) + C(i,2);                                       % Create circle path
        yyy         = R*sin(th) - C(i,1);
        pa          = patch(xxx,yyy, col);
        pa.FaceAlpha = alp;
        pa.EdgeColor = 'none';
    end
    
    if iAn == 1
        axL2 = axes('Position',[.55 .55 dim]); hold on; axis equal
    else
        axL2 = axes('Position',[.55 .1 dim]); hold on; axis equal
    end
    
    fi      = fill([0 0 20 20],[-20 0 0 -20], back, 'LineStyle','none');
    C       = unique(ccoord(:,1:2),'rows');
    R       = .25 ;                                                           	% Radius of circle
    th      = linspace(0,2*pi) ;
    axL2.YLim = [-18 -4];
    
    for i = 1:size(C,1)
        rowfind     = (ccoord(:,1:2) - C(i,:)) == 0;
        idx         = sum(rowfind,2) == 2;
        lty        	= mean(latmat(idx));                                   	% Show mean BF for overlapping coordinates
        if isnan(lty)
            continue
        end
        
        if lty > size(lcmap,1)
            col = [1 1 1];
        elseif lty <= 10
            col = [0 0 1];
        else
            col         = lcmap(round(lty),:);
        end
        
        % Define colour code
        xxx         = R*cos(th) + C(i,2);                                       % Create circle path
        yyy         = R*sin(th) - C(i,1);
        pa          = patch(xxx,yyy, col) ;
        pa.FaceAlpha = alp;
        pa.EdgeColor = 'none';
    end
    
    if iAn == 1
        axT2.YLabel.String  = {'M1';'Distance IAL [mm]'};
        axT2.XAxis.Visible  = 'off';
        axL2.XAxis.Visible  = 'off';
        axT2.Title.String   = 'Tonotopy';
        axL2.Title.String   = 'Peak Latency';
    else
        axT2.YLabel.String  = {'M2';'Distance IAL [mm]'};
    end
    
    axT2.YTick              = yax;
    axT2.YTickLabel         = ytick;
    axT2.XTick              = xax;
    axT2.FontSize           = 14;
    axT2.XLabel.String      = 'Grid position ML [mm]';
    
    axL2.XLabel.String      = 'Grid position ML [mm]';
    axL2.YTick              = yax;
    axL2.YTickLabel         = ytick;
    axL2.XTick              = xax;
    axL2.FontSize           = 14;
    axL2.YAxis.Visible      = 'off';
    
end

cm      = [back; flipud(jet(256))];
lcmap   = [linspace(0,1,50)', linspace(0,1,50)', ones(50,1)];

axCB                        = axes('Position',[.36 .73 .001 .2]);
axCB.Visible                = 'off';
colormap(axCB, cm)
cb                          = colorbar(axCB);
cb.Color                    = [.9 .9 .9];
cb.Position(3)              = .01;
cb.Label.String             = 'Best frequency [Hz]';
cb.FontSize                 = 12;
caxis([log(frex(1)) log(frex(11))])
caxis([5 log(frex(11))])
cb.Ticks                    = [log(frex(1)),log(frex(7)),log(frex(11))];
cb.TickLabels               = {num2str(frex(1)), num2str(frex(7)), ['>' num2str(frex(11))]};

axCB                        = axes('Position',[.81 .73 .001 .2]);
axCB.Visible                = 'off';
colormap(axCB, lcmap)
cb                          = colorbar(axCB);
cb.Color                    = [.9 .9 .9];
cb.Position(3)              = .01;
cb.Label.String             = 'Latency [ms]';
cb.FontSize                 = 12;
caxis([20 80])
cb.Ticks                    = [20 50 80];
cb.TickLabels               = {'20' '40' '>80'};

addpath X:\Felix\Scripts\Stuff\export_fig-master
% export_fig([dest_dir 'SFIG1'], '-r400',f);