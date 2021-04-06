% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 3 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except muaeE muaeD lfpE lfpD

% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Eric\Summary\muae.mat')
% muaeE = muae;
% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Dollar\Summary\muae.mat')
% muaeD = muae;
% clear muae

% load('/Volumes/Felix_ExtDrive/Rec/Eric/Summary/muae.mat')
muaeE = muae;
load('/Volumes/Felix_ExtDrive/Rec/Dollar/Summary/muae.mat')
muaeD = muae;
clear muae

typ         = 'muae';
% dest_dir    = 'X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\';
dest_dir    = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
indx        = 401:600; % corresponds to 201-400ms after onset
indxR       = 201:400; % corresponds to -300 - -100ms prior decision (HI) / End of figure (MI)
alph        = .01;
nMI         = 20;

% Frequency pool
freq_pool       = 440 * 2 .^((-31:97)/24);                      % SFG frequency pool
freqStart       = 180;                                          % Tuning low freq [Hz]
steps           = 14;                                           % No of desired tones
PT(1)           = freqStart;                                    % Starting frequency [Hz]

for i = 2:steps                                                 % For no of tones...
    PT(i)       =  PT(i-1)*2^(1/2);                             % 1/2 octave steps
end
frex            = round(PT);

for iAn = 1:2
    
    antM = []; posM = [];
    
    if iAn == 1
        animalID = 'Eric';
        load([dest_dir 'raw/tMap_' animalID '_' typ  '.mat']);
        AP      = find(logical(sum(~isnan(mfr_mat),2)));
        ML      = find(logical(sum(~isnan(mfr_mat))));
        [x,y] = coreBoundary(mfr_mat,AP,ML,false,animalID);
        [antM, posM] = sortUnits(x,y,muaeE);
        
    elseif iAn == 2
        animalID = 'Dollar';
        load([dest_dir 'raw/tMap_' animalID '_' typ  '.mat']);
        AP      = find(logical(sum(~isnan(mfr_mat),2)));
        ML      = find(logical(sum(~isnan(mfr_mat))));
        [x,y] = coreBoundary(mfr_mat,AP,ML,false,animalID);
        [antM, posM] = sortUnits(x,y,muaeD);
    end
    
    for iFi = 1:2
        
        if iFi == 1
            d = antM;
            str = 'ant';
        else
            d = posM;
            str = 'pos';
        end
        
        c = 0; dFGo = [];dFGd = []; dCOHo = []; dCOHd = []; dHMo = []; dHMd = []; dMCo = []; dMCd = []; in = [];
        aucFo = []; aucCo = []; aucHMo = []; aucMCo = []; aucFd = []; aucCd = []; aucHMd = []; aucMCd = []; coord = [];
        
        for ii = 1:size(d,2)
            datestr         = str2num([d{ii}.id(1:4) d{ii}.id(6:7) d{ii}.id(9:10)]);
            incl(ii)        = check20Hz(d{ii});
            if length(d{ii}.BLslope) < 200 || incl(ii) == 0 || sum(d{ii}.nTr>=10) ~= length(d{ii}.nTr) || datestr - 20190806 >= 0 || numel(find(sum(~isnan(d{ii}.on.trlOn),2)))< 200
                in(ii)      = false;
                dFGo(ii) 	= nan;
                dFGd(ii)    = nan;
                dCOHo(ii) 	= nan;
                dCOHd(ii)   = nan;
                dHMo(ii)    = nan;
                dHMd(ii)    = nan;
                dMCo(ii)    = nan;
                dMCd(ii)    = nan;
                
                aucFo(ii)  	= nan;
                aucCo(ii) 	= nan;
                aucHMo(ii) 	= nan;
                aucMCo(ii) 	= nan;
                aucFd(ii)  	= nan;
                aucCd(ii)  	= nan;
                aucHMd(ii) 	= nan;
                aucMCd(ii) 	= nan;
                continue
            end
            
            datestr     = str2num([d{ii}.id(1:4) d{ii}.id(6:7) d{ii}.id(9:10)]);
            test        = mean([d{ii}.res.HI12(:,indxR); d{ii}.res.HI8(:,indxR)],2);
            ctrl        = mean(d{ii}.res.CR(:,indxR),2);
            pp          = anova1([test; ctrl],[zeros(size(test,1),1);ones(size(ctrl,1),1)], 'off');
            %             [h,pp]= ttest2(test,ctrl)
            
            if pp < alph
                in(ii)  = true;
            else
                in(ii)  = false;
            end
            
            coord(ii,:)	= d{ii}.coord(1:2);
            
            clear FIGo FIGd HIo8 HId8 HIo12 HId12 CRo CRd MIo MId mBL
            mBL         = nanmean(nanmean(d{ii}.on.fullAvg(:,101:500),2));        	% Average BL response
            
            FIGo     	= d{ii}.on.trlOn(d{ii}.on.cat == 3 | d{ii}.on.cat == 4,indx) ./mBL;
            HIo8      	= d{ii}.on.trlOn((d{ii}.on.cat == 3 | d{ii}.on.cat == 4) & d{ii}.on.coh == 8,indx) ./mBL;
            HIo12    	= d{ii}.on.trlOn((d{ii}.on.cat == 3 | d{ii}.on.cat == 4) & d{ii}.on.coh == 12,indx) ./mBL;
            CRo        	= d{ii}.on.trlOn(d{ii}.on.cat == 6,indx) ./mBL;
            
            if sum(d{ii}.on.cat == 4) >= nMI
                MIo    	= d{ii}.on.trlOn(d{ii}.on.cat == 4,indx) ./mBL;
            else
                MIo  	= nan(2,200);
            end
            
            HId8      	= d{ii}.res.HI8(:,indxR) ./mBL; % Response aligned trials
            HId12    	= d{ii}.res.HI12(:,indxR) ./mBL;
            CRd       	= d{ii}.res.CR(:,indxR) ./mBL;
            
            if size(d{ii}.res.MI,1) >= nMI
                MId   	= (d{ii}.res.MI(:,indxR)) ./mBL;
            else
                MId   	= nan(2,200);
            end
            FIGd     	= [HId12; HId8; MId];
            
            %%% dAB = (mA-mB)/s, where mA and mB are the mean responses in stimulus
            %%% conditions A and B, and s is the pooled standard deviation.
            %%% Poort et al. (2016) Cereb Cortex
            dFGo(ii)   	= mean((nanmean(FIGo) - nanmean(CRo)) ./ nanstd([FIGo; CRo]));
            dFGd(ii)   	= mean((nanmean(FIGd) - nanmean(CRd)) ./ nanstd([FIGd; CRd]));
            
            dCOHo(ii)  	= mean((nanmean(HIo12) - nanmean(HIo8)) ./ nanstd([HIo8; HIo12]));
            dCOHd(ii)   = mean((nanmean(HId12) - nanmean(HId8)) ./ nanstd([HId8; HId12]));
            
            dHMo(ii)   	= mean((nanmean([HIo8; HIo12]) - nanmean(MIo)) ./ nanstd([HIo8; HIo12; MIo]));
            dHMd(ii)   	= mean((nanmean([HId8; HId12]) - nanmean(MId)) ./ nanstd([HId8; HId12; MId]));
            
            dMCo(ii)    = mean((nanmean(MIo) - nanmean(CRo)) ./ nanstd([MIo; CRo]));
            dMCd(ii)  	= mean((nanmean(MId) - nanmean(CRd)) ./ nanstd([MId; CRd]));
            
            %%% AUROC %%%
            aucFo(ii) 	= f_auroc(nanmean(CRo,2),nanmean(FIGo,2));
            aucMCo(ii)	= f_auroc(nanmean(CRo,2),nanmean(MIo,2));
            aucHMo(ii) 	= f_auroc(nanmean(MIo,2),nanmean([HIo8;HIo12],2));
            aucCo(ii)  	= f_auroc(nanmean(HIo8,2),nanmean(HIo12,2));
           
            aucFd(ii) 	= f_auroc(nanmean(CRd,2),nanmean(FIGd,2));
            aucMCd(ii) 	= f_auroc(nanmean(CRd,2),nanmean(MId,2));
            aucHMd(ii)	= f_auroc(nanmean(MId,2),nanmean([HId8;HId12],2));
            aucCd(ii)  	= f_auroc(nanmean(HId8,2),nanmean(HId12,2));
        end
        
        IN{iAn,iFi}     = logical(in);
        COORD{iAn,iFi}  = coord;
        FGo{iAn,iFi}    = dFGo;
        FGd{iAn,iFi}    = dFGd;
        COHo{iAn,iFi}   = dCOHo;
        COHd{iAn,iFi}   = dCOHd;
        HMo{iAn,iFi}   	= dHMo;
        HMd{iAn,iFi}    = dHMd;
        MCo{iAn,iFi}    = dMCo;
        MCd{iAn,iFi}    = dMCd;
        
        AUCFD{iAn,iFi}	= aucFd;
        AUCFO{iAn,iFi} 	= aucFo;
        AUCC{iAn,iFi}   = aucCo;
        AUCCD{iAn,iFi} 	= aucCd;
        AUCHM{iAn,iFi} 	= aucHMo;
        AUCMC{iAn,iFi}  = aucMCo;
        AUCHMD{iAn,iFi}	= aucHMd;
        AUCMCD{iAn,iFi} = aucMCd;
        
        pFd(iAn,iFi) 	= signrank(FGd{iAn,iFi}(IN{iAn,iFi}));
        pF(iAn,iFi)    	= signrank(FGo{iAn,iFi}(IN{iAn,iFi}));
        pC(iAn,iFi)    	= signrank(COHo{iAn,iFi}(IN{iAn,iFi}));
        pCd(iAn,iFi)  	= signrank(COHd{iAn,iFi}(IN{iAn,iFi}));
        pHM(iAn,iFi)   	= signrank(HMo{iAn,iFi}(IN{iAn,iFi}));
        pHMd(iAn,iFi) 	= signrank(HMd{iAn,iFi}(IN{iAn,iFi}));
        pMC(iAn,iFi)   	= signrank(MCo{iAn,iFi}(IN{iAn,iFi}));
        pMCd(iAn,iFi)  	= signrank(MCd{iAn,iFi}(IN{iAn,iFi}));
        
        npFd(iAn,iFi)  	= signrank(FGd{iAn,iFi}(~IN{iAn,iFi}));
        npF(iAn,iFi)   	= signrank(FGo{iAn,iFi}(~IN{iAn,iFi}));
        npC(iAn,iFi)   	= signrank(COHo{iAn,iFi}(~IN{iAn,iFi}));
        npCd(iAn,iFi)  	= signrank(COHd{iAn,iFi}(~IN{iAn,iFi}));
        npHM(iAn,iFi)  	= signrank(HMo{iAn,iFi}(~IN{iAn,iFi}));
        npHMd(iAn,iFi) 	= signrank(HMd{iAn,iFi}(~IN{iAn,iFi}));
        npMC(iAn,iFi) 	= signrank(MCo{iAn,iFi}(~IN{iAn,iFi}));
        npMCd(iAn,iFi) 	= signrank(MCd{iAn,iFi}(~IN{iAn,iFi}));
    end
end

%%% POOL AC %%%
COORD_pop               = [COORD{1,1}; COORD{1,2}; COORD{2,1}; COORD{2,2}];
allAUCFD_pop            = [AUCFD{1,1} AUCFD{1,2} AUCFD{2,1} AUCFD{2,2}];
allAUCFO_pop            = [AUCFO{1,1} AUCFO{1,2} AUCFO{2,1} AUCFO{2,2}];
allAUCC_pop             = [AUCC{1,1} AUCC{1,2} AUCC{2,1} AUCC{2,2}];
allAUCCD_pop            = [AUCCD{1,1} AUCCD{1,2} AUCCD{2,1} AUCCD{2,2}];
allAUCHM_pop            = [AUCHM{1,1} AUCHM{1,2} AUCHM{2,1} AUCHM{2,2}];
allAUCMC_pop            = [AUCMC{1,1} AUCMC{1,2} AUCMC{2,1} AUCMC{2,2}];
allAUCHMD_pop           = [AUCHMD{1,1} AUCHMD{1,2} AUCHMD{2,1} AUCHMD{2,2}];
allAUCMCD_pop           = [AUCMCD{1,1} AUCMCD{1,2} AUCMCD{2,1} AUCMCD{2,2}];

inclIdx                 = [IN{1,1} IN{1,2} IN{2,1} IN{2,2}];    % Only sound-responsive units with sign. FGM
PFd                     = signrank(allAUCFD_pop(inclIdx),.5);
PF                      = signrank(allAUCFO_pop(inclIdx),.5);
PC                      = signrank(allAUCC_pop(inclIdx),.5);
PCd                     = signrank(allAUCCD_pop(inclIdx),.5);
PHM                     = signrank(allAUCHM_pop(inclIdx),.5);
PHMd                    = signrank(allAUCHMD_pop(inclIdx),.5);
PMC                     = signrank(allAUCMC_pop(inclIdx),.5);
PMCd                    = signrank(allAUCMCD_pop(inclIdx),.5);

nPFd                    = signrank(allAUCFD_pop(~inclIdx),.5);
nPF                     = signrank(allAUCFO_pop(~inclIdx),.5);
nPC                     = signrank(allAUCC_pop(~inclIdx),.5);
nPCd                    = signrank(allAUCCD_pop(~inclIdx),.5);
nPHM                    = signrank(allAUCHM_pop(~inclIdx),.5);
nPHMd                   = signrank(allAUCHMD_pop(~inclIdx),.5);
nPMC                    = signrank(allAUCMC_pop(~inclIdx),.5);
nPMCd                   = signrank(allAUCMCD_pop(~inclIdx),.5);

% PFd    	= fdr(PFd);
% PF        = fdr(PF);
% PC      = fdr(PC);
% PCd     = fdr(PCd);
% PHM     = fdr(PHM);
% PHMd    = fdr(PHMd);
% PMC     = fdr(PMC);
% PMCd    = fdr(PMCd);
%
% nPFd    = fdr(nPFd);
% nPF     = fdr(nPF);
% nPC     = fdr(nPC);
% nPCd    = fdr(nPCd);
% nPHM    = fdr(nPHM);
% nPHMd   = fdr(nPHMd);
% nPMC    = fdr(nPMC);
% nPMCd   = fdr(nPMCd);

%% 

f           = figure('Units', 'normalized', 'Position', [0 0 1 1]); axis off
set(gcf,'color', [1 1 1]);
alp         = .7;
clm         = linspace(.08, .79, 4);
row         = fliplr(linspace(.05, .76, 4));
dim         = [.2 .2];
txtsz       = 14;

bin         = linspace(0, 1,50);
% inclIdx     = allIN_pop;    % Only sound-responsive units with sign. FGM

%%% COLUMN 1: Fig vs Ctr
axO                 = axes('Position',[clm(1) row(1) dim]); hold on
vec1                = allAUCFO_pop(~inclIdx);
vec2                = allAUCFO_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axO.YLabel.String   = {'Onset-aligned';'No. units'};
Yof                 = axO.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axO.YLim            = [0 maxV+(2*Yof)];
axO.XLim            = [0 1];
axO.FontSize        = txtsz;

line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPF, 1, 1, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PF, 1, 1, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)

lg                  = legend('Unresponsive','Modulated', 'Location','Northwest');
lg.FontSize         = 12;
legend boxoff

axD                 = axes('Position',[clm(1) row(3) dim]); hold on
vec1                = allAUCFD_pop(~inclIdx);
vec2                = allAUCFD_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axD.YLabel.String   = {'Response-aligned';'No. units'};
Yof                 = axD.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axD.YLim            = [0 maxV+(2*Yof)];
axD.XLim            = [0 1];
axD.FontSize        = txtsz;

line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPFd, 1, 1, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PFd, 1, 1, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)

axO.XTick = [0 .25 .5 .75 1];
axO.XLabel.String = 'AUROC';
axD.XTick = [0 .25 .5 .75 1];
axD.XLabel.String = 'AUROC';

%%% COLUMN 2: Coherence
axO                 = axes('Position',[clm(2) row(1) dim]); hold on
vec1                = allAUCC_pop(~inclIdx);
vec2                = allAUCC_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axO.YLabel.String   = {'Onset-aligned';'No. units'};
Yof                 = axO.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axO.YLim            = [0 maxV+(2*Yof)];
axO.XLim            = [0 1];
axO.FontSize        = txtsz;

line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPC, 1, 1, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PC, 1, 1, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)


axD                 = axes('Position',[clm(2) row(3) dim]); hold on
vec1                = allAUCCD_pop(~inclIdx);
vec2                = allAUCCD_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axD.YLabel.String   = {'Response-aligned';'No. units'};
Yof                 = axD.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axD.YLim            = [0 maxV+(2*Yof)];
axD.XLim            = [0 1];
axD.FontSize        = txtsz;

line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPCd, 1, 1, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PCd, 1, 1, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)

axO.XTick = [0 .25 .5 .75 1];
axO.XLabel.String = 'AUROC';
axD.XTick = [0 .25 .5 .75 1];
axD.XLabel.String = 'AUROC';

%%% COLUMN 3: HI vs MI
axO                 = axes('Position',[clm(3) row(1) dim]); hold on
vec1                = allAUCHM_pop(~inclIdx);
vec2                = allAUCHM_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axO.YLabel.String   = {'Onset-aligned';'No. units'};
Yof                 = axO.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axO.YLim            = [0 maxV+(2*Yof)];
axO.XLim            = [0 1];
axO.FontSize        = txtsz;

line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPHM, 1, 1, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PHM, 1, 1, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)

axD                 = axes('Position',[clm(3) row(3) dim]); hold on
vec1                = allAUCHMD_pop(~inclIdx);
vec2                = allAUCHMD_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axD.YLabel.String   = {'Response-aligned';'No. units'};
Yof                 = axD.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axD.YLim            = [0 maxV+(2*Yof)];
axD.XLim            = [0 1];
axD.FontSize        = txtsz;

line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPHMd, 1, 1, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PHMd, 1, 1, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)

axO.XTick = [0 .25 .5 .75 1];
axO.XLabel.String = 'AUROC';
axD.XTick = [0 .25 .5 .75 1];
axD.XLabel.String = 'AUROC';

%%% COLUMN 4: MI vs CR
axO     = axes('Position',[clm(4) row(1) dim]); hold on

vec1                = allAUCMC_pop(~inclIdx);
vec2                = allAUCMC_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axO.YLabel.String   = {'Onset-aligned';'No. units'};
Yof                 = axO.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axO.YLim            = [0 maxV+(2*Yof)];
axO.XLim            = [0 1];
axO.FontSize        = txtsz;

line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPMC, 1, 1, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PMC, 1, 1, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)


axD     = axes('Position',[clm(4) row(3) dim]); hold on

vec1                = allAUCMCD_pop(~inclIdx);
vec2                = allAUCMCD_pop(inclIdx);
h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
axD.YLabel.String   = {'Response-aligned';'No. units'};
Yof                 = axD.YLim(2)/5;
maxV                = max([h1.Values,h2.Values]);
axD.YLim            = [0 maxV+(2*Yof)];
axD.XLim            = [0 1];
axD.FontSize        = txtsz;

line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
plotStar(nPMCd, 1, 1, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
plotStar(PMCd, 1, 1, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
% testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)

axO.XTick = [0 .25 .5 .75 1];
axO.XLabel.String = 'AUROC';
axD.XTick = [0 .25 .5 .75 1];
axD.XLabel.String = 'AUROC';


%%% MAPS %%%
back = [1 1 1];
col = [linspace(0,1,40)',linspace(1,0,40)',linspace(1,0,40)'];
col = [repmat(col(1,:),[30,1]); col; repmat(col(end,:),[30,1])];
dim = [.09 .18];
falph = .7;
for iAlig = [2 4]
for iCol = 1:4
    for iAn = 1:2
        
        if iAn == 1
            animalID = 'Eric';
        else
            animalID = 'Dollar';
        end
        
        typ     = 'muae';
        load([dest_dir 'raw/tMap_' animalID '_' typ  '.mat']);
        
        AP      = find(logical(sum(~isnan(mfr_mat),2)));                                           % Boundaries of recording field
        ML      = find(logical(sum(~isnan(mfr_mat))));
        [xx,yy] = coreBoundary(mfr_mat,AP,ML,false,animalID);                                    	% Get X & Y coordinates for field boundary
        
        x       = [];
        for r = 1:length(xx)
            x   = horzcat(x, xx(r)-.5, xx(r)+.5);
        end
        
        y       = [];
        for r = 1:length(yy)
            y   = horzcat(y, yy(r), yy(r));
        end
        
        switch animalID
            case 'Eric'
                xax         = [5 10 15];
                yax         = [-17 -7];
                ytick       = {'+5' '+15'};
            case 'Dollar'
                %             xax       = [3 8 13];
                xax         = [5 10 15];
                yax         = [-15 -6];
                ytick       = {'+4', '+14'};
        end
        
        
        if iAn == 1
            axM             = axes('Position',[clm(iCol) row(iAlig)+.015 dim]); hold on; axis equal
        else
            axM             = axes('Position',[clm(iCol)+dim(1)+.02 row(iAlig)+.013 dim]); hold on; axis equal
        end
        
        imagesc([1:size(mfr_mat,1)], [-18:-1], flipud(mfr_mat));
        plot(x,-y, 'Color', [0 0 0],'LineWidth', 4)
        
        if iAn == 1
            axM.YLim            = [-18 -5];
            axM.XLim            = [5 15];

        else
            axM.YLim            = [-16 -3];
            axM.XLim            = [4 13];
        end
        
        axM.YTick           = yax;
        axM.YTickLabel      = ytick;
        axM.YTick           = yax;
        axM.XTick         	= xax;
        axM.FontSize        = 14;
        axM.XLabel.String   = 'Grid position ML [mm]';
        axM.YLabel.String   = 'Distance IAL [mm]';
        cm                  = [back; gray(256)];
        colormap(axM, cm)
        caxis([floor(log(frex(1))) (ceil(max(max(mfr_mat))*10)/10)+.3])

        if iAn == 1 && iCol == 2 && iAlig == 4
            axM.YAxis.Visible	= 'off';
        else
            axM.YAxis.Visible	= 'off';
            axM.XAxis.Visible 	= 'off';
        end
        
        crd_a = []; crd_p = []; auc = [];
        crd_a 	= COORD{iAn,1}(IN{iAn,1},:);
        
        if iAlig == 2 && iCol == 1
            auc             = round(AUCFO{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 4 && iCol == 1
            auc             = round(AUCFD{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 2 && iCol == 2
            auc             = round(AUCC{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 4 && iCol == 2
            auc             = round(AUCCD{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 2 && iCol == 3
            auc             = round(AUCHM{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 4 && iCol == 3
            auc             = round(AUCHMD{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 2 && iCol == 4
            auc             = round(AUCMC{iAn,1}(IN{iAn,1})*100);
        elseif iAlig == 4 && iCol == 4
            auc             = round(AUCMCD{iAn,1}(IN{iAn,1})*100);
        end
            
        for iUnit = 1:size(crd_a,1)
            sc = scatter(crd_a(iUnit,2),-crd_a(iUnit,1), '^', 'filled', 'MarkerFaceAlpha', falph);
            sc.MarkerEdgeColor = 'none';
            if auc(iUnit) > 79
                sc.MarkerFaceColor = col(80,:);
            elseif auc(iUnit) < 22
                sc.MarkerFaceColor = col(21,:);
            elseif isnan(auc(iUnit))
                sc.MarkerFaceColor = 'none';
            else
                sc.MarkerFaceColor = col(auc(iUnit),:);
            end
        end
        
        crd_p   = COORD{iAn,2}(IN{iAn,2},:);
        auc     = [];
        if iAlig == 2 && iCol == 1
            auc             = round(AUCFO{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 4 && iCol == 1
            auc             = round(AUCFD{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 2 && iCol == 2
            auc             = round(AUCC{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 4 && iCol == 2
            auc             = round(AUCCD{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 2 && iCol == 3
            auc             = round(AUCHM{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 4 && iCol == 3
            auc             = round(AUCHMD{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 2 && iCol == 4
            auc             = round(AUCMC{iAn,2}(IN{iAn,2})*100);
        elseif iAlig == 4 && iCol == 4
            auc             = round(AUCMCD{iAn,2}(IN{iAn,2})*100);
        end        
        
        for iUnit = 1:size(crd_p,1)
            sc = scatter(crd_p(iUnit,2),-crd_p(iUnit,1), '^', 'filled', 'MarkerFaceAlpha', falph);
            sc.MarkerEdgeColor = 'none';
            
            if auc(iUnit) > 79
                sc.MarkerFaceColor = col(80,:);
            elseif auc(iUnit) < 22
                sc.MarkerFaceColor = col(21,:);
            elseif isnan(auc(iUnit))
                sc.MarkerFaceColor = 'none';
            else
                sc.MarkerFaceColor = col(auc(iUnit),:);
            end
            
        end
        
    end
end
end

of                          = .06;
axCB                        = axes('Position',[clm(1)-of row(2)+of .001 .1]);
axCB.Visible                = 'off';
colormap(axCB, col(31:70,:))
cb                          = colorbar(axCB);
cb.Color                    = [0 0 0];
cb.Position(3)              = .01;
cb.Label.String             = 'AUROC';
cb.FontSize                 = 12;
caxis([.3 .7])
cb.Ticks                    = [.3 .4 .5 .6 .7];


for i = 2:steps                                                             % For no of tones...
    PT(i)                   =  PT(i-1)*2^(1/2);                             % 1/2 octave steps
end
frex                        = round(PT);

cm                          = gray(256);
axCB                        = axes('Position',[clm(1)-of row(4)+of .001 .1]);
axCB.Visible                = 'off';
colormap(axCB, cm)
cb                          = colorbar(axCB);
cb.Color                    = [0 0 0];
cb.Position(3)              = .01;
cb.Label.String             = 'Best frequency [Hz]';
cb.FontSize                 = 12;
caxis([floor(log(frex(1))) ceil(log(frex(11)))])
cb.Ticks                    = [log(frex(1)),log(frex(7)),log(frex(11))];
cb.TickLabels               = {num2str(frex(1)), num2str(frex(7)), ['>' num2str(frex(11))]};

ax0                         = axes('Position',[0 0 1 1],'Visible','off');
text(clm(1)-of,.99, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)-of,.99, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)-of,.99, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)-of,.99, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(1)-of,.5, 'e', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)-of,.5, 'f', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)-of,.5, 'g', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)-of,.5, 'h', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
% text(clm(1)+.04,of, 'M1', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
% text(clm(1)+.16,of, 'M2', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')

of = .015;
text(clm(1)+of, .985, 'Figure vs Control', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)+of, .985, 'Coherence8 vs Coherence12', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)+of, .985, 'Hit vs Miss', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)+of, .985, 'Miss vs Correct rejection', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
export_fig([dest_dir '/FIG3_roc_pop_pooled'], '-r400',f);
% % exportgraphics(f,[dest_dir '/FIG3_roc_pop_pooled.pdf'],'ContentType','vector','Resolution',400)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'FIG3_roc_pop_pooled'], '-dpdf', '-r400'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Field-specific AUROC

for iAn = 1:2
    for i = 1:2
        inclIdx                 = IN{iAn,i};    % Only sound-responsive units with sign. FGM
        pFd(iAn,i)              = signrank(AUCFD{iAn,i}(inclIdx),.5);
        pF(iAn,i)            	= signrank(AUCFO{iAn,i}(inclIdx),.5);
        pC(iAn,i)               = signrank(AUCC{iAn,i}(inclIdx),.5);
        pCd(iAn,i)              = signrank(AUCCD{iAn,i}(inclIdx),.5);
        pHM(iAn,i)              = signrank(AUCHM{iAn,i}(inclIdx),.5);
        pHMd(iAn,i)             = signrank(AUCHMD{iAn,i}(inclIdx),.5);
        pMC(iAn,i)              = signrank(AUCMC{iAn,i}(inclIdx),.5);
        pMCd(iAn,i)             = signrank(AUCMCD{iAn,i}(inclIdx),.5);
        
        npFd(iAn,i)              = signrank(AUCFD{iAn,i}(~inclIdx),.5);
        npF(iAn,i)               = signrank(AUCFO{iAn,i}(~inclIdx),.5);
        npC(iAn,i)               = signrank(AUCC{iAn,i}(~inclIdx),.5);
        npCd(iAn,i)              = signrank(AUCCD{iAn,i}(~inclIdx),.5);
        npHM(iAn,i)              = signrank(AUCHM{iAn,i}(~inclIdx),.5);
        npHMd(iAn,i)             = signrank(AUCHMD{iAn,i}(~inclIdx),.5);
        npMC(iAn,i)              = signrank(AUCMC{iAn,i}(~inclIdx),.5);
        npMCd(iAn,i)             = signrank(AUCMCD{iAn,i}(~inclIdx),.5);
    end
end

pFd    	= fdrP([pFd(1,:), pFd(2,:)]);
pF    	= fdrP([pF(1,:), pF(2,:)]);
pC      = fdrP([pC(1,:), pC(2,:)]);
pCd     = fdrP([pCd(1,:), pCd(2,:)]);
pHM     = fdrP([pHM(1,:), pHM(2,:)]);
pHMd    = fdrP([pHMd(1,:), pHMd(2,:)]);
pMC     = fdrP([pMC(1,:), pMC(2,:)]);
pMCd    = fdrP([pMCd(1,:), pMCd(2,:)]);

npFd    = fdrP([npFd(1,:), npFd(2,:)]);
npF     = fdrP([npF(1,:), npF(2,:)]);
npC     = fdrP([npC(1,:), npC(2,:)]);
npCd    = fdrP([npCd(1,:), npCd(2,:)]);
npHM    = fdrP([npHM(1,:), npHM(2,:)]);
npHMd   = fdrP([npHMd(1,:), npHMd(2,:)]);
npMC    = fdrP([npMC(1,:), npMC(2,:)]);
npMCd   = fdrP([npMCd(1,:), npMCd(2,:)]);

f           = figure('Units', 'normalized', 'Position', [0 0 1 1]); axis off
set(gcf,'color', [1 1 1]);
cm          = [[0 0 .9]; [0 .9 0]; [.9 0 0]];
alp         = .7;
clm         = linspace(.1, .79, 4);
row         = fliplr(linspace(.05, .86, 8));
dim         = [.2 .1];
txtsz       = 14;

%%% PLOT %%%
for iAn = 1:2
    for iFi = 1:2
        bin         = linspace(0, 1,50);
        inclIdx     = IN{iAn,iFi};    % Only sound-responsive units with sign. FGM
        
        %%% COLUMN 1: Fig vs Ctr
        if iAn == 1 && iFi == 1
            axO     = axes('Position',[clm(1) row(1) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axO     = axes('Position',[clm(1) row(3) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axO     = axes('Position',[clm(1) row(5) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axO     = axes('Position',[clm(1) row(7) dim]); hold on
        end
        
        vec1                = AUCFO{iAn,iFi}(~inclIdx);
        vec2                = AUCFO{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        axO.FontSize        = txtsz;
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npF, iAn, iFi, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pF, iAn, iFi, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 && iFi == 1
            lg = legend('Unresponsive','Modulated', 'Location','Northwest');
            lg.FontSize = 8;
            legend boxoff
        end
        
        if iAn == 1 && iFi == 1
            axD     = axes('Position',[clm(1) row(2) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axD     = axes('Position',[clm(1) row(4) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axD     = axes('Position',[clm(1) row(6) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axD     = axes('Position',[clm(1) row(8) dim]); hold on
        end
        
        vec1                = AUCFD{iAn,iFi}(~inclIdx);
        vec2                = AUCFD{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        axD.FontSize        = txtsz;
    
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npFd, iAn, iFi, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pFd, iAn, iFi, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
        
        %%% COLUMN 2: Coherence
        if iAn == 1 && iFi == 1
            axO     = axes('Position',[clm(2) row(1) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axO     = axes('Position',[clm(2) row(3) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axO     = axes('Position',[clm(2) row(5) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axO     = axes('Position',[clm(2) row(7) dim]); hold on
        end
        
        vec1                = AUCC{iAn,iFi}(~inclIdx);
        vec2                = AUCC{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        axO.FontSize        = txtsz;
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npC, iAn, iFi, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pC, iAn, iFi, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 && iFi == 1
            axD     = axes('Position',[clm(2) row(2) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axD     = axes('Position',[clm(2) row(4) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axD     = axes('Position',[clm(2) row(6) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axD     = axes('Position',[clm(2) row(8) dim]); hold on
        end
        
        vec1                = AUCCD{iAn,iFi}(~inclIdx);
        vec2                = AUCCD{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+Yof];
        axD.FontSize        = txtsz;
       
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npCd, iAn, iFi, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pCd, iAn, iFi, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
        
        %%% COLUMN 3: HI vs MI
        if iAn == 1 && iFi == 1
            axO     = axes('Position',[clm(3) row(1) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axO     = axes('Position',[clm(3) row(3) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axO     = axes('Position',[clm(3) row(5) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axO     = axes('Position',[clm(3) row(7) dim]); hold on
        end
        
        vec1                = AUCHM{iAn,iFi}(~inclIdx);
        vec2                = AUCHM{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        axO.FontSize        = txtsz;
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npHM, iAn, iFi, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pHM, iAn, iFi, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 && iFi == 1
            axD     = axes('Position',[clm(3) row(2) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axD     = axes('Position',[clm(3) row(4) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axD     = axes('Position',[clm(3) row(6) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axD     = axes('Position',[clm(3) row(8) dim]); hold on
        end
        
        vec1                = AUCHMD{iAn,iFi}(~inclIdx);
        vec2                = AUCHMD{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        axD.FontSize        = txtsz;

        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npHMd, iAn, iFi, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pHMd, iAn, iFi, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
        
        %%% COLUMN 4: MI vs CR
        if iAn == 1 && iFi == 1
            axO     = axes('Position',[clm(4) row(1) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axO     = axes('Position',[clm(4) row(3) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axO     = axes('Position',[clm(4) row(5) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axO     = axes('Position',[clm(4) row(7) dim]); hold on
        end
        
        vec1                = AUCMC{iAn,iFi}(~inclIdx);
        vec2                = AUCMC{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];        
        axO.FontSize        = txtsz;
       
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npMC, iAn, iFi, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pMC, iAn, iFi, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 && iFi == 1
            axD     = axes('Position',[clm(4) row(2) dim]); hold on
        elseif iAn == 1 && iFi == 2
            axD     = axes('Position',[clm(4) row(4) dim]); hold on
        elseif iAn == 2 && iFi == 1
            axD     = axes('Position',[clm(4) row(6) dim]); hold on
        elseif iAn == 2 && iFi == 2
            axD     = axes('Position',[clm(4) row(8) dim]); hold on
        end
        
        vec1                = AUCMCD{iAn,iFi}(~inclIdx);
        vec2                = AUCMCD{iAn,iFi}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        axD.FontSize        = txtsz;

        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(npMCd, iAn, iFi, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(pMCd, iAn, iFi, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
    end
end

of = .05;
ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(clm(1)-of,.98, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)-of,.98, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)-of,.98, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)-of,.98, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

annotation('line',[.05 .05],[row(1)+.1-.03 row(2)+.03], 'LineWidth', 2, 'Color', 'k')
annotation('line',[.05 .05],[row(3)+.1-.03 row(4)+.03], 'LineWidth', 2, 'Color', 'k')
annotation('line',[.05 .05],[row(5)+.1-.03 row(6)+.03], 'LineWidth', 2, 'Color', 'k')
annotation('line',[.05 .05],[row(7)+.1-.03 row(8)+.03], 'LineWidth', 2, 'Color', 'k')
text(.035,row(2)+.09, 'ANT', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k','Rotation',90)
text(.035,row(4)+.09, 'POS', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k','Rotation',90)
text(.02,row(3)+.09, 'M1', 'Parent', ax0, 'FontSize', txtsz+6, 'Color', 'k', 'FontWeight', 'bold','Rotation',90)
text(.035,row(6)+.09, 'ANT', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k','Rotation',90)
text(.035,row(8)+.09, 'POS', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k','Rotation',90)
text(.02,row(7)+.09, 'M2', 'Parent', ax0, 'FontSize', txtsz+6, 'Color', 'k', 'FontWeight', 'bold','Rotation',90)

of = .015;
text(clm(1)+of, .985, 'Figure vs Control', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)+of, .985, 'Coherence8 vs Coherence12', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)+of, .985, 'Hit vs Miss', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)+of, .985, 'Miss vs Correct rejection', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')

% addpath X:\Felix\Scripts\Stuff\export_fig-master
% dest_dir = 'X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\';
export_fig([dest_dir 'FIG3_roc'], '-r400',f);
% exportgraphics(f,[dest_dir '/FIG3_roc.pdf'],'ContentType','vector','Resolution',300)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'FIG3_roc'], '-dpdf', '-r400'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pooled AUROC - Aligned responses for entire AC

for iAn = 1:2
    allIN{iAn}              = [IN{iAn,1} IN{iAn,2}];
    allAUCFD{iAn}           = [AUCFD{iAn,1} AUCFD{iAn,2}];
    allAUCFO{iAn}           = [AUCFO{iAn,1} AUCFO{iAn,2}];
    allAUCC{iAn}            = [AUCC{iAn,1} AUCC{iAn,2}];
    allAUCCD{iAn}           = [AUCCD{iAn,1} AUCCD{iAn,2}];
    allAUCHM{iAn}           = [AUCHM{iAn,1} AUCHM{iAn,2}];
    allAUCMC{iAn}           = [AUCMC{iAn,1} AUCMC{iAn,2}];
    allAUCHMD{iAn}          = [AUCHMD{iAn,1} AUCHMD{iAn,2}];
    allAUCMCD{iAn}          = [AUCMCD{iAn,1} AUCMCD{iAn,2}];
    
    inclIdx                 = allIN{iAn};    % Only sound-responsive units with sign. FGM
    PFd(iAn)                = signrank(allAUCFD{iAn}(inclIdx),.5);
    PF(iAn)                 = signrank(allAUCFO{iAn}(inclIdx),.5);
    PC(iAn)                 = signrank(allAUCC{iAn}(inclIdx),.5);
    PCd(iAn)                = signrank(allAUCCD{iAn}(inclIdx),.5);
    PHM(iAn)                = signrank(allAUCHM{iAn}(inclIdx),.5);
    PHMd(iAn)               = signrank(allAUCHMD{iAn}(inclIdx),.5);
    PMC(iAn)                = signrank(allAUCMC{iAn}(inclIdx),.5);
    PMCd(iAn)               = signrank(allAUCMCD{iAn}(inclIdx),.5);
    
    nPFd(iAn)               = signrank(allAUCFD{iAn}(~inclIdx),.5);
    nPF(iAn)                = signrank(allAUCFO{iAn}(~inclIdx),.5);
    nPC(iAn)                = signrank(allAUCC{iAn}(~inclIdx),.5);
    nPCd(iAn)               = signrank(allAUCCD{iAn}(~inclIdx),.5);
    nPHM(iAn)               = signrank(allAUCHM{iAn}(~inclIdx),.5);
    nPHMd(iAn)              = signrank(allAUCHMD{iAn}(~inclIdx),.5);
    nPMC(iAn)               = signrank(allAUCMC{iAn}(~inclIdx),.5);
    nPMCd(iAn)              = signrank(allAUCMCD{iAn}(~inclIdx),.5);
end

PFd    	= fdr(PFd);
PF    	= fdr(PF);
PC      = fdr(PC);
PCd     = fdr(PCd);
PHM     = fdr(PHM);
PHMd    = fdr(PHMd);
PMC     = fdr(PMC);
PMCd    = fdr(PMCd);

nPFd    = fdr(nPFd);
nPF     = fdr(nPF);
nPC     = fdr(nPC);
nPCd    = fdr(nPCd);
nPHM    = fdr(nPHM);
nPHMd   = fdr(nPHMd);
nPMC    = fdr(nPMC);
nPMCd   = fdr(nPMCd);

f           = figure('Units', 'normalized', 'Position', [0 0 1 1]); axis off
set(gcf,'color', [1 1 1]);
alp         = .7;
clm         = linspace(.08, .79, 4);
row         = fliplr(linspace(.05, .76, 4));
dim         = [.2 .2];
txtsz       = 18;

for iAn = 1:2    
       bin         = linspace(0, 1,50);
       inclIdx     = allIN{iAn};    % Only sound-responsive units with sign. FGM
        
        %%% COLUMN 1: Fig vs Ctr
        if iAn == 1
            axO     = axes('Position',[clm(1) row(1) dim]); hold on
        elseif iAn == 2
            axO     = axes('Position',[clm(1) row(3) dim]); hold on
        end
        
        vec1                = allAUCFO{iAn}(~inclIdx);
        vec2                = allAUCFO{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPF, 1, iAn, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PF, 1, iAn, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 
            lg = legend('Unresponsive','Modulated', 'Location','Northwest');
            lg.FontSize = 8;
            legend boxoff
        end
        
        if iAn == 1
            axD     = axes('Position',[clm(1) row(2) dim]); hold on
        elseif iAn == 2
            axD     = axes('Position',[clm(1) row(4) dim]); hold on
        end
        
        vec1                = allAUCFD{iAn}(~inclIdx);
        vec2                = allAUCFD{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPFd, 1, iAn, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PFd, 1, iAn, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end    
        
        %%% COLUMN 2: Coherence
        if iAn == 1 
            axO     = axes('Position',[clm(2) row(1) dim]); hold on
        elseif iAn == 2
            axO     = axes('Position',[clm(2) row(3) dim]); hold on
        end
        
        vec1                = allAUCC{iAn}(~inclIdx);
        vec2                = allAUCC{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPC, 1, iAn, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PC, 1, iAn, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1
            axD     = axes('Position',[clm(2) row(2) dim]); hold on
        elseif iAn == 2 
            axD     = axes('Position',[clm(2) row(4) dim]); hold on
        end
        
        vec1                = allAUCCD{iAn}(~inclIdx);
        vec2                = allAUCCD{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+Yof];
        
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPCd, 1, iAn, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PCd, 1, iAn, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2 && iFi == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
        
        %%% COLUMN 3: HI vs MI
        if iAn == 1
            axO     = axes('Position',[clm(3) row(1) dim]); hold on
        elseif iAn == 2 
            axO     = axes('Position',[clm(3) row(3) dim]); hold on

        end
        
        vec1                = allAUCHM{iAn}(~inclIdx);
        vec2                = allAUCHM{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPHM, 1, iAn, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PHM, 1, iAn, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 
            axD     = axes('Position',[clm(3) row(2) dim]); hold on
        elseif iAn == 2 
            axD     = axes('Position',[clm(3) row(4) dim]); hold on
        end
        
        vec1                = allAUCHMD{iAn}(~inclIdx);
        vec2                = allAUCHMD{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPHMd, 1, iAn, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PHMd, 1, iAn, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
        
        %%% COLUMN 4: MI vs CR
        if iAn == 1 
            axO     = axes('Position',[clm(4) row(1) dim]); hold on
        elseif iAn == 2 
            axO     = axes('Position',[clm(4) row(3) dim]); hold on
        end
        
        vec1                = allAUCMC{iAn}(~inclIdx);
        vec2                = allAUCMC{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axO.YLabel.String   = {'Onset';'No. units'};
        Yof                 = axO.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axO.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axO.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axO.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axO.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axO.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPMC, 1, iAn, nanmedian(vec1), axO.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PMC, 1, iAn, nanmedian(vec2), axO.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axO.XLim(2)/50)
        
        if iAn == 1 
            axD     = axes('Position',[clm(4) row(2) dim]); hold on
        elseif iAn == 2 
            axD     = axes('Position',[clm(4) row(4) dim]); hold on
        end
        
        vec1                = allAUCMCD{iAn}(~inclIdx);
        vec2                = allAUCMCD{iAn}(inclIdx);
        h1                  = histogram(vec1, bin,'Facecolor', 'k','Edgecolor', 'k','Facealpha', alp,'Edgealpha', 0);
        h2                  = histogram(vec2, bin,'Facecolor', 'r','Edgecolor', 'r','Facealpha', alp,'Edgealpha', 0);
        axD.YLabel.String   = {'Decision';'No. units'};
        Yof                 = axD.YLim(2)/5;
        maxV                = max([h1.Values,h2.Values]);
        axD.YLim            = [0 maxV+(2*Yof)];
        
        line([.5 .5],[axD.YLim], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
        fill([nanmedian(vec1)-axD.XLim(2)/50 nanmedian(vec1) nanmedian(vec1)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [0 0 0], 'FaceAlpha', alp)
        fill([nanmedian(vec2)-axD.XLim(2)/50 nanmedian(vec2) nanmedian(vec2)+axD.XLim(2)/50],[maxV+Yof/2 maxV+Yof/4 maxV+Yof/2], [1 0 0], 'FaceAlpha', alp, 'EdgeColor', [1 0 0])
        plotStar(nPMCd, 1, iAn, nanmedian(vec1), axD.XLim(2)/50, maxV+Yof, [0 0 0])
        plotStar(PMCd, 1, iAn, nanmedian(vec2), axD.XLim(2)/50, maxV+Yof, [1 0 0])
%         testDist(vec1,vec2,  maxV+(1.6*Yof), maxV+(1.95*Yof), axD.XLim(2)/50)
        
        if iAn == 2
            axO.XAxis.Visible = 'off';
            axD.XTick = [0 .25 .5 .75 1];
            axD.XLabel.String = 'AUROC';
        else
            axO.XAxis.Visible = 'off';
            axD.XAxis.Visible = 'off';
        end
end

of = .05;
ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(clm(1)-of,.99, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)-of,.99, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)-of,.99, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)-of,.99, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

text(.02,row(2)+.2, 'M1', 'Parent', ax0, 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold','Rotation',90)
text(.02,row(4)+.2, 'M2', 'Parent', ax0, 'FontSize', 20, 'Color', 'k', 'FontWeight', 'bold','Rotation',90)

of = .015;
text(clm(1)+of, .985, 'Figure vs Control', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(2)+of, .985, 'Coherence8 vs Coherence12', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(3)+of, .985, 'Hit vs Miss', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')
text(clm(4)+of, .985, 'Miss vs Correct rejection', 'Parent', ax0, 'FontSize', txtsz+2, 'Color', 'k', 'FontWeight', 'bold')

export_fig([dest_dir 'FIG3_roc_pooled'], '-r400',f);
% exportgraphics(f,[dest_dir '/FIG3_roc_pooled.pdf'],'ContentType','vector','Resolution',300)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'FIG3_roc_pooled'], '-dpdf', '-r400'); 
%% Helper Functions

function out = fdrP(vec)
tmp = fdr(vec);
out(1,:) = tmp(1:2);
out(2,:) = tmp(3:4);
end

function plotStar(Parr, iAn, iFi, xcntr, ofst, ylev, clr)
if Parr(iAn,iFi) < .05
    if Parr(iAn,iFi) < .001
        xstar = [xcntr-ofst xcntr xcntr+ofst];
        ystar = [ylev ylev ylev];
    elseif Parr(iAn,iFi) < .01 && Parr(iAn,iFi) > .001
        xstar = [xcntr-ofst/2 xcntr+ofst/2];
        ystar = [ylev ylev];
    elseif Parr(iAn,iFi) < .05 && Parr(iAn,iFi) > .01
        xstar = xcntr;
        ystar = ylev;
    end
    star = plot(xstar, ystar, '*', 'Color', [clr .6], 'LineWidth', 1.1);
    star.MarkerSize = 6;
end
end

function testDist(vec1,vec2, ylev, ystr, xofst)

% p = anova1([vec1 vec2], [zeros(1,length(vec1)), ones(1,length(vec2))], 'off');
p = ranksum(vec1,vec2);

if p < .05
    line([nanmean(vec1) nanmean(vec2)],[ylev ylev], 'LineWidth', 2, 'Color', [.5 .5 .5])
    xcntr = mean([nanmean(vec1), nanmean(vec2)]);
    
    if p < .001
        xstar = [xcntr-xofst xcntr xcntr+xofst];
        ystar = [ystr ystr ystr];
    elseif p < .01 && p > .001
        xstar = [xcntr-xofst/2 xcntr+xofst/2];
        ystar = [ystr ystr];
    elseif p < .05 && p > .01
        xstar = xcntr;
        ystar = ystr;
    end
    star = plot(xstar, ystar, '*', 'Color', [.5 .5 .5 .6], 'LineWidth', 1.1);
    star.MarkerSize = 6;
end
end