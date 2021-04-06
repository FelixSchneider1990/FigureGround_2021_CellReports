%%% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 2  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except muaeE muaeD lfpE lfpD

% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Eric\Summary\muae.mat')
% muaeE = muae;
% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Dollar\Summary\muae.mat')
% muaeD = muae;
% clear muae

% load('/Volumes/Felix_ExtDrive/Rec/Eric/Summary/muae.mat')
% muaeE = muae;
% load('/Volumes/Felix_ExtDrive/Rec/Dollar/Summary/muae.mat')
% muaeD = muae;
% clear muae

dest_dir    = '/Users/fschneider/ownCloud/NCL_revision/Figures/raw/';
Fs          = 1000;                                                  	% Sampling frequency
T           = 1/Fs;                                                    	% Sampling period
L           = 2250;                                                   	% Length of signal
t           = (0:L-1)*T;                                             	% Time vector
frfft       = Fs*(0:(L/2))/L;
alph       	= .01;

indxR       = 201:400;                                                  % corresponds to -300:-100ms to decision
indx        = 401:600;                                                  % cooresponds to 200:400ms after figure onset
nRep        = 1000;                                                     % No. repetitions
useSEM      = true;

fullCTR = []; FIG = []; CTR = []; AP = []; ML = [];
for iAn = 1:2
    
    if iAn == 1
        dd = muaeE;
    else
        dd = muaeD;
    end
    
    fullCtrEnv = []; mfig = []; mctr = []; co = [];
    c = 0;
    
    for ii = 1:size(dd,2)
        incl = check20Hz(dd{ii});
        if length(dd{ii}.BLslope) < 200 || incl == 0
            continue
        end
        
        datestr     = str2num([dd{ii}.id(1:4) dd{ii}.id(6:7) dd{ii}.id(9:10)]);
        test        = mean([dd{ii}.res.HI12(:,indxR); dd{ii}.res.HI8(:,indxR)],2);
        ctrl        = mean(dd{ii}.res.CR(:,indxR),2);
        pp          = anova1([test; ctrl],[zeros(size(test,1),1);ones(size(ctrl,1),1)], 'off');
        
        if pp < alph && sum(dd{ii}.nTr>=10) == length(dd{ii}.nTr) && datestr - 20190806 <= 0

            %%% Norm: subtract mBL then normalise by average control activity %%%
            mBL       	= nanmean(nanmean(dd{ii}.on.fullAvg(:,101:500),2));                      % Average BL response

            c           = c+1;                                                                      % Counter
            fstim       = length(dd{ii}.figs);                                                      % Ratio of Fig to Gnd stim
            fullCtrEnv(c,:) = nanmean((dd{ii}.on.fullAvg(fstim+1:end,:)) ./mBL);                    % Control period - average across all presented stimuli
            
            mfig(c,:)  	= nanmean((dd{ii}.on.trlOn(dd{ii}.on.cat == 3,1:600)) ./mBL);              	% Figure period - average across all figure trials
            mctr(c,:)   = nanmean((dd{ii}.on.trlOn(dd{ii}.on.cat == 6,1:600)) ./mBL);
            co(c,:)     = dd{ii}.coord;
       
            if iAn == 1 && round(co(c,2)) <= 7 && round(co(c,1)) == 10
                disp('...')
                disp(dd{ii}.id)
                disp(dd{ii}.coord)
                disp(ii)
            end
        end
    end
    
    cc = 1; ap =[]; ml = [];
    for j = 1:size(co,1)
        MLr      	= co(j,2);                              % ML penetration site on grid [mm]
        dep         = co(j,3);                              % depth of recording from GT tip [mm]
        offset      = (dep/sind(90)) * sind(15);          	% calculate offset [mm]
        adj         = MLr - offset;                      	% adjusted ML value [mm]
        ap(cc)      = co(j,1);
        ml(cc)      = adj;
        cc          = cc+1;
    end
    
    fullCTR{iAn}    = fullCtrEnv;
    FIG{iAn}        = mfig;
    CTR{iAn}        = mctr;
    AP{iAn}         = ap;
    ML{iAn}         = ml;
    
    % Bootstrap significance bar
    clear diff bootstat sig
    df                  = FIG{iAn}(:,1:600) - CTR{iAn}(:,1:600);
    %     nRep                = 5000;                                             	% No. of repetitions
    %     [bootstat,~]        = bootstrp(nRep,@mean, df);                             % Bootstrap mean for each timebin
    %     idxF                = sum(bootstat > 0);                                    % Get percentage of samples larger than zero
    %     sig                 = double(idxF < nRep*.005 | idxF > nRep*.995);          % If more than 99.5%/less than 0.5% larger than zero -> significantly different
    for iBin = 1:size(df,2)
        p(iBin) = signrank(df(:,iBin));
    end
    sig                 = double(fdr(p) < alph);
    sig(sig == 0)       = nan;
    SIG{iAn}            = sig;
end

%% PLOT FIG2 %%%
close all
f                   = figure('Units', 'normalized', 'Position', [0 0 .6 1]); set(gcf,'color', [1 1 1]);
ax0                 = axes('Position',[0 0 1 1],'Visible','off');
rows                = [0.43    0.6    0.83];
lw                  = 1.5;

%%% Site coordinates %%%
typ                 = 'muae';
alp                 = .7;

freqStart           = 180;                                          % Tuning low freq [Hz]
steps               = 14;                                           % No of desired tones
PT(1)               = freqStart;                                    % Starting frequency [Hz]
for i = 2:steps                                                 % For no of tones...
    PT(i)           =  PT(i-1)*2^(1/2);                             % 1/2 octave steps
end
frex                = round(PT);
    
for iAn = 1:2
    par = []; mfr_mat = []; mlat_mat = []; cc = 0;
    
    if iAn == 1
        animalID        = 'Eric';
        axC             = axes('Position',[.05 .11 .25 .2]); hold on
        axC.Title.String= 'M1';
    else
        animalID        = 'Dollar';
        axC             = axes('Position',[.22 .11 .25 .2]); hold on
        axC.Title.String        = 'M2';
    end
    
    load([dest_dir 'tMap_' animalID '_' typ  '.mat']);
    imagesc(1:size(mfr_mat,1),-18:-1, flipud(mfr_mat));
    caxis([floor(log(frex(1))) (ceil(max(max(mfr_mat))*10)/10)+.3])

    if iAn == 1
        sc = scatter([ML{iAn}],[-AP{iAn}]);
        axC.YLim        = [-18 -5];
        axC.YLim        = [-18 -5];
    else
        sc              = scatter([ML{iAn}],[-AP{iAn}]);
        axC.YLim        = [-16 -3];
    end
    axC.YAxis.Visible   = 'off';
    axC.XAxis.Visible   = 'off';
    axC.Title.FontSize  = 12;
    
    if iAn == 1
        axC.Title.Position(1) = 8.5;
    else
        axC.Title.Position(1) = 7.5;
    end

    sc.SizeData         = 20;
    sc.Marker           = '^';
    sc.MarkerFaceColor  = [1 0 0];
    sc.MarkerEdgeColor  = [1 0 0];
    
    colormap([[1 1 1]; gray(256)])
    
    cm                  = gray(256);
    axCB                = axes('Position',[.175 .1 .15 .001]);
    axCB.Visible        = 'off';
    colormap(axCB, cm)
    cb                  = colorbar(axCB);
    cb.Location         = 'southoutside';
    cb.Color            = [0 0 0];
    cb.Label.String     = 'Best frequency [Hz]';
    cb.FontSize         = 12;
    caxis([floor(log(frex(1))) ceil(log(frex(11)))])
    cb.Ticks            = [log(frex(1)),log(frex(7)),log(frex(11))];
    cb.TickLabels       = {num2str(frex(1)), num2str(frex(7)), ['>' num2str(frex(11))]};
end

%%% FIGURE ONSET %%%
col = [[0 0 0]; [1 0 0]];
SEM = []; m = []; sd = [];
mID = {'M1','M2'};

for kk = 1:2
    if kk == 1
        axB = axes('Position',[.1 rows(end-1) .3 .15]); hold on
        ofst = .12;
    else
        axB = axes('Position',[.1 rows(end-2) .3 .15]); hold on
        ofst = .12;
    end
    
    if useSEM
        m       = mean(CTR{kk});
        sd      = std(CTR{kk});
        n       = size(CTR{kk},1);
        sem     = sd/(sqrt(n));
        CI      = [m+sem/2; m-sem/2];
    else
        CI      = bootci(nRep,@mean,CTR{kk});
    end
    l           = 1:length(CI);
    len         = [l fliplr(l)];
    ci          = [CI(1,:) fliplr(CI(2,:))];
    fill(len,ci,col(1,:)*0.7,'EdgeColor','none', 'FaceAlpha', .3);
    hold on
    pC = plot(mean(CTR{kk}), 'Color', col(1,:), 'LineWidth', lw);
    plot(1:length(SIG{kk}), SIG{kk}+ofst,  'Color', [.5 .5 .5], 'LineWidth', 3);
    
    if useSEM
        m       = mean(FIG{kk});
        sd      = std(FIG{kk});
        n       = size(FIG{kk},1);
        sem     = sd/(sqrt(n));
        CI      = [m+sem/2; m-sem/2];
    else
        CI      = bootci(nRep,@mean,FIG{kk});
    end
    l           = 1:length(CI);
    len         = [l fliplr(l)];
    ci          = [CI(1,:) fliplr(CI(2,:))];
    fill(len,ci,col(2,:)*0.7,'EdgeColor','none', 'FaceAlpha', .3);
    hold on
    pF = plot(mean(FIG{kk}), 'Color', col(2,:), 'LineWidth', lw);
    
    axB.XLim            = [1 600];
    axB.XTick           = [1 200 400 600];
    axB.XTickLabel      = {'-200','0','200','400'};
    axB.XLabel.String   = 'Time [ms]';
    axB.YLabel.String   = {mID{kk};'MUA [norm]'};
    axB.FontSize        = 14;
    line([200 200],[-1 2], 'LineStyle', '--', 'LineWidth',lw, 'Color','k')
    
    if kk == 1
        axB.YLim = [1.01 1.12];
        of = diff(axB.YLim)*.04;
        lv = axB.YLim(1)+of;
        axB.XAxis.Visible = 'off';
    else
        axB.YLim = [1.01 1.12];
        of = diff(axB.YLim)*.04;
        lv = axB.YLim(1)+of;
    end
    
    seq = [lv lv+of lv+of lv lv lv-of lv-of lv];
    ch  = [0 5 45 50];
    for ic = 1:12
        fi  = fill([ch fliplr(ch)],seq, [0 0 0], 'LineStyle','none');
        ch  = ch+50;
    end
end

lg              = legend([pC, pF],{'CTR','FIG'});
lg.FontSize     = 10;
lg.Position (2) = rows(end-1) + .11;
lg.Position (1) = .12;
lg.Box          = 'off';

%%% CONTROL PLOT
axA         = axes('Position',[.1 rows(end) .85 .15]); hold on
col         = [[0 0 0]; [.65 .65 .65]];

if useSEM
    m       = mean(fullCTR{1});
    sd      = std(fullCTR{1});
    n       = size(fullCTR{1},1);
    sem     = sd/(sqrt(n));
    CI      = [m+sem/2; m-sem/2];
else
    CI      = bootci(nRep,@mean,fullCTR{1});
end

l           = 1:length(CI);
len         = [l fliplr(l)];
ci          = [CI(1,:) fliplr(CI(2,:))];
fill(len,ci,col(1,:)*0.7,'EdgeColor','none', 'FaceAlpha', .3);
hold on
pM1 = plot(mean(fullCTR{1}), 'Color', col(1,:), 'LineWidth', lw);

if useSEM
    m       = mean(fullCTR{2});
    sd      = std(fullCTR{2});
    n       = size(fullCTR{2},1);
    sem     = sd/(sqrt(n));
    CI      = [m+sem/2; m-sem/2];
else
    CI      = bootci(nRep,@mean,fullCTR{2});
end

l           = 1:length(CI);
len         = [l fliplr(l)];
ci          = [CI(1,:) fliplr(CI(2,:))];
fill(len,ci,col(2,:)*0.7,'EdgeColor','none', 'FaceAlpha', .3);
hold on
pM2 = plot(mean(fullCTR{2}), 'Color', col(2,:), 'LineWidth', lw);

axA.YTick           = [1 1.1 1.2];
axA.YLim            = [.97 1.3];
axA.XTick           = [500 2000 3500];
axA.XTickLabel      = {'0','1500','3000'};
axA.XLabel.String   = 'Time [ms]';
axA.YLabel.String   = 'MUA [norm]';
axA.FontSize        = 14;

line([500 500],[0 2], 'LineStyle', '--', 'LineWidth',lw, 'Color','k')
line([3500 3500],[0 2], 'LineStyle', '--', 'LineWidth',lw, 'Color','k')

lg              = legend([pM1, pM2],{'M1','M2'});
lg.FontSize     = 10;
lg.Position (2) = rows(end) + .065;
lg.Position (1) = .12;
lg.Box          = 'off';


axI = axes('Position',[.29 rows(end)+.11 .1 .05]); hold on
plot(mean(fullCTR{1}(:,1900:2100)), 'Color', col(1,:), 'Linewidth', lw);
plot(mean(fullCTR{2}(:,1900:2100)), 'Color', col(2,:), 'Linewidth', lw);
axI.YLim            = [1.02 1.12];
axI.XLim            = [0 200];
axI.XTick           = [0 200];
axI.YTick           = [.1 .2];
axI.XTickLabel      = {'1400','1600'};
axI.XLabel.String   = 'Time [ms]';
axI.YAxis.Visible   = 'off';
axI.FontSize        = 12;
axI.YLim(2)         = 1.10;

of                  = diff(axI.YLim)*.1;
lv                  = axI.YLim(2)-of;
seq                 = [lv lv+of lv+of lv lv lv-of lv-of lv];
ch                  = [0 5 45 50];
for ic = 1:4
    fi      = fill([ch fliplr(ch)],seq, [0 0 0], 'LineStyle','none');
    ch      = ch+50;
end

ax0         = axes('Position',[0 0 1 1],'Visible','off');
annotation('line',[.29 .315],[.993 .993],'LineWidth', 1.5);
text(.26, .994, '50ms', 'Parent', ax0)

axI2 = axes('Position',[.45 rows(end)+.11 .1 .05]); hold on
for k = 1:2
    fftEnv = [];
    for iCell = 1:size(fullCTR{k},1)
        y = fft(fullCTR{k}(iCell,751:3000));
        P2              = abs(y/L);
        P1              = P2(1:L/2+1);
        P1(2:end-1)     = 2*P1(2:end-1);
        fftEnv(iCell,:) = P1;
    end
    plot(frfft,mean(fftEnv)./ max(mean(fftEnv)), 'Color', col(k,:), 'LineWidth', lw)
    box off
    axI2.XLim           = [1 50];
    axI2.YLim           = [0 .04];
    axI2.YTick          = [];
    axI2.YLabel.String  = 'Power';
    axI2.XLabel.String  = 'Freq [Hz]';
    axI2.FontSize       = 12;
end

%% Get data

clearvars -except muaeE muaeD lfpE lfpD indxR indx typ dest_dir alph rows f

for iAn = 1:2
    
    ant = []; pos = []; muae = [];
    
    if iAn == 1
        animalID = 'Eric';
    else
        animalID = 'Dollar';
    end
    
    load([dest_dir 'tMap_' animalID '_' typ  '.mat']);
    APmap   = find(logical(sum(~isnan(mfr_mat),2)));
    MLmap   = find(logical(sum(~isnan(mfr_mat))));
    [x,y]   = coreBoundary(mfr_mat,APmap,MLmap,false,animalID);
    
    if strcmp(animalID, 'Eric')
        muae = muaeE;
    elseif strcmp(animalID, 'Dollar')
        muae = muaeD;
    end
    
    [ant, pos] = sortUnits(x,y,muae);
    
    for iFi = 1:2
        
        if iFi == 1
            d = ant;
            str = 'ant';
        else
            d = pos;
            str = 'pos';
        end
        
        c = 0; cc = 0; mtar8 = []; mtar12 = []; mctr = []; f8 = []; f12 = []; latency12 = []; latency8 = []; dFG8 = []; dFG12 = []; dCoh = [];
        
        for ii = 1:size(d,2)
            incl = check20Hz(d{ii});
            if length(d{ii}.BLslope) < 200 || incl == 0
                continue
            end
            
            datestr     = str2num([d{ii}.id(1:4) d{ii}.id(6:7) d{ii}.id(9:10)]);
            test        = mean([d{ii}.res.HI12(:,indxR); d{ii}.res.HI8(:,indxR)],2);
            ctrl        = mean(d{ii}.res.CR(:,indxR),2);
            pp          = anova1([test; ctrl],[zeros(size(test,1),1);ones(size(ctrl,1),1)], 'off');
        
            if pp < alph && sum(d{ii}.nTr>=10) == length(d{ii}.nTr) && datestr - 20190806 <= 0
                
                %%% Norm: subtract mBL then normalise by average control activity %%%
                mBL             = nanmean(nanmean(d{ii}.on.fullAvg(:,101:500),2));                          % Average BL response    
                c               = c+1;                                                                      % Counter
                fstim           = size(d{ii}.figfreqs,3);                                                   % Ratio of Fig to Gnd stim

                HIo8            = d{ii}.on.trlOn(d{ii}.on.cat == 3 & d{ii}.on.coh == 8,indx) ./mBL;
                HIo12           = d{ii}.on.trlOn(d{ii}.on.cat == 3 & d{ii}.on.coh == 12,indx) ./mBL;
                CRo             = d{ii}.on.trlOn(d{ii}.on.cat == 6,indx) ./mBL;

                mtar8(c,:)  	= nanmean(d{ii}.on.trlOn(d{ii}.on.cat == 3 & d{ii}.on.coh == 8,1:600) ./mBL);      	% Figure period - average across all coherence8 stimuli
                mtar12(c,:)  	= nanmean(d{ii}.on.trlOn(d{ii}.on.cat == 3 & d{ii}.on.coh == 12,1:600) ./mBL);      % Figure period - average across all coherence12 stimuli
                mctr(c,:)    	= nanmean(d{ii}.on.trlOn(d{ii}.on.cat == 6,1:600) ./mBL);                          	% Control period - average across all presented stimuli

                dCoh(c)         = mean((nanmean(HIo12) - nanmean(HIo8)) ./ nanstd([HIo8; HIo12]));                  % neural d-prime
                dFG8(c)         = mean((nanmean(HIo8) - nanmean(CRo)) ./ nanstd([HIo8; CRo]));
                dFG12(c)        = mean((nanmean(HIo12) - nanmean(CRo)) ./ nanstd([HIo12; CRo]));
            
                latency8(c,:)   = d{ii}.sigT8;
                latency12(c,:)  = d{ii}.sigT12;
                
                if ~isnan(sum(sum(d{ii}.RF)))
                    cc  = cc+1;
                    no  = [];
                    no  = getFreqRF(d{ii}.RF, d{ii}.figfreqs);
                    
                    f8  = vertcat(f8, no(1:fstim/2,:));
                    f12 = vertcat(f12, no(fstim/2+1:fstim,:));
                end
            end
        end
        
        mt8{iAn,iFi}    = mtar8;
        mt12{iAn,iFi}   = mtar12;
        mc{iAn,iFi}     = mctr;
        FG8{iAn,iFi}    = dFG8;
        FG12{iAn,iFi}   = dFG12;
        COH{iAn,iFi}    = dCoh;
        lat8{iAn,iFi}   = latency8;
        lat12{iAn,iFi}  = latency12;
        F8{iAn,iFi}     = f8;
        F12{iAn,iFi}    = f12;
    end
end

for iAn = 1:2
    for iFi = 1:2
        p8(iAn,iFi)   = signrank(FG8{iAn, iFi});
        p12(iAn,iFi)  = signrank(FG12{iAn, iFi});
        pC(iAn,iFi)   = signrank(COH{iAn, iFi});
    end
end

pC  = fdr(pC);
p8  = fdr(p8);
p12 = fdr(p12);

%% PLOT ANT vs POS %%%

lw      = 1.5;
col     = [0 0 0; 0 .9 0; .9 0 0];
dist    = 'sem';
dim     = [.2 .15];
clm     = [.53 .78];
rowBox  = [.05 .21];


for iAn = 1:2
    for iFi = 1:2
        
        if iAn == 1
            animalID = 'Eric';
            sigY        = [.4 .4;.475 .475;.55 .55];
            sigX        = [1.1 1.9; 1.1 2.9; 2.1 2.9];
            ofst        = .035;
            Yof         = .035;
        else
            animalID = 'Dollar';
            sigY        = [.27 .27;.32 .32;.37 .37];
            sigX        = [1.1 1.9; 1.1 2.9; 2.1 2.9];
            ofst        = .035;
            Yof         = .025;
        end
        
        %%% TIME TRACE %%%
        if iFi == 1 && iAn == 1
            axA = axes('Position',[clm(1) rows(2) dim]); hold on
            axA.Title.String = 'ANT';
            axA.XAxis.Visible = 'off';
            str = 'ant';
        elseif iFi == 2 && iAn == 1
            axA = axes('Position',[clm(2) rows(2)  dim]); hold on
            axA.Title.String = 'POS';
            axA.XAxis.Visible = 'off';
            str = 'pos';
        elseif iFi == 1 && iAn == 2
            axA = axes('Position',[clm(1) rows(1) dim]); hold on
            str = 'ant';
        elseif iFi == 2 && iAn == 2
            axA = axes('Position',[clm(2) rows(1) dim]); hold on
            str = 'pos';
        end
        
        [pc0,pc8,pc12] = plotFRsumm(mt8{iAn,iFi}(:,1:600),mt12{iAn,iFi}(:,1:600),mc{iAn,iFi}(:,1:600),str, dist,animalID, 0, col);
        axA.XLim = [1 600];
        
        if iAn == 1
            axA.YLabel.String = {'M1','MUA [norm]'};
        else
            axA.YLabel.String = {'M2','MUA [norm]'};
        end
        
        if iFi == 2
            axA.YLabel.String = [];
        end
        
        of      = diff(axA.YLim)*.04;
        lv      = axA.YLim(1)+of;
        seq     = [lv lv+of lv+of lv lv lv-of lv-of lv];
        ch      = [0 5 45 50];
        for ic = 1:12
            fi  = fill([ch fliplr(ch)],seq, [0 0 0], 'LineStyle','none');
            ch  = ch+50;
        end
        
        
        %%% BAR PLOT %%%
        if iFi == 1 && iAn == 1
            axB1 = axes('Position',[clm(1) rowBox(2)  dim]); hold on
            axB1.XAxis.Visible = 'off';
            axB1.YLabel.String = {'M1','d-prime'};
        elseif iFi == 2 && iAn == 1
            axB1 = axes('Position',[clm(2) rowBox(2)  dim]); hold on
            axB1.XAxis.Visible = 'off';
        elseif iFi == 1 && iAn == 2
            axB1 = axes('Position',[clm(1) rowBox(1)  dim]); hold on
            axB1.YLabel.String = {'M2','d-prime'};
        elseif iFi == 2 && iAn == 2
            axB1 = axes('Position',[clm(2) rowBox(1)  dim]); hold on
        end
        
        arr = []; mat = [];
        arr = [ones(size(FG8{iAn,iFi},2),1); ones(size(FG12{iAn,iFi},2),1)+1;ones(size(COH{iAn,iFi},2),1)+2];
        mat = [FG8{iAn,iFi}';FG12{iAn,iFi}';COH{iAn,iFi}'];
        bx  = boxplot(mat, arr, 'Colors', [0 0 0]);
        set(bx(end,:),'Visible','off')
        set(bx, {'linew'},{lw})
        box off  
        
        hgt = .58;
        ofs = .1;
        for iB = 1:3
            dpvec = [];
            if iB == 1
                dpvec = FG8{iAn,iFi};
                p = p8(iAn,iFi);
            elseif iB == 2
                dpvec = FG12{iAn,iFi};
                p = p12(iAn,iFi);
            elseif iB == 3
                dpvec = COH{iAn,iFi};
                p = pC(iAn,iFi);
            end
            
            x   = [.8 1.8 2.8];
            xx  = x(iB) + ((x(iB)+.4)-x(iB)).*rand(1,sum(arr == iB));
            sc  = scatter(xx, mat(arr == iB)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
            if iB == 1
                sc.MarkerFaceColor = [0 1 0];
                sc.MarkerEdgeColor = [0 1 0];
            elseif iB == 2
                sc.MarkerFaceColor = [1 0 0];
                sc.MarkerEdgeColor = [1 0 0];
            elseif iB == 3
                sc.MarkerFaceColor = [.5 .5 .5];
                sc.MarkerEdgeColor = [.5 .5 .5];
            end
            sc.SizeData = 20;
            
            if p < .05
                if p < .001
                    xstar = [iB-ofs iB iB+ofs];
                    ystar = [repmat(hgt,1,3)];
                elseif p < .01 && p > .001
                    xstar = [iB-ofs/2 iB+ofs/2];
                    ystar = [repmat(hgt,1,2)];
                elseif p < .05 && p > .01
                    xstar = iB;
                    ystar = hgt;
                end
                star = plot(xstar, ystar, 'k*', 'LineWidth', 1.1);
                star.MarkerSize = 8;
            end
        end
 
        axB1.XTick              = [1 2 3];
        axB1.XTickLabel         = {'Coh8-CTR', 'Coh12-CTR', 'Coh12-Coh8'};
        axB1.XTickLabelRotation = 12;
        axB1.FontSize           = 14;
        axB1.YLim               = [-.3 .7];
 
        if iFi == 1 && iAn == 1
            axB1.Position = [clm(1) rowBox(2)  dim];
        elseif iFi == 2 && iAn == 1
            axB1.Position = [clm(2) rowBox(2)  dim];
        elseif iFi == 1 && iAn == 2
            axB1.Position = [clm(1) rowBox(1)  dim];
        elseif iFi == 2 && iAn == 2
            axB1.Position = [clm(2) rowBox(1)  dim];
        end
    end
end

lg              = legend([pc0,pc8,pc12],{'CTR','Coh8','Coh12'});
lg.FontSize     = 10;
lg.Position(2)  = rows(1)+.075;
lg.Position(1)  = clm(1)+.068;
lg.Box          = 'off';

ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(0,.99, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0,.765, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.44,.765, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.44,.355, 'e', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0,.355, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

% addpath X:\Felix\Scripts\Stuff\export_fig-master
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
export_fig([dest_dir 'FIG3'], '-r400',f);
% exportgraphics(f,[dest_dir 'FIG2.pdf'],'ContentType','vector','Resolution',300)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'FIG3'], '-dpdf', '-r400'); 

