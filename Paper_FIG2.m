%% SFG SU example spiking response
% addpath \\campus\rdw\ion02\02\auditory\Felix\Documents\Publications\FigGnd_Ephys\Figures
% load('Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Eric\2018-10-09_10-04-30\Data\DataStruct_2018-10-09.mat')    
% cd Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\Eric\2018-10-09_10-04-30
% load('/Volumes/Felix_ExtDrive/Rec/Eric/2018-10-09_10-04-30/Data/DataStruct_2018-10-09.mat')
% cd /Volumes/Felix_ExtDrive/Rec/Eric/2018-10-09_10-04-30

chan        = 'ch2';
clus        = 'clus1';
% data.(chan).Spks      = getSpikesSU(data, str2num(chan(end)));           % Import spike file
stID        = data.behaviour.stimID;
s           = data.behaviour.stimNrPool;
on          = data.behaviour.figOn*1000+500;
hi          = data.behaviour.resCat == 3;
cr          = data.behaviour.resCat == 6;
row         = linspace(.05,.85,6);
avgSFG      = nanmean(data.(chan).Spks.(clus).SFG.psth);
nV          = mean(avgSFG(101:500));
bl          = 0;
% bl          = mean(avgSFG(101:500));
% nV          = max(avgSFG(500:700)) - bl;

avgSFGenv  	= nanmean(data.(chan).MUAe.SFG.sig);
blE       	= 0;
nVE        	= mean(avgSFGenv(101:500));
% blE           = mean(avgSFGenv(101:500));
% nVE           = max(avgSFGenv(500:700)) - blE;

figs        = data.stimSpecs.fig;
frMat       = data.stimSpecs.freq_mat;
perc        = .88;


cm          = [[0 0 0]; [0 .9 0]; [.9 0 0]];

f = figure('Units', 'normalized', 'Position', [0 0 1 1]);
set(gcf,'color', [1 1 1]);
axis off
ax0 = axes('Position',[0 0 1 1],'Visible','off');

%%% Get control data %%%
for iCtr = 1:8
    arr = [];
    idxC = s == stID(12+iCtr) & cr == 1;
    arr = (data.(chan).Spks.(clus).SFG.psth(idxC,:)-bl) ./nV;
    cmat(iCtr,:) = nanmean(arr);
    scmat(iCtr,:) = nanstd(arr);
    cmatE(iCtr,:) = nanmean((data.(chan).MUAe.SFG.sig(idxC,:)-blE) ./nVE);
end

ind = data.behaviour.resCat < 5;
fo = unique(data.behaviour.figOn(ind));
for iFO = 1:length(fo)
    start = (fo(iFO)*1000)+500;
    tmpSU(:,iFO) = mean(cmat(:,start:start+199),2);
    tmpMU(:,iFO) = mean(cmatE(:,start:start+199),2);
end
cFR = mean(tmpSU,2);
cFRe = mean(tmpMU,2);

lv = 2.45;
of = .05;
seq = [lv lv+of lv+of lv lv lv-of lv-of lv];
        
%%% Plot stimulus specific response %%%
for iSt = 1:length(stID)*.6
    psth = []; muae = [];
    
    idx = s == stID(iSt);
    fon = unique(on(idx));
    mat = data.(chan).Spks.(clus).SFG.psth(idx,:);
    matE = data.(chan).MUAe.SFG.sig(idx,:);
    spT = data.(chan).Spks.(clus).SFG.spTimes(idx);
    chordOn = (fon-500)/50;
    figfreqs(:,:,iSt) = frMat(:,chordOn+1:chordOn+8,iSt);
    
    if iSt < 7
        axRF = axes('Position',[row(iSt) .885 .13 .1]); hold on
        axRF.XLim = [fon-199 fon+400];
        axRF.XTick = [fon-199 fon fon+200 fon+400];
        axRF.YTick = [8 16];
        axRF.YLim = [0 17];
        axRF.XColor = [0 0 0];
        axRF.YColor = [0 0 0];
        axRF.XAxis.Visible = 'off';
        axRF.FontSize = 14;
        ex = isnan(sum(mat,2));
        spT(ex) = [];
        
        for iTr = 1:length(spT)
            [x,y] = computeRaster(spT{iTr}',iTr);
            plot(x,y,'Color','k');
        end
        
        axF = axes('Position',[row(iSt)  .68 .13 .2]); hold on
        psth = (mat(:,fon-199:fon+400) - bl) ./ nV;
        mFIG = nanmean(psth);
        sdFIG = nanstd(psth);
        semFIG = sdFIG/(sqrt(size(psth,1)));
        mFR(iSt,1) = mean(mFIG(401:600));
%         mFR(iSt,2) = mean(mFIG(401:600));
        
        muae = (matE(:,fon-199:fon+400) - blE) ./ nVE;
        mFIGe = nanmean(muae);
        sdFIGe = nanstd(muae);
        semFIGe = sdFIGe/(sqrt(size(muae,1)));
        mMUA(iSt,1) = mean(mFIGe(401:600));
%         mMUA(iSt,2) = mean(mFIGe(401:600));
        
        seq = [lv lv+of lv+of lv lv lv-of lv-of lv];
        ch = [0 5 45 50];
        for ic = 1:12
            fi  = fill([ch fliplr(ch)],seq, [0 0 0], 'LineStyle','none');
            ch = ch+50;
        end
        
        yyaxis left
        fi = fill([1:600 600:-1:1],[mFIGe+semFIGe/2 fliplr(mFIGe-semFIGe/2)], [.5 .5 .5], 'LineStyle','none');
        fi.FaceAlpha = .3;
        pEnv = plot(mFIGe, 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '-');
        axF.YTick = [0 1 2];
        axF.YLim = [.5 2.5];
        axF.YLabel.String = 'MUA envelope [norm]';
        axF.YAxis(1).Color = [.5 .5 .5];

        yyaxis right
        fi1 = fill([1:600 600:-1:1],[mFIG+semFIG/2 fliplr(mFIG-semFIG/2)], [0 0 0], 'LineStyle','none');
        fi1.FaceAlpha = .3;
        pSpks = plot(mFIG, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-');
        axF.YTick = [0 5 10];
        axF.YLim = [-1 10];
        axF.YAxis(2).Color = [0 0 0];

        axF.YLabel.String = 'Firing Rate [norm]';
        axF.XLim = [1 600];
        axF.XTick = [1 200 400 600];
        axF.XTickLabel = {- 200 0 200 400};
        axF.XColor = [0 0 0];
        axF.FontSize = 14;
        axF.YAxis(1).Label.FontSize = 14;
        axF.YAxis(2).Label.FontSize = 14;
        box off
        l = line([200 200],[-2 axF.YLim(2)]);
        l.Color = cm(2,:);
        l.LineStyle = '--';
        l.LineWidth = 2;
        tt = text(45,axF.YLim(2)*perc,'Coh8', 'Parent', axF, 'FontSize', 14, 'Color', cm(2,:));
        
        if iSt ~= 1
            axRF.YAxis.Visible = 'off';
            axF.YAxis(1).Visible = 'off';
            axF.YAxis(2).Visible = 'off';
            axF.XAxis.Visible = 'off';
        else
            axRF.YLabel.String = 'Trials';
            axF.XLabel.String = 'Time [ms]';
            
            lg = legend([pEnv, pSpks],{'Envelope','Spikes'});
            lg.FontSize     = 10;
            lg.Position(2)  = .82;
            lg.Position(1)  = .12;
            lg.Box          = 'off';
        end
        
    else
        axRF = axes('Position',[row(iSt-6) .525 .13 .1]); hold on
        axRF.XLim = [fon-199 fon+400];
        axRF.XTick = [fon-199 fon fon+200 fon+400];
        axRF.YTick = [8 16];
        axRF.YLim = [0 17];
        axRF.XColor = [0 0 0];
        axRF.YColor = [0 0 0];
        axRF.XAxis.Visible = 'off';
        axRF.FontSize = 14;
        ex = isnan(sum(mat,2));
        spT(ex) = [];
        
        for iTr = 1:length(spT)
            [x,y] = computeRaster(spT{iTr}',iTr);
            plot(x,y,'Color','k');
        end
                
        axF = axes('Position',[row(iSt-6)  .32 .13 .2]); hold on
        psth = (mat(:,fon-199:fon+400) - bl) ./ nV;
        mFIG = nanmean(psth);
        sdFIG = nanstd(psth);
        semFIG = sdFIG/(sqrt(size(psth,1)));
        mFR(iSt,1) = mean(mFIG(401:600));
%         mFR(iSt,2) = mean(mFIG(401:600));
        
        muae = (matE(:,fon-199:fon+400) - blE) ./ nVE;
        mFIGe = nanmean(muae);
        sdFIGe = nanstd(muae);
        semFIGe = sdFIGe/(sqrt(size(muae,1)));     
        mMUA(iSt,1) = mean(mFIGe(401:600));
%         mMUA(iSt,2) = mean(mFIGe(401:600));
        
        ch = [0 5 45 50];
        for ic = 1:12
            fi  = fill([ch fliplr(ch)],seq, [0 0 0], 'LineStyle','none');
            ch = ch+50;
        end
        
        yyaxis left
        fi = fill([1:600 600:-1:1],[mFIGe+semFIGe/2 fliplr(mFIGe-semFIGe/2)], [.5 .5 .5], 'LineStyle','none');
        fi.FaceAlpha = .3;
        plot(mFIGe, 'Color', [.5 .5 .5], 'LineWidth', 2, 'LineStyle', '-');
        axF.YTick = [0 1 2];
        axF.YLim = [.5 2.5];

        yyaxis right
        fi1 = fill([1:600 600:-1:1],[mFIG+semFIG/2 fliplr(mFIG-semFIG/2)], [0 0 0], 'LineStyle','none');
        fi1.FaceAlpha = .3;
        plot(mFIG, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-');
        axF.YTick = [0 5 10];
        axF.YLim = [-1 10];
        
        axF.XLim = [1 600];
        axF.XTick = [1 200 400 600];
        axF.XTickLabel = {- 200 0 200 400};
        axF.XColor = [0 0 0];
        axF.FontSize = 14;
        box off
        l = line([200 200],[-2 axF.YLim(2)]);
        l.Color = cm(3,:);
        l.LineStyle = '--';
        l.LineWidth = 2;
        tt = text(40,axF.YLim(2)*perc,'Coh12', 'Parent', axF, 'FontSize', 14, 'Color', cm(3,:));
        
        axRF.YAxis.Visible = 'off';
        axF.YAxis(1).Visible = 'off';
        axF.YAxis(2).Visible = 'off';
        axF.XAxis.Visible = 'off';
    end
    pks{iSt} = findpeaks(mFIG);
end

%%% FFT of avg stimulus %%%
Fs          = 1000;                                                  	% Sampling frequency
T           = 1/Fs;                                                    	% Sampling period
L           = 2750;                                                   	% Length of signal
t           = (0:L-1)*T;                                             	% Time vector
frfft       = Fs*(0:(L/2))/L;                                           % Frequency vector
y               = fft(mean(data.(chan).Spks.(clus).SFG.alig.on.fullAvg(:,751:3500)));                        % One-sided fast Fourier transform for signal and BL
% y               = fft(mean(data.(chan).MUAe.SFG.alig.on.fullAvg(:,751:3500)));                        % One-sided fast Fourier transform for signal and BL
P2              = abs(y/L);
P1              = P2(1:L/2+1);
P1(2:end-1)     = 2*P1(2:end-1);

axFFT = axes('Position',[row(1) .06 .15 .2]); hold on
plot(frfft,P1./ max(P1), 'k', 'LineWidth', 2)
axFFT.XColor = [0 0 0];
axFFT.YColor = [0 0 0];
axFFT.XLim = [0 50];
axFFT.YLim = [0 .2];
axFFT.YTick = [0 .1 .2];
axFFT.FontSize = 14;
axFFT.YLabel.String = 'Power [norm]';
axFFT.XLabel.String = 'Frequency [Hz]';
fill([19 20 21], [.15 .14 .15], 'k')
box off

%%% Mean FR in 1st and 2nd window %%%
cc = zeros(size(cFR,1),1);
ff = ones(size(mFR,1)/2,1);
axB = axes('Position',[row(1)+.18 .07 .12 .2]);
hold on
bx = boxplot([cFR(:,1); mFR(1:6,1); mFR(7:12,1)], ...
    [cc; ff+1; ff+2;],...
    'Colors', cm);
set(bx,'MarkerEdgeColor',cm(2,:))
set(bx, {'linew'},{2})
sc1 = scatter(ones(8,1),cFR(:,1));
sc1.MarkerFaceColor = [0 0 0];
sc1.MarkerEdgeColor = [0 0 0];
sc2 = scatter(ones(6,1)+1,mFR(1:6,1));
sc2.MarkerFaceColor = cm(2,:);
sc2.MarkerEdgeColor = cm(2,:);
sc3 = scatter(ones(6,1)+2,mFR(7:12,1));
sc3.MarkerFaceColor = cm(3,:);
sc3.MarkerEdgeColor = cm(3,:);
axB.XTick = [1 2 3];
axB.XTickLabel = {'Control', 'Coh8', 'Coh12'};
xtickangle(20)
axB.XLim = [0 4];
axB.YLim = [0 5.5];
axB.YTick = [0 2.5 5];
axB.XColor = [0 0 0];
axB.YColor = [0 0 0];
axB.FontSize = 14;
axB.Title.String = 'Threshold';
axB.YLabel.String = 'MUA [norm]';
axB.XAxis.FontSize = 14;
box off

axBB = axes('Position',[row(1)+.31 .06 .12 .2]);
hold on
bx = boxplot([cFRe(:,1); mMUA(1:6,1); mMUA(7:12,1);], ...
    [cc; ff+1; ff+2;],...
    'Colors', cm);
set(bx, {'linew'},{2})
sc1 = scatter(ones(8,1),cFRe(:,1));
sc1.MarkerFaceColor = [0 0 0];
sc1.MarkerEdgeColor = [0 0 0];
sc2 = scatter(ones(6,1)+1,mMUA(1:6,1));
sc2.MarkerFaceColor = cm(2,:);
sc2.MarkerEdgeColor = cm(2,:);
sc3 = scatter(ones(6,1)+2,mMUA(7:12,1));
sc3.MarkerFaceColor = cm(3,:);
sc3.MarkerEdgeColor = cm(3,:);
axBB.XTickLabel = {'Control', 'Coh8', 'Coh12'};
xtickangle(20)
axBB.XTick = [1 2 3];
axBB.YLim = [.9 1.5];
axBB.XLim = [0 4];
axBB.YTick = [1 1.5];
axBB.XColor = [0 0 0];
axBB.YColor = [0 0 0];
axBB.FontSize = 14;
axBB.Title.String = 'Envelope';
axBB.XAxis.FontSize = 14;
% axBB.YAxis.Visible = 'off';
box off

%%% RF plot %%%
evt   	= data.trials.Tun.evt;                                      % Trial events
ev   	= data.evCodes;                                             % Event codes
stimID 	= evt{1}(evt{1} >= ev.PT(1) & evt{1} <= ev.PT(2));          % Extract stimulus sequence
s     	= unique(stimID);
mat     = data.(chan).Spks.(clus).Tun.PT.psth;
blne    = mean(nanmean(mat(:,1:200)));

for iFreq   = 1:length(s)                                         	% For each frequency (-band)...
    fidx     = stimID == s(iFreq);                                  % Get volume information
    RF(iFreq) = mean(nanmean(mat(fidx,201:300)) - blne);            % PSTH
end

rf          = RF./max(RF);
ff          = fit((1:14)', rf', 'smoothingspline', 'SmoothingParam', .95);
hm          = max(rf) - (diff([min(rf) max(rf)]) / 2);
axRF        =  axes('Position',[row(1)+.73 .06 .2 .2]); hold on
% plot(rf, 'k', 'LineWidth', 2)
plot(ff(1:14), 'k', 'LineWidth', 2)
axRF. XLim = [1 14];
axRF.XTick = [1 7 14];
axRF.XTickLabel = {'180', '1440', '16292'};
axRF.YLim = [-.2 1];
axRF.YTick = [0 .5 1];
axRF.YLabel.String = 'MUA [norm]';
axRF.XLabel.String = 'Frequency [Hz]';
axRF.XColor = [0 0 0];
axRF.YColor = [0 0 0];
axRF.FontSize = 14;
box off
l = line([1 14],[hm hm]);
l.Color = [.5 .5 .5];
l.LineStyle = '--';
l.LineWidth = 2;
vec0 = [];
x0 = [];
for i = 2:length(rf)
    vec0 = [vec0, linspace(ff(i-1), ff(i),100)];
    x0 = [x0, linspace(i-1, i,100)];
end
vec1 = vec0;
vec0(vec0 > hm) = hm;
p = patch([x0 fliplr(x0)], [vec1,fliplr(vec0)], [.5 .5 .5]);
p.FaceAlpha = .5;
p.EdgeColor = 'none';

%%% Figure chords: No frequency elements in RF %%%
no = getFreqRF(RF, figfreqs);
f8 = no(1:6, :);
f12 = no(7:12, :);
axFE1 = axes('Position',[row(1)+.48 .06 .2 .2]); hold on

for i = 1:size(f8,1)
    lf8(i,:) = polyfit(1:8,f8(i,:),1);
    plot(0:9,polyval(lf8(i,:),0:9),'Color',[.5 .5 .5 .75], 'Marker','none')
end

for i = 1:size(f12,1)
    lf12(i,:) = polyfit(1:8,f12(i,:),1);
    plot(0:9,polyval(lf12(i,:),0:9),'Color',[.5 .5 .5 .75], 'Marker','none')
end

for h = 1:8
sc = scatter(zeros(6,1)+h+.1, f12(:,h));
sc.MarkerFaceColor = cm(3,:);
sc.MarkerFaceAlpha = .5;
sc.MarkerEdgeColor = 'none';
end

for h = 1:8
sc = scatter(zeros(6,1)+h-.1, f8(:,h));
sc.MarkerFaceColor = cm(2,:);
sc.MarkerFaceAlpha = .5;
sc.MarkerEdgeColor = 'none';
end

axFE1.FontSize = 14;
axFE1.YTick = [5 10 15];
axFE1.YLim = [0 15];
axFE1.XLim = [0 9];
axFE1.XTick = [1:8];
axFE1.XColor = [0 0 0];
axFE1.YColor = [0 0 0];
axFE1.YLabel.String = 'Frequency Elements in RF';
axFE1.XLabel.String = 'Chords from Figure Onset';
box off

text(0,.98, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0,.3, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0.20,.3, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.48,.3, 'd', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.74,.3, 'e', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

% addpath X:\Felix\Scripts\Stuff\export_fig-master
% dest_dir = 'X:\Felix\Documents\Publications\FigGnd_Ephys\Figures\';
% export_fig([dest_dir 'SFIG_SU'], '-r385',f);
% print([dest_dir 'SFIG_SU'], '-dpng', '-r400')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-8b0ba13
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
export_fig([dest_dir 'SFIG_SU'], '-r400',f);
% % exportgraphics(f,[dest_dir 'SFIG_SU.pdf'],'ContentType','vector','Resolution',300)
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'SFIG_SU'], '-dpdf', '-r400'); 