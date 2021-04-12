%%% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SFIG 2 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/fschneider/Documents/GitHub/FigureGround_Ephys_Analysis
load('PATH/Channel20Data.mat')

clus            = 'clus2';
evt             = tmp.evt;
ev              = tmp.ev;
mat             = tmp.Spks.(clus).Tun.CLK.psth;
lfp             = tmp.LFP.Tun.CLK.sig;
spT             = tmp.Spks.(clus).Tun.CLK.spikeTimes;
clk_EV          = evt{3}(evt{3} >= ev.CLK(1) & evt{3} <= ev.CLK(2));       	% Event codes
cFreq           = unique(clk_EV);                                         	% [25 50 75 100] Hz
cm              = [[.9 0 0]; [0 0 0]; [0 0 .9]];

for iCl = 1:length(cFreq)
    clkSpT{iCl} = spT(clk_EV == cFreq(iCl));
    clkLFP{iCl} = lfp(clk_EV == cFreq(iCl),:);
    clkLFP{iCl} = exludeTrials(clkLFP{iCl}, [], 5, 'PEAK', 'NAN', 1, 1000);
end

f               = figure('Units', 'normalized', 'Position', [0 0 .6 .7]); set(gcf,'color', [1 1 1]);
ax0             = axes('Position',[0 0 1 1],'Visible','off');
dim             = size(clkLFP);
maxVLFP         = 400;
row             = linspace(.09,.72, 4);
freq            = [25 50 75 100];
fs              = 44100;
dur             = 0.2;                                                      % [s]
width           = .0002;                                                    % [s]
pool            = [];
Fs              = 1000;                                                  	% Sampling frequency
T               = 1/Fs;                                                    	% Sampling period
L               = 200;                                                   	% Length of signal
t               = (0:L-1)*T;                                             	% Time vector
frq             = Fs*(0:(L/2))/L;                                           % Frequency vector
cc              = 0;

for iFreq = 1:dim(2)
    
    axR         = axes('Position',[row(iFreq) .76 .2 .2]);
    axR.XLim    = [1 500];
    hold on
    
    pool        = [];
    ex          = logical(sum(isnan(clkLFP{iFreq}),2));
    spkT        = clkSpT{iFreq};
    spkT(ex)    = [];
    
    for iTr = 1:length(spkT)
        [x,y]           = computeRaster(spkT{iTr}',iTr+cc);
        p1              = plot(x,y,'Color','k', 'Parent', axR);
        pool            = [pool;spkT{iTr}];
        
        y               = fft(clkLFP{iFreq}(iTr,201:400));                	% One-sided fast Fourier transform for signal and BL
        P2              = abs(y/L);
        P1              = P2(1:L/2+1);
        P1(2:end-1)     = 2*P1(2:end-1);
        
        yy              = fft(clkLFP{iFreq}(iTr,1:200));
        PP2             = abs(yy/L);
        PP1             = PP2(1:L/2+1);
        PP1(2:end-1)	= 2*PP1(2:end-1);
        fftSig(iTr,:)  	= P1;                                              	% Trialwise FFT
        fftBL(iTr,:)   	= PP1;
    end
    
    if iFreq == 1
        axR.YLabel.String   = 'Trials';
        axR.YTick           = [20 40 60];
        axR.FontSize        = 14;
        axR.XAxis.Visible   = 'off';
        box off
    else
        axR.XAxis.Visible   = 'off';
        axR.YAxis.Visible   = 'off';
    end
    
    ax                      = axes('Position',[row(iFreq) .355 .2 .4]);
    fi                      = fill([201 201 400 400],[0 50 50 0], [.8 .8 .8], 'LineStyle','none'); hold on
    ax.XLim                 = [1 500];
    ax.YLim                 = [0 50];
    ax.YTick                = [0 12 24 36 48];
    ax.YTickLabel           = {'0' num2str((12/60)*500) num2str((24/60)*500) num2str((36/60)*500) num2str((48/60)*500)};
    hh                      = histogram(pool, 250);
    hh.FaceColor            = [0 0 0];
    hh.FaceAlpha            = 1;
    
    if iFreq == 1
        ax.XColor           = [0 0 0];
        ax.YColor           = cm(1,:);
        ax.XTick            = [200 400];
        ax.XTickLabel       = {'0' '200'};
        ax.XLabel.String    = 'Time [ms]';
        ax.YLabel.String    = 'Firing rate [Hz]';
        ax.FontSize         = 14;
        box off
    else
        ax.YAxis.Visible    = 'off';
        ax.XAxis.Visible    = 'off';
    end
    
    yyaxis right
    fp                      = nanmean(clkLFP{iFreq});
    p1                      = plot(fp, 'LineWidth', 2, 'Color', cm(1,:), 'Parent', ax);
    if iFreq == 4
        ax.YColor           = cm(1,:);
        ax.FontSize         = 14;
        ax.YTick            = [-200 0 200 400];
        ax.YLim             = [-250 400];
        ax.YLabel.String    = 'LFP amplitude [uV]';
    else
        ax.YLim           	= [-250 400];
        ax.YAxis(2).Visible = 'off';
    end
    
    if iFreq== 4
        of                  = 15;
    else
        of                  = 25;
    end
    text(row(iFreq)+of,350, [num2str(round(freq(iFreq))) ' Hz'],'FontSize', 16, 'Color', [.5 .5 .5])
    
    axF                     = axes('Position',[row(iFreq) .08 .2 .2]);
    m                       = nanmean(fftBL);                                  	% Determine threshold
    s                       = nanstd(fftBL);
    thresh                  = m + s*2;
    avg                     = nanmean(fftSig);
    nV                      = max(max(avg), max(thresh));
    p1                      = plot(frq,avg./nV, 'LineWidth', 2, 'Color', cm(1,:), 'Parent', axF); hold on
    p2                      = plot(frq,thresh./nV, 'LineWidth', 2, 'Color', 'k', 'Parent', axF);
    
    if iFreq == 1
        axF.XColor          = [0 0 0];
        axF.YColor          = [0 0 0];
        axF.XTick           = [25 50 75 100];
        axF.XLim            = [1 110];
        axF.XLabel.String   = 'Frequency [Hz]';
        axF.YLabel.String   = 'Power [norm]';
        axF.FontSize        = 14;
        box off
        fill([freq(iFreq)-2 freq(iFreq) freq(iFreq)+2], [.5 .4 .5], 'k')
    else
        fill([freq(iFreq)-2 freq(iFreq) freq(iFreq)+2], [.3 .2 .3], 'k')
        axF.XTick           = [25 50 75 100];
        axF.XLim            = [1 110];
        axF.YAxis.Visible   = 'off';
        axF.XAxis.Visible   = 'off';
    end
    
    if iFreq == 4
        text(50,.85, 'Threshold','FontSize', 16, 'Color', [0 0 0])
        text(50,.7, 'Response','FontSize', 16, 'Color', cm(1,:))
    end
end

text(0.01,.98, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0.01,.31, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
export_fig([dest_dir 'SFIG2'], '-r400',f);

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'SFIG2'], '-dpdf', '-r400'); 
