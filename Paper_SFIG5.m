clearvars -except muaeE muaeD lfpE lfpD
addpath /Users/fschneider/ownCloud/NCL_revision/BR_control
datdir = dir('/Users/fschneider/ownCloud/NCL_revision/BR_control');
mua = [];
lfp = [];

for i = 3:size(datdir,1)/2+1
    
    clear BRsig
    load([datdir(i).folder '/' datdir(i).name])
    
    c           = i-2;
    idxBL       = 501:700;
    idxAL       = 701:900;
    
    mbl         = mean(mean(BRsig.MUA(:,501:700),2));
    mat         = BRsig.MUA ./ mbl;
    mmat        = mean(mat);
    bl          = mean(mat(:,idxBL),2);
    alig        = mean(mat(:,idxAL),2);
    
    [p(c),h(c)]	= signrank(bl,alig);
    
    
    mbl_lfp   	= mean(mean(BRsig.LFP(:,501:700),2));
    mat_lfp    	= BRsig.LFP;
    mmat_lfp    = mean(mat_lfp);
    bl_lfp   	= mean(mat_lfp(:,idxBL),2);
    alig_lfp  	= mean(mat_lfp(:,idxAL),2);
    
    [pp(c),hh(c)]	= signrank(bl_lfp,alig_lfp);
    
    cIdx        = zeros(size(bl,1),1) + c;
    mua        	= [mua;[bl, alig]];
    lfp         = [lfp;[bl_lfp, alig_lfp]];
    avgMUA(c,:) = mmat;
    avgLFP(c,:) = mmat_lfp;
end


%% PLOT

idx     = 501:1050;
CI      = bootci(1000,@mean,avgMUA(:,idx));
l     	= 1:length(CI);
len   	= [l fliplr(l)];
sem_plt = [CI(1,:) fliplr(CI(2,:))];

f = figure('Units', 'normalized', 'Position', [0 0 .6 .6]); set(gcf,'color', [1 1 1]);
ax1 = subplot(6,4,[1 2 5 6]);hold on
imagesc(avgMUA(:,idx))
colormap(gray(256));
caxis(ax1,[0 3])
ax1.YLim = [1 30];
ax1.YTick = [5 10 15 20 25 30];
ax1.YLabel.String = 'No. Units';
yyaxis right
line([500 500],[-1 2], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'w')
% fill([201 201 400 400],[-1 2 2 -1],[0 0 0],'EdgeColor','none', 'FaceAlpha', .3);
fill(len,sem_plt,[1 0 0]*0.7,'EdgeColor','none', 'FaceAlpha', .3);
plot(mean(avgMUA(:,idx)), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
ax1.FontSize = 14;
ax1.XLim = [1 length(mean(avgMUA(:,idx)))];
ax1.YLim = [.9 1.8];
ax1.XTick = [1 100 200 300 400 500];
ax1.XTickLabel = {'-500','-400','-300','-200','-100', '0'};
ax1.YLabel.String = 'MUA [norm]';
ax1.XLabel.String = 'Time [ms]';
ax1.Position(1) = .08;
ax1.YAxis(2).Color = 'r';
ax1.Title.String = 'MUA';
ax1.XAxis.Color = [1 1 1];

CI      = bootci(1000,@mean,avgLFP(:,idx));
l     	= 1:length(CI);
len   	= [l fliplr(l)];
sem_plt = [CI(1,:) fliplr(CI(2,:))];

ax2 = subplot(6,4,[3 4 7 8]);
hold on
imagesc(avgLFP(:,idx))
colormap(gray(256));
caxis(ax2,[-100 400])
ax2.YLim = [1 30];
ax2.YTick = [5 10 15 20 25 30];
ax2.YLabel.String = 'No. Units';
yyaxis right
line([500 500],[-100 500], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'w')
fill(len,sem_plt,[1 0 0]*0.7,'EdgeColor','none', 'FaceAlpha', .3);
plot(mean(avgLFP(:,idx)), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
ax2.FontSize = 14;
ax2.XLim = [1 length(mean(avgLFP(:,idx)))];
ax2.YLim = [-100 500];
ax2.XTick = [1 100 200 300 400 500];
ax2.XTickLabel = {'-500','-400','-300','-200','-100', '0'};
ax2.YLabel.String = 'LFP [uV]';
ax2.XLabel.String = 'Time [ms]';
ax2.Position(1) = .565;
ax2.YAxis(2).Color = 'r';
ax2.Title.String = 'LFP';
ax2.XAxis.Color = [1 1 1];

CI      = bootci(1000,@mean,mat(:,idx));
l     	= 1:length(CI);
len   	= [l fliplr(l)];
sem_plt = [CI(1,:) fliplr(CI(2,:))];

ax3 = subplot(6,4,[9 10 13 14]);
hold on
imagesc(mat(:,idx))
colormap(gray(256));
caxis(ax3,[0 3])
ax3.YLim = [1 size(mat,1)];
% ax3.YTick = [5 10 15 20 25 30];
ax3.YLabel.String = 'No. Trials';
yyaxis right
line([500 500],[-1 2], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'w')
% fill([201 201 400 400],[-1 2 2 -1],[0 0 0],'EdgeColor','none', 'FaceAlpha', .3);
fill(len,sem_plt,[1 0 0]*0.7,'EdgeColor','none', 'FaceAlpha', .3);
plot(mean(mat(:,idx)), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
ax3.FontSize = 14;
ax3.XLim = [1 length(mean(mat(:,idx)))];
ax3.YLim = [.9 1.8];
ax3.XTick = [1 100 200 300 400 500];
ax3.XTickLabel = {'-500','-400','-300','-200','-100', '0'};
ax3.YLabel.String = 'MUA [norm]';
ax3.XLabel.String = 'Time [ms]';
ax3.Position(1) = .08;
ax3.Position(2) = .41;
ax3.YAxis(2).Color = 'r';

ax4 = subplot(6,4,[11 12 15 16]);
hold on
imagesc(mat_lfp(:,idx))
colormap(gray(256));
caxis(ax4,[-100 400])
ax4.YLim = [1 size(mat,1)];
% ax4.YTick = [5 10 15 20 25 30];
ax4.YLabel.String = 'No. Trials';
yyaxis right
line([500 500],[-100 500], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'w')
fill(len,sem_plt,[1 0 0]*0.7,'EdgeColor','none', 'FaceAlpha', .3);
plot(mean(mat_lfp(:,idx)), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
ax4.FontSize = 14;
ax4.XLim = [1 length(mean(mat_lfp(:,idx)))];
ax4.YLim = [-100 500];
ax4.XTick = [1 100 200 300 400 500];
ax4.XTickLabel = {'-500','-400','-300','-200','-100', '0'};
ax4.YLabel.String = 'LFP [uV]';
ax4.XLabel.String = 'Time [ms]';
ax4.Position(1) = .565;
ax4.Position(2) = .41;
ax4.YAxis(2).Color = 'r';

% load('/Users/fschneider/OneDrive - Newcastle University/BR_control/tfa_BR.mat')
% fr = round(linspace(100,7,100));
% ax5 = subplot(6,4,[19 20 23 24]);
% map             = 10*log10(alltfa);
% imagesc(flipud(mean(map(:,501:1000,:),3)))
% line([500 500],[-100 500], 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k')
% ax5.XTick        = [100 200 300 400 500];
% ax5.XTickLabel   = {'-400', '-300','-200','-100', '0'};
% ax5.FontSize     = 14;
% % ax5.YAxis.Visible = 'off';
% ax5.XLabel.String = 'Time [ms]';
% ax5.YLabel.String = 'Frequency [Hz]';
% ax5.Position(1) = .565;
% ax5.Position(2) = .08;
% ax5.YTick            = [1 25 50 75 100];
% ax5.YTickLabel       = [fr(1) fr(25) fr(50) fr(75) fr(100)];
% 
% r           = rectangle('Position',[201 1 200 70]);
% r.LineWidth = 2;
% r.LineStyle = '--';
% r           = rectangle('Position',[201 75 200 25]);
% r.LineWidth = 2;
% r.LineStyle = '--';

% cb              = colorbar(ax5);
% cb.Label.String = 'Power [dB]';
% cb.FontSize     = 12;
% % cb.Ticks        = [-1 -.5 0 .5 1];
% cb.Position(1)  = .35;
% cb.Position(4)  = cb.Position(4)/2;
% cb.Position(2)  = ax5.Position(2) + (ax5.Position(4)/4);
% colormap(gca,[[1 1 1];jet(256)]);
% caxis(ax5,[-1 1])
% box off 

cb2              = colorbar(ax1);
cb2.Label.String = 'MUA [norm]';
cb2.FontSize     = 12;
% cb2.Ticks        = [0 .25 .5 .75 1];
cb2.Position(1)  = .15;
cb2.Position(4)  = cb2.Position(4)/2;
cb2.Position(2)  = .15;

cb3              = colorbar(ax2);
cb3.Label.String = 'LFP [uV]';
cb3.FontSize     = 12;
% cb2.Ticks        = [0 .25 .5 .75 1];
cb3.Position(1)  = .25;
cb3.Position(4)  = cb3.Position(4)/2;
cb3.Position(2)  = .15;

ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(0.01,.97, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(0.01,.67, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
% text(.47,.36, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
export_fig([dest_dir 'SFIG_BRcontrol'], '-r400',f);
% exportgraphics(f,[dest_dir 'SFIG_BRcontrol.pdf'],'ContentType','vector','Resolution',300)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'SFIG_BRcontrol'], '-dpdf', '-r400'); 