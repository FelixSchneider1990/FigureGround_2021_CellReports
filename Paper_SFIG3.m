%%% FIGURE_GROUND EPHYS PAPER %%%
%%% FELIX SCHNEIDER, 02/2020 %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SFIG 3 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/fschneider/Documents/GitHub/FigureGround_Ephys_Analysis
clearvars -except muaeE muaeD

% Load MUA data
% load('PATH/muae_M1.mat')
% muaeE = muae;
% load('PATH/muae_M2.mat')
% muaeD = muae;
% clear muae

%% Frequency elements in receptive field

F8                  = []; 
F12                 = []; 
lat8                = []; 
lat12               = [];
dest_dir            = 'PATH/';
typ                 = 'muae';
indxR               = 201:400; % corresponds to -300:-100ms to decision
indx                = 401:600;
alph                = .01;

for iAn = 1:2
    
    ant             = []; 
    pos             = []; 
    muae            = [];
    
    if iAn == 1
        animalID    = 'M1';
    else
        animalID    = 'M2';
    end
    
    load([dest_dir 'maps/tMap_' animalID '_' typ  '.mat']);
    AP              = find(logical(sum(~isnan(mfr_mat),2)));
    ML              = find(logical(sum(~isnan(mfr_mat))));
    [x,y]           = coreBoundary(mfr_mat,AP,ML,false,animalID);
    
    if strcmp(animalID, 'M1')
        muae        = muaeE;
    elseif strcmp(animalID, 'M2')
        muae        = muaeD;
    end
    
    [ant, pos]      = sortUnits(x,y,muae);
    
    for iFi = 1:2
        
        if iFi == 1
            d       = ant;
            str 	= 'ant';
        else
            d       = pos;
            str     = 'pos';
        end
        
        c = 0; cc = 0; diff = []; wd = []; FIG = []; CTR = []; f8 = []; f12 = []; latency12 = []; latency8 = [];
        
        for ii = 1:size(d,2)
            incl    = check20Hz(d{ii});
            if length(d{ii}.BLslope) < 200 || incl == 0
                continue
            end
            
            datestr	= str2num([d{ii}.id(1:4) d{ii}.id(6:7) d{ii}.id(9:10)]);
            test    = mean([d{ii}.res.HI12(:,indxR); d{ii}.res.HI8(:,indxR)],2);
            ctrl    = mean(d{ii}.res.CR(:,indxR),2);
            pp      = anova1([test; ctrl],[zeros(size(test,1),1);ones(size(ctrl,1),1)], 'off');
            
            if pp < alph && sum(d{ii}.nTr>=10) == length(d{ii}.nTr) && datestr - 20190806 <= 0
                
                mBL             = nanmean(nanmean(d{ii}.on.fullAvg(:,101:500),2));      	% Average BL response
                c               = c+1;                                                     	% Counter
                fstim           = size(d{ii}.figfreqs,3);                                 	% Ratio of Fig to Gnd stim
                
                latency8(c,:)   = d{ii}.sigT8;
                latency12(c,:)  = d{ii}.sigT12;
                mFIG(c,:)   	= nanmean([d{ii}.res.HI8; d{ii}.res.HI12] ./mBL);           % Figure stimuli
                mCTR(c,:)    	= nanmean(d{ii}.res.CR ./mBL);                              % Control stimuli
                diff(c)         = nanmean(mFIG(c,indxR) - nanmean(mCTR(c,indxR)));
                wd(c)           = d{ii}.Width;
                
                if ~isnan(sum(sum(d{ii}.RF)))
                    cc          = cc+1;
                    no          = [];
                    no          = getFreqRF(d{ii}.RF, d{ii}.figfreqs);
                    
                    f8          = vertcat(f8, no(1:fstim/2,:));
                    f12       	= vertcat(f12, no(fstim/2+1:fstim,:));
                end
            end
        end
        
        lat8{iAn,iFi}   = latency8;
        lat12{iAn,iFi}  = latency12;
        F8{iAn,iFi}     = f8;
        F12{iAn,iFi}    = f12;
        Df{iAn,iFi}     = diff;
        Wd{iAn,iFi}     = wd;
        [rrr, ppp]      = corrcoef(Wd{iAn,iFi}, Df{iAn,iFi});
        rCorr(iAn,iFi)  = rrr(1,2);
        pCorr(iAn,iFi)  = ppp(1,2);
    end
end

pCorr = fdr(pCorr);

%% PLOT %%%

f       = figure('Units', 'normalized', 'Position', [0 0 .8 .8]); set(gcf,'color', [1 1 1]); axis off
col     = [0 .9 0; .9 0 0];
dim     = [.2 .2];
r       = linspace(.08, .74,4);

for iAn = 1:2
    for iFi = 1:2
        
        lf8 = [];
        for i           = 1:size(F8{iAn, iFi},1)
            lf8(i,:)    = polyfit(1:8,F8{iAn, iFi}(i,:),1);
            %             plot(0:9,polyval(lf8(i,:),0:9),'Color',[.5 .5 .5 .3], 'Marker','none')
        end
        [p8(iAn,iFi),h] = signrank(lf8(:,1));
        
        lf12 = [];
        for i = 1:size(F12{iAn, iFi},1)
            lf12(i,:)   = polyfit(1:8,F12{iAn, iFi}(i,:),1);
            %             plot(0:9, polyval(lf12(i,:),0:9),'Color',[.5 .5 .5 .3], 'Marker','none')
        end
        [p12(iAn,iFi),h] = signrank(lf12(:,1));
        
    end
end

tmp             = fdr([p8,p12]);
p8              = tmp(1:2,1:2);
p12             = tmp(1:2,3:4);

for iAn = 1:2
    for iFi = 1:2
        if iAn == 1
            c   = .07;
        else
            c   = .3;
        end
        
        if iFi == 1
            rr  = r(4);
            str = 'ANT: p = ';
        else
            rr  = r(3);
            str = 'POS: p = ';
        end
        
        ax                  = axes('Position',[c rr dim]); hold on
        ax.YLabel.String    = {'Coh8';'Elements in RF'};
        ax.FontSize         = 14;
        
        lf8 = [];
        for i = 1:size(F8{iAn, iFi},1)
            lf8(i,:)        = polyfit(1:8,F8{iAn, iFi}(i,:),1);
            p1              = plot(0:9,polyval(lf8(i,:),0:9),'Color',[.5 .5 .5 .3], 'Marker','none');
        end
        
        bx                  = boxplot(F8{iAn, iFi}, 'Colors', col(1,:));
        set(bx,'MarkerEdgeColor','k')
        set(bx, {'linew'},{2})
        ax.XAxis.Visible    = 'off';
        ax.YTick            = [0 5 10 15];
        ax.YLim             = [-1 16];
        ax.FontSize         = 14;
        box off
        
        if iFi == 1
            if iAn == 1
                ax.Title.String     = {'M1'; [str ' ' num2str(p8(iAn,iFi))]};
            else
                ax.Title.String     = {'M2'; [str ' ' num2str(p8(iAn,iFi))]};
            end
        else
            ax.Title.String         = [str ' ' num2str(p8(iAn,iFi))];
        end
        ax.Title.FontSize           = 10;
        
        if iAn == 2
            ax.YLabel.String        = '';
        end
        ax.Position                 = [c rr dim];
    end
end

for iAn = 1:2
    for iFi = 1:2
        if iAn == 1
            c   = .07;
        else
            c   = .3;
        end
        
        if iFi == 1
            rr  = r(2);
            str = 'ANT: p = ';
        else
            rr  = r(1);
            str = 'POS: p = ';
        end
        
        ax      = axes('Position',[c rr dim]); hold on
        lf12    = [];
        for i = 1:size(F12{iAn, iFi},1)
            lf12(i,:)   = polyfit(1:8,F12{iAn, iFi}(i,:),1);
            p1          = plot(0:9, polyval(lf12(i,:),0:9),'Color',[.5 .5 .5 .3], 'Marker','none');
        end
        
        bx              = boxplot(F12{iAn, iFi}, 'Colors', col(2,:));
        set(bx,'MarkerEdgeColor','k')
        set(bx, {'linew'},{2})
        box off
        
        ax.Title.String     = [str ' ' num2str(p12(iAn,iFi))];
        ax.YTick            = [0 5 10 15];
        ax.YLim             = [-1 16];
        ax.YLabel.String    = {'Coh12';'Elements in RF'};
        ax.FontSize         = 14;
        ax.Title.FontSize   = 10;

        if iAn == 2
            ax.YLabel.String    = '';
        end
        
        if iFi == 2
            ax.XLabel.String    = 'Chord after figure onset';
        else
            ax.XAxis.Visible    = 'off';
        end
        ax.Position             = [c rr dim];

    end
end


%% Latency

antLat          = vertcat(lat8{1,1},lat12{1,1},lat8{2,1},lat12{2,1});
antLat          = antLat(:,3);
antgroup    	= vertcat(ones(size(lat8{1,1},1),1),...
                ones(size(lat12{1,1},1),1)+1,...
                ones(size(lat8{2,1},1),1)+2,...
                ones(size(lat12{2,1},1),1)+3);
posLat          = vertcat(lat8{1,2},lat12{1,2}, lat8{2,2},lat12{2,2});
posLat          = posLat(:,3);
posgroup       	= vertcat(ones(size(lat8{1,2},1),1),ones(size(lat12{1,2},1),1)+1,ones(size(lat8{2,2},1),1)+2,ones(size(lat12{2,2},1),1)+3);
col             = [0 .9 0; .9 0 0];

for iFi = 1:2
    if iFi == 1
        ax      = axes('Position',[.58 .64 .2 .3]); hold on
        bx      = boxplot(antLat, antgroup, 'Colors', 'k');
        set(bx, {'linew'},{1.5})
        set(bx,'MarkerEdgeColor',[1 1 1])
        
        x       = [.8 1.8 2.8 3.8];
        for i = 1:4
            xx  = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(antgroup == i));
            sc  = scatter(xx, antLat(antgroup == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
            if i == 1 || i == 3
                sc.MarkerFaceColor  = col(1,:);
                sc.MarkerEdgeColor  = col(1,:);
            else
                sc.MarkerFaceColor  = col(2,:);
                sc.MarkerEdgeColor  = col(2,:);
            end
            sc.SizeData             = 20;
        end
    else
        ax      = axes('Position',[.79 .64 .2 .3]); hold on
        bx      = boxplot(posLat, posgroup, 'Colors', 'k');
        set(bx, {'linew'},{1.5})
        set(bx,'MarkerEdgeColor',[1 1 1])
        
        x       = [.8 1.8 2.8 3.8];
        for i = 1:4
            xx  = x(i) + ((x(i)+.4)-x(i)).*rand(1,sum(posgroup == i));
            sc  = scatter(xx, posLat(posgroup == i)','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0);
            if i == 1 || i == 3
                sc.MarkerFaceColor  = col(1,:);
                sc.MarkerEdgeColor  = col(1,:);
            else
                sc.MarkerFaceColor  = col(2,:);
                sc.MarkerEdgeColor  = col(2,:);
            end
            sc.SizeData             = 20;
        end
    end
    
    ax.YTick                = [100 200 300 400 500];
    ax.YLim                 = [30 500];
    ax.FontSize             = 14;
    ax.XTickLabel           = {'Coh8', 'Coh12', 'Coh8', 'Coh12'};
    ax.XTickLabelRotation   = 10;
  	box off   
    text(1.35, 475, 'M1', 'FontSize', 14)
    text(3.35, 475, 'M2', 'FontSize', 14)

    if iFi == 2
        ax.Title.String     = 'POS';
        ax.YAxis.Visible    = 'off';
    else
        ax.YLabel.String    = 'Effect Latency [ms]';
        ax.Title.String     = 'ANT';
    end
end

%% Tuning width

for iAn = 1:2
    for iFi = 1:2
        if iAn == 1
            if iFi == 1
                ax = axes('Position',[.58 .35 dim]);
            else
                ax = axes('Position',[.79 .35 dim]);
            end
        else
            if iFi == 1
                ax = axes('Position',[.58 .1 dim]);
            else
                ax = axes('Position',[.79 .1 dim]);
            end
        end
        
        sc                  = scatter(Wd{iAn,iFi}, Df{iAn,iFi});
        sc.MarkerFaceColor  = [0 0 0];
        sc.MarkerEdgeColor  = [0 0 0];
        ax.FontSize         = 14;
        ax.YLim             = [-.4 .5];
        ax.YTick            = [-.4 0 .5];
        ax.XLim             = [0 1];
        ax.XLabel.String    = 'Tuning width [%]';
        ax.XAxis.Color      = [0 0 0];
        ax.YAxis.Color      = [0 0 0];
        l                   = lsline;
        l.LineWidth         = 2;
        l.Color             = [1 0 0];
        
        text(.8,.49, ['n = ' num2str(length(Df{iAn,iFi}))], 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold')
        text(.8,.39, ['r = ' num2str(round(rCorr(iAn,iFi),4))], 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold')
        text(.8,.29, ['p = ' num2str(round(pCorr(iAn,iFi),4))], 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold')
        
        if iFi == 1
            if iAn == 1
                ax.YLabel.String = {'M1'; 'Fig - Ctr'};
                ax.XAxis.Visible = 'off';
            else
                ax.YLabel.String = {'M2'; 'Fig - Ctr'};
            end
        else
            if iAn == 1
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
            else
                ax.YAxis.Visible = 'off';
            end
        end
    end
end

ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(0.01,.98, 'a', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.52,.98, 'b', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.52,.6, 'c', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
dest_dir = '/Users/fschneider/ownCloud/NCL_revision/Figures/';
export_fig([dest_dir 'SFIG3'], '-r400',f);

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, [dest_dir 'SFIG3'], '-dpdf', '-r400'); 