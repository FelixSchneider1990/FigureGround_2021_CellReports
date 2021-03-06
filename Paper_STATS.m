%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATS FIGURE-GROUND EPHYS PAPER %%%%%%%%%%%%%%%
%%% RUN FIGURE SCRIPT FIRST TO ACCESS VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIGURE 1: Subject-wise test %%

% D-prime
d1      = dpE(:,2);	% M1
d2      = dpE(:,3); % M2

d1      = dpD(1:end-10,2);
d2      = dpD(1:end-10,3);

% Mean RT
d1   	= mRTE(:,2);
d2      = mRTE(:,3);

d1      = mRTD(1:end-10,2);
d2      = mRTD(1:end-10,3);

% Standard deviation RT
d1  	= CVE(:,2);
d2      = CVE(:,3);

d1      = CVD(1:end-10,2);
d2      = CVD(1:end-10,3);
        
%%% TEST %%%
smd = (mean(d2) - mean(d1)) / std([d1;d2])
[P,H,STATS] = signrank(d1,d2)      
    
%% FIGURE 1: Reaction times - Linear mixed effects model %%

% load('PATH/bRaw_M1.mat')
% bM1                 = bRaw;
% load('PATH/bRaw_M2.mat')
% bM2                 = bRaw;
% clear bRaw

sess                = [];
id                  = [];
RT                  = [];
coh                 = [];
res                 = [];
c                   = 0;

for mID = 1:2
    clear arr
    if mID == 1
        arr         = bM1;
    else
        arr         = bM2(1:end-10,:);
    end
    
    for iSess = 1:size(arr,1)   
        clear dmat
        dmat        = arr{iSess,1};
             
        % Content(curr_data):
        % column  1: Hit = 1; Error = 0;
        % column  2: PreFigEarly = 1; FigEarly = 2; Hit = 3; Late = 4;
        %            FalseAlarm = 5; CorrRejection = 6
        % column  3: RT [s]
        % column  4: Figure onset
        % column  5: Minimum bar holding time [s]
        % column  6: Maximum bar holding time [s]
        % column  7: Figure coherence
        % column  8: StimNo.
        
        c           = c+1;
        sess        = vertcat(sess, repmat(c,size(dmat,1),1));
        id          = vertcat(id, repmat(mID,size(dmat,1),1));
        RT          = vertcat(RT, dmat(:,3));
        res         = vertcat(res, dmat(:,2));
        coh         = vertcat(coh, dmat(:,7));
    end
end

HIidx   = res == 3;
tbl     = table(RT(HIidx),coh(HIidx),sess(HIidx), id(HIidx), 'VariableNames',{'RT','Coh','Sess','ID'});

% lme = fitlme(tbl, 'RT~Coh + ( Coh | Sess ) + ( Coh | ID)');               % random intercept - repeated measures
lme1    = fitlme(tbl, 'RT ~1+ ( 1 | Sess ) + ( 1 | ID)');                   % random model
lme2    = fitlme(tbl, 'RT~Coh + ( 1 | Sess ) + ( 1 | ID)');                 % Coherence as predictor
stats   = compare(lme1,lme2);
[~,~,lmeStat] = fixedEffects(lme2);

%% FIGURE 3 %%

% Perceptage responsive
size(FIG{1},1)/size(muaeE,2)
size(FIG{2},1)/size(muaeD,2)

% FIG vs CTRL
indx        = 401:600;
clear d1 d2
d1          = mean(FIG{1}(:,indx),2);
d2         	= mean(CTR{1}(:,indx),2);

d1          = mean(FIG{2}(:,indx),2);
d2          = mean(CTR{2}(:,indx),2);

%%% TEST %%%
smd         = (mean(d1) - mean(d2)) / std([d1;d2])
[P,H,STATS] = signrank(d1,d2)      
    

% Peak delay
on          = 0:50:550;
[pk, loc]   = findpeaks(mean(CTR{2}), 'MinPeakHeight',1.055);
avgPK       = mean(loc - on);

%% FIGURE 3: Difference ANT vs POS

% Test FGM magnitude between ANT vs POS
clear d1 d2

d1 = [FG8{1,1} FG12{1,1}]; % M1
d2 = [FG8{1,2} FG12{1,2}];

d1 = [FG8{2,1} FG12{2,1}]; % M2
d2 = [FG8{2,2} FG12{2,2}];

%%% TEST %%%
smd = (mean(d1) - mean(d2)) / std([d1 d2])
[P,H,STATS] = ranksum(d1,d2)  

% Coherence post-hoc test
median(COH{1,1})
median(COH{1,2})
median(COH{2,1})
median(COH{2,2})

[P(1),H,STATS] = signrank(COH{1,1}) % M1 - ANT
[P(2),H,STATS] = signrank(COH{1,2}) % M1 - POS
[P(3),H,STATS] = signrank(COH{2,1}) % M2 - ANT
[P(4),H,STATS] = signrank(COH{2,2}) % M2 - POS
fdr(P)

%% SFIGURE 3: Latency ANT vs POS

clear d1 d2
d1 = [lat8{1,1}(:,3); lat8{1,2}(:,3); lat8{2,1}(:,3); lat8{2,2}(:,3)];      % COH
d2 = [lat12{1,1}(:,3); lat12{1,2}(:,3); lat12{2,1}(:,3); lat12{2,2}(:,3)];

d1 = [lat8{1,1}(:,3); lat12{1,1}(:,3); lat8{2,1}(:,3); lat12{2,1}(:,3)];    % Field
d2 = [lat8{1,2}(:,3); lat12{1,2}(:,3); lat8{2,2}(:,3); lat12{2,2}(:,3)];

%%% TEST %%%
smd = (nanmean(d1) - nanmean(d2)) / nanstd([d1; d2])
[P,H,STATS] = ranksum(d1,d2) 

% Median latencies for coherence level
nanmedian([lat8{1,1}(:,3); lat8{1,2}(:,3); lat8{2,1}(:,3); lat8{2,2}(:,3)])
nanmedian([lat12{1,1}(:,3); lat12{1,2}(:,3); lat12{2,1}(:,3); lat12{2,2}(:,3)])

%% FIGURE 4: Onset- + Reponse-aligned against 0.5

clear d1 d2 d3 d4
inclIdx	= [IN{1,1} IN{1,2} IN{2,1} IN{2,2}];    % Only sound-responsive units with sign. FGM
d1 = allAUCFD_pop(inclIdx);
d2 = allAUCFO_pop(inclIdx);
d3 = allAUCFD_pop(~inclIdx);
d4 = allAUCFO_pop(~inclIdx);

nanmedian(d1)
nanmedian(d2)
nanmedian(d3)
nanmedian(d4)

%%% TEST %%%
[P(1),H,STATS] = signrank(d1,.5) 
[P(2),H,STATS] = signrank(d2,.5) 
[P(3),H,STATS] = signrank(d3,.5) 
[P(4),H,STATS] = signrank(d4,.5) 
fdr(P)

% Percentage of excitation
for iAn = 1:2      
        fracD(iAn) = sum(allAUCFD{iAn}(allIN{iAn}) > 0.5) / size(allAUCFD{iAn}(allIN{iAn}),2);
        fracO(iAn) = sum(allAUCFO{iAn}(allIN{iAn}) > 0.5) / size(allAUCFO{iAn}(allIN{iAn}),2);
end

%% FIGURE 4: Modulated vs unmodulated units

smd = (nanmean(d1) - nanmean(d3)) / nanstd([d1 d3])
smd = (nanmean(d2) - nanmean(d4)) / nanstd([d2 d4])

[P,H,STATS] = ranksum(d1,d3) 
[P,H,STATS] = ranksum(d2,d4) 

%% FIGURE 4: Repsonse vs onset aligned

smd = (nanmean(d1) - nanmean(d2)) / nanstd([d1 d2])

[P,H,STATS] = ranksum(d1,d2) 

nanmean(d1 - d2)


%% FIGURE 4: Coherence

clear d1 d2 d3 d4
inclIdx	= [IN{1,1} IN{1,2} IN{2,1} IN{2,2}];    % Only sound-responsive units with sign. FGM
d1 = allAUCC_pop(inclIdx);
d2 = allAUCCD_pop(inclIdx);
d3 = allAUCC_pop(~inclIdx);
d4 = allAUCCD_pop(~inclIdx);

nanmedian(d1)
nanmedian(d2)
nanmedian(d3)
nanmedian(d4)

%%% TEST %%%
[P(1),H,STATS] = signrank(d1,.5) 
[P(2),H,STATS] = signrank(d2,.5) 
[P(3),H,STATS] = signrank(d3,.5) 
[P(4),H,STATS] = signrank(d4,.5) 
fdr(P)

%% FIGURE 4: CR vs MI

clear d1 d2 d3 d4
inclIdx	= [IN{1,1} IN{1,2} IN{2,1} IN{2,2}];    % Only sound-responsive units with sign. FGM
d1 = allAUCMCD_pop(inclIdx);
d2 = allAUCMC_pop(inclIdx);
d3 = allAUCMCD_pop(~inclIdx);
d4 = allAUCMC_pop(~inclIdx);

nanmedian(d1)
nanmedian(d2)
nanmedian(d3)
nanmedian(d4)

%%% TEST %%%
[P(1),H,STATS] = signrank(d1,.5) 
[P(2),H,STATS] = signrank(d2,.5) 
[P(3),H,STATS] = signrank(d3,.5) 
[P(4),H,STATS] = signrank(d4,.5) 
fdr(P)

% Modulated vs unmodulated
smd = (nanmean(d1) - nanmean(d3)) / nanstd([d1 d3])
smd = (nanmean(d2) - nanmean(d4)) / nanstd([d2 d4])

[P,H,STATS] = ranksum(d1,d3) 
[P,H,STATS] = ranksum(d2,d4) 

%% FIGURE 4: HI vs MI

clear d1 d2 d3 d4
inclIdx	= [IN{1,1} IN{1,2} IN{2,1} IN{2,2}];    % Only sound-responsive units with sign. FGM
d1 = allAUCHMD_pop(inclIdx);
d2 = allAUCHM_pop(inclIdx);
d3 = allAUCHMD_pop(~inclIdx);
d4 = allAUCHM_pop(~inclIdx);

nanmedian(d1)
nanmedian(d2)
nanmedian(d3)
nanmedian(d4)

%%% TEST %%%
[P(1),H,STATS] = signrank(d1,.5) 
[P(2),H,STATS] = signrank(d2,.5) 
[P(3),H,STATS] = signrank(d3,.5) 
[P(4),H,STATS] = signrank(d4,.5) 
fdr(P)

% Modulated vs unmodulated
smd = (nanmean(d1) - nanmean(d3)) / nanstd([d1 d3])
smd = (nanmean(d2) - nanmean(d4)) / nanstd([d2 d4])

[P,H,STATS] = ranksum(d1,d3) 
[P,H,STATS] = ranksum(d2,d4) 

% Subject comparison
clear d1 d2 d3 d4
d1 = [AUCHM{1,1}(IN{1,1}) AUCHM{1,2}(IN{1,2})];
d2 = [AUCHM{2,1}(IN{2,1}) AUCHM{2,2}(IN{2,2})];
d3 = [AUCHMD{1,1}(IN{1,1}) AUCHMD{1,2}(IN{1,2})];
d4 = [AUCHMD{2,1}(IN{2,1}) AUCHMD{2,2}(IN{2,2})];

nanmedian(d1)
nanmedian(d2)
nanmedian(d3)
nanmedian(d4)

[P,H,STATS] = signrank(d1,.5) 
[P,H,STATS] = signrank(d2,.5)
[P,H,STATS] = signrank(d3,.5) 
[P,H,STATS] = signrank(d4,.5)

%% Bar release

idxBL       = 501:700;
idxAL       = 701:900;

clear d1 d2
d1 = mean(avgMUA(:,idxBL),2);
d2 = mean(avgMUA(:,idxAL),2);

smd = (nanmean(d2) - nanmean(d1)) / nanstd([d1; d2])
[P,H,STATS] = signrank(d2,d1)

Diff = mean(avgMUA(:,idxAL),2) - mean(avgMUA(:,idxBL),2);
mDiff = mean(Diff);
sd      = nanstd(Diff);
n       = size(Diff,1);
sem     = sd/(sqrt(n));