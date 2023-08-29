%%%%%%%%%%%%%%% PPI Analysis %%%%%%%%%%%%%%%
clear all

%%%%% Step 1: Outline base directories 

basedir = 'Y:\mke3\analysis_fMRI_phasic_tonic\PPI';
atlasdir = 'C:\Users\mkegan\Documents\MATLAB\Functional_ROIs';
need_imcalc = 1;

trialtype = {'cueITI' 'cueTrial'};
runtype = {'mono' 'multi'};

networks = {'anterior_Salience', 'dorsal_DMN', 'LECN', 'RECN', 'Visuospatial'};
netlabel = {'SAL', 'DMN', 'LECN', 'RECN', 'DAN'};
roicount = [7 9 6 6 11];
netlines = cumsum(roicount);

analysisTypes = {'cueITI_only', 'noCueITI_only', 'cueTrialOnly', 'noCueTrialOnly'};

%%%%% Step 2: Create full list of all necessary ROIs
% This generates a matrix with 4 columns which are as follows:
% conFile/maskFile/Network label/ROI Seed Number
roilist = {};
m = 1;
for n = 1:length(networks)
    for roi = 1:roicount(n)
        % Pad if needed
        if roi < 10
            roiSeed = ['0' num2str(roi)];
        else
            roiSeed = num2str(roi);
        end
        
        roilist{m,1}=fullfile(networks{n}, roiSeed, 'con_0001.nii');
        roilist{m,2}=fullfile(networks{n}, roiSeed, [num2str(roi) '.nii']);
        roilist{m,3}=netlabel{n};
        roilist{m,4}=roiSeed;
        
        m = m + 1;
    end
end

%%%%% Step 3: Loop for ImCalc data
for t = 2%:length(trialtype)
    for r = 2%:length(runtype)
        for aT = 3
            
            if need_imcalc
                
                outdir = fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{r}, 'roi_masks');
                
                if exist(outdir,'dir')
                    rmdir(outdir, 's')
                    mkdir(outdir)
                else
                    mkdir(outdir)
                end
                
                % Subject loop
                %errorVal={};
                for s=1:25
                    % Set current subject variable
                    if s < 10
                        thissubj = ['subject0' num2str(s)];
                    else
                        thissubj = ['subject' num2str(s)];
                    end
                    
                    subjdir = fullfile(basedir, trialtype{t}, thissubj, 'FIND_PPI', 'cueTrial', runtype{r});
                    
                    % Loop across all ROIs (con_file)
                    for i = 1:length(roilist)
                        
                        % Loop across all ROIs (mask_file)
                        parfor j = 1:length(roilist)
                            
                            try
                                inputfiles = { % Put into ImCalc format
                                    fullfile(subjdir, [roilist{i,1} ',1']) % Con File
                                    fullfile(atlasdir, [roilist{j,2} ,',1']) % Mask
                                    };
                                
                                outpath = fullfile(outdir, [thissubj '_' roilist{i,3} '_' roilist{i,4} '_' roilist{j,3} '_' roilist{j,4} '_iso.nii']);
                                
                                spm_imcalc(inputfiles,outpath,'i1.*(i2>1)',{0,1,1,4,''});
                                
                            catch
                                [s i j]
                            end
                        end
                    end
                end
            end
        end
    end
end

tic
%%%%% Step 4: Calculating matrix values
for t = 2%1:2
    for r = 2%1:2
        for aT = 3%1:2
            
            resultsdir = fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{r});
            outdir = fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{r}, 'roi_masks');
            cd(outdir)
            
            resultsMat = cell(39,39);
            nameMat = cell(39,39);
            
            for q=1:length(roilist)
                parfor p=1:length(roilist)
                    %nameMat{q,p} = [roilist{q,3} '_' roilist{q,4} '_' roilist{p,3} '_' roilist{p,4}];
                    pairname = [roilist{q,3} '_' roilist{q,4} '_' roilist{p,3} '_' roilist{p,4}];
                    filelist = dir(['*' pairname '*.nii']);
                    tmplist = [];
                    for t=1:length(filelist)
                        tmplist = [tmplist mean2(spm_read_vols(spm_vol(fullfile(outdir, filelist(t).name))))];
                    end
                    resultsMat{q,p} = tmplist;
                end
            end
            
            save(fullfile(resultsdir, 'resultsMatFull'), 'resultsMat');
            toc
        end
    end
end

%%
%%%%% Step 5: Perform the t-tests - single type
t = 2;
r = 2;
aT = 4;
analysisLabels = {'cueITIonly', 'noCueITIonly', 'cueTrialOnly', 'noCueTrialOnly'};
load(fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{r}, 'resultsMatFull.mat'))

rawsigMat = zeros(39,39);
pvalMat = zeros(39,39);
tscoreMat = zeros(39,39);
for j=1:length(resultsMat)
    for k=1:length(resultsMat)
        try
            [h p ci stats] = ttest(resultsMat{j,k});
            if j==k
                rawsigMat(j,k) = 0;
                pvalMat(j,k) = -.1;
                tscoreMat(j,k) = 255;
            else
                rawsigMat(j,k) = h;
                pvalMat(j,k) = p;
                tscoreMat(j,k) = stats.tstat;
                
            end
        catch
            %disp(j)
            %disp(k)
        end
    end
end

netlines = cumsum(roicount);

%%% Raw Scores
figure('Units','normalized'); imagesc(0,0, rawsigMat); 
colormap([ .95 .95 .95; 1 0 0; 0 0 0]); caxis([0 2]);
title(['t-test significant or not (1 or 0) - '  analysisLabels{aT} ' ' runtype{r}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end


%%% P-values
for ii = 1:length(pvalMat)
    for jj = 1:length(pvalMat)
        if pvalMat(ii,jj)==0
            pvalMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, pvalMat); 
%colormap([0 0 0; flip(hot(20)); .9 .9 .9]); caxis([0 .1]);
colormap([hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
title(['p-value of t-test - ' analysisLabels{aT} ' ' runtype{r}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%% t-scores
for ii = 1:length(tscoreMat)
    for jj = 1:length(tscoreMat)
        if tscoreMat(ii,jj)==0
            tscoreMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, tscoreMat); 
colormap([.95 .95 .95; jet(20); 0 0 0]); caxis([-5 5]);
colorbar('westoutside')
title(['t-value of t-test - ' analysisLabels{aT} ' ' runtype{r}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%
%%%%% Step 5: Perform the t-tests - paired t-test
t = 2;
aT = 4;
analysisLabels = {'cueITIonly', 'noCueITIonly', 'cueTrialOnly', 'noCueTrialOnly'};

monoMat = load(fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{1}, 'resultsMatFull.mat'));
multiMat = load(fullfile(basedir, trialtype{t}, 'FIND_PPI', analysisTypes{aT}, runtype{2}, 'resultsMatFull.mat'));

rawsigMat = zeros(39,39);
pvalMat = zeros(39,39);
tscoreMat = zeros(39,39);
for j=1:length(resultsMat)
    for k=1:length(resultsMat)
        try
            [h p ci stats] = ttest(monoMat.resultsMat{j,k},multiMat.resultsMat{j,k});
            if j==k
                rawsigMat(j,k) = 255;
                pvalMat(j,k) = -.1;
                tscoreMat(j,k) = 255;
            else
                rawsigMat(j,k) = h;
                pvalMat(j,k) = p;
                tscoreMat(j,k) = stats.tstat;
                
            end
        catch
            disp(j)
            disp(k)
        end
    end
end

%%% Raw Scores
figure('Units','normalized'); imagesc(0,0, rawsigMat); 
colormap([ .95 .95 .95; 1 0 0; 0 0 0]); caxis([0 2]);
%title(['Mono v Multi paired t-test significant or not (1 or 0) - ' trialtype{t}])
title(['Mono v Multi paired t-test significant or not (1 or 0) - ' analysisLabels{aT}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end


%%% P-values
for ii = 1:length(pvalMat)
    for jj = 1:length(pvalMat)
        if pvalMat(ii,jj)==0
            pvalMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, pvalMat); 
colormap([hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
%title(['Mono v Multi paired p-value of t-test - ' trialtype{t}])
title(['Mono v Multi paired p-value of t-test - ' analysisLabels{aT}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%% t-scores
for ii = 1:length(tscoreMat)
    for jj = 1:length(tscoreMat)
        if tscoreMat(ii,jj)==0
            tscoreMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, tscoreMat); 
colormap([.95 .95 .95; jet(20); 0 0 0]); caxis([-5 5]);
colorbar('westoutside')
%title(['Mono v Multi paired t-value of t-test - ' trialtype{t}])
title(['Mono v Multi paired t-value of t-test - ' analysisLabels{aT}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%
%%%%% Step 6: Perform 2 way ANOVAs

%%%%%
%%%%%%% Intertrial %%%%%%%
%%%%%

monoMatCueITI = load(fullfile(basedir, trialtype{1}, 'FIND_PPI', analysisTypes{1}, runtype{1}, 'resultsMatFull.mat'));
multiMatCueITI = load(fullfile(basedir, trialtype{1}, 'FIND_PPI', analysisTypes{1}, runtype{2}, 'resultsMatFull.mat'));
monoMatNoCueITI = load(fullfile(basedir, trialtype{1}, 'FIND_PPI', analysisTypes{2}, runtype{1}, 'resultsMatFull.mat'));
multiMatNoCueITI = load(fullfile(basedir, trialtype{1}, 'FIND_PPI', analysisTypes{2}, runtype{2}, 'resultsMatFull.mat'));

structMat = {struct2cell(monoMatCueITI) struct2cell(monoMatNoCueITI)
            struct2cell(multiMatCueITI) struct2cell(multiMatNoCueITI)};
anovaMat(1,1) = structMat{1,1}; anovaMat(1,2) = structMat{1,2}; anovaMat(2,1) = structMat{2,1}; anovaMat(2,2) = structMat{2,2};

rawsigMat = zeros(39,39);
pvalMat = zeros(39,39);
fscoreMatCue = zeros(39,39);
fscoreMatContent = zeros(39,39);
fscoreMatInteraction = zeros(39,39);
pMain1Mat = zeros(39,39);
pMain2Mat = zeros(39,39);
pIntMat = zeros(39,39);
for j = 1:39
    for k = 1:39
        try
            clear p tbl stats
            tmpMat = [cell2mat(anovaMat{1}(j,k))' cell2mat(anovaMat{2}(j,k))'
                cell2mat(anovaMat{3}(j,k))' cell2mat(anovaMat{4}(j,k))'];
            
            numReps = length(tmpMat)/2;
            
            [p, tbl, stats] = anova2(tmpMat,numReps,'off');
            
            if j==k
                rawsigMat(j,k) = 255;
                fscoreMatCue(j,k) = 0;
                fscoreMatContent(j,k) = 0;
                fscoreMatInteraction(j,k) = 0;
            else
                %rawsigMat(j,k) = sum(p<0.05);
                if sum(p<0.05)>1
                    %                     if p(1)<0.05 && p(2)<0.05
                    rawsigMat(j,k) = 4;
                    %                     elseif p(1)<0.05 && p(3)<0.05
                    %                         rawsigMat(j,k) = 1;
                    %                     elseif p(2)<0.05 && p(3)<0.05
                    %                         rawsigMat(j,k) = 2;
                    %                    rawsigMat(j,k) = 4;
                    %                elseif p(1)<0.05 % columns, cueing
                    %                    rawsigMat(j,k) = 1;
                elseif p(2)<0.05 % rows, content
                    rawsigMat(j,k) = 2;
                    %                                elseif p(3)<0.05 % interaction
                    %                                    rawsigMat(j,k) = 3;
                end
                
                fscoreMatCue(j,k) = tbl{2,5};
                fscoreMatContent(j,k) = tbl{3,5};
                fscoreMatInteraction(j,k) = tbl{4,5};
            end
        catch
        end
        
    end
end

%%% Raw Scores
figure('Units','normalized'); imagesc(0,0, rawsigMat); 
colormap([ .80 .80 .80; [0 .6 1]; [1 0 0]; [0 .9 .2]; [1 1 0]; 0 0 0]); caxis([0 6]);
title('ANOVA - Intertrial - Raw Significance')
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%% P-values
for ii = 1:length(pvalMat)
    for jj = 1:length(pvalMat)
        if pvalMat(ii,jj)==0
            pvalMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, pvalMat); 
colormap([hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
%title(['Mono v Multi paired p-value of t-test - ' trialtype{t}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%%%%
%%%%%%% Trial %%%%%%%
%%%%%

monoMatCueTrial = load(fullfile(basedir, trialtype{2}, 'FIND_PPI', analysisTypes{3}, runtype{1}, 'resultsMatFull.mat'));
multiMatCueTrial = load(fullfile(basedir, trialtype{2}, 'FIND_PPI', analysisTypes{3}, runtype{2}, 'resultsMatFull.mat'));
monoMatNoCueTrial = load(fullfile(basedir, trialtype{2}, 'FIND_PPI', analysisTypes{4}, runtype{1}, 'resultsMatFull.mat'));
multiMatNoCueTrial = load(fullfile(basedir, trialtype{2}, 'FIND_PPI', analysisTypes{4}, runtype{2}, 'resultsMatFull.mat'));

structMat = {struct2cell(monoMatCueTrial) struct2cell(monoMatNoCueTrial)
            struct2cell(multiMatCueTrial) struct2cell(multiMatNoCueTrial)};
anovaMat(1,1) = structMat{1,1}; anovaMat(1,2) = structMat{1,2}; anovaMat(2,1) = structMat{2,1}; anovaMat(2,2) = structMat{2,2};

rawsigMat = zeros(39,39);
pvalMat = zeros(39,39);
tscoreMat = zeros(39,39);
pMain1Mat = zeros(39,39);
pMain2Mat = zeros(39,39);
pIntMat = zeros(39,39);
for j = 1:39
    for k = 1:39
        try
            tmpMat = [cell2mat(anovaMat{1}(j,k))' cell2mat(anovaMat{2}(j,k))'
                cell2mat(anovaMat{3}(j,k))' cell2mat(anovaMat{4}(j,k))'];
            
            numReps = length(tmpMat)/2;
            
            [p, tbl, stats] = anova2(tmpMat,numReps,'off');
            
            if j==k
                rawsigMat(j,k) = 255;
                pMain1Mat(j,k)=0;
                pMain2Mat(j,k)=0;
                pIntMat(j,k)=0;
            else
                %rawsigMat(j,k) = sum(p<0.05);
                if sum(p<0.05)>1
                    if p(1)<0.05 && p(2)<0.05
                        rawsigMat(j,k) = 4;
                    elseif p(1)<0.05 && p(3)<0.05
                        rawsigMat(j,k) = 1;
                    elseif p(2)<0.05 && p(3)<0.05
                        rawsigMat(j,k) = 2;
                    end
                    %                   rawsigMat(j,k) = 4;
                    %                elseif p(1)<0.05 % columns, cueing
                    %                    rawsigMat(j,k) = 1;
                    %                elseif p(2)<0.05 % rows, content
                    %                    rawsigMat(j,k) = 2;
                    %                elseif p(3)<0.05 % interaction
                    %                    rawsigMat(j,k) = 3;
                    
                end
                pMain1Mat(j,k)=p(1);
                pMain2Mat(j,k)=p(2);
                pIntMat(j,k)=p(3);
            end
            
        catch
        end
        
    end
end

%%% Raw Scores
figure('Units','normalized'); imagesc(0,0, rawsigMat); 
colormap([ .80 .80 .80; [0 .6 1]; [1 0 0]; [0 .9 .2]; [1 1 0]; 0 0 0]); caxis([0 6]);
title('ANOVA - Trial - Raw Significance')
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

%%% P-values
for ii = 1:length(pvalMat)
    for jj = 1:length(pvalMat)
        if pvalMat(ii,jj)==0
            pvalMat(ii,jj) = -5;
        end
    end
end

figure('Units','normalized'); imagesc(0,0, pvalMat); 
colormap([hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
%title(['Mono v Multi paired p-value of t-test - ' trialtype{t}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
%imagesc(0,0, diagmat, 'AlphaData', 0.5);
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

% Post-hoc
thisTrial = {'Intertrial' 'Trial'};
t=1;

figure('Units','normalized'); imagesc(0,0, pMain1Mat); 
colormap([0 0 0; hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
title(['Full Network ANOVA p-values - Factor "Cueing" - ' thisTrial{t}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

figure('Units','normalized'); imagesc(0,0, pMain2Mat); 
colormap([.8 .8 .8; hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
title(['Full Network ANOVA p-values - Factor "Content" - ' thisTrial{t}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end

figure('Units','normalized'); imagesc(0,0, pIntMat); 
colormap([.8 .8 .8; hot(20); .9 .9 .9]); caxis([0 .1]);
colorbar('westoutside')
title(['Full Network ANOVA p-values - Factor "Interaction" - ' thisTrial{t}])
xlim([-.5 38.5]); ylim([-.5 38.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end
%%
%%%%% Step 7: Perform NBS using ANOVAS

%%%%%
%%%%%%% Intertrial %%%%%%%
%%%%%

types = {'monoCue' 'monoNoCue' 'multiCue' 'multiNoCue'};
for s = 1:25
    for t = 1:4
        tmpMat = zeros(39,39);
        try
            for j = 1:39
                try
                    for k = 1:39
                        try
                            if j == k
                                tmpMat(j,k) = 1;
                            else
                                switch t
                                    case 1
                                        tmpMat(j,k) = anovaMat{1,1}{j,k}(s);
                                    case 2
                                        tmpMat(j,k) = anovaMat{1,2}{j,k}(s);
                                    case 3
                                        tmpMat(j,k) = anovaMat{2,1}{j,k}(s);
                                    case 4
                                        tmpMat(j,k) = anovaMat{2,2}{j,k}(s);
                                end
                                
                                if abs(tmpMat(j,k)) > 100
                                    tmpMat(j,k) = 0;
                                end
                            end
                        catch
                            tmpMat(j,k) = 0;
                        end
                    end
                catch
                    tmpMat(j,k) = 0;
                end
                
                
            end
            
          tmpMat(:,37)=[];
          tmpMat(:,37)=[];
          tmpMat(37,:)=[];
          tmpMat(37,:)=[];
            
            if s<10
                dlmwrite(['C:\Users\mkegan\Documents\MATLAB\NBS_test\NoCereb\Trial\s0' num2str(s) '_' types{t} '.txt'],tmpMat,'delimiter',' ')
            else
                dlmwrite(['C:\Users\mkegan\Documents\MATLAB\NBS_test\NoCereb\Trial\s' num2str(s) '_' types{t} '.txt'],tmpMat,'delimiter',' ')
            end
        catch
        end
    end
end

thisDNBS = dlmread('C:\Users\mkegan\Documents\MATLAB\NBSDirected1.0\Result\13-Jul-2021_11_57_45\Intertrial\11.00-Treshold-ftest-Test-5000-Perms-0.026-pValue\adjacency_matrix_network_weighted_1.txt');
figure('Units','normalized'); imagesc(0,0, thisDNBS); 
colormap([.8 .8 .8; hot(20); .9 .9 .9]); %caxis([0 .1]);
colorbar('westoutside')
title('NBS Weights - Factor "Intertrial"')
xlim([-.5 36.5]); ylim([-.5 36.5]); set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca, 'XTick', [4 11 18 24 32], 'xticklabel', netlabel, 'xaxisLocation', 'top')
set(gca, 'YTick', [4 11 18 24 32], 'yticklabel', netlabel, 'yaxisLocation', 'left')
hold on
for i=1:length(netlines)
    xline(netlines(i)-.5,'LineWidth',1);
    yline(netlines(i)-.5,'LineWidth',1);
end