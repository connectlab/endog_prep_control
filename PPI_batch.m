%%%%%%%%%%%%%%% Batch PPI %%%%%%%%%%%%%%%
clear all

%%%%% Step 1: Outline base directories to do

basedir = 'Y:\mke3\analysis_fMRI_phasic_tonic';
outdir = 'Y:\mke3\analysis_fMRI_phasic_tonic\PPI';
runType = {'mono' 'multi'};
analysisType = {'cueITI' 'cueTrial'};
scanTR = 2;
m = 1; % For Matlab Batch Line
spm('defaults', 'FMRI');
spm_jobman('initcfg');
errorVal=[];
networks = {'anterior_Salience', 'dorsal_DMN', 'LECN', 'RECN', 'Visuospatial'};

tic
% Loop for analysisType
for thisAnalysis = 1:2
    
    % Loop for all subjects
    for s = 2%:25
        
        % Set current subject variable
        if s < 10
            thissubj = ['subject0' num2str(s)];
        else
            thissubj = ['subject' num2str(s)];
        end
        
        % Loop across the networks
        for n = 1:length(networks)
            
            % Selecet network dir and get ROI folder list
            netdir = fullfile('C:\Users\mkegan\Documents\MATLAB\Functional_ROIs', networks{n});
            cd(netdir)
            theseROIs=dir('*');
            theseROIs=theseROIs([theseROIs(:).isdir]);
            theseROIs = theseROIs(~ismember({theseROIs(:).name},{'.','..'}));
            
            % Loop for mono/multi
            for type = 1:2
                
                % Loop across all ROIs in that network
                for roi = 1:length(theseROIs)
                    
                    thisROI = fullfile(netdir, theseROIs(roi).name, [num2str(roi) '.nii']);
                                     
                    try
                        
                        FIND_PPI_singlesub(basedir, outdir, scanTR, thissubj, 
                        
                    catch
                        errorVal = [errorVal; s n roi type];
                    end
                end
            end
        end
    end
end

toc