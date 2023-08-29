
%%%%%%%%%%%%%%%  Individual Subject PPI %%%%%%%%%%%%%%
function PPI_singlesubj(basedir, outdir, scanTR, thissubj, type, runtype, thisAnalysis, analysistype, network, roi, thisROI)

foldername = {'cueTrial', 'noCueTrial'};

%%%%% Step 1: Initiate new matlabbatch start
m = 1;

%%%%% Step 2: Collect all necessary files

%analysisTypes = {'cueITI', 'cueTrial', 'cueITI_solo', 'noCueITI_solo'};
filedir = fullfile(basedir, thissubj, 'functional');
cd(filedir)
xdirs = dir([runtype '*']);
glmdir = fullfile(outdir, analysistype, thissubj, [runtype 'GLM']);

% Check to see if the PPI outdir exists, if so delete
% it and remake, otherwise just make it
if roi < 10
    ppidir = fullfile(outdir, analysistype, thissubj, 'FIND_PPI', foldername{thisAnalysis-4}, runtype, network, ['0' num2str(roi)]);
else
    ppidir = fullfile(outdir, analysistype, thissubj, 'FIND_PPI', foldername{thisAnalysis-4}, runtype, network, num2str(roi));
end
if exist(ppidir,'dir')
    rmdir(ppidir, 's')
    mkdir(ppidir)
else
    mkdir(ppidir)
end

imglist = {}; % Functional files
for i=1:length(xdirs)
    tmpdir = fullfile(filedir, xdirs(i).name);
    cd(tmpdir)
    tmpfiles = dir('swrf*.img');
    for f=1:length(tmpfiles)
        imglist{length(imglist)+1} = [tmpdir '\' tmpfiles(f).name ',1']; % Functional files
    end
end
imglist = imglist'; % Needs to be inverted

% clear old batch
clear matlabbatch

%%%%% Step 6: Setup VOI extraction
voiName = [network '_' num2str(roi)];

matlabbatch{m}.spm.util.voi.spmmat = {fullfile(glmdir,'SPM.mat')};
matlabbatch{m}.spm.util.voi.adjust = 1;
matlabbatch{m}.spm.util.voi.session = 1;
matlabbatch{m}.spm.util.voi.name = voiName;

% The anatomical ROI mask
matlabbatch{m}.spm.util.voi.roi{1}.mask.image = {thisROI};

% Subject bounding mask
matlabbatch{m}.spm.util.voi.roi{2}.mask.image = {fullfile(glmdir, 'mask.nii')};

% The logical expression combining images.
matlabbatch{m}.spm.util.voi.expression = 'i1 & i2';

% m = m + 1; % Update matlabbatch counter
spm_jobman('run', matlabbatch);
clear matlabbatch

%%%%% Step 7: PPI calculation
matlabbatch{m}.spm.stats.ppi.spmmat = {fullfile(glmdir,'SPM.mat')};
matlabbatch{m}.spm.stats.ppi.type.ppi.voi = {fullfile(glmdir, ['VOI_' voiName '_1.mat'])};
if thisAnalysis == 1 % cueITI (solo)
    ppiname = ['cueITI_' voiName];
    matlabbatch{m}.spm.stats.ppi.type.ppi.u = [1 1 1];
elseif thisAnalysis == 2 % noCueITI (solo)
    ppiname = ['nocueITI_' voiName];
    matlabbatch{m}.spm.stats.ppi.type.ppi.u = [2 1 1];
elseif thisAnalysis == 3 % cueTrial (solo)
    ppiname = ['cueTrial_' voiName];
    if type == 1 % Mono
        matlabbatch{m}.spm.stats.ppi.type.ppi.u = [3 1 1];
    else % Multi
        matlabbatch{m}.spm.stats.ppi.type.ppi.u = [3 1 1; 5 1 1; 7 1 1; 9 1 1];
    end
elseif thisAnalysis == 4 % noCueTrial (solo)
    ppiname = ['noCueTrial_' voiName];
    if type == 1 % Mono
        matlabbatch{m}.spm.stats.ppi.type.ppi.u = [4 1 1];
    else % Multi
        matlabbatch{m}.spm.stats.ppi.type.ppi.u = [4 1 1; 6 1 1; 8 1 1; 10 1 1];
    end
end
matlabbatch{m}.spm.stats.ppi.name = ppiname;
matlabbatch{m}.spm.stats.ppi.disp = 1;

% m = m + 1; % Update matlabbatch counter
spm_jobman('run', matlabbatch);
clear matlabbatch

%%%%% Steps 8-10: PPI model specification, estimation, and contrast

matlabbatch{m}.spm.stats.fmri_spec.dir = {ppidir};
matlabbatch{m}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{m}.spm.stats.fmri_spec.timing.RT = scanTR;
matlabbatch{m}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{m}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{m}.spm.stats.fmri_spec.sess.scans = imglist;
matlabbatch{m}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{m}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{m}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{m}.spm.stats.fmri_spec.sess.multi_reg = {
    fullfile(glmdir, ['PPI_' ppiname '.mat'])
    %fullfile(glmdir, 'PPI_nocueITI_over_cueITI_dDMN.mat')
    %fullfile(glmdir, 'PPI_nocueITImono.mat')
    %fullfile(glmdir, 'PPI_nocueTrial_over_cueTrial.mat')
    fullfile(glmdir, 'MotionRegressors.mat')
    };
matlabbatch{m}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{m}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{m}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{m}.spm.stats.fmri_spec.volt = 1;
matlabbatch{m}.spm.stats.fmri_spec.global = 'None';
matlabbatch{m}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{m}.spm.stats.fmri_spec.mask = {''};
matlabbatch{m}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{m+1}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{m+1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{m+1}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{m+2}.spm.stats.con.spmmat(m) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{m+2}.spm.stats.con.consess{1}.tcon.name = 'PPI';
matlabbatch{m+2}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
matlabbatch{m+2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{m+2}.spm.stats.con.delete = 0;

%%%%% Step 11: Run the Batch code

spm_jobman('run', matlabbatch);
clear matlabbatch

end