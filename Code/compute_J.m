
%% Apply linear DCM to attention-to-visual-motion fMRI dataset

% 1- Check if attention folder exists;
% if not, instruct the user to:
% download it from http://www.fil.ion.ucl.ac.uk/spm/data/attention/
% then, unzip and place attention folder under Data folder

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(rootDir));
cd(rootDir)

folderName = 'attention';
downloadURL = 'http://www.fil.ion.ucl.ac.uk/spm/data/attention/';

if ~exist(folderName, 'dir')
    error(['Folder "%s" is missing.\nDownload attention.zip from: %s\n' ...
        'Extract and place attention folder inside Data folder'], ...
          folderName, downloadURL);
end

% 2- Load the precomputed SPM for attention dataset:
load(fullfile(rootDir,'Data\GLM\SPM.mat'),'SPM');
% To compute your own SPM.mat for the attention dataset,
% follow sections 35.3.1&2 of SPM12 Manual (p.316-318).
% If you do so, put the resulting GLM folder under Data.

addpath(genpath(rootDir))% in case the user added their own GLM folder
cd(rootDir)

% 3- Keep functional file names only (ie remove path to avoid dir error)
% will use cd in spm_mb_ui_mod (line 175) to switch to Data\attention\functional
for i = 1:size(SPM.xY.VY,1)
    str = SPM.xY.VY(i).fname;
    tmp = split(str,'\');
    str_new = tmp{end};
    SPM.xY.VY(i).fname = str_new;
end

% 4- Set relative file path for SPM.mat and contrast images
SPM.swd = fullfile(rootDir,'\Data\GLM'); 

% 5- Run linear DCM (takes about 25 min on avg PC)
MB  = spm_mb_ui_mod('specify',SPM); % uses spm_dcm_J_mod.m inside
J = MB{1}.J; % effective connectivity (expected value)
Cp = MB{1}.VAR; % posterior variances. This Cp is n*n (not n^2 * n^2 ).
% ie Cp contains the diagonal entries of a normal Cp = posterior covariance matrix)
J = J.* (abs(J)> 1e-4);
n = size(J,1); % #nodes
pC = zeros(n,n); 
pC(J~=0) = 1/512; % prior variance (from spm_dcm_J, line 118)

% 6- Save results in Data folder
cd(fullfile(rootDir,'/Data'))
save('J_1024_user.mat','J','Cp','pC','MB')
