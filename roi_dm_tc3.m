% To calculate time courses from the seed maps
% created by Anita Barber
% adapted by Jun Hua

function mn_roi = roi_dm_tc3(data_path,rundir,outputdir,roi_path,data_filenm_mask, roi_filenm_mask)

if(nargin<4),
    error('Not enough arguements');
end;

%% load seed files
n=1;
for i_rmask=1:length(roi_filenm_mask),
    files=dir(fullfile(roi_path,[roi_filenm_mask{i_rmask},'.nii']));
    for i_roi=1:length(files),
        P{n}=[fullfile(roi_path,files(i_roi).name)];
        mn_roi.name{n}=strrep(files(i_roi).name,'.nii','');
        n=n+1;
    end;
    clear files;
end;

P=strvcat(P{:});

V=spm_vol(P);
sM=spm_read_vols(V);    % seed maps

clear files V P;

%% calculate time courses
files=dir(fullfile(data_path,rundir,[data_filenm_mask,'*.nii']));

mn_roi.tc=zeros(size(sM,4),length(files));
for i_time=1:length(files),
    P=[fullfile(data_path,rundir,files(i_time).name),',1'];
    V=spm_vol(P);
    Y=spm_read_vols(V);
    for i_roi=1:size(sM,4),
        I=find(sM>0);
        mn_roi.tc(i_time) = mean(Y(I));
    end;
end;

tc = mn_roi.tc(:);
output_path = [data_path '/' outputdir];
if ~exist(output_path,'dir')
    mkdir(output_path); 
end
save(fullfile(data_path,outputdir,[mn_roi.name{:},'_roi_tc.mat']),'tc');

