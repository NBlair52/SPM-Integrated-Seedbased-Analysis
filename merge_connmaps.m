% To merge correlation maps
% created by Anita Barber
% adapted by Jun Hua

function merge_connmaps(data_path, anatdir, seedbased_output_dir, roi_filenm_mask, net, brain_mask_filenm)

% load brain mask
vdata = spm_vol(fullfile(data_path,anatdir,[brain_mask_filenm,'.nii']));
v0 = spm_read_vols(vdata);
tseries = zeros(sum(v0(:)>0), 1);

% load correlation maps
for iimg = 1:length(roi_filenm_mask)
    fntmp = strcat('conn_',roi_filenm_mask{iimg},'_roi_tc.nii');
    imgnm = fullfile(data_path,seedbased_output_dir,fntmp{1});
    temp = spm_read_vols(spm_vol(imgnm));
    temp = temp(:);
    roits(:,iimg) = temp(v0(:)>0);
end

% average correlation maps
for ivox = 1:size(tseries) 
    tseries(ivox,1) = mean(roits(ivox,:));
end

% write average map
result = v0;
result(result~=0) = tseries;
v1 = vdata;
v1.fname = fullfile(data_path,seedbased_output_dir,['conn_mean_',net,'_roi_tc.nii']);
v1.private.dat.fname=v1.fname;
v1.dt=[16 0];
spm_write_vol(v1, result);

% write average map with NaN
result = v0;
result(v0~=0) = tseries;
result(v0==0) = NaN;
v1 = vdata;
v1.fname = fullfile(data_path,seedbased_output_dir,['conn_mean_NaN_',net,'_roi_tc.nii']);
v1.private.dat.fname=v1.fname;
v1.dt=[16 0];
spm_write_vol(v1, result);
