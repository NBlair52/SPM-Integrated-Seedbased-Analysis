%% Created by Nicholas Blair
% with help from Dr. Jun Hua
% contact at NBlair@jhu.edu

clear; clc; format compact; warning('off', 'all');
%% Inputs
%=========================================================================
%Input data for the algorithm
%=========================================================================

%Path to data
%In general this script is meant to run on data that has all subjects in
% their own folders at a common path.  Ie
%Subject one:  /path/to/subject/subject1/
%Subject two:  /path/to/subject/subject2/.
%Control one:  /path/to/control/control_1/.
%Within each subjects folder should be subdirectories for anatomical data
%and functional data (this is how SPM outputs after preprocessing)

%Path to subject data
sub_path{1} = '/g4/kirby/jhua/SCZ7T/';
%Path to control data
sub_path{2} = '/g4/kirby/jhua/SCZ7T/';

%Path to functional data (in subject directory)
run_dir = 'run1';
%Path to anatomical data (in subject directory)
anat_dir = 'anat';
%Path to output data about each subjects connectivity (in subject directory)
seedbased_output_dir = 'nick_design_seedbased_analysis';

%Anatomical brain mask name (in subject's anatomical data directory)
brain_mask = 'brain_mask';

%Functional data file names prefix (in subject's functional data
%directory).  These prefixes are affixed by SPM during preprocessing.
%Check to see which preprocessed file you would like to use.
run_prefix = 'fsnfwa';

%Subject directories.  As mentioned before these are the names of the
%directories for each subject, located under sub_path{ }.  Within each
%directory each subject should have anatomical and run directories with the
%appropriate data in each.
sub_id{1,:} = [{'7001'},{'7002'},{'7003'},{'7004'},{'7005'},{'7006'},{'7007'},{'7008'},{'7009'},{'7010'},{'7011'},{'7012'},{'7013'},{'7014'}];
%Control files
sub_id{2,:} = [{'8001'},{'8002'},{'8003'},{'8004'},{'8006'},{'8007'},{'8008'},{'8009'},{'8010'},{'8011'},{'8012'},{'8013'},{'8014'},{'8015'}];
    

%Seed data
%Seed files should be individual .nii files with only values at 
%the desired seed region (ie all 1 for desired seed, 0 or NaN for all
%else).  All seed files should be located at the same path.
%Seed path
seed_path = '/home/nblair2/matlab/fmri/connectivity/seedbased/seeds/Thalamus/';

%Seed files
seed_id = [{'Thalamus-maxprob-thr25-2mm-1p2-spmvol'},...
           {'Thalamus-maxprob-thr25-2mm-3p4-spmvol'},...
           {'Thalamus-maxprob-thr25-2mm-5p6-spmvol'},...
           {'Thalamus-maxprob-thr25-2mm-7p-spmvol'}];

%Output path
output_path = '/g4/kirby/jhua/SCZ7T/nick_design_seedbased_2ndlevel_group/';

%Seed description (ie 'thalamus' or 'striatum')
seed_description = 'nick_thalamus';


%%Atlas data
%Atlas of regions to compare seed region correlation too.  Each seed will
%individually be compared to each region in the atlas for subject > control
%or control > subject connectivity.  Atlas should be a single .nii file
%with each voxel of a region containing a single value, and each region
%determined by what value it contains,and any area to disregard either 0.
%Ie if you only want to examine two areas, all voxels in area one
%have value 1, all voxels in area two have value 2, and all other voxels
%have value 0.
%These atlasses can easily be constructed using wfu_pickatlas, with 'write
%seperate regions' option checked.  Once an atlas is constructed using
%wfu_pickatlas' 'generate table' option generates a table which will be
%read in to label each region.  As of now, this script will not work if you
%do not use wfu_pickatlas to generate the atlas AND the .tbl file
%Atlas path
atlas_path = '/home/nblair2/matlab/fmri/connectivity/seedbased/seeds/IBASPM116/';
%Atlas file (should be in above directory, .nii)
atlas_id = 'IBASPM116_all_independant';
%table file (should be in same directory, .tbl)
atlas_tbl = 'IBASPM116_all_sig_regions';


%=========================================================================
% Options
%=========================================================================
%Options for the program to run, mainly in SPM t-test.  You can find these
%options further down under 'factorial_job{1}.spm.stats.factorial_design'
%and results_XX_job{1}.spm.stats.results.  Many additional options exist
%but have been taken out of the options here.  Using SPM matlabbatch
%editior you can create custom options and edit this file
%Aditional options include:
%Grand Mean Scaling and ANCOVA for PET data, adding covariates,
%threshholding image values (either relative or absolute), implicit
%masking, global effect calculation for PET data, grand mean scaling,
%normalization

%t-test options
%Independance of groups (0 - independant, 1 - dependant)
factorial_grp_dependance = 0;
%Variance between groups (1 - unequal, 0 - equal)
factorial_grp_variance = 1;
%Explicit mask ('' - none, '/path/to/mask/mask.nii - mask)
factorial_explicit_mask = '';

%threshold options for results
%Threshold type ('none' - none, 'FWE' - family-wise error)
results_threshold_type = 'none';
%Maximum p value
results_threshold_pval = 0.05;
%Minimum number of voxels in a group
results_extent_vox = 40;

%=========================================================================
%OPTIONS AND DATA INPUT ENDS HERE.
%=========================================================================
%Helper variables
sub_num(1) = length(sub_id{1});
sub_num(2) = length(sub_id{2});
seed_num = length(seed_id);

contrast_name = [{'D>C'}, {'C>D'}];
contrast_tmap_id = [{'spmT_0001'}, {'spmT_0002'}];

% Contsruct list of atlas regions
atlas_text = textread(fullfile(atlas_path, [atlas_tbl '.tbl']), '%s');
atlas_regions = cell(1, (size(atlas_text, 1) - 18)/ 21);
i_atlas_region_cells = 1;
for i_text_cells = 26:21:size(atlas_text, 1)
    atlas_regions{i_atlas_region_cells} = atlas_text{i_text_cells};
    i_atlas_region_cells = i_atlas_region_cells + 1;
end
atlas_trimmed_regions = cell(size(atlas_regions));
atlas_sides = atlas_trimmed_regions;
num_atlas_regions = size(atlas_regions, 2);
for i_atlas_regions = 1:num_atlas_regions
    region = atlas_regions{i_atlas_regions};
    side = region(length(region));
    if ((side ~= 'L') && (side ~= 'R'))
        side = 'N/A';
    else
        region = region(1:end-2);
    end
    atlas_trimmed_regions{i_atlas_regions} = region;
    atlas_sides{i_atlas_regions} = side;
end

%% Compute Time Course
disp('Computing time course and correlation');
for i_sub_grp = 1:2 %for each group
    for i_sub_num = 1:sub_num(i_sub_grp) %for each subject
        
        disp(['-Computing correlation for subject ' sub_id{i_sub_grp}{i_sub_num}]);
        tic; elapsed_time = 0;
        
        for i_seed_num = 1:seed_num %for each seed
            disp(['--Seed ', seed_id{i_seed_num}]);
            %Compute time courses for seed map
            roi_dm_tc3(fullfile(sub_path{i_sub_grp},sub_id{i_sub_grp}{i_sub_num}),...
                       run_dir, seedbased_output_dir, seed_path, run_prefix,seed_id(i_seed_num));
            %Calculate correlation to the seeds
            conn(fullfile(sub_path{i_sub_grp},sub_id{i_sub_grp}{i_sub_num}),...
                 run_dir, anat_dir, run_prefix, seed_id(i_seed_num),seedbased_output_dir, brain_mask);

            elapsed_time = toc - elapsed_time;
            disp(['--Done in ', num2str(elapsed_time/60), 'min']);
        end
        
        disp(['--Merging Correlation Maps for ', sub_id{i_sub_grp}{i_sub_num}]);
        
        % average network maps
        merge_connmaps(fullfile(sub_path{i_sub_grp},sub_id{i_sub_grp}{i_sub_num}),...
                       anat_dir, seedbased_output_dir, seed_id, seed_description,...
                       brain_mask);
        
        elapsed_time = toc;
        disp(['-' sub_id{i_sub_grp}{i_sub_num} ' data has been processed in ', num2str(elapsed_time/60), 'min']);
    end
end

disp('Performing 2nd level statisical analysis')
%% Perform 2nd level group statistical analysis
for i_seed_num = 1:seed_num %for each seed
    disp(['-Computing analysis for seed: ' seed_id{i_seed_num}]);
    %Path to output data
    seed_output_path = fullfile(output_path, seed_id{i_seed_num}, filesep);
    %SPM file 
    SPM_file = fullfile(seed_output_path, 'SPM.mat');
    %Name of connectivity map for this seed
    seed_NaN_conn_map = ['conn_NaN_' seed_id{i_seed_num} '_roi_tc.nii'];
    
    mkdir(seed_output_path);
    cd(seed_output_path);
    
    %Outline the matlabbatch jobs.  See SPM documentation to learn more
    % Factorial job
    factorial_job = [];
    factorial_job{1}.spm.stats.factorial_design.dir = {seed_output_path};
    for j_sub_num1 = 1:sub_num(1)
        factorial_job{1}.spm.stats.factorial_design.des.t2.scans1{j_sub_num1} =...
            fullfile(sub_path{1}, sub_id{1}{j_sub_num1}, seedbased_output_dir, seed_NaN_conn_map);
    end
    for j_sub_num2 = 1:sub_num(2)
        factorial_job{1}.spm.stats.factorial_design.des.t2.scans2{j_sub_num2} =...
            fullfile(sub_path{2}, sub_id{2}{j_sub_num2}, seedbased_output_dir, seed_NaN_conn_map);
    end
    factorial_job{1}.spm.stats.factorial_design.des.t2.dept = factorial_grp_dependance;
    factorial_job{1}.spm.stats.factorial_design.des.t2.variance = factorial_grp_variance;
    factorial_job{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    factorial_job{1}.spm.stats.factorial_design.des.t2.ancova = 0;
    factorial_job{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    factorial_job{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    factorial_job{1}.spm.stats.factorial_design.masking.im = 1;
    factorial_job{1}.spm.stats.factorial_design.masking.em = {factorial_explicit_mask};
    factorial_job{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    factorial_job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    factorial_job{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    save(fullfile(seed_output_path,'factorial_job.mat'),'factorial_job');

    % Estimation job
    estimation_job{1}.spm.stats.fmri_est.spmmat = {SPM_file};
    estimation_job{1}.spm.stats.fmri_est.method.Classical = 1;
    save(fullfile(seed_output_path,'estimation_job.mat'),'estimation_job');

    % Contrast job
    contrast_job = [];
    contrast_job{1}.spm.stats.con.spmmat = {SPM_file};
    contrast_job{1}.spm.stats.con.consess{1}.tcon.name = contrast_name{1};
    contrast_job{1}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
    contrast_job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    contrast_job{1}.spm.stats.con.consess{2}.tcon.name = contrast_name{2};
    contrast_job{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
    contrast_job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    contrast_job{1}.spm.stats.con.delete = 0;
    save(fullfile(seed_output_path,'contrast_job.mat'),'contrast_job');

    % Results for contrast 1
    results_DC_job = [];
    results_DC_job{1}.spm.stats.results.spmmat = {SPM_file};
    results_DC_job{1}.spm.stats.results.conspec(1).titlestr = 'D>C_results';
    results_DC_job{1}.spm.stats.results.conspec(1).contrasts = 1;
    results_DC_job{1}.spm.stats.results.conspec(1).threshdesc = results_threshold_type;
    results_DC_job{1}.spm.stats.results.conspec(1).thresh = results_threshold_pval;
    results_DC_job{1}.spm.stats.results.conspec(1).extent = results_extent_vox;
    results_DC_job{1}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
    results_DC_job{1}.spm.stats.results.units = 1;
    results_DC_job{1}.spm.stats.results.print = false;
    save(fullfile(seed_output_path,'results_DC_job.mat'),'results_DC_job');

    % Results for contrast 2
    results_CD_job = [];
    results_CD_job{1}.spm.stats.results.spmmat = {SPM_file};
    results_CD_job{1}.spm.stats.results.conspec(1).titlestr = 'C>D_results';
    results_CD_job{1}.spm.stats.results.conspec(1).contrasts = 2;
    results_CD_job{1}.spm.stats.results.conspec(1).threshdesc = results_threshold_type;
    results_CD_job{1}.spm.stats.results.conspec(1).thresh = results_threshold_pval;
    results_CD_job{1}.spm.stats.results.conspec(1).extent = results_extent_vox;
    results_CD_job{1}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
    results_CD_job{1}.spm.stats.results.units = 1;
    results_CD_job{1}.spm.stats.results.print = false;
    save(fullfile(seed_output_path,'results_CD_job.mat'),'results_CD_job');
    
    disp('--Running SPM jobs');

    %Run jobs
    %Run factorial (ie t-test)
    spm_jobman('run', factorial_job);
    %Estimate the model
    spm_jobman('run', estimation_job);
    %Insert contrasts (ie Disease > Control, and Control > Disease)
    spm_jobman('run', contrast_job);
    %Run results (compute significant clusters
    for i_contrast_num = 1:2 % for each contrast
        %Calculate resulting significant clusters
        spm('defaults','FMRI')
        if (i_contrast_num == 1)
            disp('---D>C contrast');
            spm_jobman('run', results_DC_job);
        elseif (i_contrast_num == 2)
            disp('---C>D contrast');
            spm_jobman('run', results_CD_job);
        end
        
        xSPM = evalin('base','xSPM;');
        Z       = spm_clusters(xSPM.XYZ);
        num     = max(Z);
        [n, ni] = sort(histc(Z,1:num), 2, 'descend');
        n       = size(ni);
        n(ni)   = 1:num;
        Z       = n(Z);
        %Save n-ary file of significant clusters
        roi_nary_id = [seed_id{i_seed_num} '_' contrast_name{i_contrast_num} '_p<' num2str(results_threshold_pval) '_c>' num2str(results_extent_vox) '_nary'];
        spm_write_filtered(Z, xSPM.XYZ, xSPM.DIM, xSPM.M, '', roi_nary_id);
        
        %% Compute ROI results
        % Load cluster data
        V = load_nii(fullfile(seed_output_path, roi_nary_id));
        mask_cluster = V.img;
        mask_cluster = uint16(mask_cluster);
        % Load atals data
        V = load_nii(fullfile(atlas_path, [atlas_id '.nii']));
        mask_atlas = V.img;
        atlas_region_num = max(mask_atlas(:));
        %pick only voxels in both atlas and cluster
        mask_atlas(mask_cluster==0) = 0;
        %load tmap
        [~, xyzmm] = spm_read_vols(spm_vol(fullfile(seed_output_path, [contrast_tmap_id{i_contrast_num} '.img'])));
        V = load_nii(fullfile(seed_output_path, [contrast_tmap_id{i_contrast_num} '.img']));
        contrast_tmap = V.img;
        temp = contrast_tmap;
        temp(isnan(temp)) = 0;
        voxel_num = sum(sum(sum(temp~=0)));
        clear V;
        
        % calculate peak x,y,z and peak/mean/std t and p for each seed in each contrast
        % setup matrices to hold data
        voxel_roi_num = zeros(atlas_region_num,1); 
        xyzmm_roi_peakt = zeros(atlas_region_num,4);
        t_roi_max = voxel_roi_num;
        t_roi_mean = voxel_roi_num;
        t_roi_std = voxel_roi_num;
        P_corr_min = voxel_roi_num;
        P_corr_mean = voxel_roi_num;
        P_corr_std = voxel_roi_num;
        p_uncorr_min = voxel_roi_num;
        p_uncorr_mean = voxel_roi_num;
        p_uncorr_std = voxel_roi_num;
        
        % calculate data
        for i_atlas_region_num = 1:atlas_region_num %for each region in atlas
            tmp = contrast_tmap;
            %Pick only voxels in this region of atlas
            tmp(mask_atlas ~= i_atlas_region_num) = 0;
            %Calculate location of max t value in area
            xyzmm_roi_peakt1 = xyzmm(:,tmp==max(max(max(tmp))))';
            tmp = tmp(mask_atlas == i_atlas_region_num);  % cannot do this earlier
            tmp = tmp(~isnan(tmp)); 
            voxel_roi_num(i_atlas_region_num) = length(tmp(:));

            if voxel_roi_num(i_atlas_region_num) > 0 %If voxels exist in area
                xyzmm_roi_peakt(i_atlas_region_num,1:3) = xyzmm_roi_peakt1(1,:);  % if multiple peaks, pick the first one
                xyzmm_roi_peakt(i_atlas_region_num,4) = size(xyzmm_roi_peakt1,1);  % number of peaks

                t_roi_mean(i_atlas_region_num) = mean(tmp(:));
                t_roi_max(i_atlas_region_num) = max(tmp(:));
                t_roi_std(i_atlas_region_num) = std(tmp(:));

                tmp = tmp(:);
                for itmp = 1:length(tmp)
                    [corrP(itmp) unp(itmp)] = spm_P(1,0,tmp(itmp),[voxel_num-1 2*voxel_num-2],'T',1,1,voxel_num);  
                    % df=2n-2 for two-sample, n-1 for one-sample; how is df(2) defined?? R=1???
                end
                P_corr_mean(i_atlas_region_num) = mean(corrP(:));
                P_corr_min(i_atlas_region_num) = min(corrP(:));
                P_corr_std(i_atlas_region_num) = std(corrP(:));
                p_uncorr_mean(i_atlas_region_num) = mean(unp(:));
                p_uncorr_min(i_atlas_region_num) = min(unp(:));
                p_uncorr_std(i_atlas_region_num) = std(unp(:));
            end %if voxel_roi_num(i_atlas_region_num) > 0
        end %for i_atlas_region_num = 1:atlas_region_num
        
        misc_output_data = [voxel_roi_num xyzmm_roi_peakt t_roi_max t_roi_mean t_roi_std ...
                            P_corr_min P_corr_mean P_corr_std p_uncorr_min p_uncorr_mean p_uncorr_std];
        
        % calculate average connectivity for each subject in each atlas region to each seed region
        conn_output_data = zeros(double(sub_num(1) + sub_num(2)), double(atlas_region_num)*3); 
        for j_sub_grp = 1:2 %for each group
            conn_roi = zeros(sub_num(j_sub_grp), atlas_region_num);
            conn_roi_std = conn_roi;
            conn_roi_num = conn_roi;
            for j_sub_num = 1:sub_num(j_sub_grp)%for each subject
                V = load_nii(fullfile(sub_path{j_sub_grp}, sub_id{j_sub_grp}{j_sub_num}, seedbased_output_dir, seed_NaN_conn_map));
                map_conn = V.img;
                %Calculate average conn for each region
                for j_atlas_region_num = 1:atlas_region_num
                    tmp = map_conn(mask_atlas == j_atlas_region_num);
                    tmp = tmp(~isnan(tmp));
                    conn_roi(j_sub_num, j_atlas_region_num) = mean(tmp(:));
                    conn_roi_std(j_sub_num, j_atlas_region_num) = std(tmp(:));
                    conn_roi_num(j_sub_num, j_atlas_region_num) = length(tmp(:));
                end
            end
            clear V;
            
            %size groups
            row_start = 1;
            row_end = sub_num(1);
            if (j_sub_grp > 1);
                for jgrp = 2:j_sub_grp
                    row_start = row_start + sub_num(jgrp-1);
                    row_end = row_end + sub_num(jgrp);
                end
            end
            conn_output_data(row_start:row_end,1:3:end) = conn_roi_num;
            conn_output_data(row_start:row_end,2:3:end) = conn_roi;
            conn_output_data(row_start:row_end,3:3:end) = conn_roi_std;
        end %for j_sub_grp = 1:2
        
        %Output data to a .csv file.  
        output_file_name = [seed_id{i_seed_num} '_' contrast_name{i_contrast_num} '_p<' num2str(results_threshold_pval) '_c>' num2str(results_extent_vox) '_results.csv'];
        output_file_fid = fopen(fullfile(seed_output_path, output_file_name), 'w');
        %titles for each collumn
        fprintf(output_file_fid, 'region, hemisphere, size, disease mean, disead std, control mean, control std, p-value, size, x, y, z, num clusters, t max, t mean, t std, P min, P mean, P std, p_uncorr min, p_uncorr mean, p_uncorr std\n');
        for j_atlas_regions = 1:num_atlas_regions
            
            %Compute region conn data
            vox_size = mean(conn_output_data(1:end, (j_atlas_regions-1)*3+1));
            d_mean = mean(conn_output_data(1:sub_num(1), (j_atlas_regions-1)*3+2));
            d_std = std(conn_output_data(1:sub_num(1), (j_atlas_regions-1)*3+2));
            c_mean = mean(conn_output_data(sub_num(1)+1:end, (j_atlas_regions-1)*3+2));
            c_std = std(conn_output_data(sub_num(1)+1:end, (j_atlas_regions-1)*3+2));
            [h, p_val] = ttest2(conn_output_data(1:sub_num(1), (j_atlas_regions-1)*3+2),...
                                conn_output_data(sub_num(1)+1:end, (j_atlas_regions-1)*3+2),...
                                0.05, 'both','unequal');
            
            %Print
            if (vox_size > 0)
                %Print region name and hemisphere
                fprintf(output_file_fid, '%s, %s, ', atlas_trimmed_regions{j_atlas_regions}, atlas_sides{j_atlas_regions});
                %Print conn and misc data
                output_vals = [vox_size, d_mean, d_std, c_mean, c_std, p_val, misc_output_data(j_atlas_regions, :)];
                fprintf(output_file_fid, [repmat('%f, ', 1, length(output_vals)-1), '%f\n'], output_vals);
            end
        end
        fclose(output_file_fid);
    end %for i_contrast_num = 1:2
end %for i_seed_num = 1:seed_num