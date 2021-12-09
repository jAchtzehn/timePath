% % The matlab batch script for non-parametric permutation with SnPM.
% 
% % Authors : Qi WANG (qiqi.wang@lis-lab.fr)

conditions = {'time-vs-dist'};
fwe = 0.05;
Nperm = 1000;
k = 10;          % minimal cluster size

for condition = conditions
    resultdir = sprintf('/home/achtzehnj/data/timePath/derivatives/nilearn/group_results_ispa_std/%s', condition{1});
    snpm_dirname = sprintf('snpm_batch');
    outputdir = sprintf('%s/%s', resultdir, snpm_dirname);
    snpmfig = sprintf('%s/%s', outputdir, snpm_dirname);
    % List of score maps

    splits = [0:23];
    files=[];
    for split = splits
        input_filename = sprintf('ispa_sl_wb_decoding-%s_split-%02d', condition{1}, split);
        files = [files; cellstr(sprintf('%s/%s.nii', resultdir, input_filename))];
    end

    job_id = 1;
    % Matlab batch script
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.P = files;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.dir = {outputdir};
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.nPerm = Nperm;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.vFWHM = [6 6 6];
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.bVolm = 1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.ST.ST_later = -1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.masking.im = 1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.masking.em = {''};
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
    matlabbatch{job_id}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
    job_id = job_id + 1;
    matlabbatch{job_id}.spm.tools.snpm.cp.snpmcfg(1) = cfg_dep('MultiSub: One Sample T test on diffs/contrasts: SnPMcfg.mat configuration file', substruct('.','val', '{}',{job_id-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','SnPMcfg'));
    job_id = job_id + 1;
    matlabbatch{job_id}.spm.tools.snpm.inference.SnPMmat(1) = cfg_dep('Compute: SnPM.mat results file', substruct('.','val', '{}',{job_id-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','SnPM'));
    matlabbatch{job_id}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = fwe;
    matlabbatch{job_id}.spm.tools.snpm.inference.Tsign = 1;
    matlabbatch{job_id}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered';
    matlabbatch{job_id}.spm.tools.snpm.inference.Report = 'MIPtable';
    job_id = job_id + 1;

    matlabbatch{job_id}.spm.util.print.fname = snpmfig;
    matlabbatch{job_id}.spm.util.print.fig.fighandle = NaN;
    matlabbatch{job_id}.spm.util.print.opts = 'png';
    job_id = job_id + 1;

    spm_jobman('run',matlabbatch);


    % cluster removal from https://en.wikibooks.org/wiki/SPM/How-to#How_to_remove_clusters_under_a_certain_size_in_a_binary_mask?
    ROI  = sprintf('%s/SnPM_filtered.nii', outputdir);  % input image (binary, ie a mask)
    ROIf = sprintf('%s/SnPM_filtered_binary_mask_cluster.nii', outputdir); % output image (filtered on cluster size)

    %-Connected Component labelling
    V = spm_vol(ROI);
    dat = spm_read_vols(V);
    [l2, num] = spm_bwlabel(double(dat>0),26);
    if ~num, warning('No clusters found.'); end

    %-Extent threshold, and sort clusters according to their extent
    [n, ni] = sort(histc(l2(:),0:num), 1, 'descend');
    l  = zeros(size(l2));
    n  = n(2:end); ni = ni(2:end)-1;
    ni = ni(n>=k); n  = n(n>=k);
    for i=1:length(n), l(l2==ni(i)) = i; end

    fprintf('Selected %d clusters (out of %d) in image.\n',length(n),num);

    %-Write new image
    V.fname = ROIf;
    spm_write_vol(V,l~=0); % save as binary image. Remove '~=0' so as to
                           % have cluster labels as their size. 
                           % or use (l~=0) if output image should be binary
end
