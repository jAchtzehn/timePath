experiment_dir = fullfile('/home/achtzehnj/data/timePath/')
rsa_dir = fullfile(experiment_dir, 'derivatives', 'rsa')

% cluster removal from https://en.wikibooks.org/wiki/SPM/How-to#How_to_remove_clusters_under_a_certain_size_in_a_binary_mask?
ROI  = sprintf('%s/corrImg_rel-time_irrel-dist_p_odd_cluster_mask.nii', rsa_dir);  % input image (binary, ie a mask)

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