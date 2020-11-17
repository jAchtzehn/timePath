import traits
import os
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.interfaces.matlab import MatlabCommand
from string import Template

def createSubjectlist(subjects):

	subjList = []
	
	for subj in subjects:
		subjList.append('sub-' + str(subj).zfill(2))
	
	return subjList


def conj_stat_maps(map1, map2, out_file=None):
	"""
	Runs a conjunction on stat maps (returning a map with the min t-stat)

	Creates a conjunction of two statistical maps (typically the result of
	EstimateContrast or Threshold from the spm interface).
	
	Args:
		map1 (str): filename of the first stat map in the conjunction
		map2 (str): filename of the second stat map in the conjunction
		
	Optional:
		out_file (str): output filename. If None (default), creates
			'conj_map.nii' in current directory
		
	Returns:
		out_file (str): output filename (absolute path)
	
	"""
	
	# Imports all required packages so the func can be run as a nipype
	# Function Node
	import nilearn.image as nli
	import os.path

	if out_file is None:
		out_file = 'conj_map.nii'
		
	conj_tmap = nli.math_img('np.minimum(img1*(img1>0), img2*(img2>0)) + ' +
							 'np.maximum(img1*(img1<0), img2*(img2<0))',
							 img1=map1,
							 img2=map2)
	
	conj_tmap.to_filename(out_file)
	return os.path.abspath(out_file)


class ClusterThreshInputSpec(BaseInterfaceInputSpec):

	in_file = File(exists=True, copyfile=True, mandatory=True)
	extent_threshold = traits.Int(mandatory=True)
	matlab_path = traits.Str('', usedefault=True)
	out_fname = traits.Str('conj_map_thr.nii', usedefault=True)


class ClusterThreshOutputSpec(TraitedSpec):
	
	out_file = File(exists=True)
	cluster_label_file = File()
	
	
class ClusterThresh(BaseInterface):

	'''
	Intended to be run as a Node after conjunction of thresholded maps (that
	already have voxel-level thresholding, assumes nan in all voxels that
	do not pass threshold)
	
	Examples
	--------
	>>> clusterThresh = ClusterThresh()
	>>> clusterThresh.inputs.in_file = '/in_path/conj_map.nii'
	>>> clusterThresh.inputs.extent_threshold = 10
	>>> clusterThresh.run()
	outputs: out_file (path to output .nii file with cluster-thresholded map)
	'''
	
	input_spec = ClusterThreshInputSpec
	output_spec = ClusterThreshOutputSpec
	
	def _run_interface(self, runtime):
		
		# Sets output files (in same dir)
		out_file = os.path.join(os.path.dirname(self.inputs.in_file),
						self.inputs.out_fname)
		cluster_label_file = os.path.join(os.path.dirname(self.inputs.in_file),
		                                  'cluster_labels.nii')
		
		d = dict(in_file=self.inputs.in_file, out_file=out_file,
		         cluster_label_file=cluster_label_file,
		         extent_threshold=self.inputs.extent_threshold,
		         matlab_path=self.inputs.matlab_path)
		
		# Matlab code that runs the cluster thresholding, based on spm
		# Threshold interface
		script = Template("""in_file='$in_file';
						  out_file = '$out_file';
						  cluster_label_file = '$cluster_label_file';
						  matlab_path = '$matlab_path';
						  extent_threshold = $extent_threshold;
						  
						  % If matlab SPM path is given - adds it to path
						  if (~isempty(matlab_path))
							  addpath(matlab_path);
						  end
						  
						  stat_map_vol = spm_vol(in_file);
						  [stat_map_data, stat_map_XYZmm] = ...
							  spm_read_vols(stat_map_vol);
						  
						  Z = stat_map_data(:)';
						  [x,y,z] = ind2sub(size(stat_map_data), ...
							  (1:numel(stat_map_data))');
						  XYZ = cat(1, x', y', z');
							
						  XYZth = XYZ(:, ~isnan(Z));
						  Zth = Z(~isnan(Z));
						  
						  max_size = 0;
						  max_size_index = 0;
						  th_nclusters = 0;
						  nclusters = 0;
						  if isempty(XYZth)
							  thresholded_XYZ = [];
							  thresholded_Z = [];
						  else
							  voxel_labels = spm_clusters(XYZth);
							  nclusters = max(voxel_labels);
							  
							  thresholded_XYZ = [];
							  thresholded_Z = [];
							  new_vox_labels = [];
							  
							  for i = 1:nclusters
								  cluster_size = sum(voxel_labels==i);
								  if cluster_size > extent_threshold
									  thresholded_XYZ = cat(2, ...
										thresholded_XYZ, ...
										XYZth(:,voxel_labels == i));
									  thresholded_Z = cat(2, ...
										thresholded_Z, Zth(voxel_labels == i));
									  th_nclusters = th_nclusters + 1;
									  new_vox_labels = cat(2, new_vox_labels,...
										th_nclusters*ones(1,cluster_size));
								  end
							  end
						  end
						  
						  if isempty(thresholded_XYZ)
							  thresholded_Z = [0];
							  thresholded_XYZ = [1 1 1]';
							  th_nclusters = 0;
						  end
						  
						  spm_write_filtered(thresholded_Z,thresholded_XYZ,...
							stat_map_vol.dim',stat_map_vol.mat,...
							'thresholded map', out_file);
						  
						  if (th_nclusters>0)
							  spm_write_filtered(new_vox_labels,...
								thresholded_XYZ,stat_map_vol.dim',...
								stat_map_vol.mat,'cluster_labels map',...
								cluster_label_file);
						  end
					  
						  exit;
						  """).substitute(d)
						
		mlab = MatlabCommand(script=script, mfile=True)
		result = mlab.run()
		return result.runtime
	
	def _list_outputs(self):
		outputs = self._outputs().get()
		outputs['out_file'] = os.path.join(os.path.dirname(
				self.inputs.in_file), self.inputs.out_fname)
		outputs['cluster_label_file'] = os.path.join(os.path.dirname(
				self.inputs.in_file), 'cluster_labels.nii')
		return outputs