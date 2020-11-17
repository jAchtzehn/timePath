from os.path import join as opj
from os.path import abspath
from sys import platform
import os
from mvpa2.suite import fmri_dataset, zscore, mean_group_sample, poly_detrend
from mvpa2.datasets.sources.native import load_tutorial_data, load_example_fmri_dataset
from mvpa2.measures.searchlight import sphere_searchlight
from mvpa2.base.learner import ChainLearner
from mvpa2.mappers.shape import TransposeMapper
import numpy as np
from mvpa2.measures import rsa
import pandas as pd
import pylab as pl
from scipy.spatial.distance import squareform
from os.path import join as pjoin
import nipy
from mvpa2 import pymvpa_dataroot


# plot datamatrices
def plot_mtx(mtx, labels, title):
    pl.figure(figsize=(8, 8))
    pl.imshow(mtx, interpolation='nearest')
    pl.xticks(range(len(mtx)), labels, rotation=-45)
    pl.yticks(range(len(mtx)), labels)
    pl.title(title)
    pl.clim((0, 1))
    pl.colorbar()
    pl.show()
    


datapath = pymvpa_dataroot
ds = load_tutorial_data(path=datapath, roi=(15, 16, 23, 24, 36, 38, 39, 40, 48))

poly_detrend(ds, polyord=1, chunks_attr='chunks')
zscore(ds, chunks_attr='chunks', param_est=('targets', 'rest'))
ds = ds[ds.sa.targets != 'rest']


mtgs = mean_group_sample(['targets'])
mtds = mtgs(ds)

dsm = rsa.PDist(square=True)
res = dsm(mtds)
plot_mtx(res, mtds.sa.targets, 'ROI pattern correlation distances')

# mtcgs = mean_group_sample(['targets', 'chunks'])
# mtcds = mtcgs(ds)
#
# dscm = rsa.PDistConsistency()
# sl_cons = sphere_searchlight(dscm, 2)
# slres_cons = sl_cons(mtcds)
# mean_consistency = np.mean(slres_cons, axis=0)
# print('Most stable pattern: %s' % mtcds.fa.voxel_indices[mean_consistency.argmax()])
#
# plot_mtx(squareform(slres_cons.samples[:, mean_consistency.argmax()]), mtcds.sa.targets,
#          'Most consistent searchlight pattern correlation distances')

# tdsm = rsa.PDistTargetSimilarity(slres_cons.samples[:, mean_consistency.argmax()])
#
# sl_tdsm = sphere_searchlight(ChainLearner([tdsm, TransposeMapper()]), 2)
# slres_tdsm = sl_tdsm(mtcds)


# # plot the spatial distribution using NiPy
# vol = ds.a.mapper.reverse1(slres_tdsm.samples[0])
# import nibabel as nb
# anat = nb.load(pjoin(datapath, 'sub001', 'anatomy', 'highres001.nii.gz'))
#
# from nipy.labs.viz_tools.activation_maps import plot_map
# pl.figure(figsize=(15, 4))
# sp = pl.subplot(121)
# pl.title('Distribution of target similarity structure correlation')
# slices = plot_map(
#             vol,
#             ds.a.imgaffine,
#             cut_coords=np.array((12, -42, -20)),
#             threshold=.5,
#             cmap="bwr",
#             vmin=0,
#             vmax=1.,
#             axes=sp,
#             anat=anat.get_data(),
#             anat_affine=anat.affine,
#          )
# img = pl.gca().get_images()[1]
# cax = pl.axes([.05, .05, .05, .9])
# pl.colorbar(img, cax=cax)
#
# sp = pl.subplot(122)
# pl.hist(slres_tdsm.samples[0],
#         normed=False,
#         bins=30,
#         color='0.6')
# pl.ylabel("Number of voxels")
# pl.xlabel("Target similarity structure correlation")
