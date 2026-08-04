# PVH SCENIC+ pipeline
# Assume you already built your cisTopic object
# This script generates topic/DAR region sets, motif enrichment results, and the SCENIC+ object

import os
import pickle
import re
import traceback

import dill
import numpy as np
import pandas as pd
import pyranges as pr
import scanpy as sc

from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycistarget.utils import region_names_to_coordinates
from pycistarget.motif_enrichment_cistarget import *
from scenicplus.wrappers.run_pycistarget import *
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus


# ---------------------------
# Define working directory
# ---------------------------

# Root working directory for SCENIC+ outputs
work_dir = "/path/to/PVH_SCENIC/"


# ---------------------------
# Topic binarization
# ---------------------------

# Load cistopic object

cistopic_obj = pickle.load(
    open(
        os.path.join(
            '/path/to/pvh_updated_cistopic_obj.pkl'
        ),
        'rb'
    )
)


# Binarize topics using Otsu thresholding
# Converts continuous region-topic probabilities into binary region sets
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')

# Binarize topics by selecting the top 3000 regions per topic
# Keeps highest-scoring regions for each topic
region_bin_topics_top3k = binarize_topics(
    cistopic_obj,
    method='ntop',
    ntop=3000
)


# ---------------------------
# Differential accessibility workflow
# ---------------------------

# Impute chromatin accessibility per cell from topic model
# selected_cells=None -> use all cells
# selected_regions=None -> use all regions
# scale_factor=10^6 -> normalize to 1,000,000 counts
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

# Normalize imputed accessibility scores
# scale_factor=10^4 -> rescale to 10,000 (analogous to RNA normalization)
normalized_imputed_acc_obj = normalize_scores(
    imputed_acc_obj,
    scale_factor=10**4
)

# Identify highly variable regions across cells
# plot=False -> suppress diagnostic variability plot
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    plot=False
)

# Identify differentially accessible regions between clusters
# variable='cluster' -> group comparison variable
# var_features=variable_regions -> restrict to variable regions
# split_pattern='-' -> parsing rule for cluster names if needed
markers_dict = find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='cluster',
    var_features=variable_regions,
    split_pattern='-'
)



# ---------------------------
# Ensure output directory exists
# ---------------------------

# If candidate enhancer directory does not exist, create it
if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))


# ---------------------------
# Save results to disk
# ---------------------------

# Save Otsu-binarized topic regions
pickle.dump(
    region_bin_topics_otsu,
    open(
        os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'),
        'wb'
    )
)

# Save top-3000 binarized topic regions
pickle.dump(
    region_bin_topics_top3k,
    open(
        os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'),
        'wb'
    )
)

# Save differential accessibility markers dictionary
pickle.dump(
    markers_dict,
    open(
        os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'),
        'wb'
    )
)


# ---------------------------
# Motif enrichment
# ---------------------------

# Use results generated in the binarization step
binarized_topic_region = region_bin_topics_otsu
DARs_dict = markers_dict

# ---------------------------
# Prepare region sets for pycistarget
# ---------------------------

# Initialize dictionary that will contain region sets
region_sets = {}


# ---------------------------
# Format region sets
# ---------------------------

# Convert topic-based regions into PyRanges objects
# Each topic becomes a named region set
region_sets['Topics'] = {key: pr.PyRanges(region_names_to_coordinates(binarized_topic_region[key].index.tolist())) for key in binarized_topic_region.keys()}


# Convert DARs into PyRanges objects
# Clean keys to ensure safe naming (replace non-alphanumeric characters with "_")
region_sets['DARs'] = {
    re.sub('[^A-Za-z0-9]+', '_', str(key)): pr.PyRanges(region_names_to_coordinates(DARs_dict[key].index.tolist()))
    for key in DARs_dict.keys()
}

# Print available region set names for verification
for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

# ---------------------------
# Run motif enrichment
# ---------------------------

motif_enrichment_dir = os.path.join(work_dir, "motif_enrichment_results")

run_pycistarget(region_sets,
                ctx_db_path = '/path/to/pycisTopic/databases/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather',
                species = 'mus_musculus',
                save_path = motif_enrichment_dir,
                dem_db_path = '/path/to/pycisTopic/databases/mm10_screen_v10_clust.regions_vs_motifs.scores.feather',
                run_without_promoters = True, # Also generate promoter-filtered motif enrichment results; ; this was required for this dataset to produce usable cistromes for SCENIC+
                biomart_host = 'http://www.ensembl.org',
                promoter_space = 500, #All default parameters
                ctx_auc_threshold = 0.005,
                ctx_nes_threshold = 3.0,
                ctx_rank_threshold = 0.05,
                dem_log2fc_thr = 0.5,
                dem_motif_hit_thr = 3.0,
                dem_max_bg_regions = 500,
                path_to_motif_annotations = '/path/to/pycisTopic/databases/motifs-v10-nr.mgi-m0.00001-o0.0.tbl',
                annotation_version = 'v10nr_clust',
                annotation = ['Direct_annot', 'Orthology_annot'],
                n_cpu = 1)


# ---------------------------
# Convert regions to bed
# ---------------------------

# Define output base directory
out_dir = "/path/to/new_region_sets"

# Create Topics directory and write .bed files
os.makedirs(f"{out_dir}/Topics", exist_ok=True)
for topic, df in binarized_topic_region.items():
    coords = region_names_to_coordinates(df.index.tolist())
    bed = pd.DataFrame(coords)
    bed = bed[['Chromosome', 'Start', 'End']]
    bed.to_csv(f"{out_dir}/Topics/{topic}.bed", sep='\t', header=False, index=False)

# Create DARs directory and write .bed files
os.makedirs(f"{out_dir}/DARs", exist_ok=True)
for key, df in markers_dict.items():
    name = re.sub('[^A-Za-z0-9]+', '_', str(key))
    coords = region_names_to_coordinates(df.index.tolist())
    bed = pd.DataFrame(coords)
    bed = bed[['Chromosome', 'Start', 'End']]
    bed.to_csv(f"{out_dir}/DARs/{name}.bed", sep='\t', header=False, index=False)


# ---------------------------
# Modify scanpy object (AnnData)
# ---------------------------

adata_path = "/path/to/new_pvh_20250218.h5ad"

# Load AnnData
adata = sc.read_h5ad(adata_path)

# Create split_clusters = cluster_condition 
adata.obs["split_clusters"] = (
    adata.obs["cluster"].astype(str) + "_" + adata.obs["condition"].astype(str)
)

# Filter out rows where split_clusters contains "NA"
adata = adata[~adata.obs["split_clusters"].str.contains("NA")].copy()

# If adata.raw exists, remove a problematic '_index' column from raw.var (if present)
if adata.raw is not None:
    raw_var = adata.raw.var.copy()
    if "_index" in raw_var.columns:
        raw_var = raw_var.drop(columns=["_index"])
        adata.raw._var = raw_var  # overwrite raw.var

# Save once (after all modifications)
adata.write(adata_path)


# ---------------------------
# Modify cisTopic object
# ---------------------------

cistopic_path = "/path/to/pvh_updated_cistopic_obj.pkl"

# Add split_clusters = cluster_group
cistopic_obj.cell_data["split_clusters"] = (
    cistopic_obj.cell_data["cluster"].astype(str) + "_" + cistopic_obj.cell_data["group"].astype(str)
)

# Save
with open(cistopic_path, "wb") as f:
    pickle.dump(cistopic_obj, f)


# ---------------------------
# Create SCENICPLUS object
# ---------------------------

with open(os.path.join(motif_enrichment_dir, "menr.pkl"), "rb") as f:
    menr = dill.load(f)
    
    
# Harmonize cell barcodes with cisTopic/scATAC side if needed
adata.obs_names = (adata.obs_names + "-PVH").str.replace("_1-PVH", "-PVH")

# Build GEX anndata from raw
test_data = adata.raw.to_adata()

# Normalize var naming into a clean index
# (This assumes the current var has a single column that represents gene names)
test_data.var.columns = ["genes"]
test_data.var.set_index("genes", inplace=True)

# Create SCENICPLUS object
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata=test_data,
    cisTopic_obj=cistopic_obj,
    menr=menr
)

# Convert expression matrix to dense numpy array (memory-heavy; only do if required)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

# Keep only first 2 dims for SCope visualization
scplus_obj.dr_cell["GEX_X_pca"] = scplus_obj.dr_cell["GEX_X_pca"].iloc[:, 0:2]
scplus_obj.dr_cell["GEX_X_umap"] = scplus_obj.dr_cell["GEX_X_umap"].iloc[:, 0:2]

# Ensure split_clusters exists in SCENICPLUS metadata
scplus_obj.metadata_cell["split_clusters"] = (
    scplus_obj.metadata_cell["GEX_cluster"].astype(str) + "_" + scplus_obj.metadata_cell["ACC_group"].astype(str)
)

# Save scplus object once
output_pkl_path = os.path.join(work_dir, "scplus_obj.pkl")
with open(output_pkl_path, "wb") as f:
    dill.dump(scplus_obj, f)


# ---------------------------
# Run SCENIC+
# ---------------------------

biomart_host = "nov2020.archive.ensembl.org"
tmp_dir = "/path/to/temp"

try:
    run_scenicplus(
        scplus_obj=scplus_obj, # SCENIC+ object
        variable=["split_clusters"], # Metadata variable for grouping 
        species="mmusculus", # Species
        assembly="mm10", # Genome assembly
        tf_file="/path/to/mouse_tf_gerg.gsc.riken.jp.txt", # TF list file
        save_path=os.path.join(work_dir, "scenicplus"), # Output directory
        biomart_host=biomart_host, # Ensembl host
        upstream=[1000, 150000], # Upstream regulatory window
        downstream=[1000, 150000], # Downstream regulatory window
        calculate_TF_eGRN_correlation=False, # Skip TF-expression vs eRegulon-activity correlation because RNA and ATAC are not measured in the same cells
        calculate_DEGs_DARs=True,  # Compute DEGs and DARs
        export_to_loom_file=True, # Export loom file
        export_to_UCSC_file=True,  # Export UCSC browser files
        path_bedToBigBed="/path/to/dir/with/binary", # Path to bedToBigBed binary
        n_cpu=12, # Number of CPUs
        _temp_dir=tmp_dir # Temporary directory
    )
except Exception as e:
    with open("error_log.txt", "w") as f:  # Write traceback to file if error occurs
        traceback.print_exc(file=f)
    print("An error occurred. Check 'error_log.txt' for details.")  # Print error message
    raise  # Re-raise exception
