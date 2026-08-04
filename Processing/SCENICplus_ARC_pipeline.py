# ARC SCENIC+ pipeline
# Assume you already built your cisTopic and motif enrichment objects

# Load libraries
import os
import pickle
import traceback
import dill
import numpy as np
import scanpy as sc

from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus


work_dir = "path/to/ARC_SCENIC/"  # Define working/output directory

# ---------------------------
# Binarize topics
# ---------------------------

cistopic_obj = pickle.load(
    open(
        os.path.join(work_dir, "ARC_updated_cistopic_obj.pkl"), # Load file
        "rb"
    )
)

region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')  # Binarize topics using Otsu thresholding
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop=3000)  # Binarize topics by selecting top 3000 regions per topic

imputed_acc_obj = impute_accessibility(  # Impute chromatin accessibility from topic model
    cistopic_obj,
    selected_cells=None,      # Use all cells
    selected_regions=None,    # Use all regions
    scale_factor=10**6        # Scale accessibility values to 1,000,000
)
normalized_imputed_acc_obj = normalize_scores(  # Normalize imputed accessibility scores
    imputed_acc_obj,
    scale_factor=10**4        # Scale normalized values to 10,000
)
variable_regions = find_highly_variable_features(  # Identify highly variable regions
    normalized_imputed_acc_obj,
    plot=False                # Do not generate diagnostic plot
)
markers_dict = find_diff_features(  # Identify differentially accessible regions per cluster
    cistopic_obj,
    imputed_acc_obj,
    variable='cluster',       # Group comparison variable
    var_features=variable_regions,  # Restrict analysis to variable regions
    split_pattern='-'         # Define pattern for splitting cluster names
)

if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):  # Check if output directory exists
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))  # Create directory if it does not exist

pickle.dump(  # Save Otsu-binarized topic regions to disk
    region_bin_topics_otsu,
    open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb')  # Open file in binary write mode
)
pickle.dump(  # Save top-3000 binarized topic regions to disk
    region_bin_topics_top3k,
    open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb')
)
pickle.dump(  # Save differential accessibility markers dictionary to disk
    markers_dict,
    open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb')
)


# ---------------------------
# Run SCENIC+
# ---------------------------

adata = sc.read_h5ad(  # Load AnnData object from disk
    os.path.join(work_dir, 'new_arc_20250218.h5ad')
)

menr = dill.load(  # Load motif enrichment results object
    open(os.path.join(work_dir, 'motif_enrichment_results/menr.pkl'), 'rb')
)

adata.obs_names = adata.obs_names + "-ARC"  # Append "-ARC" suffix to cell barcodes
adata.obs_names = adata.obs_names.str.replace("_1-ARC", "-ARC")  # Clean barcode formatting

test_data = adata.raw.to_adata()  # Extract raw expression layer as new AnnData

test_data.var.columns = ['genes']  # Rename variable columns to 'genes'

test_data.var.set_index('genes', inplace=True)  # Set gene names as index

scplus_obj = create_SCENICPLUS_object(  # Create SCENICPLUS object
    GEX_anndata=test_data,              # Provide gene expression AnnData
    cisTopic_obj=cistopic_obj,          # Provide cisTopic object
    menr=menr                           # Provide motif enrichment results
)

scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())  # Convert sparse expression matrix to dense NumPy array

# Keep only first two PCA dimensions for visualization compatibility
scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

# Keep only first two UMAP dimensions for visualization compatibility
scplus_obj.dr_cell['GEX_X_umap'] = scplus_obj.dr_cell['GEX_X_umap'].iloc[:, 0:2]

tmp_dir = '/path/to/temp'  # Define temporary directory for SCENIC+

biomart_host = "nov2020.archive.ensembl.org"  # Define Ensembl archive host for gene annotation

# Create combined grouping label from RNA cluster and ATAC group metadata
scplus_obj.metadata_cell['split_clusters'] = (
    scplus_obj.metadata_cell['GEX_cluster'].astype(str) + '_' +
    scplus_obj.metadata_cell['ACC_group']
)

try:
    run_scenicplus(  # Execute SCENIC+ pipeline
        scplus_obj=scplus_obj,  # SCENICPLUS object
        variable=['split_clusters'],  # Metadata variable for grouping 
        species='mmusculus',  # Species name
        assembly='mm10',  # Genome assembly
        tf_file='/path/to/mouse_tf_gerg.gsc.riken.jp.txt',  # TF list file
        save_path=os.path.join(work_dir, 'scenicplus'),  # Output directory
        biomart_host=biomart_host,  # Ensembl host
        upstream=[1000, 150000],  # Upstream regulatory window
        downstream=[1000, 150000],  # Downstream regulatory window
        calculate_TF_eGRN_correlation=False,  # Skip TF-expression vs eRegulon-activity correlation because RNA and ATAC are not measured in the same cells
        calculate_DEGs_DARs=True,  # Compute DEGs and DARs
        export_to_loom_file=True,  # Export loom file
        export_to_UCSC_file=True,  # Export UCSC browser files
        path_bedToBigBed='/path/to/dir/with/binary',  # Path to bedToBigBed binary
        n_cpu=4,  # Number of CPUs
        _temp_dir=tmp_dir  # Temporary directory
    )
except Exception as e:
    with open("error_log.txt", "w") as f:  # Write traceback to file if error occurs
        traceback.print_exc(file=f)
    print("An error occurred. Check 'error_log.txt' for details.")  # Print error message
    raise  # Re-raise exception
