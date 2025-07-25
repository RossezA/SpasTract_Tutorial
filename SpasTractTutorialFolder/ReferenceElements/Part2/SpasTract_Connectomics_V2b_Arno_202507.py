# Full exploratory pipeline with threading, connectome generation, and visualization
import os
import datetime
import subprocess
import tkinter as tk
from tkinter import Tk, filedialog, messagebox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
# from collections import defaultdict
import nibabel as nib
# import colorcet as cc
import cmcrameri.cm as cmc
import glob
from scipy.stats import ttest_ind
import fnmatch
from matplotlib.colors import Normalize
import math
import matplotlib.colors as mcolors  # For color blending
from sys import exit

def get_mrtrix_version():
    try:
        result = subprocess.run(['mrconvert', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = result.stdout.strip()
        # Example output: "== mrconvert 3.0.6 =="
        if 'mrconvert' in output:
            parts = output.split()
            for part in parts:
                if part.count('.') == 2:  # e.g., '3.0.6'
                    return part
        return None
    except FileNotFoundError:
        return None

def run_command(command, cwd=None):
    """Function to run command lines in terminal but allowing to not write the whole script each time with checkings"""

    print(f"\nRunning: {' '.join(command)}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(command)}")
    return result


def load_LUT(path):
    """
    Loads a LUT text file and returns:
    - `labels`: list of acronyms ordered by index (e.g. 'FRP_L')
    - `colors`: list of colors in hex format (e.g. '#268F45')
    - `fullnames`: list of full region names
    """
    labels = []
    colors = []
    fullnames = []

    lut = {}
    with open(path) as f:
        for line in f:
            parts = re.findall(r'"([^"]+)"|(\d+)', line)
            flat = [x for tup in parts for x in tup if x]
            if len(flat) < 4:
                continue
            idx = int(flat[0])
            acronym = flat[1]
            fullname = flat[2]
            color = flat[3]
            lut[idx] = (acronym, fullname, f"#{color.strip('#')}")

    # Ensure order by region index
    for idx in sorted(lut.keys()):
        acronym, fullname, hex_color = lut[idx]
        labels.append(acronym)
        colors.append(hex_color)
        fullnames.append(fullname)

    return labels, colors, fullnames


def plot_connectome(csv_file, output_dir, labels=None, colors=None):
    """Plots the connectome matrices based on MRtrix3 tck2connectome outputs and showcase global, inter- and intra-hemispheric structural connectivity.
       Also takes a look at GM and WM connectivity."""    
       
       
    matrix = pd.read_csv(csv_file, header=None).values
    n = matrix.shape[0]
    if n != matrix.shape[1]:
        print(f"Skipping {csv_file} (not square).")
        return
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    
    if labels is None :
        labels = ["No" for i in range(n//50)]
    if colors is None :
        colors = ["#999999"] * (n//50)
        
     
    def is_left(label):
        return label.endswith("_L")

    def is_right(label):
        return label.endswith("_R")
    
    def strip_suffix(label):
        return re.sub(r'_[LR]$', '', label)

    # Define CSF-related structures
    csf_acronyms = {
        "VL", "SEZ", "V3", "AQ", "V4", "V4r", "c"}
    
    # White matter structures that break lowercase convention
    wm_exceptions = {
        "moV", "sV", "sptV", "II", "gVIIn", "vVIIIn"
    }
    
    # def strip_suffix(label):
    #     """Remove suffix like '_L' or '_R' from label."""
    #     return label.split("_")[0]
    
    def is_csf(label):
        base = strip_suffix(label)
        return base in csf_acronyms
    
    def is_wm(label):
        base = strip_suffix(label)
        return (
            (base.islower() or base in wm_exceptions) and
            not is_csf(label)
        )
    
    def is_gm(label):
        base = strip_suffix(label)
        return (
            any(c.isupper() for c in base) and
            base not in wm_exceptions and
            not is_csf(label)
        )

    
    def extract_subset(labels, colors, condition_func):
        """
        Filters labels and corresponding colors based on a condition function.
        Returns:
            - filtered_labels: list of labels that match
            - filtered_indices: indices in the original list
            - filtered_colors: colors corresponding to filtered_labels
        """
        filtered = [
            (i, lbl, colors[i])
            for i, lbl in enumerate(labels)
            if i < len(colors) and condition_func(lbl)
        ]
        if not filtered:
            return [], [], []
        idxs, lbls, cols = zip(*filtered)
        return list(lbls), list(idxs), list(cols)
    
    # Creating index groups for L, R, GM and WM based on suffix and upper/lowercase characters in acronyms (a.k.a labels here)
    L_labels, L_idx, L_colors = extract_subset(labels, colors, is_left)
    R_labels, R_idx, R_colors = extract_subset(labels, colors, is_right)
    GM_labels, GM_idx, GM_colors = extract_subset(labels, colors, is_gm)
    WM_labels, WM_idx, WM_colors = extract_subset(labels, colors, is_wm)
    
    GM_L_labels, GM_L_idx, GM_L_colors = extract_subset(GM_labels, GM_colors, is_left)
    GM_R_labels, GM_R_idx, GM_R_colors = extract_subset(GM_labels, GM_colors, is_right)
    
    
    def valid_tick(idx):
     """Return True if the label is NOT a 'No' label with grey color."""
     return not (labels[idx] == "No" and colors[idx].lower() == "#999999")
    
 
    def label_mapper(ax, matrixsize, labelsleft, labelsright, colorsleft, colorsright):
        """
        Downsampled label visualization for margins,
        like set_colored_ticks but without blur of superimposed labels.
        Skips labels with color #999999 and text 'No'.
        """
        # Filter labels and colors based on condition
        filteredleft = [
            (i, label, color)
            for i, (label, color) in enumerate(zip(labelsleft, colorsleft))
            if not (label == "No" and color.lower() == "#999999")
        ]
        filteredright = [
            (i, label, color)
            for i, (label, color) in enumerate(zip(labelsright, colorsright))
            if not (label == "No" and color.lower() == "#999999")
        ]
    
        if not filteredleft and not filteredright:
            return  # Nothing to draw
    
        # Unpack filtered values
        positionsleft, filtered_labelsleft, filtered_colorsleft = zip(*filteredleft)
        positionsleft = np.array(positionsleft)
        positionsright, filtered_labelsright, filtered_colorsright = zip(*filteredright)
        positionsright = np.array(positionsright)

        
        # Set ticks and labels
        ax.set_xticks(positionsright)
        ax.set_xticklabels(filtered_labelsright, rotation=90)
        for tick, color in zip(ax.get_xticklabels(), filtered_colorsright):
            tick.set_color(color)
    
        ax.set_yticks(positionsleft)
        ax.set_yticklabels(filtered_labelsleft, rotation=0)
        for tick, color in zip(ax.get_yticklabels(), filtered_colorsleft):
            tick.set_color(color)
        
    def extract_submatrix(matrix, rows, cols):
        return matrix[np.ix_(rows, cols)]
    
        
    #Manual handling of colormaps; diverging or sequential colormap for difference and p-values : avoid diverging colormap with 0 white (too bright), have sequential colormap with brightest at 0
    cmaps_list = ["inferno"] #default for classical matrices
    is_pval_matrix = "p_values" in base_name.lower()
    is_diff_matrix = "minus" in base_name.lower()
    if is_pval_matrix:
        # # Threshold p-values: set any value > 0.05 to 0, keep the rest
        # matrix = np.where(matrix > 0.05, 0, matrix)
        cmaps_list = ['turbo'] #Wistia too bright
    elif is_diff_matrix :
        # Difference matrices use managua colormap
        cmaps_list = [cmc.managua] # otherwise "berlin" or "vanimo" but this one is used for symmetry analysis... 
        
        
    # Handle % matrices: Remove ±100% values that come from 0-division
    if "percent" in base_name.lower():
        matrix[np.abs(matrix) >= 100] = 0  # remove noisy % effects
    # # Also Fill holes: convert all NaN or inf/-inf to 0
    # matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    
    #Manual handling of truncation factor, empirically 10 is best
    truncation_factor = 10
    
    
    #Defining boundaries of colorscales
    abs_max = np.nanmax(np.abs(matrix)) #Better than np.max and -np.max because : 
                                        #Ensures symmetric scaling around 0;
                                        #Captures the most extreme absolute value (whether it's negative or positive); 
                                        #Ideal for visualizing difference and percentage connectomes using diverging colormaps like seismic or coolwarm.
    if is_pval_matrix:
        abs_min = 0
        center_value = None
        value_trunc_min = 0
        value_trunc_max = abs_max / 5 #To go from p<0.05 to p<0.01
    elif is_diff_matrix :
        abs_min =-abs_max
        center_value = 0
        value_trunc_max = abs_max / truncation_factor
        value_trunc_min = -value_trunc_max 
    else : 
        abs_min = 0 
        center_value = None
        value_trunc_max = abs_max / truncation_factor
        value_trunc_min = 0
        
        
    for cmaps in cmaps_list :
    
        # Full Matrix
        plt.figure(figsize=(14, 12))
        ax = sns.heatmap(matrix, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
        label_mapper(ax, n , labels, labels, colors, colors)
        plt.title(f"{base_name}: Full Matrix")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{base_name}_full.png"),dpi=300, bbox_inches='tight')
        plt.close()
        
        if "percent" not in base_name.lower():
            # Full Matrix truncated in value by truncation_factor
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(matrix, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
            label_mapper(ax, n , labels, labels, colors, colors)
            plt.title(f"{base_name}: Full Matrix truncated in value by a factor {truncation_factor}")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_full_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
            plt.close()
        
        
        if len(L_idx)>0 :
            # Intra-hemispheric L
            intra_L = extract_submatrix(matrix, L_idx, L_idx)
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(intra_L, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
            label_mapper(ax, len(L_idx) , L_labels, L_labels, L_colors, L_colors)
            plt.title(f"{base_name}: Intra-Hemisphere (Left)")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_intra_L.png"),dpi=300, bbox_inches='tight')
            plt.close()
            
            if "percent" not in base_name.lower():
                # Intra-hemispheric L truncated in value by truncation_factor
                plt.figure(figsize=(14, 12))
                ax = sns.heatmap(intra_L, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
                label_mapper(ax, len(L_idx) , L_labels, L_labels, L_colors, L_colors)
                plt.title(f"{base_name}: Intra-Hemisphere (Left) truncated in value by a factor {truncation_factor}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_intra_L_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
                plt.close()
            
        if len(R_idx)>0 :
            # Intra-hemispheric R
            intra_R = extract_submatrix(matrix, R_idx, R_idx)
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(intra_R, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
            label_mapper(ax, len(R_idx) , R_labels, R_labels, R_colors, R_colors)
            plt.title(f"{base_name}: Intra-Hemisphere (Right)")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_intra_R.png"),dpi=300, bbox_inches='tight')
            plt.close()
            
            if "percent" not in base_name.lower():
                # Intra-hemispheric R truncated in value by truncation_factor
                plt.figure(figsize=(14, 12))
                ax = sns.heatmap(intra_R, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
                label_mapper(ax, len(R_idx) , R_labels, R_labels, R_colors, R_colors)
                plt.title(f"{base_name}: Intra-Hemisphere (Right) truncated in value by a factor {truncation_factor}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_intra_R_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
                plt.close()
            
        if len(L_idx + R_idx)>0 :
            # Inter-hemispheric L vs R
            inter_LR = extract_submatrix(matrix, L_idx, R_idx) 
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(inter_LR, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
            
            label_mapper(ax, max(len(L_idx),len(R_idx)) , L_labels, R_labels, L_colors, R_colors)

            plt.title(f"{base_name}: Inter-Hemisphere (L vs R)")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_inter_LR.png"),dpi=300, bbox_inches='tight')
            plt.close()
            
            if "percent" not in base_name.lower():
            # Inter-hemispheric L vs R truncated in value by truncation_factor
                plt.figure(figsize=(14, 12))
                ax = sns.heatmap(inter_LR, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
                
                label_mapper(ax, max(len(L_idx),len(R_idx)) , L_labels, R_labels, L_colors, R_colors)
                
                plt.title(f"{base_name}: Inter-Hemisphere (L vs R) truncated in value by a factor {truncation_factor}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_inter_LR_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
                plt.close()
            
            
        if len(GM_idx)>0 :
            # Grey Matter Matrix
            GM_matrix = extract_submatrix(matrix, GM_idx, GM_idx)
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(GM_matrix, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
            label_mapper(ax, len(GM_idx) , GM_labels, GM_labels, GM_colors, GM_colors)
            plt.title(f"{base_name}: Grey Matter")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_GM.png"),dpi=300, bbox_inches='tight')
            plt.close()
            
            
            if "percent" not in base_name.lower():
                # GM Matrix truncated in value by truncation_factor
                plt.figure(figsize=(14, 12))
                ax = sns.heatmap(GM_matrix, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
                label_mapper(ax, len(GM_idx) , GM_labels, GM_labels, GM_colors, GM_colors)
                plt.title(f"{base_name}: Grey Matter truncated in value by a factor {truncation_factor}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_GM_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
                plt.close()
        
        if len(WM_idx)>0 :
            # White Matter Matrix
            WM_matrix = extract_submatrix(matrix, WM_idx, WM_idx)
            
            #for WhiteMatter having lower connectivity ?
            # abs_max_WM = np.nanmax(np.abs(WM_matrix))
            # value_min_WM = -abs_max_WM / truncation_factor
            # value_max_WM = abs_max_WM / truncation_factor
            
            plt.figure(figsize=(14, 12))
            ax = sns.heatmap(WM_matrix, cmap=cmaps, center=center_value, vmin=abs_min, vmax=abs_max)
            label_mapper(ax, len(WM_idx) , WM_labels, WM_labels, WM_colors, WM_colors)
            plt.title(f"{base_name}: White Matter")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base_name}_WM.png"),dpi=300, bbox_inches='tight')
            plt.close()
            
            if "percent" not in base_name.lower():
                # WM Matrix truncated in value by truncation_factor
                plt.figure(figsize=(14, 12))
                ax = sns.heatmap(WM_matrix, cmap=cmaps, center=center_value, vmin=value_trunc_min, vmax=value_trunc_max)
                label_mapper(ax, len(WM_idx) , WM_labels, WM_labels, WM_colors, WM_colors)
                plt.title(f"{base_name}: White Matter truncated in value by a factor {truncation_factor}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f"{base_name}_WM_trunc{truncation_factor}.png"),dpi=300, bbox_inches='tight')
                plt.close()

        print(f"Plots saved for {base_name}")
        

#Function crop_to_mask and extract_and_crop_voxels_by_values for VOXELS EXTRACTION
def crop_to_mask(data, affine):
    coords = np.array(np.nonzero(data))
    if coords.size == 0:
        return None, None, None  # No nonzero voxels
    min_coords = coords.min(axis=1)
    max_coords = coords.max(axis=1) + 1  # +1 for inclusive slice
    cropped_data = data[
        min_coords[0]:max_coords[0],
        min_coords[1]:max_coords[1],
        min_coords[2]:max_coords[2]
    ]

    new_origin_voxel = min_coords
    new_affine = affine.copy()
    new_affine[:3, 3] = (
        affine[:3, :3] @ new_origin_voxel + affine[:3, 3]
    )
    return cropped_data, new_affine, min_coords

def extract_and_crop_voxels_by_values(input_file, labelslist, target_values, output_dir):
    # Load the NIfTI file
    img = nib.load(input_file)
    data = img.get_fdata()
    affine = img.affine
    header = img.header
    
    # Clean filename base
    base_name = os.path.basename(input_file)
    root_name, ext = os.path.splitext(base_name)
    if ext == '.gz': # Handle .nii.gz
        root_name, _ = os.path.splitext(root_name)
        ext = '.nii.gz'
    else:
        ext = '.nii'
    
    # Process each target value
    for value in target_values:
        print(f"Extracting and cropping value: {value}; matching region : {labelslist[value-1]}")
        mask = (data == value)
        masked_data = np.zeros_like(data, dtype=np.uint16)
        masked_data[mask] = value

        cropped_data, new_affine, origin = crop_to_mask(masked_data, affine)

        if cropped_data is None:
            print(f"No voxels found for region {labelslist[value-1]} of value {value}. Skipping.")
            continue

        new_img = nib.Nifti1Image(cropped_data, affine=new_affine, header=header)
        output_file = f"{root_name}_{labelslist[value-1]}_{value}_cropped{ext}"
        output_filepath = os.path.join(output_dir,output_file)
        nib.save(new_img, output_filepath)
        print(f"Saved cropped file: {output_file}")
##End functions voxel extraction       
 
#Functions for Avg, Diff, Pval WT & Spastin-KO conditions
def select_folder(title="Select folder", initialDir=None):
    root = Tk()
    root.withdraw()
    return filedialog.askdirectory(title=title, initialdir=initialDir)

def find_connectome_files(base_folder, connectome_type):
    pattern = os.path.join(base_folder, "*", "ConnectomeProject_*", connectome_type)
    return glob.glob(pattern, recursive=True)

def load_matrices(csv_paths):
    matrices = []
    for path in csv_paths:
        try:
            mat = pd.read_csv(path, header=None).values
            if mat.shape[0] == mat.shape[1]:
                mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
                matrices.append(mat)
        except Exception as e:
            print(f"Error loading {path}: {e}")
    return np.array(matrices)

def save_matrix(matrix, path, fmt="%.4f"):
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    pd.DataFrame(matrix).to_csv(path, header=False, index=False, float_format=fmt, na_rep="0")

def process_connectome_type(connectome_type, wt_folder, spastin_folder, output_folder):
    print(f"\nProcessing: {connectome_type}")

    wt_paths = find_connectome_files(wt_folder, connectome_type)
    spastin_paths = find_connectome_files(spastin_folder, connectome_type)
    
    print(f"search pattern: {'*/*ConnectomeProject*/' + connectome_type}")
    print(f"Found WT files: {wt_paths}")
    print(f"Found SPASTIN files: {spastin_paths}")

    wt_matrices = load_matrices(wt_paths)
    spastin_matrices = load_matrices(spastin_paths)  
    
    if len(wt_matrices) == 0 or len(spastin_matrices) == 0:
        print(f"Skipping {connectome_type}: missing data in one or both groups.")
        return

    assert wt_matrices.shape[1:] == spastin_matrices.shape[1:], "Matrix dimensions must match."

    avg_wt = np.mean(wt_matrices, axis=0)
    avg_spastin = np.mean(spastin_matrices, axis=0)

    base_name = connectome_type.replace(".csv", "")

    # Save average matrices
    os.makedirs(os.path.join(output_folder,"AvgWT"), exist_ok=True)
    save_matrix(avg_wt, os.path.join(output_folder,"AvgWT", f"{base_name}_WT_average.csv"))
    os.makedirs(os.path.join(output_folder, "AvgSpastin-KO"), exist_ok=True)
    save_matrix(avg_spastin, os.path.join(output_folder, "AvgSpastin-KO", f"{base_name}_Spastin-KO_average.csv"))

    # Difference matrix
    diff = avg_spastin - avg_wt

    # Safer percent difference
    safe_denominator = np.where((avg_wt + avg_spastin) == 0, np.nan, avg_wt + avg_spastin)
    diff_percent = 100 * diff / safe_denominator
    diff_percent = np.nan_to_num(diff_percent, nan=0.0, posinf=0.0, neginf=0.0)
    
    os.makedirs(os.path.join(output_folder, "Diff_Spastin-KOvsWT"), exist_ok=True)
    save_matrix(diff, os.path.join(output_folder,"Diff_Spastin-KOvsWT", f"{base_name}_Spastin-KOvsWT.csv"))
    save_matrix(diff_percent, os.path.join(output_folder, "Diff_Spastin-KOvsWT", f"{base_name}_Spastin-KOvsWT_percent.csv"))

    # T-test
    t_stat, p_val = ttest_ind(spastin_matrices, wt_matrices, axis=0, equal_var=False, nan_policy="omit")
    p_val = np.nan_to_num(p_val, nan=1.0)  # fallback: not significant

    N_WT = len(wt_matrices)
    N_SPASTIN = len(spastin_matrices)
    p_threshold = 0.05
    significant_mask = p_val < p_threshold

    filtered_diff = np.where(significant_mask, diff, 0)
    filtered_diff_percent = np.where(significant_mask, diff_percent, 0)

    # Apply mask to p-values: non-significant entries set to threshold
    masked_p_val = np.where(significant_mask, p_val, p_threshold)

    # Save matrices
    os.makedirs(os.path.join(output_folder, "Pval_Spastin-KOvsWT"), exist_ok=True)
    save_matrix(masked_p_val, os.path.join(output_folder,"Pval_Spastin-KOvsWT", f"{base_name}_p_values.csv"))
    save_matrix(filtered_diff, os.path.join(output_folder, "Pval_Spastin-KOvsWT", f"{base_name}_p<{p_threshold}_diff_{N_SPASTIN}vs{N_WT}.csv"))
    save_matrix(filtered_diff_percent, os.path.join(output_folder, "Pval_Spastin-KOvsWT", f"{base_name}_p<{p_threshold}_diff_{N_SPASTIN}vs{N_WT}_percent.csv"))

    print(f"Done: {connectome_type}")        
##End functions for Avg, Diff, Pval WT & Spastin-KO conditions  
 
#Functions for symmetry analysis
def compute_symmetry_matrix(connectome):
    N = connectome.shape[0]
    half = N // 2
    RR = connectome[:half, :half]
    LL = connectome[half:, half:]
    # LL_mirrored = np.flipud(np.fliplr(LL))
    return RR - LL

def extract_clean_R_labels(lut_labels):
    return [lbl.replace('_R', '') for lbl in lut_labels if lbl.endswith('_R')]

def strip_suffix(label):
    return label.rsplit("_", 1)[0]

# White matter structures that break lowercase convention
wm_exceptions = {
    "moV", "sV", "sptV", "II", "gVIIn", "vVIIIn"
}
csf_acronyms = {
    "VL", "SEZ", "V3", "AQ", "V4", "V4r", "c"
}

def is_csf(label):
    return strip_suffix(label) in csf_acronyms

def is_wm(label):
    base = strip_suffix(label)
    return (base.islower() or base in wm_exceptions) and not is_csf(label)

def is_gm(label):
    base = strip_suffix(label)
    return any(c.isupper() for c in base) and base not in wm_exceptions and not is_csf(label)

def get_file_map(folder, keyword):
    mapping = {}
    for f in os.listdir(folder):
        if f.endswith(".csv"):
            base = re.sub(rf'_{keyword}_.*\.csv$', '', f, flags=re.IGNORECASE)
            mapping[base] = os.path.join(folder, f)
    return mapping

def plot_symmetry_pair(sym1, sym2, labels, lut_colors, title, filename, vmax=None):
    vmax = vmax or max(np.abs(sym1).max(), np.abs(sym2).max())
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    cmap = cmc.vanimo

    for ax, mat, t in zip(axs, [sym1, sym2], ['Spastin-KO (Right - Left)', 'WT (Right - Left)']):
        im = ax.imshow(mat, cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_title(t)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticklabels(labels, fontsize=6)
        for lbl, color in zip(ax.get_xticklabels(), lut_colors):
            lbl.set_color(color)
        for lbl, color in zip(ax.get_yticklabels(), lut_colors):
            lbl.set_color(color)

    fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.75, label="Difference (R - L)")
    fig.suptitle(title, fontsize=13)
    plt.savefig(filename, dpi=300)
    plt.close()
##End functions symmetry analysis  

#Functions for connectograms
def plot_connectogram_circle(matrix, labels, colors, output_path, threshold=0.001, top_n=20):
    """
    Plots a circular connectogram:
    - Only top N connections are labeled with larger font and higher visibility.
    - Node labels are hidden for clarity.
    """
    assert matrix.shape[0] == len(labels) == len(colors), "Matrix/Labels/Colors size mismatch"

    n = len(labels)
    angle_step = 2 * math.pi / n
    angles = [i * angle_step for i in range(n)]
    positions = [(math.cos(a), math.sin(a)) for a in angles]

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw nodes (dots only)
    for i, (x, y) in enumerate(positions):
        ax.plot(x, y, 'o', markersize=12, color=colors[i], markeredgecolor='k')

    # Normalize for consistent linewidth/alpha scaling
    norm = Normalize(vmin=0, vmax=np.max(matrix))

    # Extract all valid connections above threshold
    connections = [
        (i, j, matrix[i, j])
        for i in range(n) for j in range(i + 1, n)
        if matrix[i, j] > threshold
    ]

    # Get top N connections by strength
    top_connections = sorted(connections, key=lambda x: x[2], reverse=True)[:top_n]
    top_pairs = set((i, j) for i, j, _ in top_connections)

    # Draw background (non-top) connections
    for i, j, val in connections:
        if (i, j) in top_pairs:
            continue
        x1, y1 = positions[i]
        x2, y2 = positions[j]
        c1 = np.array(mcolors.to_rgb(colors[i]))
        c2 = np.array(mcolors.to_rgb(colors[j]))
        blended = (c1 + c2) / 2
        lw = 0.5 + 2 * norm(val)
        alpha = 0.2 + 0.4 * norm(val)
        ax.plot([x1, x2], [y1, y2], color=blended, lw=lw, alpha=alpha)

    # Draw top N connections in foreground
    for i, j, val in top_connections:
        x1, y1 = positions[i]
        x2, y2 = positions[j]
        c1 = np.array(mcolors.to_rgb(colors[i]))
        c2 = np.array(mcolors.to_rgb(colors[j]))
        blended = (c1 + c2) / 2
        lw = 1 + 5 * norm(val)
        ax.plot([x1, x2], [y1, y2], color=blended, lw=lw, alpha=1.0)

        # Label at midpoint
        xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(
            xm, ym, f"{labels[i]}-{labels[j]}",
            fontsize=9, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, lw=0)
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved connectogram to: {output_path}")

##End function for connectograms 

########################
#MAIN SCRIPT STARTS HERE
########################

#0-INITIALIZATION OF THE SCRIPT AND CHOOSING MODE----------------------------------------------------------------------------------
print("ConnectomeProject Pipeline (Exploratory + Visuals Modes)")
num_threads = os.cpu_count()
print(f"Using {num_threads} CPU threads")
os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(num_threads)

root = tk.Tk()
root.withdraw()


print("You just launched the connectomics script, please pick the mode you want:\n")
print("1) Generate connectomes, plot connectivity matrices and extract specific tracts (streamlines & voxels) for one sample")
print("2) Average WT and Spastin-KO connectomes and compute differences and significant differences connectomes")
print("3) Plot only connectivity matrices based on connectomes in a folder (to plot average matrices, difference matrices, p-values matrices...)")
print("4) Plot analysis of symmetry based on WT and Spastin-KO average connectomes in respective folders")
print("5) Plot connectograms based on connectomes in a folder")


while True:
    choice = input("\nEnter your choice (1, 2, 3, 4, 5) or abort script (X): ")
    if choice in {"1", "2", "3", "4", "5"}:
        print(f"\nYou selected mode {choice}.")
        if choice in {"1"}:
            print("\nGeneration of connectomes for one sample starting...")
        elif choice in {"2"}:
            print("\nGeneration of WT Spastin-KO average connectomes and their differences starting...")
        elif choice in {"3"}:
            print("\nPlotting structural connectivity matrices based on connectomes in future chosen folder...")
        elif choice in {"4"}:
            print("\nPlotting analysis of symmetry based on WT and Spastin-KO average connectomes...")
        elif choice in {"5"}:
            print("\nPlotting connectograms based on connectomes in future chosen folder...")
        break
    elif choice in {"X"}:
        print(f"\nYou choose to abort the script ({choice}), process will be terminated...")
        exit()
    else:
        print("Invalid input. Please enter 1, 2, 3, 4, 5.")
#0-END----------------------------------------------------------------------------------------------------------------------------------------------------
####MODE1####
#1-1READING FILES FOR ANTs ATLAS COREGISTRATION--------------------------------------------------------------------------------------------------------------------------

if choice in {"1"} :
    print("Selecting Atlas Template NIfTI file...")
    atlas_path = filedialog.askopenfilename(title="Select Atlas reference template image (e.g. P56_Atlas_Template_REOR.nii)",filetypes=[("Atlas reference", "*.nii")]) #.nii.gz not compatible with MacOS + Nibabel
    print("Selecting Atlas remapped Parcellation NIfTI file wanted for connectomics (can be grouped, lateralized or grouped&lateralized)...")
    atlas_parcellation_path = filedialog.askopenfilename(title="Select Atlas remapped Parcellation/Annotation image (e.g. P56_Atlas_Parcellation_REOR_remapped*.nii, can be grouped,lateralized or both)",filetypes=[("Atlas Parcellation/Annotation", "*.nii")]) #.nii.gz
    print("Selecting Atlas LUT text file matching the chosen Parcellation (continuous grouped/lateralized/grouped&lateralized)...")
    labels_txt_path = filedialog.askopenfilename(title="Select the Atlas LUT (e.g. AllenBrain_Metadata_labelsshortfullhexacolor_indexcontinuous*.txt file, continuous grouped/lateralized/grouped&lateralized)",filetypes=[("Atlas LUT", ".txt")])
    
    print("Selecting preprocessed DWI file...")
    dwi_path = filedialog.askopenfilename(title="Select preprocessed DWI file (e.g. 'Sample'_dwi_den_unr_N4cat.mif)",filetypes=[("DWI files", "*.mif")])
    print("Selecting tractogram .tck file matching the chosen sample...")
    tck_path = filedialog.askopenfilename(title="Select tractogram .tck file matching the sample chosen (e.g. tracts_*.tck)",filetypes=[("Tractogram files", "*.tck")])
    
    base_output_dir = os.path.dirname(tck_path)
    print(f"By default the following directory has been chosen to host the new ConnectomeProject : {base_output_dir}")
    
    #Here we load the LookUpTable or LUT to get the regions' acronyms, fullnames and colors based on their ID...
    labels, colors, _ = load_LUT(labels_txt_path)
    
    project_dir = os.path.join(base_output_dir, "ConnectomeProject_")
    now = datetime.datetime.now()
    project_dir_date = project_dir + (now.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(project_dir_date, exist_ok=True)
    os.chdir(project_dir_date)
    
    b0_mean_name = os.path.splitext(os.path.split(dwi_path)[1])[0] + "_b0_mean.nii" #.nii.gz
    b0_mean_path = os.path.join(project_dir_date, b0_mean_name)
    
    # First process: dwiextract
    extract = subprocess.Popen(["dwiextract", dwi_path, "-", "-bzero", "-nthreads", str(num_threads)],stdout=subprocess.PIPE)
    # Second process: mrmath
    avg = subprocess.Popen(["mrmath", "-", "mean", b0_mean_path, "-axis", "3", "-nthreads", str(num_threads)],stdin=extract.stdout)
    # Ensure extract doesn't hang
    extract.stdout.close()
    # Wait for the processes to finish
    avg.communicate()
    
    
    use_cs = messagebox.askyesno("Use Compressed Sensing image?", "Do you want to use a CS high-resolution anatomical image as registration reference?")
    
    
    if use_cs:
            print("Selecting Compressed Sensing of the chosen sample...")
            cs_path = filedialog.askopenfilename(initialdir=base_output_dir, title="Select Compressed Sensing of the chosen sample (e.g. 3DCS_FLASH_50um_'Sample'_REOR_aligned2dwi.nii)",filetypes=[("CS Images", "*.nii")]) #.nii.gz
            triggerquestion=False
            if triggerquestion :
                correct_cs = messagebox.askyesno("CS image warping to DWI space (resolution altering) ?", "Does the CS image need to be warped to the DWI space (keep in mind that this will alter the resolution back to DWI one) ?")
                
                if correct_cs : 
                    # Register CS to DWI (rigid only)
                    run_command(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", b0_mean_path, "-m", cs_path, "-t", "a", "-o", "cs2dwi_"])
                    
                    # Register Atlas to CS (nonlinear)
                    run_command(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", cs_path, "-m", atlas_path, "-o", "atlas2cs_"])
    
                    # Apply chained transform to atlas labels
                    warped_labels_path = os.path.join(project_dir_date, "atlas_labels_in_dwi.nii") #.nii.gz
                    run_command([
                        "antsApplyTransforms", "-d", "3",
                        "-i", atlas_parcellation_path,
                        "-r", b0_mean_path,
                        "-t", "cs2dwi_0GenericAffine.mat",
                        "-t", "atlas2cs_1Warp.nii.gz", #.nii.gz
                        "-t", "atlas2cs_0GenericAffine.mat",
                        "-n", "NearestNeighbor",
                        "-o", warped_labels_path
                    ])
                else :
                    # Register Atlas to CS (nonlinear)
                    run_command(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", cs_path, "-m", atlas_path, "-o", "atlas2cs_"])
    
                    # Apply Atlas2CS transform to atlas labels
                    warped_labels_path = os.path.join(project_dir_date, "atlas_labels_in_CS.nii") #.nii.gz
                    run_command([
                        "antsApplyTransforms", "-d", "3",
                        "-i", atlas_parcellation_path,
                        "-r", cs_path,
                        "-t", "atlas2cs_1Warp.nii.gz", #.nii.gz
                        "-t", "atlas2cs_0GenericAffine.mat",
                        "-n", "NearestNeighbor",
                        "-o", warped_labels_path
                    ])
            
            if not(triggerquestion) :
                # Register Atlas to CS (nonlinear)
                run_command(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", cs_path, "-m", atlas_path, "-o", "atlas2cs_"])
    
                # Apply Atlas2CS transform to atlas labels
                warped_labels_path = os.path.join(project_dir_date, "atlas_labels_in_CS.nii") #.nii.gz
                run_command([
                    "antsApplyTransforms", "-d", "3",
                    "-i", atlas_parcellation_path,
                    "-r", cs_path,
                    "-t", "atlas2cs_1Warp.nii.gz", #.nii.gz
                    "-t", "atlas2cs_0GenericAffine.mat",
                    "-n", "NearestNeighbor",
                    "-o", warped_labels_path
                ])
            
    else:
        # No CS, register Atlas directly to DWI
        run_command(["antsRegistrationSyNQuick.sh", "-d", "3", "-f", b0_mean_path, "-m", atlas_path, "-o", "atlas2dwi_"])
        warped_labels_path = os.path.join(project_dir_date, "atlas_labels_in_dwi.nii") #.nii.gz
        run_command([
            "antsApplyTransforms", "-d", "3",
            "-i", atlas_parcellation_path,
            "-r", b0_mean_path,
            "-t", "atlas2dwi_1Warp.nii.gz", #.nii.gz
            "-t", "atlas2dwi_0GenericAffine.mat",
            "-n", "NearestNeighbor",
            "-o", warped_labels_path
        ])
       
#1-1-END----------------------------------------------------------------------------------------------------------------------------------------------------  
#2-1-MRTRIX3 TCK2CONNECTOME USE TO GENERATE INITIAL CONNECTOMES-----------------------------------------------------------------------------------------------------------------------------
    connectome_outputs = []
    
    out_csv_default = os.path.join(project_dir_date, "connectome_default.csv")
    assignmentsname = "connectome_default_assignments.txt"
    run_command(["tck2connectome", "–symmetric","–zero_diagonal", tck_path, warped_labels_path, out_csv_default, "-out_assignments", assignmentsname, "-nthreads", str(num_threads)])
    connectome_outputs.append(out_csv_default)
    
    out_csv_invnodevol = os.path.join(project_dir_date, "connectome_invnodevol.csv")
    assignmentsname_invnodevol = "connectome_invnodevol_assignments.txt"
    run_command(["tck2connectome","–symmetric","–zero_diagonal","-scale_invnodevol", tck_path, warped_labels_path, out_csv_invnodevol, "-out_assignments", assignmentsname_invnodevol, "-nthreads", str(num_threads)])
    connectome_outputs.append(out_csv_invnodevol)
    
    # #Defining radius for radial search in mm:
    # radius = 6
    # out_csv_invnodevol_RADIAL = os.path.join(project_dir_date, f"connectome_invnodevol_RADIAL{radius}.csv")
    # assignmentsname_invnodevol_RADIAL = f"connectome_invnodevol_RADIAL{radius}_assignments.txt"
    # run_command(["tck2connectome","–symmetric","–zero_diagonal","-scale_invnodevol", tck_path, warped_labels_path, out_csv_invnodevol_RADIAL, "-out_assignments", assignmentsname_invnodevol_RADIAL, "-assignment_radial_search", str(radius), "-nthreads", str(num_threads)])
    # connectome_outputs.append(out_csv_invnodevol_RADIAL)
    
    out_csv_length = os.path.join(project_dir_date, "connectome_scalelength.csv")
    assignmentsname_length = "connectome_scalelength_assignments.txt"
    run_command(["tck2connectome","-symmetric","–zero_diagonal","-scale_length", tck_path, warped_labels_path, out_csv_length, "-out_assignments", assignmentsname_length, "-nthreads", str(num_threads)])
    connectome_outputs.append(out_csv_length)
    
    out_csv_invlength = os.path.join(project_dir_date, "connectome_scaleinvlength.csv")
    assignmentsname_invlength = "connectome_scaleinvlength_assignments.txt"
    run_command(["tck2connectome","-symmetric","–zero_diagonal","-scale_invlength", tck_path, warped_labels_path, out_csv_invlength, "-out_assignments", assignmentsname_invlength, "-nthreads", str(num_threads)])
    connectome_outputs.append(out_csv_invlength)
    
#2-1-END-------------------------------------------------------------------------------------------------------------------------------------------------------------
#3-1-INITIAL CONNECTOMES PLOTTING------------------------------------------------------------------------------------------------------------------------------------
    for csv in connectome_outputs:
        plot_connectome(csv, project_dir_date, labels=labels, colors=colors)
#3-1-END------------------------------------------------------------------------------------------------------------------------------------
#4-1-EXTRACTING TCK AND VOXELS FROM CONNECTOMES, AND STATS(TODO)----------------------------------------------------------------------------------------------------------------------------
    for csv in connectome_outputs:
        
        [csvdirpath, csvname] = os.path.split(csv)
        os.chdir(csvdirpath)
        connectome2tckdir = os.path.join(csvdirpath,"Connectome2Tck")
        os.makedirs(connectome2tckdir, exist_ok=True)
        
        csvassignmentsname = os.path.splitext(csv)[0] + "_assignments.txt"
    
        # csvedgesname = csvname + "_edges.tck"
        # csvstatsname = csvname + "_stats.txt"
        # edge_tracks_path = os.path.join(project_dir_date, csvedgesname)
        # stats_output_path = os.path.join(project_dir_date, csvstatsname)
        
        os.chdir(connectome2tckdir)
        if csvname in {"connectome_default.csv"} :
        
        #WARNING : in MRtrix3.0.6 it seems to require .tck at the end of the prefix out for connectome2tck however in MRtrix3.0.4 and MRtrix3.0.5 it wasn't the case...    
            version = get_mrtrix_version()
            if version == "3.0.6":
                print("MRtrix 3.0.6 is installed, proceeding with .tck extension in prefix out for connectome2tck.")
                
                #Example here (in Left-Right Parcellation) : STR striatum - MB midbrain connections (26,150-38,162)
                STRMBname = "STR-MB.tck" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "26,150,38,162",
                             tck_path, csvassignmentsname, STRMBname, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here : Anterior Commissures aco 46,170
                aconame = "aco.tck" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "46,170",
                             tck_path, csvassignmentsname, aconame, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here : Corpus Callosum cc with : 
                #76,200 scwm supra-callosal_cerebral_white_matter; 
                #77,201 fa corpus_callosum_anterior_forceps;
                #78,202 ec corpus_callosum_external_capsule;
                #79,203 ee corpus_callosum_extreme_capsule;
                #80,204 ccg corpus_callosum_genu;
                #81,205 fp corpus_callosum_posterior_forceps;
                #82,206 ccb corpus_callosum_body;
                #83,207 ccs corpus_callosum_splenium;
                ccname = "cc.tck" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "76,77,78,79,80,81,82,83,200,201,202,203,204,205,206,207",
                             tck_path, csvassignmentsname, ccname, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here :  Cortico-spinal tract cst with Pyramidal tract py, Pyramidal decussation pyd, Pons P, Medulla MY (84,208-87,211-88,212-39,163-40,164)
                cstname = "Cortico-spinal.tck" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "84,208,87,211,88,212,39,163,40,164",
                             tck_path, csvassignmentsname, cstname, "-files", "single", "-nthreads", str(num_threads)])
                
            else:
                print(f"MRtrix version is {version} — not 3.0.6, no need for .tck in prefix out for connectome2tck")
                
                #Example here (in Left-Right Parcellation) : STR striatum - MB midbrain connections (26,150-38,162)
                STRMBname = "STR-MB" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "26,150,38,162",
                             tck_path, csvassignmentsname, STRMBname, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here : Anterior Commissures aco 46,170
                aconame = "aco" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "46,170",
                             tck_path, csvassignmentsname, aconame, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here : Corpus Callosum cc with : 
                #76,200 scwm supra-callosal_cerebral_white_matter; 
                #77,201 fa corpus_callosum_anterior_forceps;
                #78,202 ec corpus_callosum_external_capsule;
                #79,203 ee corpus_callosum_extreme_capsule;
                #80,204 ccg corpus_callosum_genu;
                #81,205 fp corpus_callosum_posterior_forceps;
                #82,206 ccb corpus_callosum_body;
                #83,207 ccs corpus_callosum_splenium;
                ccname = "cc" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "76,77,78,79,80,81,82,83,200,201,202,203,204,205,206,207",
                             tck_path, csvassignmentsname, ccname, "-files", "single", "-nthreads", str(num_threads)])
                
                #Example here :  Cortico-spinal tract cst with Pyramidal tract py, Pyramidal decussation pyd, Pons P, Medulla MY (84,208-87,211-88,212-39,163-40,164)
                cstname = "Cortico-spinal" # + "_" + os.path.basename(csv).replace(".csv", "")
                run_command(["connectome2tck", "-nodes", "84,208,87,211,88,212,39,163,40,164",
                             tck_path, csvassignmentsname, cstname, "-files", "single", "-nthreads", str(num_threads)])

            
           
            #Extract the voxels of the nodes used and also get the names right on each voxels extracted
        
            voxelextract_dir = os.path.join(connectome2tckdir,"VoxelExtraction")    
            nowvoxel = datetime.datetime.now()
            voxelextract_dir_date = voxelextract_dir + (nowvoxel.strftime("%Y%m%d_%H%M%S"))
            os.makedirs(voxelextract_dir_date, exist_ok=True)
            extract_and_crop_voxels_by_values(atlas_parcellation_path, labels, [26,150,38,162,46,170,76,77,78,79,80,81,82,83,200,201,202,203,204,205,206,207,84,208,87,211,88,212,39,163,40,164], voxelextract_dir_date)
        
    
        #TODO: for connectome visualization in mrview , below an example extracted from the BATMAN Tutorial    
        # label2mesh hcpmmp1_parcels_coreg.mif hcpmmp1_mesh.obj
        """connectome2tck sift_1mio.tck assignments_hcpmmp1.csv exemplar –files
        single –exemplars hcpmmp1_parcels_coreg.mif
        “exemplar” is the prefix name for your output file, such that the resulting file will be called
        “exemplar.tck”
        With “-files single” you specify that your output should be merged into one file. This is necessary
        for the downstream visualization in MRview."""
        # Generate a single track file containing edge exemplar trajectories:
        # $ connectome2tck tracks.tck assignments.txt exemplars.tck -files single -exemplars nodes.mif
        # This produces the track file that is required as input when attempting to display connectome edges using the streamlines or streamtubes geometries within the meview connectome tool.
        
        #TODO: connectome2stats (no example given by MRtrix)
        # Connectome group-wise statistics at the edge level using non-parametric permutation testing
        # connectomestats [ options ]  input algorithm design contrast output
        
        #     input: a text file listing the file names of the input connectomes
        #     algorithm: the algorithm to use in network-based clustering/enhancement. Options are: nbs, tfnbs, none
        #     design: the design matrix
        #     contrast: the contrast matrix
        #     output: the filename prefix for all output.

#4-1-END------------------------------------------------------------------------------------------------------------------------------------
    print(f"\n Mode {choice} done! Check your outputs in:")
    print(f"{project_dir_date}")
    
####ENDMODE1####
####MODE2####

if choice in {"2"} :
    
    print("Selecting WT folder containing all WT samples...")
    wt_folder = select_folder("Select WT folder")
    print("Selecting SPASTIN-KO folder containing all Spastin-KO samples...")
    spastin_folder = select_folder("Select SPASTIN-KO folder", os.path.dirname(wt_folder))
    output_folder = os.path.join(os.path.dirname(wt_folder),"ConnectomesAvgDiffpval")
    nowavgdiffpval = datetime.datetime.now()
    output_folder_date = output_folder + (nowavgdiffpval.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(output_folder_date, exist_ok=True)
    print(f"Created output folder for average, difference, significant differences and p-values connectomes at : {output_folder_date}")
    
    connectome_types = [
        "connectome_default.csv",
        "connectome_invnodevol.csv",
        "connectome_scalelength.csv",
        "connectome_scaleinvlength.csv"
    ]
    
    for connectome_type in connectome_types:
        process_connectome_type(connectome_type, wt_folder, spastin_folder, output_folder_date)
    
    print("\nAll connectome types processed successfully.")
    print(f"\n Mode {choice} done! Check your outputs in: {output_folder_date}")
    
###ENDMODE2####
####MODE3####

if choice in {"3"} :

#0-3-INITIALIZATION OF THE SCRIPT AND READING FILES FOR ATLAS COREGISTRATION-------------------------------------------
    
    print("Selecting folder containing connectomes to plot into structural connectivity matrices...")
    base_output_dir = filedialog.askdirectory(title="Select directory containing connectome csv files for plotting structural connectivity matrices")
    os.chdir(base_output_dir)
    print("Selecting LUT used for the connectomes present in the chosen folder...")
    label_txt_path = filedialog.askopenfilename(title="Select the LUT used for said connectomes (e.g. Allen Brain Metadata .txt file (grouped, lateralized,...))",filetypes=[("Atlas LUT", ".txt")])
    
    
    #Here we load the LookUpTable or LUT to get the regions' acronyms, fullnames and colors based on their ID...
    labels, colors, _ = load_LUT(label_txt_path)
    
    plotting_dir = os.path.join(base_output_dir, "ConnectomePlottingOnly_")
    now = datetime.datetime.now()
    plotting_dir_date = plotting_dir + (now.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(plotting_dir_date, exist_ok=True)
    print(f"Created output folder for structural connectivity matrices at : {plotting_dir_date}")
    
#0-3-END----------------------------------------------------------------------------------------------------------------------------------------------------
#1-3-CONNECTOME PLOTTING------------------------------------------------------------------------------------------------------------------------------------
    
    #Retrieving csv files' names to then plot them
    connectomeVector = []
    for files in os.listdir(base_output_dir):
        if fnmatch.fnmatch(files.lower(), "*connectome*.csv"):
            connectomeVector.append(files)
            
    for csv in connectomeVector:
        plot_connectome(csv, plotting_dir_date, labels=labels, colors=colors)
    
#1-3-END------------------------------------------------------------------------------------------------------------------------------------
    print(f"\n Mode {choice} done! Check your outputs in:")
    print(f"{plotting_dir_date}")

###ENDMODE3####
####MODE4####

if choice in {"4"} :

    # GUI: folder selection
    print("Selecting AvgWT folder containing the average connectomes of the WT condition...")
    wt_folder = filedialog.askdirectory(title="Select AvgWT folder containing average connectomes of the WT condition (usually in 'ConnectomesAvgDiffpval' folder)")
    print("Selecting AvgSpastin-KO folder containing the average connectomes of the Spastin-KO condition...")
    spastin_folder = filedialog.askdirectory(title="Select AvgSpastin-KO folder containing average connectomes of the Spastin-KO condition (usually in 'ConnectomesAvgDiffpval' folder)")
    output_folder = os.path.join(os.path.dirname(wt_folder),"Spastin-KOvsWT_SymmetryAnalysis_")
    nowsymmetry = datetime.datetime.now()
    output_folder_date = output_folder + (nowsymmetry.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(output_folder_date, exist_ok=True)
    print(f"Created output folder for symmetry analysis at : {output_folder_date}")

    print("Selecting LUT text file that was used for the connectomes selected...")
    lut_path = filedialog.askopenfilename(title="Select LUT file used for these connectomes (e.g. AllenBrain Metadata .txt file (grouped, lateralized,...))", filetypes=[("Atlas LUT", "*.txt")])
    
    # Load LUT and indices
    lut_labels, lut_colors, lut_fullnames = load_LUT(lut_path)
    cropped_labels = extract_clean_R_labels(lut_labels)
    right_indices = [i for i, lbl in enumerate(lut_labels) if lbl.endswith('_R')]
    
    gm_indices = [i for i in right_indices if is_gm(lut_labels[i])]
    wm_indices = [i for i in right_indices if is_wm(lut_labels[i])]
    gm_labels = [cropped_labels[i] for i in range(len(cropped_labels)) if right_indices[i] in gm_indices]
    wm_labels = [cropped_labels[i] for i in range(len(cropped_labels)) if right_indices[i] in wm_indices]
    gm_colors = [lut_colors[right_indices.index(i)] for i in gm_indices]
    wm_colors = [lut_colors[right_indices.index(i)] for i in wm_indices]
    
    # File mapping
    wt_files = get_file_map(wt_folder, "WT")
    spastin_files = get_file_map(spastin_folder, "Spastin-KO")
    common_keys = sorted(set(wt_files.keys()) & set(spastin_files.keys()))
    
    # Main processing loop
    for base in common_keys:
        print(f"Processing: {base}")
        wt_mat = pd.read_csv(wt_files[base], header=None).values
        spastin_mat = pd.read_csv(spastin_files[base], header=None).values
    
        wt_sym = compute_symmetry_matrix(wt_mat)
        spastin_sym = compute_symmetry_matrix(spastin_mat)
    
        vmax = max(np.abs(wt_sym).max(), np.abs(spastin_sym).max())
    
        # Full symmetry matrix
        plot_symmetry_pair(spastin_sym, wt_sym, cropped_labels,
                           lut_colors[:len(cropped_labels)],
                           f"Symmetry Analysis: {base}",
                           os.path.join(output_folder_date, f"{base}_Symmetry_Full.png"),
                           vmax)
    
        # GM-only
        spastin_gm = spastin_sym[np.ix_(gm_indices, gm_indices)]
        wt_gm = wt_sym[np.ix_(gm_indices, gm_indices)]
        plot_symmetry_pair(spastin_gm, wt_gm, gm_labels, gm_colors,
                           f"GM Symmetry: {base}",
                           os.path.join(output_folder_date, f"{base}_Symmetry_GM.png"),
                           vmax)
    
        # WM-only
        spastin_wm = spastin_sym[np.ix_(wm_indices, wm_indices)]
        wt_wm = wt_sym[np.ix_(wm_indices, wm_indices)]
        plot_symmetry_pair(spastin_wm, wt_wm, wm_labels, wm_colors,
                           f"WM Symmetry: {base}",
                           os.path.join(output_folder_date, f"{base}_Symmetry_WM.png"),
                           vmax)
    
        print(f"Saved all plots for {base}")
        
    print(f"\n Mode {choice} done! Check your outputs in:")
    print(f"{output_folder_date}")

###ENDMODE4####
####MODE5####

if choice in {"5"} :

    print("Selecting folder containing connectomes to plot into connectograms...")
    output_dir = filedialog.askdirectory(title="Select directory containing connectome csv files for plotting connectograms")
    
    print("Selecting LUT text file that was used for the connectomes selected...")
    lut_path = filedialog.askopenfilename(title="Select LUT file used for these connectomes (e.g. AllenBrain Metadata .txt file (grouped, lateralized,...))", filetypes=[("Atlas LUT", "*.txt")])

    labels, colors, _ = load_LUT(lut_path)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_subdir = os.path.join(output_dir, f"ConnectogramBatch_{timestamp}")
    os.makedirs(output_subdir, exist_ok=True)
    print(f"Created output folder for connectograms at : {output_subdir}")


    for fname in os.listdir(output_dir):
        if fname.lower().endswith(".csv"):
            full_path = os.path.join(output_dir, fname)
            print(f"Processing: {fname}")
            try:
                matrix = pd.read_csv(full_path, header=None).values
                if matrix.shape[0] != len(labels):
                    print(f"Skipped {fname}: size mismatch with LUT ({matrix.shape[0]} vs {len(labels)})")
                    continue
                output_img = os.path.join(output_subdir, fname.replace(".csv", "_connectogram.png"))
                plot_connectogram_circle(matrix, labels, colors, output_img)
            except Exception as e:
                print(f"Error processing {fname}: {e}")
                
    print(f"\n Mode {choice} done! Check your outputs in:")
    print(f"{output_subdir}")

###ENDMODE5####