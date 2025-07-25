import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog
from scipy.stats import ttest_ind
import re
from datetime import datetime

def load_LUT(path):
    labels, colors, acronyms = [], [], []
    label_map, color_map = {}, {}

    with open(path) as f:
        for line in f:
            parts = re.findall(r'"([^"]+)"|(\d+)', line)
            flat = [x for tup in parts for x in tup if x]
            if len(flat) < 4:
                continue
            idx = int(flat[0])
            acronym = flat[1]
            fullname = flat[2]
            hex_color = f"#{flat[3].strip('#')}"
            labels.append(idx)
            acronyms.append(acronym)
            colors.append(hex_color)
            label_map[idx] = acronym
            color_map[idx] = hex_color

    return labels, acronyms, colors, label_map, color_map

def load_parcellation(file_path):
    img = nib.load(file_path)
    data = img.get_fdata()
    voxel_volume = np.prod(img.header.get_zooms())  # in mm³
    return data, voxel_volume

def compute_volumes(parcellation_data, voxel_volume):
    labels, counts = np.unique(parcellation_data, return_counts=True)
    volumes = {int(label): count * voxel_volume for label, count in zip(labels, counts) if label != 0}
    return volumes

def extract_parcellation_volumes(group_folder):
    sample_names = sorted([
        d for d in os.listdir(group_folder)
        if os.path.isdir(os.path.join(group_folder, d))
    ])
    
    volumes_list = []
    sample_ids = []

    for sample in sample_names:
        sample_path = os.path.join(group_folder, sample)
        connectome_dirs = [
            d for d in os.listdir(sample_path)
            if d.startswith("ConnectomeProject_") and os.path.isdir(os.path.join(sample_path, d))
        ]
        if not connectome_dirs:
            print(f"Warning: No ConnectomeProject_* folder in {sample}")
            continue
        
        # Take the first (or only) connectome directory
        connectome_path = os.path.join(sample_path, connectome_dirs[0])
        atlas_path = os.path.join(connectome_path, "atlas_labels_in_CS.nii")
        if not os.path.exists(atlas_path):
            print(f"Warning: atlas_labels_in_CS.nii not found in {connectome_path}")
            continue

        data, voxel_volume = load_parcellation(atlas_path)
        vols = compute_volumes(data, voxel_volume)
        volumes_list.append(vols)
        sample_ids.append(sample)

    # Collect all labels and build dataframe
    all_labels = sorted(set(k for d in volumes_list for k in d.keys()))
    df = pd.DataFrame([{label: v.get(label, 0) for label in all_labels} for v in volumes_list])
    df['subject'] = sample_ids
    return df.set_index('subject')


def perform_ttests(wt_df, sp_df, label_map):
    results = []
    for col in wt_df.columns:
        stat, pval = ttest_ind(wt_df[col], sp_df[col], equal_var=False)
        results.append((col, label_map.get(col, str(col)), pval))
    df = pd.DataFrame(results, columns=['label', 'region', 'p_value'])
    return df

def extract_clean_R_labels(lut_labels):
    """Extract only right hemisphere labels and crop '_R'."""
    return [lbl[:-2] for lbl in lut_labels if lbl.endswith('_R')]  # safer than replace


def sum_left_right_volumes(df, label_map):
    """
    Sum left and right volumes (df columns are integer labels) for lateralized regions.
    Uses only right labels (ending with '_R' in acronyms) as keys, cropping '_R'.
    Non-lateralized regions (no _L/_R suffix) are kept as is.
    Returns a DataFrame with integer label columns (right label IDs for lateralized, or original label IDs for others).
    """
    # Invert label_map: acronym -> int_label
    acronym_to_label = {v: k for k, v in label_map.items()}

    # Extract right side cleaned labels (without '_R')
    clean_right_labels = extract_clean_R_labels(acronym_to_label.keys())

    summed_data = {}

    for base_acronym in clean_right_labels:
        left_ac = base_acronym + '_L'
        right_ac = base_acronym + '_R'

        left_label = acronym_to_label.get(left_ac, None)
        right_label = acronym_to_label.get(right_ac, None)

        # Average left and right if both exist, else use whichever exists
        if left_label in df.columns and right_label in df.columns:
            summed_vol = (df[left_label] + df[right_label])
            summed_data[right_label] = summed_vol
        elif right_label in df.columns:
            summed_data[right_label] = df[right_label]
        elif left_label in df.columns:
            summed_data[left_label] = df[left_label]
        # else: no data for this region; skip

    # Add non-lateralized regions (acronyms without _L or _R)
    non_lateralized = [ac for ac in acronym_to_label.keys() if not (ac.endswith('_L') or ac.endswith('_R'))]
    for ac in non_lateralized:
        label = acronym_to_label[ac]
        if label in df.columns:
            summed_data[label] = df[label]

    summed_df = pd.DataFrame(summed_data, index=df.index)
    # Sort columns for consistent ordering
    summed_df = summed_df.reindex(sorted(summed_df.columns), axis=1)
    return summed_df



def create_combined_label_and_color_maps(label_map, color_map):
    acronym_to_label = {v: k for k, v in label_map.items()}
    combined_label_map = {}
    combined_color_map = {}

    clean_right_labels = extract_clean_R_labels(acronym_to_label.keys())

    # For lateralized regions, keep the base acronym and use right side color
    counter_idx = 1 #start at 1 no 0 value in the atlas for connectomes
    for base_ac in clean_right_labels:
        right_ac = base_ac + '_R'
        right_label = acronym_to_label.get(right_ac)
        combined_label_map[counter_idx] = base_ac
        combined_color_map[counter_idx] = color_map.get(right_label, '#000000')
        counter_idx+=1

    return combined_label_map, combined_color_map




def plot_volumes(wt_df, sp_df, label_map, color_map, ttest_df, acronyms, colors, output_dir):
    x = np.arange(len(wt_df.columns))
    labels = wt_df.columns
    wt_means = wt_df.mean()
    sp_means = sp_df.mean()

    fig, ax = plt.subplots(figsize=(18, 6))

    # Scatter
    for i, label in enumerate(labels):
        ax.scatter([x[i]] * len(wt_df), wt_df[label], color='blue', alpha=0.6)
        ax.scatter([x[i]] * len(sp_df), sp_df[label], color='red', alpha=0.6)

    # Cumulative mean lines
    ax.plot(x, wt_means.cumsum(), color='blue', linewidth=2)
    ax.fill_between(x, 0, wt_means.cumsum(), color='blue', alpha=0.1)

    ax.plot(x, sp_means.cumsum(), color='red', linewidth=2)
    ax.fill_between(x, 0, sp_means.cumsum(), color='red', alpha=0.1)

    # Significance stars
    pvals = ttest_df.set_index('label')['p_value']
    for i, label in enumerate(labels):
        p = pvals.get(label, 1.0)
        if p <= 0.001:
            star = '***'
        elif p <= 0.01:
            star = '**'
        elif p <= 0.05:
            star = '*'
        else:
            continue
        max_val = max(wt_df[label].max(), sp_df[label].max())
        ax.text(x[i], max_val + 5, star, ha='center', va='bottom', fontsize=12, color='black')

    ax.set_ylabel("Volume (mm³)")
    ax.set_title("Region Volumes: WT vs Spastin-KO")
    ax.set_xticks(x)

    for i, (acro, color) in enumerate(zip(acronyms, colors)):
        ax.text(x[i], -10, acro, rotation=90, ha='center', va='top', color=color, fontsize=8)

    ax.set_xticklabels([''] * len(labels))  # hide tick labels, replaced with text
    ax.legend(['WT (mean)', 'SPASTIN-KO (mean)'])
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "volumetry_comparison_plot.png"), dpi=300)
    plt.show()


def plot_volumesV2(wt_df, spastin_df, label_map, color_map, ttest_results, output_dir=None, suffix=None):
    labels = list(wt_df.columns)
    label_names = [label_map.get(l, str(l)) for l in labels]
    x = np.arange(len(labels))
    
    wt_means = wt_df.mean()
    sp_means = sp_df.mean()

    fig, ax1 = plt.subplots(figsize=(20, 6))

    # Scatter: blue for WT, red for KO
    ax1.scatter([], [], color='blue', label=f'WT N={len(wt_df)}')  # dummy for legend
    ax1.scatter([], [], color='red', label=f'Spastin-KO N={len(sp_df)}')
    for i, label in enumerate(labels):
        ax1.scatter([x[i]], wt_means[label], color='blue', alpha=0.5, s=30)
        ax1.scatter([x[i]], sp_means[label], color='red', alpha=0.5, s=30)
        # ax1.scatter([x[i]] * len(wt_df), wt_df[label], color='blue', alpha=0.5, s=30)
        # ax1.scatter([x[i]] * len(spastin_df), spastin_df[label], color='red', alpha=0.5, s=30)

    ax1.set_ylabel('Region Volume (mm³)', color='black')

    # Twin axis for cumulative volume
    ax2 = ax1.twinx()
    wt_cum_mean = wt_df.mean(axis=0).cumsum()
    sp_cum_mean = spastin_df.mean(axis=0).cumsum()

    ax2.plot(x, wt_cum_mean, color='blue', linewidth=2, alpha=0.7, label='WT cumulative')
    ax2.fill_between(x, 0, wt_cum_mean, color='blue', alpha=0.1)

    ax2.plot(x, sp_cum_mean, color='red', linewidth=2, alpha=0.7, label='Spastin-KO cumulative')
    ax2.fill_between(x, 0, sp_cum_mean, color='red', alpha=0.1)

    ax2.set_ylabel('Cumulative Volume (mm³)', color='black')

    # Significance stars (on primary y-axis)
    pval_map = ttest_results.set_index('label')['p_value'].to_dict()
    for i, label in enumerate(labels):
        p = pval_map.get(label, 1.0)
        if p <= 0.001:
            star = '***'
        elif p <= 0.01:
            star = '**'
        elif p <= 0.05:
            star = '*'
        else:
            continue
        max_val = max(wt_df[label].max(), spastin_df[label].max())
        ax1.text(x[i], max_val + 5, star, ha='center', va='bottom', fontsize=12, color='black')

    # X-axis labels using LUT colors
    ax1.set_xticks(x)
    for i, label in enumerate(labels):
        acronym = label_map.get(label, str(label))
        color = color_map.get(label, '#000000')
        ax1.text(
            x[i], ax1.get_ylim()[0] - 0.03 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]),
            acronym, color=color, ha='right', va='top', rotation=45, fontsize=7, clip_on=False
        )

    ax1.set_xticklabels([''] * len(labels))  # Hide default ticks
    ax2.set_xticks([])  # Prevent xticks from appearing on the twin axis

    # Legend
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

    ax1.set_title('Region Volumes & Cumulative Volumes (WT vs Spastin-KO)')
    plt.tight_layout()
    if output_dir:
        plt.savefig(os.path.join(output_dir, f"volumetry_comparison_plot_Spastin-KOvsWT_{suffix}.jpg"), dpi=300)
    plt.show()


########## MAIN SCRIPT ##########

# File/folder selection GUI
root = Tk()
root.withdraw()
print("Select WT folder")
wt_folder = filedialog.askdirectory(title="Select WT Folder")
print("Select SPASTIN-KO folder")
spastin_folder = filedialog.askdirectory(title="Select SPASTIN-KO Folder")
print("Select label LUT matching parcellation used (Allen Metadata .txt file, grouped and/or lateralized LR)")
lut_path = filedialog.askopenfilename(title="Select LUT (Allen Metadata .txt file, grouped and/or lateralized LR)")
print("Default output directory located in common root of WT and Spastin-KO folders")
output_base = os.path.dirname(wt_folder)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = os.path.join(output_base, f"Volumetry_{timestamp}")
os.makedirs(output_dir, exist_ok=True)

# Load LUT
labels, acronyms, colors, label_map, color_map = load_LUT(lut_path)

# Process groups
print("Processing WT group...")
wt_df = extract_parcellation_volumes(wt_folder)
print("Processing SPASTIN group...")
sp_df = extract_parcellation_volumes(spastin_folder)

print(f"WT samples processed: {len(wt_df)}")
print(f"Spastin-KO samples processed: {len(sp_df)}")
# Align columns
wt_df = wt_df[labels]
sp_df = sp_df[labels]

# Sum left/right volumes, keep only right labels combined with left:
wt_df_sum = sum_left_right_volumes(wt_df, label_map)
sp_df_sum = sum_left_right_volumes(sp_df, label_map)
combined_label_map, combined_color_map = create_combined_label_and_color_maps(label_map, color_map)


# T-tests
print("Performing t-tests...")
ttest_df = perform_ttests(wt_df, sp_df, label_map)
ttest_df.to_csv(os.path.join(output_dir, "ttest_results.csv"), index=False)
ttest_df_grouped = perform_ttests(wt_df_sum, sp_df_sum, combined_label_map)
ttest_df_grouped.to_csv(os.path.join(output_dir, "ttest_results.csv"), index=False)

# Plot
print("Plotting volumetry...")
suffixLR = "Left-Right"
suffixSym = "Symmetrical"
plot_volumesV2(wt_df, sp_df, label_map, color_map, ttest_df, output_dir, suffixLR)
plot_volumesV2(wt_df_sum, sp_df_sum, combined_label_map, combined_color_map, ttest_df_grouped, output_dir, suffixSym)

# Save data
wt_df.to_csv(os.path.join(output_dir, "wt_volumes_LeftRight.csv"), index=False)
sp_df.to_csv(os.path.join(output_dir, "spastin_volumes_LeftRight.csv"), index=False)
wt_df_sum.to_csv(os.path.join(output_dir, "wt_volumes_symmetrical.csv"), index=False)
sp_df_sum.to_csv(os.path.join(output_dir, "spastin_volumes_symmetrical.csv"), index=False)

print(f"Analysis complete. Results saved in: {output_dir}")
