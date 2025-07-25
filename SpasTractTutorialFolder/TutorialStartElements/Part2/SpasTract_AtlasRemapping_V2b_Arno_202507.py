import nibabel as nib
import numpy as np
import os
import re
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

# Regions to exclude from the *indexed* file (acronym + fullname)
EXCLUDED_LABELS = {
    ('root', 'root'),
    ('fiber-tracts', 'fiber tracts'),
    ('tspd', 'direct tectospinal pathway'),
    ('c', 'central canal, spinal cord/medulla')
}
# ('VL', 'lateral ventricle'),
# ('V3', 'third ventricle'),
# ('V4', 'fourth ventricle'),
# ('V4r', 'lateral recess'),
# ('SEZ', 'subependymal zone'),
# ('AQ', 'cerebral aqueduct'),
# 

def parse_allen_label_file(txt_path):
    """Parses Allen-style atlas label file into label and color dictionaries."""
    labels = {}
    colors = {}

    with open(txt_path) as f:
        for line in f:
            match = re.match(r"\{'(\d+)',\s*'([0-9A-Fa-f]{6})',\s*'([^']+)',\s*'([^']+)'\}", line.strip())
            if match:
                region_id = int(match.group(1))
                hex_color = match.group(2).upper()
                acronym = match.group(3)
                fullname = match.group(4)

                labels[region_id] = (acronym, fullname)
                colors[region_id] = hex_color

    if not labels:
        raise ValueError("No valid entries parsed from label file.")

    return labels, colors

def get_unique_values_from_nifti(nifti_path):
    img = nib.load(nifti_path)
    data = img.get_fdata()
    unique_vals = np.unique(data.astype(int))
    return unique_vals.tolist()

def write_reformatted_label_file(original_labels, colors, output_path):
    """Write all label entries with original region IDs."""
    with open(output_path, "w") as f:
        for region_id in sorted(original_labels.keys()):
            acronym, fullname = original_labels[region_id]
            fullname = fullname.replace(" ", "_")
            color = colors.get(region_id, "999999")
            f.write(f'{region_id} "{acronym}" "{fullname}" "{color}"\n')

def write_indexed_label_file(unique_values, original_labels, colors, output_path):
    """
    Write only regions present in the NIfTI image, excluding the unwanted ones,
    and assign continuous indices starting from 1.
    """
    sorted_values = sorted(val for val in unique_values if val in original_labels)

    with open(output_path, "w") as f:
        idx = 1  # Start at 1
        for val in sorted_values:
            acronym, fullname = original_labels[val]
            if (acronym, fullname) in EXCLUDED_LABELS:
                continue  # Skip excluded regions

            fullname_clean = fullname.replace(" ", "_")
            color = colors.get(val, "999999")
            f.write(f'{idx} "{acronym}" "{fullname_clean}" "{color}"\n')
            idx += 1

def load_label_table(filepath):
    """
    Load a label file formatted as:
    value acronym fullname hexacolor
    Returns a dict: value -> (acronym, fullname)
    """
    label_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            value = int(parts[0])
            acronym = parts[1]
            fullname = parts[2]
            label_dict[value] = (acronym, fullname)
    return label_dict

def build_remap_dict(file1_path, file2_path):
    """
    Returns a dict: original_value (in file1) -> new_index (in file2)
    based on matching acronym + fullname
    """
    file1 = load_label_table(file1_path)
    file2 = load_label_table(file2_path)

    # Reverse lookup: (acronym, fullname) -> new_index
    reverse_lookup = {v: k for k, v in file2.items()}

    remap_dict = {}
    for orig_value, label_tuple in file1.items():
        if label_tuple in reverse_lookup:
            remap_dict[orig_value] = reverse_lookup[label_tuple]
    return remap_dict

def remap_nifti_values(nifti_path, remap_dict, output_path):
    img = nib.load(nifti_path)
    data = img.get_fdata().astype(np.int32)

    remapped = np.zeros_like(data, dtype=np.int16)

    for orig_val, new_val in remap_dict.items():
        remapped[data == orig_val] = new_val

    new_img = nib.Nifti1Image(remapped, affine=img.affine, header=img.header)
    nib.save(new_img, output_path)
    print(f"Remapped NIfTI saved to: {output_path}")


def load_indexed_table(path):
    d = {}
    with open(path) as f:
        for line in f:
            parts = re.findall(r'"([^"]+)"|(\d+)', line)
            flat = [x for tup in parts for x in tup if x]
            if len(flat) < 4: continue
            idx, acronym, fullname, color = int(flat[0]), flat[1], flat[2], flat[3]
            d[idx] = (acronym, fullname, color)
    return d


def extract_stem(acronym, fullname):
    # Manual exceptions
    cc_exceptions = {"fa", "ec", "ee", "ccg", "fp", "ccb", "ccs"}
    thalamus_group = {"TH", "VAL", "VM", "VPL", "VPM", "PoT", "SPF", "SPA", "PP", "MG", "LG", "LP", "PO", "POL", "SGN", "Eth", "AV", "AM", "AD", "IAM",
                      "IAD", "LD", "IMD", "MD", "SMT", "PR", "PVT", "PT", "RE", "Xi", "RH", "CM", "PCN", "CL", "PF", "PIL", "RT", "IGL", "IntG", "SubG",
                      "MH", "LH"}
    hypothalamus_group = {"HY", "SO", "ASO", "PVH", "PV", "ARH", "ADP", "AVP", "AVPV", "DMH", "MEPO", "MPO", "OV", "PD", "PS", "SBPV", "SCH", "SFO", 
                          "VMPO", "VLPO", "AHN", "LM", "Mmme", "Mml", "Mmm", "Mmp", "Mmd", "SUM", "TM", "MPN", "PM", "VMH", "PH", "LHA", "LPO", 
                          "PST", "PSTN", "PeF", "RCH", "STN", "TU", "ZI", "FF", "ME"} 
    olfactory_group = {"OLF", "MOB", "AOB", "AON", "TT", "DP", "PIR", "NLOT", "BA", "OT"}
    striatum_group = {"STR", "CP", "ACB", "FS", "SH"}
    BST_group = {"BST", "BAC"}
    amygdala_group = {"COA", "PAA", "LA", "BLA", "BMA", "PA", "AAA", "CEA", "IA", "MEA"} #to AMG later with a posteriori modification of COA
    medial_septum_group = {"MS", "NDB"}
    hippocampalformation_group = {"HPF", "CA", "DG", "FC", "IG", "ENT", "PAR", "POST", "PRE", "SUB", "ProS"}
    midbrain_group = {"MB", "SC", "IC", "NB", "SAG", "PBG", "MEV", "SCO", "SN", "VTA", "PN", "RR", "MRN", "PAG", "PRC", "INC", "ND", "Su3", "APN", "MPT",
                      "NOT", "NPC", "OP", "PPT", "RPF", "CUN", "RN", "III", "EW", "IV", "Pa4", "VTN", "AT", "LT", "DT", "MT", "PPN", "IF",
                      "IPN", "IPR", "IPC", "IPA", "IPL", "IPI", "IPDM", "IPDL", "IPRL", "RL", "CLI", "DR"}
    pons_group = {"P", "NLL", "PSV", "PB", "KF", "POR", "SOC", "DTN", "PDT", "PCG", "PG", "PRN", "SG", "SUT", "TRN", "V", "P5", "Acs5", "PC",
                  "I5", "CS", "LC", "LDT", "NI", "RPO", "SLC", "SLD"}
    medulla_group = {"MY", "DCO", "VCO", "CU", "GR", "ECU", "NTB", "NTS", "SPVC", "SPVI", "SPVO", "Pa5", "VI", "VII", "ACVII", "AMB", "DMX", "GRN",
                     "ICB", "IO", "IRN", "ISN", "LIN", "LRN", "MARN", "MDRN", "PARN", "PAS", "PGRN", "NR", "PRP", "PPY", "LAV", "MV", "SPIV", "SUV",
                     "x", "XII", "y", "RM", "RPA", "RO"} 
    cerebellum_group = {"CB","LING", "CENT", "CUL", "DEC", "FOTU", "PYR", "UVU", "NOD", "SIM", "AN", "PRM", "COPY", "PFL", "FL", "FN", "IP", "DN", "VeCB"}
    
    if acronym in olfactory_group:
        return "OLF"                  
    if acronym in thalamus_group:
        return "TH"
    if acronym in hypothalamus_group:
        return "HY"
    if acronym in amygdala_group:
        return "AMG"
    if acronym in hippocampalformation_group:
        return "HPF"
    if acronym in medial_septum_group:
        return "MS"
    if acronym in BST_group:
        return "BST"
    if acronym in striatum_group:
        return "STR"
    if acronym in midbrain_group:
        return "MB"
    if acronym in pons_group:
        return "P"
    if acronym in medulla_group:
        return "MY"
    if acronym in cerebellum_group:
        return "CB"
    if acronym in cc_exceptions:
        return acronym

    # General case: extract longest prefix of uppercase letters (at least 3 chars if possible)
    m = re.match(r'^([A-Z]{3,})', acronym)
    if m:
        return m.group(1)
    
    # Fallback: use first 2 uppercase letters if available
    m = re.match(r'^([A-Z]{2})', acronym)
    if m:
        return m.group(1)

    return acronym


def get_groups(label_dict):
    # First pass: build initial mapping from index to stem
    first_pass = {}
    for idx, (acr, full, col) in label_dict.items():
        if re.match(r'^[a-z0-9]', acr):  # lowercase/WM: keep as-is
            first_pass[idx] = acr
        else:
            first_pass[idx] = extract_stem(acr, full)

    # Second pass: reprocess stems to group subregions under shared higher-order group
    second_pass = {}
    for idx, stem in first_pass.items():
        second_pass[idx] = extract_stem(stem, "")  # fullname not needed here

    # Now build groups based on second pass
    groups = {}
    order = []
    for idx, (acr, full, col) in label_dict.items():
        key = second_pass[idx]
        if key not in groups:
            groups[key] = {
                "members": [],
                "color": col,
                "fullname": full.split(",")[0]
            }
            order.append(key)
        groups[key]["members"].append(idx)

    return order, groups


def build_grouped_table(order, groups):
    lut = []
    mapping = {}
    idx = 1
    for key in order:
        full = groups[key]["fullname"].replace(" ", "_")
        col = groups[key]["color"]
        lut.append((idx, key, full, col))
        mapping[key] = idx
        idx += 1
    return lut, mapping

def remap_to_grouped(img, label_dict, groups, mapping):
    data = img.get_fdata().astype(int)
    new_data = np.zeros_like(data, dtype=int)

    group_map = {}
    for group_key, info in groups.items():
        for member_id in info["members"]:
            group_map[member_id] = mapping[group_key]

    for old_val, new_val in group_map.items():
        new_data[data == old_val] = new_val

    return new_data


def load_grouped_lut(path):
    lut = {}
    with open(path) as f:
        for line in f:
            parts = re.findall(r'"([^"]+)"|(\d+)', line)
            flat = [x for tup in parts for x in tup if x]
            if len(flat) < 4:
                continue
            idx = int(flat[0])
            acronym, fullname, color = flat[1], flat[2], flat[3]
            lut[idx] = {'acronym': acronym, 'fullname': fullname, 'color': color}
    return lut


def compute_symmetrical_hemisphere(annotation_img, label_info, midline_offset=0):
    data = annotation_img.get_fdata().astype(int)
    new_data = np.zeros_like(data, dtype=int)

    # Determine left-right axis from affine matrix
    lr_axis = int(np.argmax(np.abs(annotation_img.affine[:3, 0])))
    midline = data.shape[lr_axis] // 2 + midline_offset

    new_labels = {}
    missing_regions = []

    # Map from old region id to new left and right ids
    # Right side keeps original ID
    # Left side gets new ID starting after max existing ID
    max_id = max(label_info.keys())
    next_left_id = max_id + 1

    # Prepare masks for left and right hemispheres
    left_mask = np.zeros_like(data, dtype=bool)
    right_mask = np.zeros_like(data, dtype=bool)

    # Hemisphere masks by slicing the data array
    slicer_left = [slice(None)] * data.ndim
    slicer_right = [slice(None)] * data.ndim
    slicer_left[lr_axis] = slice(0, midline)
    slicer_right[lr_axis] = slice(midline, data.shape[lr_axis])

    left_mask[tuple(slicer_left)] = True
    right_mask[tuple(slicer_right)] = True

    # For each region in label_info
    for region_id, info in label_info.items():
        region_mask = (data == region_id)

        # Voxels on each hemisphere
        region_left = region_mask & left_mask
        region_right = region_mask & right_mask

        acronym = info['acronym']
        fullname = info['fullname']
        color = info['color']

        # Right hemisphere keeps original ID
        if np.any(region_right):
            new_data[region_right] = region_id
            new_labels[region_id] = {
                'acronym': acronym + '_R',
                'color': color,
                'fullname': fullname + '_Right'
            }
        else:
            # Still create label for right hemisphere with no voxels
            new_labels[region_id] = {
                'acronym': acronym + '_R',
                'color': color,
                'fullname': fullname + '_Right'
            }

        # Left hemisphere gets new ID regardless of presence
        new_data[region_left] = next_left_id
        new_labels[next_left_id] = {
            'acronym': acronym + '_L',
            'color': color,
            'fullname': fullname + '_Left'
        }
        next_left_id += 1

    return new_data, new_labels, missing_regions

# def close_window():
#     win.quit()
#     win.destroy()


########################
########################
# MAIN SCRIPT
########################
########################
#0-INITIALIZATION OF THE SCRIPT AND READING FILES FOR NEW LUTs WITH ACTUAL VALUES IN PARCELLATION (PARCELLATION & LUT FILES)-------------------------------------------

print("SpasTract Connectomics : Atlas Remapping, Creating new LUTs based on actual values in parcellation")
num_threads = os.cpu_count()
print(f"Using {num_threads} CPU threads")
os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(num_threads)


root = tk.Tk()
root.withdraw()

#TODO : functional menu that allows to exit without bugging
# win = tk.Toplevel()
# win.title("Choose Remapping Type")
# win.geometry("400x400")

# ttk.Label(win, text="Select Remapping of parcellation :").pack(pady=10)

# options = [
#     "Continuous Parcellation",
#     "Continuous Grouped Parcellation",
#     "Continuous Lateralized Parcellation",
#     "Continuous Grouped Lateralized Parcellation"
# ]

# combo = ttk.Combobox(win, values=options, state="readonly")
# combo.current(0)
# combo.pack()

# ttk.Label(win, text="Then click OK to continue:").pack(pady=10)

# ok_button = ttk.Button(win, text="OK", command=close_window)
# ok_button.pack(pady=5)

# win.mainloop()

# # Get the selection AFTER the window is closed
# selection = combo.get()
# print(f"User selected: {selection}")

print("Selecting original parcellation NIfTI (e.g. P56_Atlas_Parcellation_REOR.nii)...")
nifti_path = filedialog.askopenfilename(title="Select original parcellation NIfTI (e.g. P56_Atlas_Parcellation_REOR.nii)", filetypes=[("NIfTI", "*.nii")]) #.nii.gz not compatible with MacOS + Nibabel
print("Selecting original LUT (e.g. AllenBrain_Metadata_LabelHexColorShortFullname.txt file)...")
label_txt_path = filedialog.askopenfilename(title="Select original LUT (e.g. AllenBrain_Metadata_LabelHexColorShortFullname.txt file)", filetypes=[("Text files", "*.txt")])

if nifti_path and label_txt_path:
    output_dir = os.path.dirname(nifti_path)

    output_file1 = os.path.join(output_dir, "AllenBrain_Metadata_labelsshortfullhexacolor_reformatted.txt")
    print("Finished extracting existing labels of selected parcellation (AllenBrain_Metadata_labelsshortfullhexacolor_reformatted.txt)")
    output_file2 = os.path.join(output_dir, "AllenBrain_Metadata_labelsshortfullhexacolor_indexcontinuous.txt")
    print("Finished sorting existing labels of selected parcellation so that they are continuous, starting from 1 (AllenBrain_Metadata_labelsshortfullhexacolor_indexcontinuous.txt)")


    labels, colors = parse_allen_label_file(label_txt_path)
    unique_vals = get_unique_values_from_nifti(nifti_path)

    write_reformatted_label_file(labels, colors, output_file1)  # All regions, no filtering
    write_indexed_label_file(unique_vals, labels, colors, output_file2)  # Filtered + reindexed starting at 1

    print(f"Reformatted label file (based on exisiting NIfTI parcellation values) written to: {output_file1}")
    print(f"Indexed label file (based on NIfTI values, sorted into continuous values) written to: {output_file2}")

#0-END-----------------------------------------------------------------------------------------------------------------------------------
#1-REMAPPING ATLAS PARCELLATION WITH CONTINUOUS VALUES FOR MRTRIX------------------------------------------------------------------------

print("SpasTract Connectomics : Atlas Remapping, Remapping parcellation with continuous values for MRtrix")

print("Selecting Reformatted LUT file (AllenBrain_Metadata_labelsshortfullhexacolor_reformatted.txt)...")
file1_path = filedialog.askopenfilename(title="Select Reformatted LUT file (AllenBrain_Metadata_labelsshortfullhexacolor_reformatted.txt)", filetypes=[("text files", "*.txt")])

print("Selecting Indexed Continuous LUT file (AllenBrain_Metadata_labelsshortfullhexacolor_indexcontinuous.txt)...")
file2_path = filedialog.askopenfilename(title="Select Indexed Continuous LUT file (AllenBrain_Metadata_labelsshortfullhexacolor_indexcontinuous.txt)", filetypes=[("text files", "*.txt")])

# === Remap and Save ===
remap_dict = build_remap_dict(file1_path, file2_path)

output_dir = os.path.dirname(nifti_path)
basename = os.path.basename(nifti_path).replace(".nii.gz", "").replace(".nii", "")
output_continuous_path = os.path.join(output_dir, f"{basename}_remapped_continuous.nii")

remap_nifti_values(nifti_path, remap_dict, output_continuous_path)

#1-END-----------------------------------------------------------------------------------------------------------------------------------
#2-REMAPPING ATLAS PARCELLATION BY GROUPING SMALL REGIONS (OPTIONAL)---------------------------------------------------------------------

# if selection == "Continuous Grouped Parcellation" or selection == "Continuous Grouped Lateralized Parcellation" : 
group_parcellation = messagebox.askyesno("Group small regions in parcellation?", "Do you want to group small regions in the parcellation for better matrix visualization ?")
if group_parcellation : 
           
    print("SpasTract Connectomics : Atlas Remapping, Remapping parcellation by grouping small regions")
    print("You chose to indeed group small regions together, now processing...")
    
    nii_in = output_continuous_path
    lut_in = output_file2
    
    nii_prefix = os.path.splitext(os.path.basename(nii_in))[0]
    nii_out_path = os.path.join(os.path.split(nii_in)[0], nii_prefix)
    
    lut_prefix = os.path.splitext(os.path.basename(lut_in))[0]
    lut_out_path = os.path.join(os.path.split(nii_in)[0], lut_prefix)
    
    # Load and process
    label_dict = load_indexed_table(lut_in)
    order, groups = get_groups(label_dict)
    lut, mapping = build_grouped_table(order, groups)
    
    img = nib.load(nii_in)
    new_data = remap_to_grouped(img, label_dict, groups, mapping)
    
    # Save parcellation
    new_img = nib.Nifti1Image(new_data.astype(np.int16), img.affine, img.header)
    nib.save(new_img, f"{nii_out_path}_grouped.nii")
    
    # === Patch grouped LUT entries a posteriori ===
    lut_patched = []
    for idx, acr, full, col in lut:
        if acr == "AMG" and full == "Cortical_amygdalar_area":
            acr, full = "AMG", "Amygdala"
        elif acr == "AUD" and full == "Dorsal_auditory_area":
            full = "Auditory_area"
        elif acr == "VIS" and full == "Anterolateral_visual_area":
            full = "Visual_area"
        elif acr == "fa" : 
            full = "corpus_callosum_anterior_forceps"
        elif acr == "ec" : 
            full = "corpus_callosum_external_capsule"
        elif acr == "ee" : 
            full = "corpus_callosum_extreme_capsule"
        elif acr == "ccg" : 
            full = "corpus_callosum_genu"
        elif acr == "fp" : 
            full = "corpus_callosum_posterior_forceps"
        elif acr == "ccb" : 
            full = "corpus_callosum_body"
        elif acr == "ccs" : 
            full = "corpus_callosum_splenium"
        lut_patched.append((idx, acr, full, col))
    
    # Save grouped LUT
    with open(f"{lut_out_path}_grouped.txt", "w") as f:
        for entry in lut_patched:
            f.write(f'{entry[0]} "{entry[1]}" "{entry[2]}" "{entry[3]}"\n')
    
    print("Grouped parcellation and LUT saved.")

#2-END-----------------------------------------------------------------------------------------------------------------------------------
#3-REMAPPING ATLAS PARCELLATION BY LATERALIZING IT FOR LEFT-RIGHT COMPARISON (OPTIONAL)--------------------------------------------

lateralize_parcellation = messagebox.askyesno("Lateralize regions in parcellation?", "Do you want to lateralize regions in the parcellation for left-right comparison ?")
if group_parcellation and lateralize_parcellation :
    
# elif selection == "Continuous Lateralized Parcellation" or selection == "Continuous Grouped Lateralized Parcellation" :
    print("SpasTract Connectomics : Atlas Remapping, Remapping parcellation by lateralizing it for left-right comparison")
    print("You chose to follow the grouping with a left-right lateralization of the parcellation, now processing...")
    
    nii_path = nii_out_path + "_grouped.nii"
    lut_path = os.path.splitext(output_file2)[0] + "_grouped.txt"
    
    if nii_path and lut_path:
        img = nib.load(nii_path)
        # Dict Format : {'acronym' : 'FRP', 'fullname' : 'Frontal pole', 'color' : '268F45'}
        label_dict = load_grouped_lut(lut_path)
        new_data, new_labels, missing = compute_symmetrical_hemisphere(img, label_dict)
    
        # Save new NIfTI
        out_img_path = os.path.splitext(nii_path)[0] + "_lateralizedLeftRight.nii"
        nib.save(nib.Nifti1Image(new_data, img.affine, img.header), out_img_path)
        print(f"Saved split NIfTI to: {out_img_path}")
    
        # Save new LUT
        out_lut_path = os.path.splitext(lut_path)[0] + "_lateralizedLeftRight.txt"
        with open(out_lut_path, "w") as f:
            for idx in sorted(new_labels):
                info = new_labels[idx]
                f.write(f'{idx} "{info["acronym"]}" "{info["fullname"]}" "{info["color"]}"\n')
        print(f"Saved split LUT to: {out_lut_path}")

# elif selection == "Continuous Grouped Lateralized Parcellation" : 
elif lateralize_parcellation and not(group_parcellation) : 
    print("SpasTract Connectomics : Atlas Remapping, Remapping parcellation by lateralizing it for left-right comparison, without grouping")
    print("You chose to not group small regions yet apply a left-right lateralization of the parcellation, now processing...")
    
    nii_path = output_continuous_path
    lut_path = output_file2
    
    if nii_path and lut_path:
        img = nib.load(nii_path)
        # Dict Format : {'acronym' : 'FRP', 'fullname' : 'Frontal pole', 'color' : '268F45'}
        label_dict = load_grouped_lut(lut_path)
        new_data, new_labels, missing = compute_symmetrical_hemisphere(img, label_dict)
    
        # Save new NIfTI
        out_img_path = os.path.splitext(nii_path)[0] + "_lateralizedLeftRight.nii"
        nib.save(nib.Nifti1Image(new_data, img.affine, img.header), out_img_path)
        print(f"Saved split NIfTI to: {out_img_path}")
    
        # Save new LUT
        out_lut_path = os.path.splitext(lut_path)[0] + "_lateralizedLeftRight.txt"
        with open(out_lut_path, "w") as f:
            for idx in sorted(new_labels):
                info = new_labels[idx]
                f.write(f'{idx} "{info["acronym"]}" "{info["fullname"]}" "{info["color"]}"\n')
        print(f"Saved split LUT to: {out_lut_path}")
    
#3-END-----------------------------------------------------------------------------------------------------------------------------------
