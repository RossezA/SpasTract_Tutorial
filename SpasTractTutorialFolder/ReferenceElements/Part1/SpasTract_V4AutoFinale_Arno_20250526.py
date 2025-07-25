# -*- coding: utf-8 -*

    
#MODIFIED SCRIPT BY ARNO ROSSEZ 11/03/25 (Annotations, Modification des dialog boxes, Modification pour les intervalles de l_max)
#0-INITIALISATION & REQUIRED PACKAGES IMPORT-------------------------------------------------------------------------- 

import os
import re
import subprocess
import math
import datetime


def str2bool(v):
  return v.lower() in ("yes", "true", "1","ok","oui")

def removesuffix(s, suffix): # Recoding of removesuffix function available in Python 3.9+ so that Python 3.8 users can have access to it
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s

import fnmatch

#Function to find the parameters files through the Working Directory            
def find_parameters_file(workPath):
    parameters_filepath = None
    for file in os.listdir(workPath):
        if fnmatch.fnmatch(file, "*parameters-tractograms*.txt"):
            parameters_filepath = os.path.join(workPath, file)
            break
    return parameters_filepath


#Function to extract the numeric portion of the subject name and sort by that numerically, while maintaining prefix grouping (F, M, etc.)...
#Caution : only works on list of strings like ['F1006', 'F756', 'M974'] to sort it like : ['F756', 'F1006', 'M974']
def alphanum_key(s):
    match = re.match(r"([A-Za-z]+)(\d+)", s)
    if match:
        prefix, number = match.groups()
        return (prefix, int(number))
    return (s, 0)  # fallback

#Function to extract **from a path**, the numeric portion of the subject name and sort by that numerically, while maintaining prefix grouping (F, M, etc.)...
def extract_subject_id_key(path):
    basename = os.path.basename(path)
    match = re.search(r"([FMfm]\d{3,4})", basename)
    if match:
        return alphanum_key(match.group(1))
    return (basename, 0)  # fallback

#Function to Search Recursively through the Working Directory for Masks,BvecsBvals,DWI NIfTI files
def find_tractofiles(workPath, niiFilePathVector, maskvector, bvecsbvalsvector):
    # Iterate through subdirectories in workPath
    for root, dirs, files in os.walk(workPath):
        # Sort directories naturally (e.g., F756 before F1006 and F1006 before M974)
        # & Remove "TractographyResults" from the list of directories to be processed to avoid descending into it...
        dirs[:] = sorted([d for d in dirs if d != "TractographyResults"], key=extract_subject_id_key)

        # Search for files in the subdirectory        
        for subdir in dirs: #except "TractographyResults"
            subdir_path = os.path.join(root, subdir)

            sub_files = sorted(os.listdir(subdir_path), key=extract_subject_id_key)

            for file in sub_files:
                file_path = os.path.join(subdir_path, file)
                
                # Check for mask file
                if fnmatch.fnmatch(file, "*mask*.nii") or fnmatch.fnmatch(file, "*Mask*.nii"):
                    maskvector.append(file_path)
                
                # Check for nii file with subdir name
                if fnmatch.fnmatch(file, f"{subdir}*.nii"):
                    niiFilePathVector.append(file_path)
                
                # Check for bvecs_bvals file
                if any(keyword in file for keyword in ["bvecs_bvals", "bvecsbvals", "BvecsBvals", "Bvecs_Bvals"]):
                    if file.endswith(".txt"):
                        bvecsbvalsvector.append(file_path)

    return niiFilePathVector, bvecsbvalsvector, maskvector

# # Example usage
# workPath = "/path/to/your/directory"
# niiFilePathVector, maskvector, bvecsbvalsvector = [], [], []
# niiFilePathVector, maskvector, bvecsbvalsvector = find_tractofiles(workPath, niiFilePathVector, maskvector, bvecsbvalsvector)

# # Print the vectors
# print("\nnii file vector:", niiFilePathVector)
# print("mask vector:", maskvector)
# print("bvecs_bvals vector:", bvecsbvalsvector)

#Function to find files from SIFT analysis in the Working Directory to avoid redoing them, to overwrite them however it is possible to override this with triggerOverride            
def find_SIFT_files(workPath, SIFTname, triggerOverride):
    SIFTname_lowercase = SIFTname.lower() #Switching to lowercase to prevent issue with SIFT and sift or Sift...
    SIFT_filepath = None
    for file in os.listdir(workPath):
        if fnmatch.fnmatch(file.lower(), SIFTname_lowercase):
            SIFT_filepath = os.path.join(workPath, file)
            print("SIFT file already made at :", SIFT_filepath)
            if triggerOverride==True :
                print("Overriding skip triggered, redoing SIFT")
                return False
            else : 
                return True
    else :
        return False


import shutil
import glob

def Split4Dto3D_N4BiasCat(input_file, gradient_table, voxelsize ):
    
    """Splits a 4D series MIF file into separate 3D volumes in a temporary directory, 
        averages the b=0 volume(s) together, then apply N4BiasFieldCorrection on said average, outputs the estimated bias field,
        applies the bias field to all volumes (b=0 and diffusion) of the 4D series by dividing them by the bias field, 
        retrieves the result and removes said temp. directory."""
         
    # Handling (only) MIF file (using MRtrix)
    if input_file.endswith('.mif'):
        
        # Step 0 : Retrieving the input filename and a copy of the file in a newly created temporary directory
        input_name = removesuffix(input_file, ".mif")
        #Python 3.9+ line : #input_name = input_file.removesuffix(".mif") # Retrieving the input filename 
        temp_dir = 'TEMPN4Bias_' + input_name
        CWD = os.getcwd()
        N4BiasCatPath = os.path.join(CWD, temp_dir)
        subprocess.run(['mkdir', N4BiasCatPath], stdout=subprocess.PIPE)
        subprocess.run(['cp', input_file, os.path.join(N4BiasCatPath,input_file)])
            #Going into temporary directory (before switching back again to main directory later)
        os.chdir(N4BiasCatPath)
                
        # Step 1: Getting number of volumes in input file 
            # Run mrinfo to capture the number of frames involved in directions and also the number of b0 frames
        mrinfoglobal = subprocess.run(['mrinfo', input_file, '-size'], capture_output=True, text=True)
            # Extract the output and split it into a list of integers to retrieve the number of frames from it
        size_values = list(map(int, mrinfoglobal.stdout.strip().split()))
        num_global = size_values[3]
        
        # Step 2a : Extracting b=0 volumes and Computing average of said b=0 volumes in one command
        # First process: dwiextract
        extract = subprocess.Popen(['dwiextract', input_file, '-', '-bzero'],stdout=subprocess.PIPE)
        # Second process: mrmath
        avg = subprocess.Popen(['mrmath', '-', 'mean', 'mean_bzero.mif', '-axis', '3'],stdin=extract.stdout)
        # Ensure extract doesn't hang
        extract.stdout.close()
        # Wait for the processes to finish
        avg.communicate()
        
        # Step 2b : Converting b=0 average file from MIF format to NIfTI for N4BiasFieldCorrection
        subprocess.run(['mrconvert','mean_bzero.mif', 'mean_bzero.nii'], stdout=subprocess.PIPE)
                
        # Step 2c : Computing N4BiasFieldCorrection function on NIfTI b=0 average and output bias field estimation, converting it then to MIF format
        N4Biasop = subprocess.run("N4BiasFieldCorrection -i mean_bzero.nii -o [mean_bzero_N4.nii, N4biasfield.nii]",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)  # optional: get strings instead of bytes
        print("STDOUT:\n", N4Biasop.stdout)
        print("STDERR:\n", N4Biasop.stderr)
        subprocess.run(['mrconvert','N4biasfield.nii','-vox', voxelsize, 'N4biasfield.mif'], stdout=subprocess.PIPE) 
        
        # Step 3 : Looping through and extract each volume separately in temp directory
        MIFfile_names = []
        for i in range(num_global):
           MIFfile_name = f"volume_{(i):04d}.mif"
                      
           # Extract each volume using mrconvert
           subprocess.run(['mrconvert', input_file, MIFfile_name, '-coord', '3', str(i)], stdout=subprocess.PIPE)
                      
           MIFfile_names.append(MIFfile_name)
         
        # Step 4: Looping through the 3D files extracted to temp directory and perform dividing operation hence applying estimated bias field to all volumes.
        for MIFfile in MIFfile_names:
            MIFprocessedfile = 'N4Bias' + MIFfile
            print(f"Applying estimated N4 bias field to {MIFfile} ...")
            
            subprocess.run(['mrcalc', MIFfile, 'N4biasfield.mif', '-div', MIFprocessedfile], stdout=subprocess.PIPE)

        dwistrayCat = input_name + '_N4cat' + '.mif' # Re-adding .mif extension and specifying that it's a concatenation of N4Bias results 
        dwiN4Catpathout = os.path.join(CWD, dwistrayCat) 
        
        # Step 5 : Concatenating processed volumes back into a proper 4D series in MIF format 
        # Expand the wildcard manually because subprocess doesn't compute well 'N4Bias*.mif'
        MIF_files = glob.glob("N4Bias*.mif")
        MIF_files.sort() #sort the files in lexical order (and here since they are written the same except for number then they get sorted in increasing order)
        
        # Ensure there are enough files
        if len(MIF_files) < 2:
            raise ValueError("Mr.cat requires at least two input files in its bowl :/ ...")
        
        # Run Mr.cat with sorted list of volumes and proper wildcard
        subprocess.run(['mrcat'] + MIF_files + [dwistrayCat], stdout=subprocess.PIPE) # by default concatenation is done following -axis 3 a.k.a the time axis
        # Extracting Mr.cat and retrieving voxel size and gradient table information at the same time to be on the safe side...
        subprocess.run(['mrconvert', dwistrayCat,'-grad', gradient_table, '-vox', voxelsize, dwiN4Catpathout], stdout=subprocess.PIPE)

        # Cleanup: Remove the temporary directory after that mrcat produced dwistrayCat outside of it thanks to dwiN4Catfilepath
        os.chdir(CWD)
        shutil.rmtree(temp_dir)
        print("Temporary directory removed.")
        
        return dwistrayCat
        
    else:
        raise ValueError("Unsupported file format. Use .mif ! ^(;..;)^ !")
    
    return dwistrayCat


def transform2singleshell(lines):
    """
    Transforms the input lines by modifying the shell value in the header, 
    adjusting the command history, and removing the third line (index 2).
    """
    transformed_lines = []

    for i, line in enumerate(lines):
        if i == 2:  # Remove the third line (index 2 in Python)
            continue

        if line.startswith("# Shells:"):
            original_shells = line.strip().split(": ")[1].split(",")
            new_shell = original_shells[-1]  # Keep only the last shell
            transformed_lines.append(f"# Shells: {new_shell}\n")
            continue

        if line.startswith("# command_history:"):
            # Extract and modify the command arguments
            parts = line.split("-shells '")
            if len(parts) > 1:
                cmd_prefix = parts[0]
                shell_lmax_part = parts[1].split("' -lmax '")
                new_shells = shell_lmax_part[0].split(",")[-1]  # Keep only last shell
                new_lmax = shell_lmax_part[1].split("'")[0].split(",")[-1]  # Extract last lmax
                new_command = f"{cmd_prefix}-shells '{new_shells}' -lmax '{new_lmax}'  (version=3.0.4-180-g3af4fe87)\n"
                transformed_lines.append(new_command)
            else:
                transformed_lines.append(line)
            continue

        transformed_lines.append(line)  # Add all other lines unchanged

    return transformed_lines


result = subprocess.run(['reset'], stdout=subprocess.PIPE)
print("\n\n#####################################################################################################")
print("\n\nData processing autoscript for MRtrix3 - Coded by Diego ALVES RODRIGUES DE SOUZA - Modified by Arno ROSSEZ")
print("\n\n#####################################################################################################")

#Initializing the Vectors that will have each
niiFilePathVector, bvecBvalFilePathVector, externalMaskPathVector=[], [], []

#0-END----------------------------------------------------------------------------------------------------------------

#1-ASKING FOR PARAMETERS IN INPUT FILE THEN COLLECTING DATA IN WORKING FOLDER-----------------------------------------

print("\n\n#####################################################################################################")
print("\n\nPlease provide the required parameters, then everything will unfold as planned...")
print("\n\n#####################################################################################################")


#Read file with parameters
workPath =os.getcwd()
inputFilePath = find_parameters_file(workPath)
if not(inputFilePath):
    print("No parameters to feed on... bye-bye...T_T.")
    exit()
    
    
print("Input file path:")
print(inputFilePath)

important = []
keep_phrases = ["QuantityOfSamples:","VoxelDimensions(mm):","NumberOfDirections:","NumberOfSeedsperVoxel:","FractionalAnisotropyCutoff:","StepSize(mm):","Curvature:","PerformN4BiasonDWIVolumes:","ViewPreprocessingResultsduringRun:","SphericalDeconvolutionResponseAlgorithm:","TractographyAlgorithm:","ROIsToBeAnalysed:","PerformROIFiltering:","PerformROIStatistics:","PerformWholeBrainStatistics:","PerformSIFTOnly:","SIFTmode:","SIFT1TerminationOptions:","SIFT1TerminationValues:","WorkPath:","DataPath:","ProvideListOfImages:","ProvideExternalMask:","SameMaskforAllImages:","SameBvecsBvalsforAllImages:","SubFolderComment:","NumberOfThreads:"]

with open(inputFilePath, 'r', encoding='utf-8', errors='replace') as f:
    f = f.readlines()

for line in f:
    for phrase in keep_phrases:
        if phrase in line:
            important.append(line)
            break

print(important)
important=[w.replace('QuantityOfSamples: ', '').replace('VoxelDimensions(mm): ', '').replace('NumberOfDirections: ', '').replace('NumberOfSeedsperVoxel: ', '').replace('FractionalAnisotropyCutoff: ', '').replace('StepSize(mm): ', '').replace('Curvature: ', '').replace('PerformN4BiasonDWIVolumes: ','').replace('ViewPreprocessingResultsduringRun: ','').replace('SphericalDeconvolutionResponseAlgorithm: ','').replace('TractographyAlgorithm: ','').replace('ROIsToBeAnalysed: ','').replace('PerformROIFiltering: ','').replace('PerformROIStatistics: ','').replace('PerformWholeBrainStatistics: ','').replace('PerformSIFTOnly: ','').replace('SIFTmode: ','').replace('SIFT1TerminationOptions: ','').replace('SIFT1TerminationValues: ','').replace('WorkPath: ','').replace('DataPath: ','').replace('ProvideListOfImages: ','').replace('ProvideExternalMask: ','').replace('SameMaskforAllImages: ','').replace('SameBvecsBvalsforAllImages: ','').replace('SubFolderComment: ','').replace('NumberOfThreads: ','').replace('\n', '') for w in important]
important=[w.replace('QuantityOfSamples:', '').replace('VoxelDimensions(mm):', '').replace('NumberOfDirections:', '').replace('NumberOfSeedsperVoxel:', '').replace('FractionalAnisotropyCutoff:', '').replace('StepSize(mm):', '').replace('Curvature:', '').replace('PerformN4BiasonDWIVolumes:','').replace('ViewPreprocessingResultsduringRun:','').replace('SphericalDeconvolutionResponseAlgorithm:','').replace('TractographyAlgorithm:','').replace('ROIsToBeAnalysed:','').replace('PerformROIFiltering:','').replace('PerformROIStatistics:','').replace('PerformWholeBrainStatistics:','').replace('PerformSIFTOnly:','').replace('SIFTmode:','').replace('SIFT1TerminationOptions:','').replace('SIFT1TerminationValues:','').replace('WorkPath:','').replace('DataPath:','').replace('ProvideListOfImages:','').replace('ProvideExternalMask:','').replace('SameMaskforAllImages:','').replace('SameBvecsBvalsforAllImages:','').replace('SubFolderComment:','').replace('NumberOfThreads:','').replace('\n', '') for w in important]
important=[w.replace('\n', '').replace(':', '').replace(' ', '') for w in important]
print(important)

####################################

quantityofSamples= int(important[0])
voxelSizeString = important[1]
numberOfDirections = int(important[2])
numberofSeedsVector= list(map(int,important[3].split(",")))
cutoffVector = list(map(float,important[4].split(",")))
stepSizeVector = list(map(float,important[5].split(",")))
curvature = float(important[6])

performN4BiasonVolumes = str2bool(important[7])
viewPreprocessingResults = str2bool(important[8])
sphericalRFAlgo = list(map(str,important[9].split(",")))
tractAlgo = list(map(str,important[10].split(",")))

whichROI = important[11].split(",")
performROIFiltering = str2bool(important[12])
performROIStatistics = str2bool(important[13])

performWholeBrainStatistics = str2bool(important[14])

performSIFTOnly = str2bool(important[15])
#performAdditionalSIFT = str2bool(important[14]) > Old parameter forsaken as of 19.05.2025 since it was redundant with siftMode
siftMode = int(important[16])
siftTerminationOptionsArray = important[17].split(",")
if important[18]=='SIFT1TerminationValues' or important[18]=='':
    siftTerminationValuesArray = []
else :
    siftTerminationValuesArray = list(map(float,important[18].split(",")))

workPath = important[19]
dataPath = important[20]
provideListOfImages = str2bool(important[21])
externalMask = str2bool(important[22])
sameExternalMask = str2bool(important[23])
sameBvecsBvals = str2bool(important[24])
folderComment = important[25]

nThreads = int(important[26])

####################################

print("The following parameters have been retrieved from the parameter file :")
print(quantityofSamples) #Number of samples to be analyzed
print(voxelSizeString) #MRI Voxels dimensions (i.e '0.08 0.08 0.08' in mm) 
print(numberOfDirections) #Number of directions done in DWI
print(numberofSeedsVector) #Number of Seeds per voxels for tractography
print(cutoffVector) #List of different used cutoff thresholds of FA between 0 and 1 (usually 0.23,0.24,...,0.28)
print(stepSizeVector) #Parameter used to indicate the step size used in the tractography functions (to define the size of a step, usually 1/10 of a voxel length)
#TODO:change curvature to angle to avoid having to compute angle ? (Check with Hervé first)
print(curvature) #Parameter used to calculate the angle by doing angle=math.round(math.degrees(2*asin(step/(2*curvature))),4), From MRtrix webpage : "set the maximum angle in degrees between successive steps (defaults: 60 for deterministic algorithms; 15 for iFOD1 / nulldist1; 45 for iFOD2 / nulldist2)"


print(performN4BiasonVolumes) #Boolean parameter to ask whether or not N4Bias should be done on all DWI Volumes rather than just for the mask (yes or no)
print(viewPreprocessingResults) #Boolean parameter to ask to view the results during the data handling (yes or no) 
if sphericalRFAlgo[0]=='' or sphericalRFAlgo[0]=='SphericalDeconvolutionResponseAlgorithm' : #when a blank was left for SphericalDeconvolutionResponseAlgorithm parameter (either with space (empty space) or not (then  TractographyAlgorithm is read)) : defaulting to iFOD2
    print('no algorithm detected, defaulting to tournier algorithm for spherical deconvolution response function estimation') #Algorithm used by default for dwi2response
    TriggerDefaultTournier = True
    sphericalRFAlgo[0] = 'tournier'
else :
    print(sphericalRFAlgo) #Algorithm(s) chosen in the parameters for spherical deconvolution response function estimation (i.e tournier, dhollander, ...)
    TriggerDefaultTournier = False
    
if tractAlgo[0]=='' or tractAlgo[0]=='TractographyAlgorithm' : #when a blank was left for TractographyAlgorithm parameter (either with space (empty space) or not (then  TractographyAlgorithm is read)) : defaulting to iFOD2
    print('no algorithm detected, defaulting to iFOD2 (Probabilistic)') #Algorithm used by default in tckgen
    TriggerDefaultiFOD2 = True
    tractAlgo[0] = 'iFOD2'
else :
    print(tractAlgo) #Algorithm(s) chosen in the parameters (probabilistic or determinist or both)
    TriggerDefaultiFOD2 = False
    

print(whichROI) #List of strings giving which ROIs will be used in ROI Filtering
print(performROIFiltering) #Boolean parameter to ask if there will be ROI Filtering done on pre-existing data (yes or no)
print(performROIStatistics) #Boolean parameter to ask if ROI statistics will be computed on pre-existing data (yes or no)
print(performWholeBrainStatistics) #Boolean parameter to ask if whole-brain statistics will be computed on pre-existing data (yes or no)
print(siftMode) #Trinary logical parameter : [0 means no SIFT], [1 means classical SIFT will be done ], [2 means SIFT2 will be done] 
#print(performAdditionalSIFT) #Boolean parameter to ask if SIFT will be applied on generated data (tracts, ROIs, ...) (yes or no) > Old parameter forsaken as of 19.05.2025 since it was redundant with siftMode
print(performSIFTOnly) #Boolean parameter to ask if there will be SIFT done on pre-existing data (yes or no)
if siftTerminationOptionsArray=='None' or siftTerminationOptionsArray[0]=='' or siftTerminationOptionsArray[0]=='SiftTerminationOptions':
    print("No SIFT chosen for this run... \nIf this was a mistake, you can interrupt this script and change parameter SiftTerminationOptions to : 'standard' or 'mu' or 'ratio' or several of them in one comma-separated list...")
    SIFTAborted = True
else : 
    print(siftTerminationOptionsArray) #List of strings or single strings indicating which type of option will be used in SIFT to control when it terminates filtering, options are 'standard', 'mu' for term_mu , 'ratio' for term_ratio
    print(siftTerminationValuesArray) #List of values given that will be used to define the term_mu or term_ratio SIFT options.


print(workPath) #Path to the working directory a.k.a the directory where the code will work from (to select file for example)
print(dataPath) #Path to the directory containing pre-existing data (a.k.a whole-brain tractograms) to be processed further (with SIFT, ROIFiltering, Statistics...)
print(provideListOfImages) #Boolean (yes or no) parameter if there is or not a provided list of images with their respective paths
print(externalMask) #Boolean (yes or no) parameter if there is prepared masks for the data processing 
print(sameExternalMask) #Boolean (yes or no) parameter if there is or not the same mask used for every image
print(sameBvecsBvals) #Boolean (yes or no) parameter if there is or not the same BvecsBvals used for every image
print(folderComment) #String parameter used to name results folders/sub-folders


print(nThreads) #Number of computer threads that will be dedicated to data processing

#1-END-----------------------------------------------------------------------------------------------------------------------------

#2-RETRIEVING AUTOMATICALLY IMAGES, BVECS-BVALS & MASKS-------------------------------------------------------------------------------------------
if ((not(performROIFiltering)) and (not(performSIFTOnly)) and not(performROIStatistics) and not(performWholeBrainStatistics)):
    
    niiFilePathVector, bvecBvalFilePathVector, externalMaskPathVector= find_tractofiles(workPath,niiFilePathVector, bvecBvalFilePathVector, externalMaskPathVector)

    print("\nNii Path vector:")
    print(len(niiFilePathVector))
    print(niiFilePathVector)
    print("\nBvec bval path vector:")
    print(len(bvecBvalFilePathVector))
    print(bvecBvalFilePathVector)
    print("\nExternal mask path vector:")
    print(len(externalMaskPathVector))
    print(externalMaskPathVector)

#2-END------------------------------------------------------------------------------------------------------------------------------------

#3-CREATING OUTPUT FOLDER & COPYING IN USED PARAMETERS-----------------------------------------------------------------------------------

    print("Making an output folder for the tractograpy...")
    outputPath = workPath + '/TractographyResults'
    if not(os.path.isdir(outputPath)) :
        proc = subprocess.run(['mkdir', outputPath], stdout=subprocess.PIPE)

    now = datetime.datetime.now()

    newOutputFolderName = (now.strftime("%Y%m%d_%H%M%S_")) + 'MRtrix-output_' + folderComment
    newOutputPath = outputPath + '/' + newOutputFolderName
    os.chdir(outputPath)

    proc = subprocess.run(['rm','-R', newOutputFolderName], stdout=subprocess.PIPE)
    proc = subprocess.run(['mkdir', newOutputFolderName], stdout=subprocess.PIPE)
    os.chdir(newOutputPath)

    # CHANGE - 2020.03.12 - Copy parameter file to the output directory
    inputFileFolder, inputFileName = os.path.split(inputFilePath)
    newInputFileName = (now.strftime("%Y%m%d_%H%M%S_")) + inputFileName
    result = subprocess.run(['cp', inputFilePath, os.getcwd() + '/' + newInputFileName], stdout=subprocess.PIPE)
    # CHANGE - END ----------------------------------------------------


    print("\n\n#####################################################################################################")
    print("\nInformation about paths:")
    print("\n######################################################################################################")

    print("\nOutput directories for each sample:")

    print("\nNifTI file(s) to be processed:")
    print(niiFilePathVector)

    # CHANGE - 2021.05.27 - Keeping track of the image(s) used in a .txt file -------------------------------------------
    with open("inputImages.txt", "w") as txt_file:
        for line in niiFilePathVector:
            txt_file.write(line + "\n")
    # CHANGE - END ------------------------------------------------------------------------------------------------------


    print("\nBvec and Bval file(s) to be considered:")
    print(bvecBvalFilePathVector)

    # CHANGE - 2021.05.27 - Keeping track of the bvecs bvals used in a .txt file ----------------------------------------
    with open("inputBvecsBvals.txt", "w") as txt_file:
        for line in bvecBvalFilePathVector:
            txt_file.write(line + "\n")
    # CHANGE - END ------------------------------------------------------------------------------------------------------


    if (externalMask):
        print("\nExternal mask file(s) to be considered:")
        print(externalMaskPathVector)

        # CHANGE - 2021.05.27 - Keeping track of the mask(s) used in a .txt file ----------------------------------------
        with open("inputMasks.txt", "w") as txt_file:
            for line in externalMaskPathVector:
                txt_file.write(line + "\n")
        # CHANGE - END --------------------------------------------------------------------------------------------------



    print("\nOld output directory")
    print(outputPath)

    print("\nNew output directory:")
    print(newOutputPath)
    print("#####################################################################################################")

#3-END------------------------------------------------------------------------------------------------------------------------------------

#4-CREATING, IN THE OUPUT FOLDER, OUTPUT SUB-FOLDERS FOR EACH SAMPLE----------------------------------------------------------------------

    individualFolderNameVector=[]

    pathIndividualFolderNameVector = []

    quantityFlag=quantityofSamples
    while (quantityFlag > 0):
        individualFolderNameVector.append('Sample-' + str(quantityofSamples-quantityFlag+1))
        proc = subprocess.run(['mkdir', individualFolderNameVector[quantityofSamples-quantityFlag]], stdout=subprocess.PIPE)

        # CHANGE - 2021.05.27 - Prepping for sub-folder making ----------------------------------------------------
        pathIndividualFolderNameVector.append(newOutputPath + '/' + individualFolderNameVector[quantityofSamples-quantityFlag])
        # CHANGE - END --------------------------------------------------------------------------------------------

        quantityFlag -= 1

    print("\n\n#####################################################################################################")
    print("\n\nOutput directories for each sample:")
    print("\n######################################################################################################")
    proc = subprocess.run(['ls'])
    print("#####################################################################################################")

    # CHANGE - 2021.05.27 -  Keeping track of the sub-folders created in a .txt file ------------------------------------
    with open("outputDataFolders.txt", "w") as txt_file:
        for line in pathIndividualFolderNameVector:
            txt_file.write(line + "\n")
    # CHANGE - END -----------------------------------------------------------------------------------------------------

#4-END-------------------------------------------------------------------------------------------------------------------------

#5-DATA HANDLING START-------------------------------------------------------------------------------------------------------
#-5a)-DATA INFO-----------------------------------------------------------------------------------------------------------------
    print("\n\n#####################################################################################################")
    print("\nInformation about NifTIs to be processed:")
    print("\n######################################################################################################")
    quantityFlag=quantityofSamples
    while (quantityFlag > 0):

        print("\n----------------------------------------------------------------------------------------------")
        print(individualFolderNameVector[quantityofSamples-quantityFlag])
        print("----------------------------------------------------------------------------------------------")

        result = subprocess.run(['mrinfo', niiFilePathVector[quantityofSamples-quantityFlag]], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
        quantityFlag -= 1

        print("------------------------------------------------------------------------------------------------------\n")

    print("#####################################################################################################")

#-5b)-DATA CONVERSION NIfTI TO MIF-----------------------------------------------------------------------------------------------
    print("\n\n#####################################################################################################")
    print("\nConversion to .mif files:")
    print("\n######################################################################################################")

    quantityFlag=quantityofSamples

    while (quantityFlag > 0):

        os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples-quantityFlag])

        print("\n------------------------------------------------------------------------------------------------------")
        print(individualFolderNameVector[quantityofSamples-quantityFlag])
        print("----------------------------------------------------------------------------------------------")

        result = subprocess.run(['mrconvert', niiFilePathVector[quantityofSamples-quantityFlag],'-grad', bvecBvalFilePathVector[quantityofSamples-quantityFlag], '-vox', voxelSizeString, 'dwi_raw.mif'], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))

        result = subprocess.run(['mrinfo','dwi_raw.mif'], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))

        quantityFlag -= 1
        print("------------------------------------------------------------------------------------------------------\n")

    print("#####################################################################################################")


#-5c)-DATA PROCESSING (DENOISING, UNRINGING, LOADING MASK(S))--------------------------------------------------------------------------------------------------
    print("\n\n#####################################################################################################")
    print("\nData processing:")
    print("\n######################################################################################################")

    quantityFlag=quantityofSamples

    while (quantityFlag > 0):

        os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples-quantityFlag])

        print("\n------------------------------------------------------------------------------------------------------")
        print(individualFolderNameVector[quantityofSamples-quantityFlag])
        print("------------------------------------------------------------------------------------------------------")

	#Denoising
        result = subprocess.run(['dwidenoise','dwi_raw.mif','dwi_den.mif','-noise','noise.mif'],stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
	
	#Obtaining residual noise
        result = subprocess.run(['mrcalc','dwi_raw.mif','dwi_den.mif','-subtract','residual.mif'],stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))

        if viewPreprocessingResults:
            result = subprocess.run(['mrview','noise.mif'],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

            result = subprocess.run(['mrview','residual.mif'],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

	#Unringing using deGibbs3D
        result = subprocess.run(['deGibbs3D','dwi_den.mif','dwi_den_unr.mif'],stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
    #Assigning unringed volumes as the volumes ready for later processing
        PreProcName = 'dwi_den_unr.mif'
	
	#Obtaining residual noise after unringing
        result = subprocess.run(['mrcalc','dwi_den.mif','dwi_den_unr.mif','-subtract','residualUnringed.mif'],stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))


        if viewPreprocessingResults:
            result = subprocess.run(['mrview','dwi_den_unr.mif'],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

            result = subprocess.run(['mrview','residualUnringed.mif'],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

    #Applying N4BiasFieldCorrection on unringed DWI series 
        if performN4BiasonVolumes:
            N4CatName = Split4Dto3D_N4BiasCat('dwi_den_unr.mif', bvecBvalFilePathVector[quantityofSamples-quantityFlag],voxelSizeString) #equivalent to this but without the crash at the end : result = subprocess.run(['dwibiascorrect ants', 'dwi_den_unr.mif', 'N4Bias_den_unr.mif'],stdout=subprocess.PIPE)
            #Re-assigning newly-made N4Bias volumes as the volumes ready for later processing
            PreProcName = N4CatName

	#Converting external mask to .mif format
        if (externalMask):
            print('External Mask used :')
            print(externalMaskPathVector[quantityofSamples - quantityFlag])
            result = subprocess.run(['cp', externalMaskPathVector[quantityofSamples - quantityFlag], './'],stdout=subprocess.PIPE)
            externalMaskActualPath, externalMaskName = os.path.split(externalMaskPathVector[quantityofSamples - quantityFlag])
            MaskName = 'externalmask.mif'
            result = subprocess.run(['mrconvert', externalMaskName , MaskName],stdout=subprocess.PIPE)
	
	#or (if there's no external mask) Creating from scratch a mask using MRtrix built-in function and ANTS N4BiasFieldCorrection (in a in-house function)
        else: 
            print ('No External Mask provided > Creating Mask from scratch')
            if not(performN4BiasonVolumes):
                N4CatName = Split4Dto3D_N4BiasCat('dwi_den_unr.mif', bvecBvalFilePathVector[quantityofSamples-quantityFlag],voxelSizeString) #equivalent to this but without the crash at the end : result = subprocess.run(['dwibiascorrect ants', 'dwi_den_unr.mif', 'N4Bias_den_unr.mif'],stdout=subprocess.PIPE)
            MaskName = 'mask_den_unr_N4Bias.mif'
            result = subprocess.run(['dwi2mask', N4CatName, MaskName],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))
            
            



#-5d)-MD(Mean Diffusivity)==ADC ((Mean) Apparent Diffusion Coefficient), RD(Radial Diffusivity), AD(Axial Diffusivity) & FA(Fractional Anisotropy) MAPS MAKING 
#+ METRICS (cl, cs, cp, eigenvalue/eigenvector)--------------------------------------------------------------------------------------------------------------- 
#[code snippet added on 30.01.2021]------------------------- -------------------------------------------------------------------------------------------------
        result = subprocess.run(['dwi2tensor', '-mask', MaskName, PreProcName, 'dwi_tensor.mif'], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))

        result = subprocess.run(['tensor2metric', 'dwi_tensor.mif', '-mask', MaskName, '-adc', 'ADC_map.nii', '-fa', 'FA_map.nii', '-rd', 'RD_map.nii', '-ad', 'AD_map.nii', '-cl', 'CL_map.nii', '-cp', 'CP_map.nii', '-cs', 'CS_map.nii', '-value', 'EigenValue_map.nii', '-vector', 'EigenVector_map.nii'], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
	
	#
        result = subprocess.run(['tensor2metric', 'dwi_tensor.mif', '-mask', MaskName, '-adc', 'ADC_map.mif', '-fa', 'FA_map.mif','-rd', 'RD_map.mif', '-ad', 'AD_map.mif', '-cl', 'CL_map.mif', '-cp', 'CP_map.mif', '-cs', 'CS_map.mif','-value', 'EigenValue_map.mif', '-vector', 'EigenVector_map.mif'], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))

        if viewPreprocessingResults:
            result = subprocess.run(['mrview', PreProcName, '-overlay.load', MaskName],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

            result = subprocess.run(['mrview', 'FA_map.mif'],stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

            result = subprocess.run(['mrview', 'ADC_map.mif'], stdout=subprocess.PIPE)
            print(result.stdout.decode('utf-8'))

#[END code snippet added]-------------------------------------------------------------------------------------------------------------------------------------

#-5e)-SPHERICAL DECONVOLUTION RESPONSE FUNCTION ESTIMATION ALGORITHM(S) & FOD (Fibre Orientation Distributions)
        for RFalgorithms in sphericalRFAlgo : 
            if RFalgorithms in ['tournier']:
            #TOURNIER ALGORITHM
            # CHANGE - 10.03.2025 - Arno ROSSEZ - Allowing intervals for lmax values with respect to the numberOfDirections--------------
                if numberOfDirections>=6 and numberOfDirections<15:
                    result = subprocess.run(['dwi2response', 'tournier', '-mask',MaskName,'-lmax','2',PreProcName,'wm_response_tournier.txt','-voxel','voxelstournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=15 and numberOfDirections<28:
                    result = subprocess.run(['dwi2response', 'tournier', '-mask',MaskName,'-lmax','4',PreProcName,'wm_response_tournier.txt','-voxel','voxelstournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=28 and numberOfDirections<66:
                    result = subprocess.run(['dwi2response', 'tournier', '-mask', MaskName, '-lmax', '6', PreProcName,'wm_response_tournier.txt', '-voxel', 'voxelstournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=66 and numberOfDirections<91:
                    result = subprocess.run(['dwi2response', 'tournier', '-mask', MaskName, '-lmax', '8', PreProcName,'wm_response_tournier.txt', '-voxel', 'voxelstournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                else:
                    break
                
                if viewPreprocessingResults:
                    result = subprocess.run(['mrview',PreProcName,'-overlay.load','voxelstournier.mif'],stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                    
                if numberOfDirections>=6 and numberOfDirections<15:
                    result = subprocess.run(['dwi2fod', 'csd', PreProcName,'-mask',MaskName,'-lmax','2','wm_response_tournier.txt','wmfod_tournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=15 and numberOfDirections<28:
                    result = subprocess.run(['dwi2fod', 'csd', PreProcName,'-mask',MaskName,'-lmax','4','wm_response_tournier.txt','wmfod_tournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=28 and numberOfDirections<66:
                    result = subprocess.run(['dwi2fod', 'csd', PreProcName, '-mask',MaskName,'-lmax','6','wm_response_tournier.txt','wmfod_tournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                elif numberOfDirections>=66 and numberOfDirections<91:
                    result = subprocess.run(['dwi2fod', 'csd', PreProcName, '-mask',MaskName,'-lmax','8','wm_response_tournier.txt','wmfod_tournier.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                else:
                    break
                # CHANGE - END --------------------------------------------------------------------------------------------------------------
        
                #Generating image out of non-normalized FOD
                result = subprocess.run(['mrconvert','-coord','3','0','wmfod_tournier.mif','vf.mif'], stdout=subprocess.PIPE)
                print(result.stdout.decode('utf-8'))
        
                if viewPreprocessingResults:
                    result = subprocess.run(['mrview','vf.mif','-odf.load_sh','wmfod_tournier.mif'],stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
        
                #Normalizing FOD generated
                result = subprocess.run(['mtnormalise','wmfod_tournier.mif','wmfod_norm_tournier.mif','-mask',MaskName], stdout=subprocess.PIPE)
                print(result.stdout.decode('utf-8'))
                #Generating image out of normalized FOD
                result = subprocess.run(['mrconvert','-coord','3','0','wmfod_norm_tournier.mif','vf_norm.mif'], stdout=subprocess.PIPE)
                print(result.stdout.decode('utf-8'))
        
                if viewPreprocessingResults:
                    result = subprocess.run(['mrview','vf_norm.mif','-odf.load_sh','wmfod_norm_tournier.mif'],stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
            
            
            elif RFalgorithms in ['dhollander']: 
                #DHOLLANDER ALGORITHM (Recommended for Single-Shell)
                #dwi2response dhollander input out_sfwm out_gm out_csf [ options ]
                    # input: Input DWI dataset
                    # out_sfwm: Output single-fibre WM response function text file
                    # out_gm: Output GM response function text file
                    # out_csf: Output CSF response function text file
    
                # CHANGE - 10.03.2025 - Arno ROSSEZ - Allowing intervals for lmax values with respect to the numberOfDirections--------------
                    if numberOfDirections>=6 and numberOfDirections<15:
                        result = subprocess.run(['dwi2response', 'dhollander', PreProcName, 'sfwm_response_dhollander.txt', 'gm_response_dhollander.txt', 'csf_response_dhollander.txt', '-mask',MaskName, '-lmax','0,2', '-voxel','voxelsdhollander.mif'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=15 and numberOfDirections<28:
                        result = subprocess.run(['dwi2response', 'dhollander', PreProcName, 'sfwm_response_dhollander.txt', 'gm_response_dhollander.txt', 'csf_response_dhollander.txt', '-mask',MaskName, '-lmax','0,4', '-voxel','voxelsdhollander.mif'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=28 and numberOfDirections<66:
                        result = subprocess.run(['dwi2response', 'dhollander', PreProcName, 'sfwm_response_dhollander.txt', 'gm_response_dhollander.txt', 'csf_response_dhollander.txt', '-mask',MaskName, '-lmax','0,6', '-voxel','voxelsdhollander.mif'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=66 and numberOfDirections<91:
                        result = subprocess.run(['dwi2response', 'dhollander', PreProcName, 'sfwm_response_dhollander.txt', 'gm_response_dhollander.txt', 'csf_response_dhollander.txt', '-mask',MaskName, '-lmax','0,8', '-voxel','voxelsdhollander.mif'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    else:
                        break
    
                    #Directly modifying sfwm_response.txt into a vector
                    # Read the input file
                    sfwm_file_path = os.getcwd()+"/sfwm_response_dhollander.txt"  
                    with open(sfwm_file_path, "r", encoding="utf-8", errors='replace') as file:
                        file_contents = file.readlines()
                    # Apply the transformation
                    transformed_contents = transform2singleshell(file_contents)
                    # Write the transformed lines to a new file
                    with open(sfwm_file_path, "w", encoding="utf-8") as file:
                        file.writelines(transformed_contents)
    
                    print(f"File '{sfwm_file_path}' has been overwritten to keep only single-shell parameters.")
                    
                    if viewPreprocessingResults:
                        result = subprocess.run(['mrview',PreProcName,'-overlay.load','voxelsdhollander.mif'],stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    
                    if numberOfDirections>=6 and numberOfDirections<15:
                        result = subprocess.run(['dwi2fod', 'csd', PreProcName,'sfwm_response_dhollander.txt','wmfod_dhollander.mif', '-mask',MaskName, '-lmax','2'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=15 and numberOfDirections<28:
                        result = subprocess.run(['dwi2fod', 'csd', PreProcName,'sfwm_response_dhollander.txt','wmfod_dhollander.mif', '-mask',MaskName,'-lmax','4'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=28 and numberOfDirections<66:
                        result = subprocess.run(['dwi2fod', 'csd', PreProcName,'sfwm_response_dhollander.txt','wmfod_dhollander.mif', '-mask',MaskName,'-lmax','6'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    elif numberOfDirections>=66 and numberOfDirections<91:
                        result = subprocess.run(['dwi2fod', 'csd', PreProcName,'sfwm_response_dhollander.txt','wmfod_dhollander.mif', '-mask',MaskName,'-lmax','8'], stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
                    else:
                        break
                    
                    #Generating image out of non-normalized FOD
                    result = subprocess.run(['mrconvert','-coord','3','0','wmfod_dhollander.mif','vf_dhollander.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
            
                    if viewPreprocessingResults:
                        result = subprocess.run(['mrview','vf_dhollander.mif','-odf.load_sh','wmfod_dhollander.mif'],stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))
            
                    #Normalizing FOD generated
                    result = subprocess.run(['mtnormalise','wmfod_dhollander.mif','wmfod_norm_dhollander.mif','-mask',MaskName], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
                    #Generating image out of normalized FOD
                    result = subprocess.run(['mrconvert','-coord','3','0','wmfod_norm_dhollander.mif','vf_dhollander_norm.mif'], stdout=subprocess.PIPE)
                    print(result.stdout.decode('utf-8'))
            
                    if viewPreprocessingResults:
                        result = subprocess.run(['mrview','vf_dhollander_norm.mif','-odf.load_sh','wmfod_norm_dhollander.mif'],stdout=subprocess.PIPE)
                        print(result.stdout.decode('utf-8'))

#5-END----------------------------------------------------------------------------------------------------------------------------------------------------

#6-TRACTOGRAPHY START-------------------------------------------------------------------------------------------------------------------------------------
        
        siftTerminii=[]
        if siftMode==1 :
            for siftType in siftTerminationOptionsArray:
                for siftValue in siftTerminationValuesArray:
    
                    if (siftType == 'standard'):
                        siftTerminus = '_SIFT-standard'
                        siftTerminii.append(siftTerminus)
    
                    elif (siftType == 'mu'):
                        siftTerminus = '_SIFT-term-mu-' + str(siftValue)
                        siftTerminii.append(siftTerminus)
    
                    elif (siftType == 'ratio'):
                        siftTerminus = '_SIFT-term-ratio-' + str(siftValue)
                        siftTerminii.append(siftTerminus)
        
        for numberOfSeeds in numberofSeedsVector:
            for cutoff in cutoffVector:
                for stepSize in stepSizeVector:

                    strCutoff = ('%.4f' % cutoff)
                    angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                    
                    for RFalgorithms in sphericalRFAlgo :
                        if tractAlgo :                        
                            for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFOD1,...)
                            
                                if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2']):
                                    result = subprocess.run(
                                        ['tckgen', '-algorithm' , algorithms, '-seed_random_per_voxel', MaskName, str(numberOfSeeds),
                                         '-cutoff', strCutoff, 'wmfod_norm_' + RFalgorithms + '.mif', '-minlength', '0.1', '-maxlength', '50', '-step',
                                         str(stepSize), '-angle', str(angle), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                    print(result.stdout.decode('utf-8'))
                                    
                                    if siftMode==1 :
                                        print('\nSIFT mode 1 chosen, tcksift in progress...')
                                        #tcksift [ options ]  in_tracks in_fod out_tracks
                                        for siftTerminus in siftTerminii :
                                                
                                                # SIFT_name = 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                # + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck' 
                                                # BooleanSIFT = find_SIFT_files(os.getcwd(), SIFT_name , True) #Overriding because the aim is to produce SIFT here, but not overriding could also be used for SIFT checking
                    
                                                if ('standard' in siftTerminus) :
                                                    print('Performing standard sift...')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                                                    
                    
                                                elif ('mu' in siftTerminus) :
                                                    print('Performing mu sift...')
                                                    #Retrieving mu value
                                                    SIFTmuStringvalue = siftTerminus.split('-')[-1]  
                                   
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_mu', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                    
                    
                                                elif ('ratio' in siftTerminus):
                                                    print('Performing ratio sift...')
                                                    #Retrieving ratio value
                                                    SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                     
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_ratio', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")  
                                    
                                ##########################
                                    elif siftMode==2 :
                                        print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                        #tcksift2 [ options ]  in_tracks in_fod out_weights
                                        SIFT2 = subprocess.run(
                                            ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                        print(SIFT2.stdout.decode('utf-8'))
                                        print("SIFT2 Success !")
                                        
                                        
            
                                elif (algorithms in ['Tensor_Det','Tensor_Prob']):
                                    result = subprocess.run(
                                        ['tckgen', '-algorithm', algorithms, '-seed_random_per_voxel',
                                         MaskName, str(numberOfSeeds),
                                         '-cutoff', strCutoff, PreProcName, '-minlength', '0.1', '-maxlength', '50',
                                         '-step',
                                         str(stepSize), '-angle', str(angle), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads)],
                                        stdout=subprocess.PIPE)
                                    print(result.stdout.decode('utf-8'))
                                    
                                    if siftMode==1 :
                                        print('\nSIFT mode 1 chosen, tcksift in progress...')
                                        #tcksift [ options ]  in_tracks in_fod out_tracks
                                        for siftTerminus in siftTerminii :
                    
                                                if ('standard' in siftTerminus) :
                                                    print('Performing standard sift...')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                                                    
                    
                                                elif ('mu' in siftTerminus) :
                                                    print('Performing mu sift...')
                                                    #Retrieving mu value
                                                    SIFTmuStringvalue = siftTerminus.split('-')[-1]  
                                   
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_mu', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                    
                    
                                                elif ('ratio' in siftTerminus):
                                                    print('Performing ratio sift...')
                                                    #Retrieving ratio value
                                                    SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                     
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_ratio', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                                    
                                ##########################
                                    elif siftMode==2 :
                                        print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                        #tcksift2 [ options ]  in_tracks in_fod out_weights
                                        SIFT2 = subprocess.run(
                                            ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                        print(SIFT2.stdout.decode('utf-8'))
                                        print("SIFT2 Success !")            
                                    
            
                                elif (algorithms in ['FACT']):
                                    result = subprocess.run(
                                        ['tckgen', '-algorithm', algorithms, '-seed_random_per_voxel',
                                         MaskName, str(numberOfSeeds),
                                         '-cutoff', strCutoff, 'EigenVector_map.mif', '-minlength', '0.1', '-maxlength', '50',
                                         '-step',
                                         str(stepSize), '-angle', str(angle), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads)],
                                        stdout=subprocess.PIPE)
                                    print(result.stdout.decode('utf-8'))
                                    
                                    if siftMode==1 :
                                        print('\nSIFT mode 1 chosen, tcksift in progress...')
                                        #tcksift [ options ]  in_tracks in_fod out_tracks
                                        for siftTerminus in siftTerminii :
                    
                                                if ('standard' in siftTerminus) :
                                                    print('Performing standard sift...')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                                                    
                    
                                                elif ('mu' in siftTerminus) :
                                                    print('Performing mu sift...')
                                                    #Retrieving mu value
                                                    SIFTmuStringvalue = siftTerminus.split('-')[-1]  
                                   
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_mu', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                    
                    
                                                elif ('ratio' in siftTerminus):
                                                    print('Performing ratio sift...')
                                                    #Retrieving ratio value
                                                    SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                     
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_ratio', str(siftValue), 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    print("SIFT Succeeded with :", siftTerminus, " termination option.")
                                    
                                #########################
                                    elif siftMode==2 :
                                        print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                        #tcksift2 [ options ]  in_tracks in_fod out_weights
                                        SIFT2 = subprocess.run(
                                            ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                        print(SIFT2.stdout.decode('utf-8'))
                                        print("SIFT2 Success !")  
                                
                                # # Useless ADDITION because there is a -mask option used previously for tckgen... 
                                # # ADDITION - 05.05.2025 + 15.05.2025 - Arno ROSSEZ - Cropping tractography using mask to refine the results and remove extra fibers (-mask in tckgen is used to plant seeds in voxels contained in the mask, not to crop, hence leading to some spurious fibers sometimes in the end result, especially with agarose casing)--------------            
                                # tractoName = 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                # cropName = 'Crop_' + tractoName 
                                # tractoNamenoRF= 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                # cropNamenoRF = 'Crop_' + tractoNamenoRF 
                                # tractoNamenoRFnoAlgo= 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                # cropNamenoRFnoAlgo = 'Crop_' + tractoNamenoRFnoAlgo 
                                
                                # cropping = subprocess.run(['tckedit', tractoName, cropName,'-mask',MaskName, '-nthreads', str(nThreads), '-force'],stdout=subprocess.PIPE)
                                # if cropping.returncode==0 : #Process success
                                #     print(cropping.stdout.decode('utf-8'))
                                # else :
                                #     cropping = subprocess.run(['tckedit', tractoNamenoRF, cropNamenoRF,'-mask',MaskName, '-nthreads', str(nThreads), '-force'],stdout=subprocess.PIPE)
                                #     if cropping.returncode==0 : #Process success
                                #         print(cropping.stdout.decode('utf-8'))
                                #     else :
                                #         cropping = subprocess.run(['tckedit', tractoNamenoRFnoAlgo, cropNamenoRFnoAlgo,'-mask',MaskName, '-nthreads', str(nThreads), '-force'],stdout=subprocess.PIPE)
                                #         if cropping.returncode==0 : #Process success
                                #             print(cropping.stdout.decode('utf-8'))
                                        
                                # # END ADDITION --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
                                
                                # ADDITION - 05.05.2025 + 15.05.2025 + 16.05.2025 - Arno ROSSEZ - Adding use of tckmap function to generate 'high-resolution images' thanks to tracts data being used as a form of contrast--------------------------------            
                                tractoName = 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                tractoNamenoRF= 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                tractoNamenoRFnoAlgo= 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck'
                                tckweightstxtName = 'SIFT2weights' + removesuffix(tractoName, ".tck") + ".txt"
                                tckweightstxtNamenoRF = 'SIFT2weights' + removesuffix(tractoNamenoRF, ".tck") + ".txt"
                                tckweightstxtNamenoRFnoAlgo = 'SIFT2weights' + removesuffix(tractoNamenoRFnoAlgo, ".tck") + ".txt"
                                tckmapName = 'Map_' + removesuffix(tractoName, ".tck") + ".mif"
                                tckmapNamenoRF = 'Map_' + removesuffix(tractoNamenoRF, ".tck") + ".mif"
                                tckmapNamenoRFnoAlgo = 'Map_' + removesuffix(tractoNamenoRFnoAlgo, ".tck") + ".mif"
                                                                
                                mapping = subprocess.run(['tckmap', tractoName, tckmapName,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                if mapping.returncode==0 : #Process success
                                    print(mapping.stdout.decode('utf-8'))    
                                    if siftMode==2 :
                                        mappingSIFT2 = subprocess.run(['tckmap', tractoName, tckmapName,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-tck_weights_in',tckweightstxtName, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                        if mappingSIFT2.returncode==0 : #Process success
                                            print(mappingSIFT2.stdout.decode('utf-8'))
                                    
                                else :
                                    mapping = subprocess.run(['tckmap', tractoNamenoRF, tckmapNamenoRF,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                    if mapping.returncode==0 : #Process success
                                        print(mapping.stdout.decode('utf-8'))
                                        if siftMode==2 :
                                            mappingSIFT2 = subprocess.run(['tckmap', tractoNamenoRF, tckmapNamenoRF,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-tck_weights_in',tckweightstxtNamenoRF, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                            if mappingSIFT2.returncode==0 : #Process success
                                                print(mappingSIFT2.stdout.decode('utf-8'))
                                    else :
                                        mapping = subprocess.run(['tckmap', tractoNamenoRFnoAlgo, tckmapNamenoRFnoAlgo,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                        if mapping.returncode==0 : #Process success
                                            print(mapping.stdout.decode('utf-8'))
                                            if siftMode==2 :
                                                mappingSIFT2 = subprocess.run(['tckmap', tractoNamenoRFnoAlgo, tckmapNamenoRFnoAlgo,'-dec','-template',PreProcName,'-vox',voxelSizeString, '-tck_weights_in',tckweightstxtNamenoRFnoAlgo, '-nthreads', str(nThreads),'-force'],stdout=subprocess.PIPE)
                                                if mappingSIFT2.returncode==0 : #Process success
                                                    print(mappingSIFT2.stdout.decode('utf-8'))
                                # END ADDITION ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   
        quantityFlag -= 1

        print("------------------------------------------------------------------------------------------------------\n")
    print("\nTractography Complete !")
    print("#####################################################################################################")

#6-END-----------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------
###The following sections of this script are designed for analysis and processing of already existing tractography data (namely SiftOnly 7, StatsOnly 9 and ROI Filtering 8)
#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

#7-PERFORMING SIFT ONLY -------------------------------------------------------------------------------------------------------------------------------------------
if (performSIFTOnly and not(SIFTAborted)):
    print("\n#######")
    print("\nPerforming SIFT only on pre-existing data...")
    print("\n#######")
    newOutputPath = dataPath
    
    nowSIFT = datetime.datetime.now()
    SIFTOutputPath = newOutputPath + "/" + (nowSIFT.strftime("%Y%m%d_%H%M%S_")) + "SiftOnlyReprocessing_" + folderComment + "/"
    subprocess.run(['mkdir', SIFTOutputPath], stdout=subprocess.PIPE)
    
    # Copy parameter file to the SIFT subdirectory
    inputFileFolder, inputFileName = os.path.split(inputFilePath)
    SIFTFileName = (nowSIFT.strftime("%Y%m%d_%H%M%S_")) + inputFileName
    result = subprocess.run(['cp', inputFilePath, SIFTOutputPath + SIFTFileName], stdout=subprocess.PIPE)

    # Creating output folder vector for each sample
    # Creating output folders for each sample
    individualSIFTFolderNameVector=[]
    individualFolderNameVector=[]
    quantityFlag=quantityofSamples
    while (quantityFlag > 0):
        individualFolderNameVector.append('Sample-' + str(quantityofSamples-quantityFlag+1)) 
        individualSIFTFolderNameVector.append(SIFTOutputPath + 'SIFTSample-' + str(quantityofSamples-quantityFlag+1))
        proc = subprocess.run(['mkdir', individualSIFTFolderNameVector[quantityofSamples-quantityFlag]], stdout=subprocess.PIPE)
        quantityFlag -= 1
    print("\n\n#####################################################################################################")
    print("\n\nOutput folders for each sample's SIFT:")
    print("\n######################################################################################################")
    os.chdir(SIFTOutputPath)
    # Showing in terminal the different folders created using ls command
    proc = subprocess.run(['ls'])
    
    
    siftTerminii=[]
    if siftMode==1 :    
        for siftType in siftTerminationOptionsArray:
            for siftValue in siftTerminationValuesArray:
    
                if (siftType == 'standard'):
                    siftTerminus = '_SIFT-standard'
                    siftTerminii.append(siftTerminus)
    
                elif (siftType == 'mu'):
                    siftTerminus = '_SIFT-term-mu-' + str(siftValue)
                    siftTerminii.append(siftTerminus)
    
                elif (siftType == 'ratio'):
                    siftTerminus = '_SIFT-term-ratio-' + str(siftValue)
                    siftTerminii.append(siftTerminus)

    os.chdir(newOutputPath)

    quantityFlag = quantityofSamples
    
    while (quantityFlag > 0):
        
        os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples - quantityFlag])

        print("\n------------------------------------------------------------------------------------------------------")
        print(individualFolderNameVector[quantityofSamples - quantityFlag])
        print("------------------------------------------------------------------------------------------------------")

        for numberOfSeeds in numberofSeedsVector:
            for cutoff in cutoffVector:
                for stepSize in stepSizeVector:
                    
                    angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                    strCutoff = ('%.4f' % cutoff)
                                                 
                           
                    for RFalgorithms in sphericalRFAlgo : 
                        for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFod1,...)
                                    
                            if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2','Tensor_Det','Tensor_Prob','FACT']): #Making sure we don't take other cases than the algorithm names by accident
                                
                                if siftMode==1 :
                                    print('\nSIFT mode 1 chosen, tcksift in progress...')
                                    #tcksift [ options ]  in_tracks in_fod out_tracks
                                    for siftTerminus in siftTerminii :
                                    
                                        if ('standard' in siftTerminus):
                                            print('Performing standard sift...\r\n')
                                            SIFT = subprocess.run(
                                                ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            if SIFT.returncode==0 : #Process success
                                                print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                print(SIFT.stdout.decode('utf-8'))
                                                
                                            else :
                                                SIFT = subprocess.run(
                                                    ['tcksift', 'tracts_' + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                if SIFT.returncode==0 : #Process success
                                                    print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    
                                                else : 
                                                    SIFT = subprocess.run(['tcksift', 'tracts_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    
                                                        
                                        elif ('mu' in siftTerminus):
                                            
                                            #Retrieving mu value
                                            SIFTmuStringvalue = siftTerminus.split('-')[-1]
                                                                                
                                            print('Performing mu sift...\r\n')
                                            SIFT = subprocess.run(
                                                ['tcksift', '-term_mu', SIFTmuStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            if SIFT.returncode==0 : #Process success
                                                print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                print(SIFT.stdout.decode('utf-8'))
                                                
                                            else :
                                                 SIFT = subprocess.run(
                                                     ['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + algorithms + '_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                 if SIFT.returncode==0 : #Process success
                                                     print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                     print(SIFT.stdout.decode('utf-8'))
                                                     
                                                 else : 
                                                     SIFT = subprocess.run(['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                     print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                     print(SIFT.stdout.decode('utf-8'))
                                                                
                                        elif ('ratio' in siftTerminus):
                                            
                                            #Retrieving ratio value
                                            SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                            
                                            print('Performing ratio sift...\r\n')
                                            SIFT = subprocess.run(
                                                ['tcksift', '-term_ratio', SIFTratioStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            if SIFT.returncode==0 : #Process success
                                                print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                print(SIFT.stdout.decode('utf-8'))
                                                
                                            else :
                                                 SIFT = subprocess.run(
                                                     ['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + algorithms + '_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                 if SIFT.returncode==0 : #Process success
                                                     print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                     print(SIFT.stdout.decode('utf-8'))
                                                     
                                                 else : 
                                                     SIFT = subprocess.run(['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                         numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                         stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                     print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                     print(SIFT.stdout.decode('utf-8'))
                                                            
                            ##########################
                                elif siftMode==2 :
                                    print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                    #tcksift2 [ options ]  in_tracks in_fod out_weights
                                    
                                    SIFT2 = subprocess.run(
                                        ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                    if SIFT2.returncode==0 : #Process success
                                        print("New name convention detected, SIFT2 Success !")
                                        print(SIFT2.stdout.decode('utf-8'))
                                        copy = subprocess.run(['cp','SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                        stepSize) + '_angle-' + str(angle) + '.txt', 
                                                        individualSIFTFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                        # print(copy.stdout.decode('utf-8'))
                                    else :
                                        SIFT2 = subprocess.run(
                                            ['tcksift2', 'tracts_' + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                        if SIFT2.returncode==0 : #Process success
                                            print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT2 Success nonetheless !")
                                            print(SIFT2.stdout.decode('utf-8'))
                                            copy = subprocess.run(['cp','SIFT2weightstracts_' + algorithms + '_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.txt', 
                                                            individualSIFTFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                            # print(copy.stdout.decode('utf-8'))
                                        else : 
                                            SIFT2 = subprocess.run(
                                                ['tcksift2', 'tracts_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            print("Older name convention detected, SIFT2 Success nonetheless !")
                                            print(SIFT2.stdout.decode('utf-8'))
                                            copy = subprocess.run(['cp','SIFT2weightstracts_' + str(
                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                            stepSize) + '_angle-' + str(angle) + '.txt', 
                                                            individualSIFTFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                            # print(copy.stdout.decode('utf-8'))
                                                                                                                    
                        
        quantityFlag -= 1
    os.chdir(newOutputPath)
    
    print("------------------------------------------------------------------------------------------------------\n")
    print("\nSIFT Complete !")
    print("#####################################################################################################")
#7-END-------------------------------------------------------------------------------------------------------------------------------------------------------------

#8-ROI FILTERING AND COMPUTING ROI STATISTICS (on data with and without SIFT) -------------------------------------------------------------------------------------------

# #- ADD USELESS - 2025.03-04 - working countdown function that can be used to give time for user input...--------
# #Importing necessary modules for countdown
# import time
# import threading
# #Definition of countdown function
# def countdown(t, timeout_event): 
#     while t: 
#         if timeout_event.is_set():
#             return  # Exit early if timeout was already triggered elsewhere
#         mins, secs = divmod(t, 60) 
#         timer = '{:02d}:{:02d}'.format(mins, secs) 
#         print(timer, end="\r") 
#         time.sleep(1) 
#         t -= 1
#     print('\nCountdown over!')
#     timeout_event.set()  # Notify the main thread that time is up

# timeout_event = threading.Event()
# # Start countdown in a separate thread
# timer_thread = threading.Thread(target=countdown, args=(60, timeout_event))  # Set your timer duration here in seconds, here for 1 minute
# timer_thread.start()
# #- END ADD USELESS - 2025.03-04 ----------------------------------------------------------------------------------


# #- ADD USELESS - 2025.03-04 - working countdown function that can be used to give time for user input...--------
#  # Check result
# if timeout_event.is_set() :
#     # Timeout happened, no response
#     print("\n\n#####################################################################################################")
#     print("\n No user response: Timeout, Skipping ROI Filtering step")
#     print("\n######################################################################################################")
# #- END ADD USELESS - 2025.03-04 ----------------------------------------------------------------------------------

     
if (performROIFiltering or performROIStatistics):
    print("\n\n####################################################################################################")
    print("\n user chose : Proceeding with ROI Filtering on existing data")
    print("\n######################################################################################################")
    
    print("\nFor ROI full names matching their abbreviations please refer to the table below :")
    print("\nAnterior Commissure == ac")
    print("\nCorpus Callosum == cc")
    print("\nFasciculus Retroflexus == fr")
    print("\nright Fasciculus Retroflexus == frR")
    print("\nleft Fasciculus Retroflexus == frL")
    print("\nright Internal Capsule and Cerebral Peduncle == icR-cpR")
    print("\nleft Internal Capsule and Cerebral Peduncle == icL-cpL")
    print("\nright MammilloThalamic tract == mtR")
    print("\nleft MammilloThalamic tract == mtL")
    print("\nOptical tract == opt")
    print("\nFornix == f")
    print("\nright Fornix == fR")
    print("\nright Postcommissural Fornix == pfR")
    print("\nleft Postcommissural Fornix == pfL")
    print("\nPyramidal tract == py")
    print("\nright Stria Terminalis == stR")
    print("\nleft Stria Terminalis == stL")
    print("\nright Stria Medularis == smR")
    print("\nleft Stria Medularis == smL")
    print("\nVentral Hippocampal Commissure == vhc")
    print("\nOptic Chiasm == och")
    print("\n...")
    print("\n#####################################################################################################")

    #Choose folder with data where all ROIs are already made
    newOutputPath = dataPath        
    
    nowROI = datetime.datetime.now()
    ROIOutputPath = newOutputPath + "/" + (nowROI.strftime("%Y%m%d_%H%M%S_")) + "ROIFiltering_" + folderComment + "/"
    subprocess.run(['mkdir', ROIOutputPath], stdout=subprocess.PIPE)
    
    # Copy parameter file to the ROI subdirectory
    inputFileFolder, inputFileName = os.path.split(inputFilePath)
    ROIFileName = (nowROI.strftime("%Y%m%d_%H%M%S_")) + inputFileName
    result = subprocess.run(['cp', inputFilePath, ROIOutputPath + ROIFileName], stdout=subprocess.PIPE)
    
    
    #Create output folder vector for each sample
    #Create output folders for each sample
    individualFolderNameVector=[]
    individualROIFolderNameVector=[]
    quantityFlag=quantityofSamples
    while (quantityFlag > 0):
        individualFolderNameVector.append('Sample-' + str(quantityofSamples-quantityFlag+1))
        individualROIFolderNameVector.append(ROIOutputPath + 'ROIResultsSample-' + str(quantityofSamples-quantityFlag+1))
        proc = subprocess.run(['mkdir', individualROIFolderNameVector[quantityofSamples-quantityFlag]], stdout=subprocess.PIPE)
        quantityFlag -= 1
    print("\n\n####################################################################################################")
    print("\n\nOutput folders for each sample's ROI analysis:")
    print("\n######################################################################################################")
    os.chdir(newOutputPath)
    # Showing in terminal the different folders created using ls command
    proc = subprocess.run(['ls'])

    print("\n#####################################################################################################")
    print("\nStatistics and ROI Filtering on data...")
    print("\n#####################################################################################################")
    
   
    
    for regionSelected in whichROI:

        os.chdir(newOutputPath)
        f = open("Results" + "_" + regionSelected + ".txt","w+")
        fTitle = "Compilation of " +regionSelected+ " ROI statistics\r\n"
        f.write(fTitle)


        quantityFlag = quantityofSamples
        
        while (quantityFlag > 0):
            f.write("\n%s\r\n" % (individualFolderNameVector[quantityofSamples-quantityFlag]))
            f.write("Region RFAlgo TractoAlgo Cutoff Step-size StreamlinesCount Mean Median Std Min Max\r\n")

            os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples-quantityFlag])

            print("\n------------------------------------------------------------------------------------------------------")
            print(individualFolderNameVector[quantityofSamples-quantityFlag])
            print("------------------------------------------------------------------------------------------------------")


            for numberOfSeeds in numberofSeedsVector:
                for cutoff in cutoffVector:
                    for stepSize in stepSizeVector:
                        
                        angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                        strCutoff = ('%.4f' % cutoff)
                        
                        
                        for RFalgorithms in sphericalRFAlgo : 
                            for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFOD1,...)
                            
                                if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2','Tensor_Det','Tensor_Prob','FACT']): #Making sure we don't take other cases than the algorithm names by accident
                                    
                                
                                    ####
                                    if (regionSelected in ['frR','frL'] and performROIFiltering): #"ROI1&not(ROI3)" scheme
                                        if regionSelected=='frR':
                                            print('Analysis of the right Fasciculus Retroflexus')
                                        elif regionSelected=='frL':
                                            print('Analysis of the left Fasciculus Retroflexus')
                                        
                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                 + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                        if result.returncode==0 : #Process success
                                            print("ROI of :", regionSelected, "done.")
                                            subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                            + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                            print(result.stdout.decode('utf-8'))
                                        
                                        else : 
                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                     'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                     + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                     + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                     + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                            if result.returncode==0 : #Process success
                                                print("ROI of :", regionSelected, "done.")
                                                subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                print(result.stdout.decode('utf-8'))
                                            
                                            else :
                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                         'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                         + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + str(numberOfSeeds) 
                                                                         + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                         + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                print(result.stdout.decode('utf-8'))
                                            
                                            
                                    #####
                                    elif (regionSelected in ['och'] and performROIFiltering): 
                                    #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                        if regionSelected=='och':
                                            print('Analysis of the Optic Chiasm')
                                        #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                        #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                        ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                           ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                           ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                        #Removed from cmd line : '-exclude','ROI-'+regionSelected+ochCombinations[k][2]+'.mif', '-exclude','ROI-'+regionSelected+ochCombinations[k][3]+'.mif',
                                        
                                        for k in range(len(ochCombinations)):
                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                     'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                     + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                     + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                     + str(stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
    
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            if result.returncode==0 : #Process success
                                                print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                print(result.stdout.decode('utf-8'))
                                            
                                            else : 
                                                
                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                         'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                         + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                         + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                         + str(stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                                         
                                                if result.returncode==0 : #Process success
                                                    print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                    subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                    + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                    individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    print(result.stdout.decode('utf-8'))
                                                
                                                else :
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                             'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                             + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                             + '_'  + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                             + str(stepSize) + '_angle-' + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                    
                                                    print("Older name convention detected, ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                    subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                    + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                    individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    print(result.stdout.decode('utf-8'))
                                    
                                    #####
                                    elif (regionSelected not in ['frR', 'frL', 'och'] and performROIFiltering): #"ROI1&ROI2&not(ROI3)" default scheme
                                        if regionSelected=='ac':
                                            print('Analysis of the Anterior Commissure')
                                        elif regionSelected=='cc':
                                            print('Analysis of the Corpus Callosum')
                                        elif regionSelected=='f':
                                            print('Analysis of the Fornix')
                                        elif regionSelected=='fr':
                                            print('Analysis of the Fasciculus Retroflexus')
                                        elif regionSelected=='icR-cpR':
                                            print('Analysis of the right Internal Capsule and Cerebral Peduncle')
                                        elif regionSelected=='icL-cpL':
                                            print('Analysis of the left Internal Capsule and Cerebral Peduncle')
                                        elif regionSelected=='mtR':
                                            print('Analysis of the right MammilloThalamic tract')
                                        elif regionSelected=='mtL':
                                            print('Analysis of the left MammilloThalamic tract')
                                        elif regionSelected=='opt':
                                            print('Analysis of the Optical tract')
                                        elif regionSelected=='pfR':
                                            print('Analysis of the right Postcommissural Fornix')
                                        elif regionSelected=='pfL':
                                            print('Analysis of the left Postcommissural Fornix')
                                        elif regionSelected=='py':
                                            print('Analysis of the Pyramidal tract')
                                        elif regionSelected=='smR':
                                            print('Analysis of the right Stria Medularis')
                                        elif regionSelected=='smL':
                                            print('Analysis of the left Stria Medularis')
                                        elif regionSelected=='stR':
                                            print('Analysis of the right Stria Terminalis')
                                        elif regionSelected=='stL':
                                            print('Analysis of the left Stria Terminalis')
                                        elif regionSelected=='vhc':
                                            print('Analysis of the Ventral Hippocampal Commissure')
                                        elif regionSelected=='fR':
                                            print('Analysis of the right Fornix')
                                            
                                        #TODO: [OPTIONAL] define function able to find ROI names based on the string 'ROI', the number in and the regionSelected... 
                                        #hence the subprocess lines would be simplified with just ROI1, ROI2, ROI3not... but drawback would be that one loses the ability to check directly the name of the ROI used + might lead to issues with files having ROI in their names or else...
                                            
                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                 + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                        if result.returncode==0 : #Process success
                                            print("ROI of :", regionSelected, "done.")
                                            subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                            + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                            print(result.stdout.decode('utf-8'))
                                        
                                        else : 
                                            
                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                     'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                     + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                     + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                     + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                            if result.returncode==0 : #Process success
                                                print("ROI of :", regionSelected, "done.")
                                                subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                print(result.stdout.decode('utf-8'))
                                            
                                            else :
                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                         'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                         + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',regionSelected + '_' + str(numberOfSeeds) 
                                                                         + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                         + str(angle) + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                print(result.stdout.decode('utf-8'))
                                            
                                    ############
                                    if performROIStatistics :
                                        if (regionSelected in ['och']): 
                                            #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                            #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                            #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                            ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                               ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                               ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                            for k in range(len(ochCombinations)):
                                                # Computing Statistics of current Decussation in the loop, getting : count, mean, median, std, min, max 
                                                
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                        regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                if statsROI.returncode==0 : #Process success
                                                    statsROI_results = statsROI.stdout.decode('utf-8')
                                                
                                                else :
                                                    statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                            regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                            '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                    if statsROI.returncode==0 : #Process success
                                                        statsROI_results = statsROI.stdout.decode('utf-8')
                                                    
                                                    else : 
                                                        statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                                '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                        statsROI_results = statsROI.stdout.decode('utf-8')
                                                    
                                                # Writing all the computed statistics into the results text file 'f'
                                                notStatsProperties = regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + " " + RFalgorithms + " " + algorithms + " " + strCutoff + " " + str(stepSize) + " "
                                                f.write(notStatsProperties)
                                                f.write(statsROI_results + '\r\n')
                                                
                                                # + Creating a histogram of fibers count as a function of fibers length
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                histogramROI = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                if histogramROI.returncode==0 : #Process success
                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                          + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                          + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                          + str(angle) + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                    
                                                    subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
            
                                                else :
                                                    histogramROI = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                    '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if histogramROI.returncode==0 : #Process success
                                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                              + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                              + str(angle) + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                        
                                                        subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                    
                                                    else :                                        
                                                        histogramROI = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                        '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        if histogramROI.returncode==0 : #Process success
                                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                  + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                  + str(angle) + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                            
                                                            subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                            + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                
                                        else :  
                                            # Computing Statistics of current ROI in the loop, getting : count, mean, median, std, min, max 
                                            
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                    regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                    '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                            if statsROI.returncode==0 : #Process success
                                                statsROI_results = statsROI.stdout.decode('utf-8')
                                            
                                            else :
                                                statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                        regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                if statsROI.returncode==0 : #Process success
                                                    statsROI_results = statsROI.stdout.decode('utf-8')
                                                
                                                else : 
                                                    statsROI = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                            regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                            '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                    statsROI_results = statsROI.stdout.decode('utf-8')
                                                
                                            # Writing all the computed statistics into the results text file 'f'
                                            notStatsProperties = regionSelected + " " + RFalgorithms + " " + algorithms + " " + strCutoff + " " + str(stepSize) + " "
                                            f.write(notStatsProperties)
                                            f.write(statsROI_results + '\r\n')
                                            
                                            # + Creating a histogram of fibers count as a function of fibers length
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                            histogramROI = subprocess.run(['tckstats', regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                            '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            if histogramROI.returncode==0 : #Process success
                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + '.csv', "for :", regionSelected, " ROI.")
                                                
                                                subprocess.run(['cp', 'histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
        
                                            else :
                                                histogramROI = subprocess.run(['tckstats', regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                if histogramROI.returncode==0 : #Process success
                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + '.csv', "for :", regionSelected, " ROI.")
                                                    
                                                    subprocess.run(['cp', 'histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                
                                                else :                                        
                                                    histogramROI = subprocess.run(['tckstats', regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                    '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if histogramROI.returncode==0 : #Process success
                                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + '.csv', "for :", regionSelected, " ROI.")
                                                        
                                                        subprocess.run(['cp', 'histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                             
                 
            quantityFlag -= 1
        
        f.close()
        copy = subprocess.run(['cp', newOutputPath + "/Results_" + regionSelected + ".txt", ROIOutputPath], stdout=subprocess.PIPE)
        print("------------------------------------------------------------------------------------------------------\n")
        
    
    
    #Checking if ROI Analysis is wanted with SIFT too, if so then proceeding...
    if siftMode!=0 :     

        print("#####################################################################################################")
        print("Statistics and ROI Filtering done on data without SIFT, switching to Statistics and ROI Filtering on data with SIFT...")
        print("#####################################################################################################")
        
        
        #Here the goal is to get the names of different SIFT done on tractograms so that we can keep this information in the SIFT ROI names 
        #& Create a txt file to compile ROI statistics
        for regionSelected in whichROI:
            os.chdir(newOutputPath)
            siftTerminii=[]
            if siftMode==1 :
                fSIFT = open("Results_SIFT_" + regionSelected + ".txt", "w+")
                fSIFTTitle = "Compilation of " +regionSelected+ " ROI statistics after SIFT\r\n"
                fSIFT.write(fSIFTTitle)
                for siftType in siftTerminationOptionsArray:
                    for siftValue in siftTerminationValuesArray:
            
                        if (siftType == 'standard'):
                            siftTerminus = '_SIFT-standard'
                            siftTerminii.append(siftTerminus)
            
                        elif (siftType == 'mu'):
                            siftTerminus = '_SIFT-term-mu-' + str(siftValue)
                            siftTerminii.append(siftTerminus)
            
                        elif (siftType == 'ratio'):
                            siftTerminus = '_SIFT-term-ratio-' + str(siftValue)
                            siftTerminii.append(siftTerminus)
            elif siftMode==2 : 
                fSIFT = open("Results_SIFT2_" + regionSelected + ".txt", "w+")
                fSIFTTitle = "Compilation of " +regionSelected+ " ROI statistics after SIFT2\r\n"
                fSIFT.write(fSIFTTitle)
        
            quantityFlag = quantityofSamples
        
            while (quantityFlag > 0):
                fSIFT.write("\n%s\r\n" % (individualFolderNameVector[quantityofSamples - quantityFlag]))
                fSIFT.write("Region RFAlgo TractoAlgo SIFT Cutoff Step-size StreamlinesCount Mean Median Std Min Max\r\n")
        
                os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples - quantityFlag])
        
                print("\n------------------------------------------------------------------------------------------------------")
                print(individualFolderNameVector[quantityofSamples - quantityFlag])
                print("------------------------------------------------------------------------------------------------------")
    
    
                for numberOfSeeds in numberofSeedsVector:
                    for cutoff in cutoffVector:
                        for stepSize in stepSizeVector:
                            
                            angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                            strCutoff = ('%.4f' % cutoff)
                            
                            
                            for RFalgorithms in sphericalRFAlgo : 
                                for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFod1,...)
                                            
                                    if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2','Tensor_Det','Tensor_Prob','FACT']): #Making sure we don't take other cases than the algorithm names by accident
     
                                        if siftMode==1 :
                                            print('\nSIFT mode 1 chosen, tcksift in progress...')
                                            #tcksift [ options ]  in_tracks in_fod out_tracks
                                            for siftTerminus in siftTerminii :
                        
                                                #Searching if SIFT has already been done in the subfolder, in which case there is no need to redo SIFT...
                                                SIFT_nameRFalgo = 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'  
                                                BooleanSIFT_RFalgo = find_SIFT_files(os.getcwd(), SIFT_nameRFalgo , False) 
                                                SIFT_nameAlgo = 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'
                                                BooleanSIFT_Algo = find_SIFT_files(os.getcwd(), SIFT_nameAlgo , False) 
                                                SIFT_name = 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'
                                                BooleanSIFT = find_SIFT_files(os.getcwd(), SIFT_name , False)
                                                
                                                if (BooleanSIFT==False and BooleanSIFT_Algo==False and BooleanSIFT_RFalgo==False) :
                                                    triggerSIFT=True #Nothing was found so triggering SIFT 
                                                else : 
                                                    triggerSIFT=False #SIFT was found so no need to do it again...
                                                    
                                                                     
                                                if ('standard' in siftTerminus) and triggerSIFT==True :
                                                    print('Performing standard sift...\r\n')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if SIFT.returncode==0 : #Process success
                                                        print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                        print(SIFT.stdout.decode('utf-8'))
                                                        
                                                    else :
                                                        SIFT = subprocess.run(
                                                            ['tcksift', 'tracts_' + algorithms + '_' + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        if SIFT.returncode==0 : #Process success
                                                            print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                            print(SIFT.stdout.decode('utf-8'))
                                                            
                                                        else : 
                                                            SIFT = subprocess.run(['tcksift', 'tracts_' + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                            print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                            print(SIFT.stdout.decode('utf-8'))
                                                            
                                                elif ('mu' in siftTerminus) and triggerSIFT==True:
                                                    
                                                    #Retrieving mu value
                                                    SIFTmuStringvalue = siftTerminus.split('-')[-1]
                                                                                        
                                                    print('Performing mu sift...\r\n')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_mu', SIFTmuStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if SIFT.returncode==0 : #Process success
                                                        print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                        print(SIFT.stdout.decode('utf-8'))
                                                        
                                                    else :
                                                         SIFT = subprocess.run(
                                                             ['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + algorithms + '_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                         if SIFT.returncode==0 : #Process success
                                                             print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                             print(SIFT.stdout.decode('utf-8'))
                                                             
                                                         else : 
                                                             SIFT = subprocess.run(['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                             print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                             print(SIFT.stdout.decode('utf-8'))
                                                             
                                                elif ('ratio' in siftTerminus) and triggerSIFT==True:
                                                    
                                                    #Retrieving ratio value
                                                    SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                                    
                                                    print('Performing ratio sift...\r\n')
                                                    SIFT = subprocess.run(
                                                        ['tcksift', '-term_ratio', SIFTratioStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if SIFT.returncode==0 : #Process success
                                                        print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                        print(SIFT.stdout.decode('utf-8'))
                                                        
                                                    else :
                                                         SIFT = subprocess.run(
                                                             ['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + algorithms + '_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                         if SIFT.returncode==0 : #Process success
                                                             print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                             print(SIFT.stdout.decode('utf-8'))
                                                             
                                                         else : 
                                                             SIFT = subprocess.run(['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                                 numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                 stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                             print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                             print(SIFT.stdout.decode('utf-8'))
                                                             
                                                ####
                                                if (regionSelected in ['frR','frL'] and performROIFiltering): #"ROI1&not(ROI3)" scheme
                                                    if regionSelected=='frR':
                                                        print('Analysis of the right Fasciculus Retroflexus after SIFT')
                                                    elif regionSelected=='frL':
                                                        print('Analysis of the left Fasciculus Retroflexus after SIFT')
                                                    
                                                        
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                             'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                             + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                             + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                             + str(angle) + siftTerminus +'.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if result.returncode==0 : #Process success
                                                        print("ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))

                                                    else : 
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                                 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                 + str(angle) + siftTerminus +'.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                        if result.returncode==0 : #Process success
                                                            print("ROI of :", regionSelected, "done.")
                                                            subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                        
                                                        else : 
                                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                                     'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                     + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck',regionSelected + '_' + str(numberOfSeeds) 
                                                                                     + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                     + str(angle) + siftTerminus +'.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                            print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                            subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                            
                                                
                                                ####
                                                elif (regionSelected in ['och'] and performROIFiltering): 
                                                #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                                    if regionSelected=='och':
                                                        print('Analysis of the Optic Chiasm')
                                                    #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                                    #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                                    ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                                       ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                                       ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                                    #Removed from cmd line : '-exclude','ROI-'+regionSelected+ochCombinations[k][2]+'.mif', '-exclude','ROI-'+regionSelected+ochCombinations[k][3]+'.mif',
                                                    
                                                    for k in range(len(ochCombinations)):
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                                 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                 + str(stepSize) + '_angle-' + str(angle) +siftTerminus + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                                 + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                 + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                
                                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                        if result.returncode==0 : #Process success
                                                            print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                            subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                        
                                                        else : 
                                                            
                                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                                     'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                     + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                                     + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                     + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                                                     
                                                            if result.returncode==0 : #Process success
                                                                print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                                subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                                print(result.stdout.decode('utf-8'))
                                                            
                                                            else :
                                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                                         'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                         + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                                         + '_'  + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                         + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                                
                                                                print("Older name convention detected, ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                                subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                                individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                                print(result.stdout.decode('utf-8'))     
                                                
                                                #####
                                                elif (regionSelected not in ['frR', 'frL', 'och'] and performROIFiltering): #"ROI1&ROI2&not(ROI3)" default scheme
                                                    if regionSelected=='ac':
                                                        print('Analysis of the Anterior Commissure after SIFT')
                                                    elif regionSelected=='cc':
                                                        print('Analysis of the Corpus Callosum after SIFT')
                                                    elif regionSelected=='f':
                                                        print('Analysis of the Fornix after SIFT')
                                                    elif regionSelected=='fr':
                                                        print('Analysis of the Fasciculus Retroflexus after SIFT')
                                                    elif regionSelected=='icR-cpR':
                                                        print('Analysis of the right Internal Capsule and Cerebral Peduncle after SIFT')
                                                    elif regionSelected=='icL-cpL':
                                                        print('Analysis of the left Internal Capsule and Cerebral Peduncle after SIFT')
                                                    elif regionSelected=='mtR':
                                                        print('Analysis of the right MammilloThalamic tract after SIFT')
                                                    elif regionSelected=='mtL':
                                                        print('Analysis of the left MammilloThalamic tract after SIFT')
                                                    elif regionSelected=='opt':
                                                        print('Analysis of the Optical tract after SIFT')
                                                    elif regionSelected=='pfR':
                                                        print('Analysis of the right Postcommissural Fornix after SIFT')
                                                    elif regionSelected=='pfL':
                                                        print('Analysis of the left Postcommissural Fornix after SIFT')
                                                    elif regionSelected=='py':
                                                        print('Analysis of the Pyramidal tract after SIFT')
                                                    elif regionSelected=='smR':
                                                        print('Analysis of the right Stria Medularis after SIFT')
                                                    elif regionSelected=='smL':
                                                        print('Analysis of the left Stria Medularis after SIFT')
                                                    elif regionSelected=='stR':
                                                        print('Analysis of the right Stria Terminalis after SIFT')
                                                    elif regionSelected=='stL':
                                                        print('Analysis of the left Stria Terminalis after SIFT')
                                                    elif regionSelected=='vhc':
                                                        print('Analysis of the Ventral Hippocampal Commissure after SIFT')
                                                        
                                                        
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                             'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                             + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                             + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                             + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if result.returncode==0 : #Process success
                                                        print("ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                                    
                                                    else :
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                                 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus +'.tck',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                 + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                        if result.returncode==0 : #Process success
                                                            print("ROI of :", regionSelected, "done.")
                                                            subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                        
                                                        else :
                                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                                     'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                     + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',regionSelected + '_' + str(numberOfSeeds) 
                                                                                     + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                     + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                            print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                            subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                                         
                                                ########
                                                if performROIStatistics :
                                                    if (regionSelected in ['och']): 
                                                        #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                                        #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                                        #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                                        ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                                           ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                                           ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                                        for k in range(len(ochCombinations)):
                                                            # Computing Statistics of current Decussation in the loop, getting : count, mean, median, std, min, max 
                                                            
                                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                            statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                    regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                    '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                            if statsROISIFT.returncode==0 : #Process success
                                                                statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                                
                                                            else :
                                                                statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                        regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                if statsROISIFT.returncode==0 : #Process success
                                                                    statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                                
                                                                else :
                                                                    statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                            regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                            '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                    statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                            
                                                            # Writing all the computed statistics into the results text file 'fSIFT'
                                                            notStatsProperties = regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + " " + RFalgorithms + " " + algorithms + " " + siftTerminus.lstrip('_') + " " + strCutoff + " " + str(stepSize)+ " "
                                                            fSIFT.write(notStatsProperties)
                                                            fSIFT.write(statsROISIFT_results + '\r\n')
                                                            
                                                            # + Creating a histogram of fibers count as a function of fibers length
                                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                            histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                            '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                            + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                            if histogramROISIFT.returncode==0 : #Process success
                                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + siftTerminus + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                                
                                                                subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                                                                    
                                                            else :
                                                                histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                                if histogramROISIFT.returncode==0 : #Process success
                                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + siftTerminus + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                                    
                                                                    subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                                
                                                                else : 
                                                                    histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                    '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                    + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                                    if histogramROISIFT.returncode==0: #Process success
                                                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + siftTerminus + '.csv', "for :", regionSelected, " ROI.")
                                                                        
                                                                        subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                                        
                                                    else :
                                                        # Computing Statistics of current 'SIFTed' ROI in the loop, getting : count, mean, median, std, min, max 
                                                        
                                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                        statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                        if statsROISIFT.returncode==0 : #Process success
                                                            statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                            
                                                        else :
                                                            statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                    regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                    '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                            if statsROISIFT.returncode==0 : #Process success
                                                                statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                            
                                                            else :
                                                                statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                        regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                        
                                                        # Writing all the computed statistics into the results text file 'fSIFT'
                                                        notStatsProperties = regionSelected + " " + RFalgorithms + " " + algorithms + " " + siftTerminus.lstrip('_') + " " + strCutoff + " " + str(stepSize)+ " "
                                                        fSIFT.write(notStatsProperties)
                                                        fSIFT.write(statsROISIFT_results + '\r\n')
                                                        
                                                        # + Creating a histogram of fibers count as a function of fibers length
                                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                        histogramROISIFT = subprocess.run(['tckstats', regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                        '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        if histogramROISIFT.returncode==0 : #Process success
                                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + siftTerminus + '.csv', "for :", regionSelected, " ROI.")
                                                            
                                                            subprocess.run(['cp', 'histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                                                                
                                                        else :
                                                            histogramROISIFT = subprocess.run(['tckstats', regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                            '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                            + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                            if histogramROISIFT.returncode==0 : #Process success
                                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + siftTerminus + '.csv', "for :", regionSelected, " ROI.")
                                                                
                                                                subprocess.run(['cp', 'histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            
                                                            else : 
                                                                histogramROISIFT = subprocess.run(['tckstats', regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                                '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                                if histogramROISIFT.returncode==0: #Process success
                                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + siftTerminus + '.csv', "for :", regionSelected, " ROI.")
                                                                    
                                                                    subprocess.run(['cp', 'histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + siftTerminus + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                            
                                    ###########################
                                        elif siftMode==2 :
                                            print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                            #tcksift2 [ options ]  in_tracks in_fod out_weights
                                            
                                            #Defining names here to simplify command lines 
                                            SIFT2RFAlgoName = 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            SIFT2noRFAlgoName = 'SIFT2weightstracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            SIFT2noRFnoAlgoName = 'SIFT2weightstracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            
                                            SIFT2_ROI_RFAlgoName = 'SIFT2weights_' + regionSelected + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            SIFT2_ROI_noRFAlgoName = 'SIFT2weights_' + regionSelected + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            SIFT2_ROI_noRFnoAlgoName = 'SIFT2weights_' + regionSelected + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            
                                            
                                            #Searching if SIFT2 has already been done in the subfolder, in which case there is no need to redo SIFT2...
                                            SIFT2_nameRFalgo = 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'  
                                            BooleanSIFT2_RFalgo = find_SIFT_files(os.getcwd(), SIFT2_nameRFalgo , False) 
                                            SIFT2_nameAlgo = 'SIFT2weightstracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            BooleanSIFT2_Algo = find_SIFT_files(os.getcwd(), SIFT2_nameAlgo , False) 
                                            SIFT2_name = 'SIFT2weightstracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                            BooleanSIFT2 = find_SIFT_files(os.getcwd(), SIFT2_name , False)
                                            
                                            if (BooleanSIFT2==False and BooleanSIFT2_Algo==False and BooleanSIFT2_RFalgo==False) :
                                                triggerSIFT2=True #Nothing was found so triggering SIFT2 
                                            else : 
                                                triggerSIFT2=False #SIFT2 was found so no need to do it again...
                                            
                                            if triggerSIFT2==True :
                                            
                                                SIFT2 = subprocess.run(
                                                    ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                                if SIFT2.returncode==0 : #Process success
                                                    print("New name convention detected, SIFT2 Success !")
                                                    print(SIFT2.stdout.decode('utf-8'))
                                                    copy = subprocess.run(['cp','SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                    individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    # print(copy.stdout.decode('utf-8'))
                                                else :
                                                    SIFT2 = subprocess.run(
                                                        ['tcksift2', 'tracts_' + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if SIFT2.returncode==0 : #Process success
                                                        print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT2 Success nonetheless !")
                                                        print(SIFT2.stdout.decode('utf-8'))
                                                        copy = subprocess.run(['cp','SIFT2weightstracts_' + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        # print(copy.stdout.decode('utf-8'))
                                                    else : 
                                                        SIFT2 = subprocess.run(
                                                            ['tcksift2', 'tracts_' + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + str(
                                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                                stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        print("Older name convention detected, SIFT2 Success nonetheless !")
                                                        print(SIFT2.stdout.decode('utf-8'))
                                                        copy = subprocess.run(['cp','SIFT2weightstracts_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        # print(copy.stdout.decode('utf-8'))
                                                                 
                                                    
                                            ####
                                            if (regionSelected in ['frR','frL'] and performROIFiltering): #"ROI1&not(ROI3)" scheme
                                                if regionSelected=='frR':
                                                    print('Analysis of the right Fasciculus Retroflexus after SIFT2')
                                                elif regionSelected=='frL':
                                                    print('Analysis of the left Fasciculus Retroflexus after SIFT2')
                                                
                                                    
                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                         'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                         + '_step-' + str(stepSize) + '_angle-' + str(angle) +'.tck', regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                         + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                         + str(angle) + '.tck', '-tck_weights_in',SIFT2RFAlgoName,'-tck_weights_out',SIFT2_ROI_RFAlgoName,
                                                                         '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                if result.returncode==0 : #Process success
                                                    print("ROI of :", regionSelected, "done.")
                                                    subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                    + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                    individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    subprocess.run(['cp',SIFT2_ROI_RFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    print(result.stdout.decode('utf-8'))

                                                else : 
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                             'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                             + '_step-' + str(stepSize) + '_angle-' + str(angle) +'.tck', regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                             + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                             + str(angle) +'.tck', '-tck_weights_in',SIFT2noRFAlgoName,'-tck_weights_out',SIFT2_ROI_noRFAlgoName,
                                                                             '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                    if result.returncode==0 : #Process success
                                                        print("ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        subprocess.run(['cp',SIFT2_ROI_RFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                                    
                                                    else : 
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-exclude','ROI-fr_3_not.mif',
                                                                                 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) +'.tck', regionSelected + '_' + str(numberOfSeeds) 
                                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                 + str(angle) +'.tck', '-tck_weights_in',SIFT2noRFnoAlgoName,'-tck_weights_out',SIFT2_ROI_noRFnoAlgoName,
                                                                                 '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                        print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) +'.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        subprocess.run(['cp',SIFT2_ROI_RFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                            
                                                        
                                            #####
                                            elif (regionSelected in ['och'] and performROIFiltering): 
                                            #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                                if regionSelected=='och':
                                                    print('Analysis of the Optic Chiasm')
                                                #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                                #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                                ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                                   ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                                   ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                                #Removed from cmd line : '-exclude','ROI-'+regionSelected+ochCombinations[k][2]+'.mif', '-exclude','ROI-'+regionSelected+ochCombinations[k][3]+'.mif',
                                                
                                                for k in range(len(ochCombinations)):
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                             'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                             + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                             + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                             + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2RFAlgoName,'-tck_weights_out',SIFT2_ROI_RFAlgoName,
                                                                             '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
            
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    if result.returncode==0 : #Process success
                                                        print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                        subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        subprocess.run(['cp',SIFT2_ROI_RFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                                    
                                                    else : 
                                                        
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                                 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                 + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                                 + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                 + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFAlgoName,'-tck_weights_out',SIFT2_ROI_noRFAlgoName,
                                                                                 '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                                                 
                                                        if result.returncode==0 : #Process success
                                                            print("ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                            subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            subprocess.run(['cp',SIFT2_ROI_noRFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))
                                                        
                                                        else :
                                                            result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+ochCombinations[k][0]+'.mif','-include','ROI-'+regionSelected+ochCombinations[k][1]+'.mif',
                                                                                     'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                     + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] 
                                                                                     + '_'  + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                                     + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFnoAlgoName,'-tck_weights_out',SIFT2_ROI_noRFnoAlgoName,
                                                                                     '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                            
                                                            print("Older name convention detected, ROI of :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], "done.")
                                                            subprocess.run(['cp',regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                            + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                            individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            subprocess.run(['cp',SIFT2_ROI_noRFnoAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            print(result.stdout.decode('utf-8'))      
                                            
                                            
                                            #####
                                            elif (regionSelected not in ['frR', 'frL', 'och'] and performROIFiltering): #"ROI1&ROI2&not(ROI3)" default scheme
                                                if regionSelected=='ac':
                                                    print('Analysis of the Anterior Commissure after SIFT2')
                                                elif regionSelected=='cc':
                                                    print('Analysis of the Corpus Callosum after SIFT2')
                                                elif regionSelected=='f':
                                                    print('Analysis of the Fornix after SIFT2')
                                                elif regionSelected=='fr':
                                                    print('Analysis of the Fasciculus Retroflexus after SIFT2')
                                                elif regionSelected=='icR-cpR':
                                                    print('Analysis of the right Internal Capsule and Cerebral Peduncle after SIFT2')
                                                elif regionSelected=='icL-cpL':
                                                    print('Analysis of the left Internal Capsule and Cerebral Peduncle after SIFT2')
                                                elif regionSelected=='mtR':
                                                    print('Analysis of the right MammilloThalamic tract after SIFT2')
                                                elif regionSelected=='mtL':
                                                    print('Analysis of the left MammilloThalamic tract after SIFT2')
                                                elif regionSelected=='opt':
                                                    print('Analysis of the Optical tract after SIFT2')
                                                elif regionSelected=='pfR':
                                                    print('Analysis of the right Postcommissural Fornix after SIFT2')
                                                elif regionSelected=='pfL':
                                                    print('Analysis of the left Postcommissural Fornix after SIFT2')
                                                elif regionSelected=='py':
                                                    print('Analysis of the Pyramidal tract after SIFT2')
                                                elif regionSelected=='smR':
                                                    print('Analysis of the right Stria Medularis after SIFT2')
                                                elif regionSelected=='smL':
                                                    print('Analysis of the left Stria Medularis after SIFT2')
                                                elif regionSelected=='stR':
                                                    print('Analysis of the right Stria Terminalis after SIFT2')
                                                elif regionSelected=='stL':
                                                    print('Analysis of the left Stria Terminalis after SIFT2')
                                                elif regionSelected=='vhc':
                                                    print('Analysis of the Ventral Hippocampal Commissure after SIFT2')
                                                    
                                                    
                                                result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                         'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                         + '_step-' + str(stepSize) + '_angle-' + str(angle) +'.tck', regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                         + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                         + str(angle) + '.tck', '-tck_weights_in',SIFT2RFAlgoName,'-tck_weights_out',SIFT2_ROI_RFAlgoName,
                                                                         '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                if result.returncode==0 : #Process success
                                                    print("ROI of :", regionSelected, "done.")
                                                    subprocess.run(['cp',regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                    + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                    individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    subprocess.run(['cp',SIFT2_ROI_RFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    print(result.stdout.decode('utf-8'))
                                                
                                                else :
                                                    result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                             'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                             + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                             + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                             + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFAlgoName,'-tck_weights_out',SIFT2_ROI_noRFAlgoName,
                                                                             '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                    if result.returncode==0 : #Process success
                                                        print("ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        subprocess.run(['cp',SIFT2_ROI_noRFAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                                    
                                                    else :
                                                        result = subprocess.run(['tckedit','-include','ROI-'+regionSelected+'_1_and.mif','-include','ROI-'+regionSelected+'_2_and.mif','-exclude','ROI-'+regionSelected+'_3_not.mif',
                                                                                 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                                                 + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', regionSelected + '_' + str(numberOfSeeds) 
                                                                                 + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                 + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFnoAlgoName,'-tck_weights_out',SIFT2_ROI_noRFnoAlgoName,
                                                                                 '-nthreads', str(nThreads),'-force'], stdout=subprocess.PIPE)
                                                        print("Older name convention detected, ROI of :", regionSelected, "done.")
                                                        subprocess.run(['cp',regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' 
                                                                        + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                        individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        subprocess.run(['cp',SIFT2_ROI_noRFnoAlgoName,individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        print(result.stdout.decode('utf-8'))
                                                                         
                                            ########
                                            if performROIStatistics :
                                                if (regionSelected in ['och']): 
                                                    #Optic Chiasm Scheme (and more generally Decussation scheme) : R1&R3(ipsiR),R1&R4(RcontraL),R2&R4(ipsiL),R2&R3(LcontraR),R1&R2(controlEyes),R3&R4(controlBrain)
                                                    #TODO: For future decussation studies, implement either another set of combinations or make it so the combinations are general (not righteye but rightext for example), but beware of the retrocompatibility aspect !
                                                    #Defining in an array the different combinations used to explore the och decussation, according to Decussation Scheme
                                                    ochCombinations = [['-1righteye','-3rightbrain','-2lefteye','-4leftbrain'],['-1righteye','-4leftbrain','-2lefteye','-3rightbrain'],
                                                                       ['-2lefteye','-4leftbrain','-1righteye','-3rightbrain'],['-2lefteye','-3rightbrain','-1righteye','-4leftbrain'],
                                                                       ['-1righteye','-2lefteye','-3rightbrain','-4leftbrain'],['-3rightbrain','-4leftbrain','-1righteye','-2lefteye']]
                                                    for k in range(len(ochCombinations)):
                                                        # Computing Statistics of current Decussation in the loop, getting : count, mean, median, std, min, max 
                                                        
                                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                        statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_RFAlgoName,
                                                                                '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                        if statsROISIFT.returncode==0 : #Process success
                                                            statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                            
                                                        else :
                                                            statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                    regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_noRFAlgoName,
                                                                                    '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                            if statsROISIFT.returncode==0 : #Process success
                                                                statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                            
                                                            else :
                                                                statsROISIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                        regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_noRFnoAlgoName,
                                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                statsROISIFT_results = statsROISIFT.stdout.decode('utf-8')
                                                        
                                                        # Writing all the computed statistics into the results text file 'fSIFT'
                                                        notStatsProperties = regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + " " + RFalgorithms + " " + algorithms + " " + 'SIFT2'+ " " + strCutoff + " " + str(stepSize)+ " "
                                                        fSIFT.write(notStatsProperties)
                                                        fSIFT.write(statsROISIFT_results + '\r\n')
                                                        
                                                        # + Creating a histogram of fibers count as a function of fibers length
                                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                        histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_RFAlgoName,
                                                                        '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        if histogramROISIFT.returncode==0 : #Process success
                                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                            
                                                            subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                                                                
                                                        else :
                                                            histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_noRFAlgoName,
                                                                            '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                            + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                            if histogramROISIFT.returncode==0 : #Process success
                                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected+ochCombinations[k][0]+ochCombinations[k][1], " ROI.")
                                                                
                                                                subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                            
                                                            else : 
                                                                histogramROISIFT = subprocess.run(['tckstats', regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2_ROI_noRFnoAlgoName,
                                                                                '-histogram', os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                                + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                                if histogramROISIFT.returncode==0: #Process success
                                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected, " ROI.")
                                                                    
                                                                    subprocess.run(['cp', 'histogram_' + regionSelected+ochCombinations[k][0]+ochCombinations[k][1] + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                
                                                else :
                                                    # Computing Statistics of current 'SIFTed' ROI in the loop, getting : count, mean, median, std, min, max 
                                                    
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    statsROISIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                            regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                            '-tck_weights_in',SIFT2_ROI_RFAlgoName, '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                            
                                                    if statsROISIFT2.returncode==0 : #Process success
                                                        statsROISIFT2_results = statsROISIFT2.stdout.decode('utf-8')
                                                        
                                                    else :
                                                        statsROISIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                                '-tck_weights_in',SIFT2_ROI_noRFAlgoName, '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                        
                                                        if statsROISIFT2.returncode==0 : #Process success
                                                            statsROISIFT2_results = statsROISIFT2.stdout.decode('utf-8')
                                                        
                                                        else :
                                                            statsROISIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                                    regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                                    '-tck_weights_in',SIFT2_ROI_noRFnoAlgoName, '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                            statsROISIFT2_results = statsROISIFT2.stdout.decode('utf-8')
                                                    
                                                    # Writing all the computed statistics into the results text file 'fSIFT'
                                                    notStatsProperties = regionSelected + " " + RFalgorithms + " " + algorithms + " " + "SIFT2" + " " + strCutoff + " " + str(stepSize)+ " "
                                                    fSIFT.write(notStatsProperties)
                                                    fSIFT.write(statsROISIFT2_results + '\r\n')
                                                    
                                                    
                                                    # + Creating a histogram of fibers count as a function of fibers length
                                                    ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                    histogramROISIFT2 = subprocess.run(['tckstats', regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', 
                                                                    '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + '_SIFT2' + '.csv', '-tck_weights_in',SIFT2_ROI_RFAlgoName, '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if histogramROISIFT2.returncode==0 : #Process success
                                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected, " ROI.")
                                                        
                                                        subprocess.run(['cp', 'histogram_' + regionSelected + '_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                                                            
                                                    else :
                                                        histogramROISIFT2 = subprocess.run(['tckstats', regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                        '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                        + str(angle) + '_SIFT2' + '.csv', '-tck_weights_in',SIFT2_ROI_noRFAlgoName, '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        if histogramROISIFT2.returncode==0 : #Process success
                                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected, " ROI.")
                                                            
                                                            subprocess.run(['cp', 'histogram_' + regionSelected + '_' + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                        
                                                        else : 
                                                            histogramROISIFT2 = subprocess.run(['tckstats', regionSelected + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                            '-histogram', os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                            + str(angle) + '_SIFT2' + '.csv', '-tck_weights_in',SIFT2_ROI_noRFnoAlgoName, '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                            if histogramROISIFT2.returncode==0: #Process success
                                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '_SIFT2' + '.csv', "for :", regionSelected, " ROI.")
                                                                
                                                                subprocess.run(['cp', 'histogram_' + regionSelected + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '_SIFT2' + '.csv', individualROIFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                                                                                           
                            
                quantityFlag -= 1
            
            fSIFT.close()
            if siftMode==1:
                copy = subprocess.run(['cp', newOutputPath + "/Results_SIFT_" + regionSelected + ".txt", ROIOutputPath], stdout=subprocess.PIPE) 
            elif siftMode==2:
                copy = subprocess.run(['cp', newOutputPath + "/Results_SIFT2_" + regionSelected + ".txt", ROIOutputPath], stdout=subprocess.PIPE) 

            print("------------------------------------------------------------------------------------------------------\n")
                                                                                                   
        
        print("\n\n####################################################################################################")
        print("\nROI Filtering and Statistics done !")
        print("\n######################################################################################################")
                
#8-END-------------------------------------------------------------------------------------------------------------------------

#9-COMPUTING GLOBAL STATISTICS (on data with and without SIFT) ----------------------------------------------------------------

if performWholeBrainStatistics : 
    newOutputPath = dataPath

    nowStats = datetime.datetime.now()
    statsOutputPath = newOutputPath + "/" + (nowStats.strftime("%Y%m%d_%H%M%S_")) + "StatsOnly_" + folderComment + "/"
    subprocess.run(['mkdir', statsOutputPath], stdout=subprocess.PIPE)
    
    # Copy parameter file to the statistics subdirectory
    inputFileFolder, inputFileName = os.path.split(inputFilePath)
    statsFileName = (nowStats.strftime("%Y%m%d_%H%M%S_")) + inputFileName
    result = subprocess.run(['cp', inputFilePath, statsOutputPath + statsFileName], stdout=subprocess.PIPE)
    
    
    # Creating output folder vector for each sample
    # Creating output folders for each sample
    individualStatsFolderNameVector=[]
    individualFolderNameVector=[]
    quantityFlag=quantityofSamples
    while (quantityFlag > 0):
        individualFolderNameVector.append('Sample-' + str(quantityofSamples-quantityFlag+1)) 
        individualStatsFolderNameVector.append(statsOutputPath + 'StatsSample-' + str(quantityofSamples-quantityFlag+1))
        proc = subprocess.run(['mkdir', individualStatsFolderNameVector[quantityofSamples-quantityFlag]], stdout=subprocess.PIPE)
        quantityFlag -= 1
    print("\n\n#####################################################################################################")
    print("\n\nOutput folders for each sample's statistics:")
    print("\n######################################################################################################")
    os.chdir(statsOutputPath)
    # Showing in terminal the different folders created using ls command
    proc = subprocess.run(['ls'])
    
    os.chdir(newOutputPath)
    
    fstats = open("Results" + "_" + "WholeBrainStatistics.txt","w+")
    fstats.write("Compilation of whole-brain tractograms statistics\r\n")

    quantityFlag = quantityofSamples    
    while (quantityFlag > 0):

        fstats.write("\n%s\r\n" % (individualFolderNameVector[quantityofSamples-quantityFlag]))
        fstats.write("Region RFAlgo TractoAlgo Cutoff Step-size StreamlinesCount Mean Median Std Min Max\r\n")

        os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples-quantityFlag])

        print("\n------------------------------------------------------------------------------------------------------")
        print(individualFolderNameVector[quantityofSamples-quantityFlag])
        print("------------------------------------------------------------------------------------------------------")

        for numberOfSeeds in numberofSeedsVector:
            for cutoff in cutoffVector:
                for stepSize in stepSizeVector:
                    
                    angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                    strCutoff = ('%.4f' % cutoff)
                    
                    
                    for RFalgorithms in sphericalRFAlgo : 
                        for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFod1,...)
                                    
                            if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2','Tensor_Det','Tensor_Prob','FACT']): #Making sure we don't take other cases than the algorithm names by accident
                             
    
                                ########
                                ### Computing Statistics of whole-brain tractogram in the loop, getting : count, mean, median, std, min, max 
                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)                                
                                stats = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                        'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                stats_results = stats.stdout.decode('utf-8')
                                
                                if stats.returncode==0 : #Process success
                                    stats_results = stats.stdout.decode('utf-8')
                                    
                                else :
                                    stats = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                            'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                            '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                    stats_results = stats.stdout.decode('utf-8')
                                    
                                    if stats.returncode==0 : #Process success
                                        stats_results = stats.stdout.decode('utf-8')
                                    
                                    else :
                                        stats = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                        stats_results = stats.stdout.decode('utf-8')
                                
                                # Writing all the computed statistics into the results text file 'fstats'
                                notStatsProperties ="Whole-brain" + " " + RFalgorithms + " " + algorithms + " " + strCutoff + " " + str(stepSize)+ " "
                                fstats.write(notStatsProperties)
                                fstats.write(stats_results + '\r\n')
                                
                                # + Creating a histogram of fibers count as a function of fibers length
                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                histogram = subprocess.run(['tckstats', 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                               
                                if histogram.returncode==0 : #Process success
                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                    + str(angle) + '.csv')
                                    
                                    subprocess.run(['cp', 'histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                    + str(angle) + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                           
                                else :
                                    histogram = subprocess.run(['tckstats', 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                    '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                                   
                                    if histogram.returncode==0 : #Process success
                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                        + str(angle) + '.csv')
                                        
                                        subprocess.run(['cp', 'histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                        + str(angle) + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                        
                                    else : 
                                        histogram = subprocess.run(['tckstats', 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                        '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                        if histogram.returncode==0: #Process success
                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                            + str(angle) + '.csv')
                                            
                                            subprocess.run(['cp', 'histogram_' + 'tracts_' + str(numberOfSeeds) 
                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                            + str(angle) + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
        
        quantityFlag -=1
    
    fstats.close()
    copy = subprocess.run(['cp', newOutputPath + "/Results_WholeBrainStatistics.txt", statsOutputPath], stdout=subprocess.PIPE) 
    
    print("\n\n####################################################################################################")
    print("\nComputing whole-brain tractograms Statistics done !")
    print("\n######################################################################################################")
        
###############
    if siftMode!=0 :
    
        print("\n\n####################################################################################################")
        print("\nSwitching to computing whole-brain tractograms Statistics after SIFT...")
        print("\n######################################################################################################")
        
        os.chdir(newOutputPath)
        
        siftTerminii=[]
        if siftMode==1 :
            fstatsSIFT = open("Results_SIFT_WholeBrainStatistics.txt", "w+")
            fstatsSIFT.write("Compilation of whole-brain tractograms statistics after SIFT\r\n")
            for siftType in siftTerminationOptionsArray:
                for siftValue in siftTerminationValuesArray:
        
                    if (siftType == 'standard'):
                        siftTerminus = '_SIFT-standard'
                        siftTerminii.append(siftTerminus)
        
                    elif (siftType == 'mu'):
                        siftTerminus = '_SIFT-term-mu-' + str(siftValue)
                        siftTerminii.append(siftTerminus)
        
                    elif (siftType == 'ratio'):
                        siftTerminus = '_SIFT-term-ratio-' + str(siftValue)
                        siftTerminii.append(siftTerminus)
        elif siftMode==2 :
            fstatsSIFT = open("Results_SIFT2_WholeBrainStatistics.txt", "w+")
            fstatsSIFT.write("Compilation of whole-brain tractograms statistics after SIFT2\r\n")
    
        quantityFlag = quantityofSamples
        while (quantityFlag > 0):
            fstatsSIFT.write("\n%s\r\n" % (individualFolderNameVector[quantityofSamples - quantityFlag]))
            fstatsSIFT.write("Region RFAlgo TractoAlgo SIFT Cutoff Step-size StreamlinesCount Mean Median Std Min Max\r\n")
    
            os.chdir(newOutputPath + '/' + individualFolderNameVector[quantityofSamples - quantityFlag])
    
            print("\n------------------------------------------------------------------------------------------------------")
            print(individualFolderNameVector[quantityofSamples - quantityFlag])
            print("------------------------------------------------------------------------------------------------------")
            
            for numberOfSeeds in numberofSeedsVector:
                for cutoff in cutoffVector:
                    for stepSize in stepSizeVector:
                        
                        angle = round(math.degrees(2 * math.asin(stepSize / (2 * curvature))), 4)
                        strCutoff = ('%.4f' % cutoff)
                        
                        
                        for RFalgorithms in sphericalRFAlgo : 
                            for algorithms in tractAlgo : #Looping through algorithms specified in tractAlgo parameter (Tensor_Det, iFod1,...)
                                        
                                if (algorithms in ['iFOD1','iFOD2','SD_STREAM','Nulldist1','Nulldist2','Tensor_Det','Tensor_Prob','FACT']): #Making sure we don't take other cases than the algorithm names by accident
                                    
                                    if siftMode==1 :
                                        print('\nSIFT mode 1 chosen, tcksift in progress...')
                                        #tcksift [ options ]  in_tracks in_fod out_tracks
                                        for siftTerminus in siftTerminii :
                    
                                            #Searching if SIFT has already been done in the subfolder, in which case there is no need to redo SIFT...
                                            SIFT_nameRFalgo = 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'  
                                            BooleanSIFT_RFalgo = find_SIFT_files(os.getcwd(), SIFT_nameRFalgo , False) 
                                            SIFT_nameAlgo = 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'
                                            BooleanSIFT_Algo = find_SIFT_files(os.getcwd(), SIFT_nameAlgo , False) 
                                            SIFT_name = 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck'
                                            BooleanSIFT = find_SIFT_files(os.getcwd(), SIFT_name , False)
                                            
                                            if (BooleanSIFT==False and BooleanSIFT_Algo==False and BooleanSIFT_RFalgo==False) :
                                                triggerSIFT=True #Nothing was found so triggering SIFT 
                                            else : 
                                                triggerSIFT=False #SIFT was found so no need to do it again...
                                                
                                                                 
                                            if ('standard' in siftTerminus) and triggerSIFT==True :
                                                print('Performing standard sift...\r\n')
                                                SIFT = subprocess.run(
                                                    ['tcksift', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                if SIFT.returncode==0 : #Process success
                                                    print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    
                                                else :
                                                    SIFT = subprocess.run(
                                                        ['tcksift', 'tracts_' + algorithms + '_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if SIFT.returncode==0 : #Process success
                                                        print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                        print(SIFT.stdout.decode('utf-8'))
                                                        
                                                    else : 
                                                        SIFT = subprocess.run(['tcksift', 'tracts_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                        print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                        print(SIFT.stdout.decode('utf-8'))
                                                        
                                            elif ('mu' in siftTerminus) and triggerSIFT==True:
                                                
                                                #Retrieving mu value
                                                SIFTmuStringvalue = siftTerminus.split('-')[-1]
                                                                                    
                                                print('Performing mu sift...\r\n')
                                                SIFT = subprocess.run(
                                                    ['tcksift', '-term_mu', SIFTmuStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                if SIFT.returncode==0 : #Process success
                                                    print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    
                                                else :
                                                     SIFT = subprocess.run(
                                                         ['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + algorithms + '_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                     if SIFT.returncode==0 : #Process success
                                                         print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                         print(SIFT.stdout.decode('utf-8'))
                                                         
                                                     else : 
                                                         SIFT = subprocess.run(['tcksift', '-term_mu', SIFTmuStringvalue, 'tracts_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                         print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                         print(SIFT.stdout.decode('utf-8'))
                                                         
                                            elif ('ratio' in siftTerminus) and triggerSIFT==True:
                                                
                                                #Retrieving ratio value
                                                SIFTratioStringvalue = siftTerminus.split('-')[-1]
                                                
                                                print('Performing ratio sift...\r\n')
                                                SIFT = subprocess.run(
                                                    ['tcksift', '-term_ratio', SIFTratioStringvalue,'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo) 
                                                if SIFT.returncode==0 : #Process success
                                                    print("New name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                    print(SIFT.stdout.decode('utf-8'))
                                                    
                                                else :
                                                     SIFT = subprocess.run(
                                                         ['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + algorithms + '_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + algorithms + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                     if SIFT.returncode==0 : #Process success
                                                         print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                         print(SIFT.stdout.decode('utf-8'))
                                                         
                                                     else : 
                                                         SIFT = subprocess.run(['tcksift', '-term_ratio', SIFTratioStringvalue, 'tracts_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'tracts_' + str(
                                                             numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                             stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                         print("Older name convention detected, SIFT Succeeded with :", siftTerminus, "termination option.")
                                                         print(SIFT.stdout.decode('utf-8'))
                                                
                                            ########
                                            # Computing Statistics of current 'SIFTed' tractograms in the loop, getting : count, mean, median, std, min, max 
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                            statsSIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                    'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                    '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                            statsSIFT_results = statsSIFT.stdout.decode('utf-8')
                                            
                                            if statsSIFT.returncode==0 : #Process success
                                                statsSIFT_results = statsSIFT.stdout.decode('utf-8')
                                                
                                            else : 
                                                statsSIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                        'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                        '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                statsSIFT_results = statsSIFT.stdout.decode('utf-8')
                                                
                                                if statsSIFT.returncode==0 : #Process success
                                                    statsSIFT_results = statsSIFT.stdout.decode('utf-8')
                                                
                                                else : 
                                                    statsSIFT = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                            'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                            '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                    statsSIFT_results = statsSIFT.stdout.decode('utf-8')
                                    
                                            # Writing all the computed statistics into the results text file 'fstatsSIFT'
                                            notStatsProperties = "Whole-brain" + " " + RFalgorithms + " " + algorithms + " " + siftTerminus.lstrip('_') + " " + strCutoff + " " + str(stepSize)+ " "
                                            fstatsSIFT.write(notStatsProperties)
                                            fstatsSIFT.write(statsSIFT_results + '\r\n')
                                            
                                            
                                            # + Creating a histogram of fibers count as a function of fibers length
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                            histogramSIFT = subprocess.run(['tckstats', 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck', 
                                                            '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            
                                            if histogramSIFT.returncode==0 : #Process success
                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + siftTerminus +'.csv')
                                                
                                                subprocess.run(['cp', 'histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + siftTerminus + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
            
                                            else :
                                                histogramSIFT = subprocess.run(['tckstats', 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                
                                                if histogramSIFT.returncode==0 : #Process success
                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + siftTerminus +'.csv')
                                                    
                                                    subprocess.run(['cp', 'histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + siftTerminus + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                
                                                else : 
                                                    histogramSIFT = subprocess.run(['tckstats', 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + siftTerminus + '.tck',
                                                                    '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                    + str(angle) + siftTerminus + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    if histogramSIFT.returncode==0: #Process success
                                                        print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + siftTerminus + '.csv')
                                                        
                                                        subprocess.run(['cp', 'histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                        + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                        + str(angle) + siftTerminus +'.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                
                                #########################
                                    elif siftMode==2 :
                                        print('\nSIFT mode 2 chosen, tcksift2 in progress...')
                                        #tcksift2 [ options ]  in_tracks in_fod out_weights
                                        
                                        #Searching if SIFT2 has already been done in the subfolder, in which case there is no need to redo SIFT2...
                                        SIFT2_nameRFalgo = 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'  
                                        BooleanSIFT2_RFalgo = find_SIFT_files(os.getcwd(), SIFT2_nameRFalgo , False) 
                                        SIFT2_nameAlgo = 'SIFT2weightstracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                        BooleanSIFT2_Algo = find_SIFT_files(os.getcwd(), SIFT2_nameAlgo , False) 
                                        SIFT2_name = 'SIFT2weightstracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                        BooleanSIFT2 = find_SIFT_files(os.getcwd(), SIFT2_name , False)
                                        
                                        if (BooleanSIFT2==False and BooleanSIFT2_Algo==False and BooleanSIFT2_RFalgo==False) :
                                            triggerSIFT2=True #Nothing was found so triggering SIFT2 
                                        else : 
                                            triggerSIFT2=False #SIFT2 was found so no need to do it again...
                                        
                                        if triggerSIFT2==True :
                                        
                                            SIFT2 = subprocess.run(
                                                ['tcksift2', 'tracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm_' + RFalgorithms + '.mif', 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                            if SIFT2.returncode==0 : #Process success
                                                print("New name convention detected, SIFT2 Success !")
                                                print(SIFT2.stdout.decode('utf-8'))
                                                copy = subprocess.run(['cp','SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(
                                                numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                # print(copy.stdout.decode('utf-8'))
                                            else :
                                                SIFT2 = subprocess.run(
                                                    ['tcksift2', 'tracts_' + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + algorithms + '_' + str(
                                                        numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                        stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                if SIFT2.returncode==0 : #Process success
                                                    print("New name convention detected but without", RFalgorithms,"algorithm name for wmfod files, SIFT2 Success nonetheless !")
                                                    print(SIFT2.stdout.decode('utf-8'))
                                                    copy = subprocess.run(['cp','SIFT2weightstracts_' + algorithms + '_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                    individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    # print(copy.stdout.decode('utf-8'))
                                                else : 
                                                    SIFT2 = subprocess.run(
                                                        ['tcksift2', 'tracts_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.tck', 'wmfod_norm.mif', 'SIFT2weightstracts_' + str(
                                                            numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                            stepSize) + '_angle-' + str(angle) + '.txt', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                    print("Older name convention detected, SIFT2 Success nonetheless !")
                                                    print(SIFT2.stdout.decode('utf-8'))
                                                    copy = subprocess.run(['cp','SIFT2weightstracts_' + str(
                                                    numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(
                                                    stepSize) + '_angle-' + str(angle) + '.txt', 
                                                                    individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                                    # print(copy.stdout.decode('utf-8'))
                                                                                                   
                                        ########   
                                        # Computing Statistics of current 'SIFTed' tractograms in the loop, getting : count, mean, median, std, min, max 
                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                        SIFT2RFAlgoName = 'SIFT2weightstracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                        SIFT2noRFAlgoName = 'SIFT2weightstracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'
                                        SIFT2noRFnoAlgoName = 'SIFT2weightstracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.txt'

                                        statsSIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                '-tck_weights_in', SIFT2RFAlgoName ,'-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                        statsSIFT2_results = statsSIFT2.stdout.decode('utf-8')
                                        
                                        if statsSIFT2.returncode==0 : #Process success
                                            statsSIFT2_results = statsSIFT2.stdout.decode('utf-8')
                                            
                                        else : 
                                            
                                            statsSIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                    'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                    + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                    '-tck_weights_in', SIFT2noRFAlgoName, '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                                    
                                            statsSIFT2_results = statsSIFT2.stdout.decode('utf-8')
                                            
                                            if statsSIFT2.returncode==0 : #Process success
                                                statsSIFT2_results = statsSIFT2.stdout.decode('utf-8')
                                            
                                            else : 

                                                statsSIFT2 = subprocess.run(['tckstats', '-output', 'count', '-output', 'mean', '-output', 'median', '-output', 'std', '-output', 'min', '-output', 'max',
                                                                        'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck',
                                                                        '-tck_weights_in', SIFT2noRFnoAlgoName, '-nthreads', str(nThreads)], stdout=subprocess.PIPE)
                                                statsSIFT2_results = statsSIFT2.stdout.decode('utf-8')
                                
                                        # Writing all the computed statistics into the results text file 'fstatsSIFT'
                                        notStatsProperties = "Whole-brain" + " " + RFalgorithms + " " + algorithms + " " + "SIFT2" + " " + strCutoff + " " + str(stepSize)+ " "
                                        fstatsSIFT.write(notStatsProperties)
                                        fstatsSIFT.write(statsSIFT2_results + '\r\n')
                                        
                                        # + Creating a histogram of fibers count as a function of fibers length
                                        ##Followed check scheme : Algo&RFAlgo ? > Algo&not(RFAlgo) ? > not(Algo)&not(RFAlgo)
                                        histogramSIFT2 = subprocess.run(['tckstats', 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                        + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2RFAlgoName,
                                                        '-histogram', os.getcwd() + '/histogram_' + 'tracts_' 
                                                        + RFalgorithms + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' 
                                                        + '_step-' + str(stepSize) + '_angle-' + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), 
                                                        '-force'], stdout=subprocess.PIPE)
                                                        
                                        
                                        if histogramSIFT2.returncode==0 : #Process success
                                            print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                            + str(angle) + '_SIFT2' +'.csv')
                                            
                                            subprocess.run(['cp', 'histogram_' + 'tracts_' + RFalgorithms + algorithms + '_' + str(numberOfSeeds) 
                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                            + str(angle) + '_SIFT2' + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
        
                                        else :
                                            histogramSIFT2 = subprocess.run(['tckstats', 'tracts_' + algorithms + '_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                            + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFAlgoName,
                                                            '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                            + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                            + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                            
                                            if histogramSIFT2.returncode==0 : #Process success
                                                print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + '_SIFT2' +'.csv')
                                                
                                                subprocess.run(['cp', 'histogram_' + 'tracts_' + algorithms + '_' + str(numberOfSeeds) 
                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                + str(angle) + '_SIFT2' + '.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE)
                                            
                                            else : 
                                                histogramSIFT2 = subprocess.run(['tckstats', 'tracts_' + str(numberOfSeeds) + '_cutoff-' + strCutoff 
                                                                + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' + str(angle) + '.tck', '-tck_weights_in',SIFT2noRFnoAlgoName,
                                                                '-histogram', os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                                + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                                + str(angle) + '_SIFT2' + '.csv', '-nthreads', str(nThreads), '-force'], stdout=subprocess.PIPE)
                                                if histogramSIFT2.returncode==0: #Process success
                                                    print("Distribution of fibers count over length computed and saved at :", os.getcwd() + '/histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + '_SIFT2' + '.csv')
                                                    
                                                    subprocess.run(['cp', 'histogram_' + 'tracts_' + str(numberOfSeeds) 
                                                    + '_cutoff-' + strCutoff + '_minl-0.1_maxl-50' + '_step-' + str(stepSize) + '_angle-' 
                                                    + str(angle) + '_SIFT2' +'.csv', individualStatsFolderNameVector[quantityofSamples-quantityFlag] + '/' ], stdout=subprocess.PIPE) 
                                            
                                            
            quantityFlag -=1
            
        fstatsSIFT.close()
        if siftMode==1:
            copy = subprocess.run(['cp', newOutputPath + "/Results_SIFT_WholeBrainStatistics.txt", statsOutputPath], stdout=subprocess.PIPE)
        elif siftMode==2:
            copy = subprocess.run(['cp', newOutputPath + "/Results_SIFT2_WholeBrainStatistics.txt", statsOutputPath], stdout=subprocess.PIPE)
        
        print("\n\n####################################################################################################")
        print("\nComputing whole-brain tractograms Statistics after SIFT done !")
        print("\n######################################################################################################")    

#9-END-------------------------------------------------------------------------------------------------------------------------

print("\n\n####################################################################################################")
print("\nEnd of Script !")
print("\n######################################################################################################")
