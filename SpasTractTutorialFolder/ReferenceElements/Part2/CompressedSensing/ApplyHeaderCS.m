%% This Script transfers .nii header properties from a reference file onto another .nii file's header while preserves voxel dimensions etc. %%

%% Script Start and Checkpoint
clear all;

%Define the path to NIfTI(s) folder
%from where other inputs asked will start from this directory
[WorkingFolder] = uigetdir('/','Select the folder containing NIfTI CS file(s)...');
if isequal(WorkingFolder, 0)
    error('No directory selected.');
end
%Example of WorkingFolder : '/SUMO/DATA_WORK'; (no end slash !)

%Asking user how many NIfTI to process : 
prompt = "How many NIfTI file(s) do you want to process ?";
niiNumber = input(prompt);

for k = 1:1:niiNumber

    %Define the path to NIfTI(s) to modify
    [nii_filename, nii_path] = uigetfile({'*.nii', 'Supported Files (*.nii)';'*.*','All Files (*.*)'}, ...
        'Select the correctly oriented NIfTI file', fullfile(WorkingFolder, 'No File selected yet'));
    if isequal(nii_filename, 0)
        error('No reference file selected.');
    end
    niifilePath = fullfile(nii_path, nii_filename);
    %Example of niifilePath : '/SUMO/DATA_WORK/Mjo_M71_3DCS_FLASH.nii';
    
    % Changing NIfTI header [DONE FOR 50µm SETUP]
    % mrconvert step
    % Constructing the shell command
    niifilePathREOR = erase(niifilePath,".nii") + "_REOR.nii";
    cmda = sprintf('mrconvert %s %s -axes 2,0,1 -vox 0.05,0.05,0.05 -strides 3,-2,1 -force', niifilePath, niifilePathREOR);
    
    % Running the command in the terminal
    [status, result] = system(cmda);
    
    % Optional: checking for errors
    if status ~= 0
        error('mrconvert failed: %s', result);
    else
        disp('mrconvert completed successfully.');
    end
    
    % mrtransform step
    %Generating transform txt file for mrtransform
    % Define the 3x4 matrix (no offset for dwi otherwise replace -3 by -2.2)
    T = [1  -0  0  -5.8;
         0  -1  0   8.9;
         0  -0  1  -3];              
    % Write the matrix to a text file
    filename = 'TransformMRIFilesManager_3DCS-1+8.9-AxAP.txt';
    transformFilePath = fullfile(WorkingFolder, filename);
    writematrix(T, transformFilePath, 'Delimiter', ' ');                 
    disp(['File "', filename, '" written successfully.']);
    
    % Constructing the shell command 
    cmdb = sprintf('mrtransform %s %s -replace %s -force', niifilePathREOR, niifilePathREOR, transformFilePath);
    
    % Running the command in the terminal
    [status, result] = system(cmdb);
    
    % Optional: checking for errors
    if status ~= 0
        error('mrtransform failed: %s', result);
    else
        disp('mrtransform completed successfully.');
    end

end

disp('ApplyHeaderCS script done !')

