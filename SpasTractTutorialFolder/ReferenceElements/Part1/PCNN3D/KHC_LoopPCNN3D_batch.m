%% PCNN3D batch auto brain extraction WORK IN PROGRESS
% saves as *_mask**parameters**.nii.gz
% requires nifti toolbox
%
% 2017/07/27 ver1.0
% 2017/08/29 ver1.1 bug fix; add ZoomFactor; use mean image for 4D data
% 2017/09/12 ver1.2 save as 8-bit mask
% 2018/01/23 ver1.3 use with command line bash script
% Kai-Hsiang Chuang, QBI/UQ
% 2025/02 Slight changes by Arno ROSSEZ
% 2025/03 More changes by Arno ROSSEZ to loop and have .nii in working folder

%% This Script creates a compressed copy of every .nii files contained in the working folder into .nii.gz files & loops with them for different parameters. %%
%% Masks are then saved in a subfolders with a name indicating which parameters were used on which file. %%
%% Script Start and Checkpoint
clear all;

%Define the directory from where you want to start from, like a checkpoint
%where other input asked will be asked starting from this directory so
%that it easier to navigate around upstream or downstream.
[MaskFolder] = uigetdir(pwd,'Select the directory where the .nii files to be processed for maskmaking are...');
if isequal(MaskFolder, 0)
    error('No directory selected.');
end
%MaskFolder = '/SUMO/DATA_WORK/';

%% Get all .nii files in WorkingFolder to then compress them to .nii.gz format used by PCNN3D
niiFiles = dir(fullfile(MaskFolder, '*.nii'));

% Process each .nii file and create a compressed .nii.gz copy
for i = 1:length(niiFiles)
    % Full path of .nii file
    niiFilePath = fullfile(MaskFolder, niiFiles(i).name);
    
    % Define the output .nii.gz filename
    gzFilePath = strcat(niiFilePath, '.gz');
    
    % Compress the .nii file (creates a .nii.gz copy)
    gzip(niiFilePath);
    
    % Display confirmation
    fprintf('Created copy: %s -> %s\n', niiFilePath, gzFilePath);
end

%% Getting .nii.gz files created for PCNN3D
niigzFiles = dir(fullfile(MaskFolder, '*.nii.gz'));

% Construct output folder for all the masks
[output_path] = uigetdir(MaskFolder,'Select the output folder for masks');
if isequal(output_path, 0)
    error('No directory selected.');
end

% BrSizeArray=[[350,550];[475,575]]; % brain size array of optimal MOUSE and classical MOUSE(mm3)
BrSizeArray=[50,150]; % brain size array test for eyes+optic chiasm MOUSE (mm3)
%BrSizeArray=[475,575]; 
[BrSizenumberconfig,~]=size(BrSizeArray); % getting the number of configurations used (number of lines)
%StrucRadiusArray=[3,4,5,6,7]; % use =3 for low resolution, use 5 or 7 for highres data
StrucRadiusArray=[3]; 
ZoomFactor=1; % resolution magnification factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run PCNN loop with batch

for j=1:1:length(niigzFiles)
    for k=1:1:length(StrucRadiusArray)
        for l=1:1:BrSizenumberconfig
            StrucRadius = StrucRadiusArray(k); %Getting the StrucRadius to be used
            BrSize = BrSizeArray(l,:); %Getting the BrSize configuration

            % Reload NIfTI data at every loop to ensure no data corruption/overwriting
            [nii] = load_untouch_nii(fullfile(MaskFolder,niigzFiles(j).name));
            mtx=size(nii.img);
            if length(mtx)==4
                disp('Data is 4D, use the average image to generate mask')
                nii.img=mean(nii.img,4);
            end
            voxdim=nii.hdr.dime.pixdim(2:4);

            %Action of PCNN3D below
            [I_border, G_I, optG] = PCNN3D(single(nii.img), StrucRadius, voxdim, BrSize*ZoomFactor^3);
            V=zeros(mtx);
            for n=1:mtx(3)
                V(:,:,n)=I_border{optG}{n};
            end

            %% save data
            disp(['Saving mask at ',output_path,'/',niigzFiles(j).name(1:end-7),'_mask.nii.gz....'])
            nii.img=V;
            nii.hdr.dime.dim(1)=3; nii.hdr.dime.dim(5)=1;
            nii.hdr.dime.datatype=2; nii.hdr.dime.bitpix=8; % save as unsigned char
            p=strfind(niigzFiles(j).name,'.nii'); % Find the position of '.nii' to remove the extension correctly

    %         % Extract filename without extension
    %         [~, name, ~] = fileparts(filename); 
    %         % Construct the output filename (WITHOUT .gz)
    %         mask_filename = sprintf('%s_mask%d_%d_%d-%d.nii', ...
    %             name, StrucRadius, ZoomFactor, BrSize(1), BrSize(2));

            mask_filename = [niigzFiles(j).name(1:p-1), '_mask', num2str(StrucRadius), '_', num2str(ZoomFactor), '_', ...
                num2str(BrSize(1)), '-', num2str(BrSize(2)), '.nii.gz']; % Construct new filename with the required suffix

            % Full path for the new .nii file
            mask_filepath = fullfile(output_path, mask_filename);

            % Save the modified NIfTI file as .nii
            save_untouch_nii(nii, mask_filepath);
            disp(['Mask file saved as: ', mask_filepath]);

        end
    end
end
disp('Done with the loop')

