%% PCNN3D auto brain extraction
% saves as *_mask**parameters**.nii.gz
% requires nifti toolbox
%
% 2017/07/27 ver1.0
% 2017/08/29 ver1.1 bug fix; add ZoomFactor; use mean image for 4D data
% 2017/09/12 ver1.2 save as 8-bit mask
% 2018/01/23 ver1.3 use with command line bash script
% Kai-Hsiang Chuang, QBI/UQ
% 2025/02 Slight changes by Arno ROSSEZ
% 2025/02 More changes by Arno ROSSEZ
%% This Script creates a compressed copy of selected .nii file into a .nii.gz file & runs PCNN3D on it with set parameters. %%
%% The resulting mask is then saved in the working folder with a name indicating which parameters were used on which file. %%
%% Script Start and Checkpoint
clear all;

%Define the directory from where you want to start from, like a checkpoint
%where other input asked will be asked starting from this directory so
%that it easier to navigate around upstream or downstream.
[MaskFolder] = uigetdir(pwd,'Select the directory where you want to start working from (it will be your checkpoint)');
if isequal(MaskFolder, 0)
    error('No directory selected.');
end
%MaskFolder = '/SUMO/DATA_WORK/';

%% Selecting .nii file for PCNN3D
[filename,path]=uigetfile({'*.nii', 'Supported Files (*.nii)';'*.*','All Files (*.*)'}, ...
    'Select a .nii File for mask-making', fullfile(MaskFolder, 'No File selected yet')); % data path 
% Check if the user selected a file or canceled
if isequal(filename, 0)
    disp('User canceled file selection.');
else
    datpath = fullfile(path, filename); %datpath is the path leading to the exact selected .nii.gz file
    disp(['User selected: ', datpath]);
end

%% Compressing  a copy of selected .nii file into .nii.gz format used by PCNN3D
% Define the output .nii.gz filename
gzdatpath = strcat(datpath, '.gz');

% Compress the .nii file (creates a .nii.gz copy)
gzip(datpath);

% Display confirmation
fprintf('Created copy: %s -> %s\n', datpath, gzdatpath);

%% Setting parameters for PCNN3D 
%BrSize=[350,550]; % brain size range for MOUSE (mm3).
%BrSize=[1500,3000]; % brain size range for RAT > infinite iterations on MOUSE and mask is full picture (mm3)
%BrSize=[1,1500]; % brain size range Olivier (mm3) > NOPE
%BrSize=[350,2000]; % brain size range after measurement mrview/mri-conv(mm3) > NOPE
%BrSize=[425,475]; % brain size short range mrview+prior testing (mm3)
%BrSize=[20,500]; % brain size range David (mm3)
BrSize=[475,575]; % brain size optimal MOUSE ?(mm3)
StrucRadius=3; % use =3 for low resolution, use 5 or 7 for highres data
ZoomFactor=1; % resolution magnification factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run PCNN
[nii] = load_untouch_nii(gzdatpath);
mtx=size(nii.img);
if length(mtx)==4
    disp('Data is 4D, use the average image to generate mask')
    nii.img=mean(nii.img,4);
end
voxdim=nii.hdr.dime.pixdim(2:4);
[I_border, G_I, optG] = PCNN3D(single(nii.img), StrucRadius, voxdim, BrSize*ZoomFactor^3);
V=zeros(mtx);
for n=1:mtx(3)
    V(:,:,n)=I_border{optG}{n};
end

%% save data
disp(['Saving mask at ',datpath(1:end-7),'_mask.nii.gz....'])
nii.img=V;
nii.hdr.dime.dim(1)=3; nii.hdr.dime.dim(5)=1;
nii.hdr.dime.datatype=2; nii.hdr.dime.bitpix=8; % save as unsigned char
p=strfind(datpath,'.nii'); % Find the position of '.nii' to remove the extension correctly

new_filename = [datpath(1:p-1), '_mask', num2str(StrucRadius), '_', num2str(ZoomFactor), '_', ...
    num2str(BrSize(1)), '-', num2str(BrSize(2)), '.nii.gz']; % Construct new filename with the required suffix

% Save the modified NIfTI file
save_untouch_nii(nii, new_filename);
%save_untouch_nii(nii,[datpath(1:p-1),'_mask.nii.gz'])

disp('Done')