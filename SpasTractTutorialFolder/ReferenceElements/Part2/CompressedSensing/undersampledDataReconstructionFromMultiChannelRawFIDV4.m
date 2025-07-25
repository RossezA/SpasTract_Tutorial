function [undersampledKSpace, undersampledKSpaceView, undersampledZPImage, undersampledZPImageView, undersamplingPatterRetrieved] = undersampledDataReconstructionFromMultiChannelRawFIDV4(dataPath, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, correctImage, fftShiftCorrection, offsetCorrection, showResults, saveResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, rotAngleSliceX, rotAngleSliceY, rotAngleSliceZ, savePathAndFileNameundersampled, pathTxtDoc, imgstoragePath)

% -------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------
% Fully sampled -----------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------

% Redifining whichSliceX for first view (wihtout fftshift correction) ---
whichSliceXTemp = whichSliceX;
whichSliceX = 1;
% -----------------------------------------------------------------------

undersampledKSpace = zeros(finalReadoutDim,finalPhase1Dim,finalPhase2Dim,nChannels,nVolumes);
undersamplingPatterRetrieved = zeros(finalReadoutDim,finalPhase1Dim,finalPhase2Dim,nChannels,nVolumes);

fprintf("Number of channels to reconstruct: %d\n",nChannels);
for channel=0:(nChannels-1)    
    fprintf("Reconstructing channel %d\n", channel);
    fileID = strcat(dataPath,num2str(channel));
    [undersampledKSpacePerChannel, ~, csMaskRetrieved] =  readParavisionUndersampledRawFIDV3(fileID, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, saveResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, imgstoragePath);           
    undersampledKSpace(:,:,:,channel+1,:) = squeeze(undersampledKSpacePerChannel);    
    undersamplingPatterRetrieved(:,:,:,channel+1,:) = squeeze(csMaskRetrieved);
    
end

%[undersampledKSpace0, undersampledImage0] =  readParavisionUndersampledRawFIDV3(fileID, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, dataTypeFS, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, imgstoragePath);


%fileID = strcat(dataPath,'1');
%[undersampledKSpace1, ~, csMaskRetrieved1] =  readParavisionUndersampledRawFIDV3(fileID, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, imgstoragePath);




%undersampledKSpace(:,:,:,1,:) = squeeze(undersampledKSpace0);
%undersampledKSpace(:,:,:,2,:) = squeeze(undersampledKSpace1);
%undersampledKSpace(:,:,:,3,:) = squeeze(undersampledKSpace2);
%undersampledKSpace(:,:,:,4,:) = squeeze(undersampledKSpace3);


% Redifining whichSliceX for correct view ---
whichSliceX = whichSliceXTemp;
% -------------------------------------------

% Performing phase correction -------------------------------------------------------------------------------------------------------------------------------------------------
%[undersampledKSpace, recoImageCorrected] =  PhaseCorrection (undersampledKSpace, 2, phaseCorrectionValue, 1, whichSliceX, whichSliceY, whichSliceZ, whichChannel, whichVolume);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Making k-space for visualization ------------------------
undersampledKSpaceView = squeeze(bart('rss 8', undersampledKSpace));
% ---------------------------------------------------------

% Creating image ----------------------------------------
undersampledZPImage = bart('fft -i 7', undersampledKSpace);
undersampledZPImageView = bart('rss 8', undersampledZPImage);
undersampledZPImageView = squeeze(undersampledZPImageView);
% -------------------------------------------------------


if (correctImage)
    % Correcting image with FFTshift --------------------------------------------
    if (fftShiftCorrection(1))
        undersampledZPImageView = fftshift(undersampledZPImageView,1);
    end
    if (fftShiftCorrection(2))
        undersampledZPImageView = fftshift(undersampledZPImageView,2);
    end
    if (fftShiftCorrection(3))
        undersampledZPImageView = fftshift(undersampledZPImageView,3);
    end
    % ---------------------------------------------------------------------------
    
    
    % Correct image offset ----------------------------------------------------
    undersampledZPImageView = circshift(undersampledZPImageView, offsetCorrection);
    % -------------------------------------------------------------------------

end



if (showResults && saveResults)
    close all
    figure;
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))))
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))))
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))))
    saveas(gcf,strcat(imgstoragePath, '/UndersampledDATAReco_Kspaceview','.jpg'));

    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))
    saveas(gcf,strcat(imgstoragePath, '/UndersampledDATAReco_ZProtatedview','.jpg'));

elseif (showResults && not(saveResults))
    close all
    figure;
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))))
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))))
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))))

    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

elseif (not(showResults) && saveResults)
    close all
    figure("Visible","off");
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))))
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))))
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))))
    saveas(gcf,strcat(imgstoragePath, '/UndersampledDATAReco_Kspaceview','.jpg'));

    figure("Visible","off");
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))
    saveas(gcf,strcat(imgstoragePath, '/UndersampledDATAReco_ZProtatedview','.jpg'));

end

disp(' ')
disp("Saving undersampled image...");
% Save .nii image in a specific folder -------------------------------------------
niftiwrite(single(abs(undersampledZPImageView)),strcat(savePathAndFileNameundersampled,'_ZP-reconstruction.nii'));
%Commented line for 4-ch because it seems not good and also shfit entailing aliasing
%niftiwrite(single(abs(fftshift(undersampledZPImage,1))),strcat(savePathAndFileNameundersampled,'_ZP-reconstruction_4-ch.nii'))
disp(strcat("Done!!! Undersampled image saved as NifTI in: ",savePathAndFileNameundersampled));
disp(' ')
%%%%%%%%%%%%%%
% Changing NIfTI header [DONE FOR 50µm SETUP]
% mrconvert step
% Constructing the shell command
pathAndfilename = strcat(savePathAndFileNameundersampled,'_ZP-reconstruction.nii');
pathAndfilenameREOR = savePathAndFileNameundersampled + '_ZP-reconstruction_REOR.nii';
cmda = sprintf('mrconvert %s %s -axes 2,0,1 -vox 0.05,0.05,0.05 -strides 3,-2,1 -force', pathAndfilename, pathAndfilenameREOR);
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
[filepath, ~, ~] = fileparts(pathAndfilename);
transformFilePath = fullfile(filepath, filename);
writematrix(T, transformFilePath, 'Delimiter', ' ');                  
disp(['Textfile "', filename, '" written successfully.']);
% Constructing the shell command 
cmdb = sprintf('mrtransform %s %s -replace %s -force', pathAndfilenameREOR, pathAndfilenameREOR, transformFilePath);
% Running the command in the terminal
[status, result] = system(cmdb);
% Optional: checking for errors
if status ~= 0
    error('mrtransform failed: %s', result);
else
    disp('mrtransform completed successfully.');
end
%%%%%%%%%%%
% --------------------------------------------------------------------------------

fileID = fopen(pathTxtDoc,'a');
fprintf(fileID,'%s',"#Reconstructed ZP image from undersampledDataReconstructionFromMultiChannelRawFID,RAW");
fprintf(fileID,'%s\r\n',strcat(savePathAndFileNameundersampled,'_ZP-reconstruction.nii'));
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Transform used with mrtransform -replace");
fprintf(fileID,'%s\r\n',transformFilePath);
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Copy, REORiented (and with proper header)");
fprintf(fileID,'%s\r\n',pathAndfilenameREOR);
fprintf(fileID,'\n');
fclose(fileID);
fprintf("\nPath of reconstructed undersampled image added in a txt file (as well as related files' paths).\n")

% Pause step 1 ---
pause(2)
% ----------------