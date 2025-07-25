function [undersampledImageView] = zeroPaddingReconstruction_V3 (undersampledKSpace, correctImage, fftShiftCorrection, offsetCorrection, showResults, saveResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, rotAngleSliceX, rotAngleSliceY, rotAngleSliceZ, savePathAndFileNameUndersampled, simplerSavePath, pathTxtDoc, imgstoragePath)

% Making k-space for visualization ------------------------
undersampledKSpaceView = squeeze(bart('rss 8', undersampledKSpace));
% ---------------------------------------------------------

% Creating image ----------------------------------------
undersampledImage = bart('fft -i 7', undersampledKSpace);
undersampledImageView = bart('rss 8', undersampledImage);
undersampledImageView = squeeze(undersampledImageView);
% -------------------------------------------------------


if (correctImage)
    % Correcting image with FFTshift --------------------------------------------
    if (fftShiftCorrection(1))
        undersampledImageView = fftshift(undersampledImageView,1);
    end
    if (fftShiftCorrection(2))
        undersampledImageView = fftshift(undersampledImageView,2);
    end
    if (fftShiftCorrection(3))
        undersampledImageView = fftshift(undersampledImageView,3);
    end
    % ---------------------------------------------------------------------------
    
    
    % Correct image offset ----------------------------------------------------
    undersampledImageView = circshift(undersampledImageView, offsetCorrection);
    % -------------------------------------------------------------------------

end


if (showResults && saveResults)
    
    figure;
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))))
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))))
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))))
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);

    saveas(gcf,strcat(imgstoragePath, '/UndersampledZeroPaddingReco_Kspaceview', '.jpg'));

    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

    saveas(gcf,strcat(imgstoragePath, '/UndersampledZeroPaddingReco_rotatedview','.jpg'));

elseif (showResults)

    figure;
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))))
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))))
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))))
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);

    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

elseif (saveResults)

    figure("Visible","off");
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);
    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))))
    title(['XY plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceZ)]);
    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))))
    title(['XZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceY)]);
    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))))
    title(['YZ plane - Ch:  - Slice: ', num2str(whichVolume), ' - Slice: ', num2str(whichSliceX)]);

    saveas(gcf,strcat(imgstoragePath, '/UndersampledZeroPaddingReco_Kspaceview', '.jpg'));

    figure("Visible","off");
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))
    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))
    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

    saveas(gcf,strcat(imgstoragePath, '/UndersampledZeroPaddingReco_rotatedview','.jpg'));

end

disp(' ')
disp("Saving Zero-Padding image...");
% Save .nii image in a specific folder -------------------------------------------
% +Producing a copy that can go directly to Reco-CS-Only folder
niftiwrite((single(abs(undersampledImageView))),simplerSavePath);

niftiwrite(single(abs(undersampledImageView)),strcat(savePathAndFileNameUndersampled,'_ZeroPadding.nii'));
%ZeroPadding 4-ch could be commented but left active for further study as it's not well understood why Diego did that
niftiwrite(single(abs(fftshift(undersampledImage,1))),strcat(savePathAndFileNameUndersampled,'_ZeroPadding_4-ch.nii'))
disp(strcat("Done!!! Zero-Padding image saved as NifTI in: ",savePathAndFileNameUndersampled));
disp(' ')
%%%%%%%%%%%%%%
% Changing NIfTI header [DONE FOR 50µm SETUP]
% mrconvert step
% Constructing the shell command
pathAndfilename = strcat(savePathAndFileNameUndersampled,'_ZeroPadding.nii');
pathAndfilenameREOR = savePathAndFileNameUndersampled + 'ZeroPadding_REOR.nii';
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
disp(['File "', filename, '" written successfully.']);
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
% Copying modified NIfTI into simpler name
simplerSavePathREOR = erase(simplerSavePath,".nii") + "_ZP_REOR.nii";
cmdc = sprintf('cp %s %s', pathAndfilenameREOR, simplerSavePathREOR);
[status, result] = system(cmdc);
if status ~= 0
    error('cp failed: %s', result);
else
    disp('cp completed successfully.');
end
% --------------------------------------------------------------------------------

fileID = fopen(pathTxtDoc,'a');
fprintf(fileID,'\r%s',"#Reconstructed ZP image RAW");
fprintf(fileID,'%s\r\n\n',strcat(savePathAndFileNameUndersampled,'_ZeroPadding.nii'));
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Simpler-named copy");
fprintf(fileID,'%s\r\n',simplerSavePath);
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Transform used with mrtransform -replace");
fprintf(fileID,'%s\r\n',transformFilePath);
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Copy, REORiented (and with proper header)");
fprintf(fileID,'%s\r\n',pathAndfilenameREOR);
fprintf(fileID,'\n');
fprintf(fileID,'%s',"#Simpler-named copy, REORiented (and with proper header)");
fprintf(fileID,'%s\r\n',simplerSavePathREOR);
fprintf(fileID,'\n');
fclose(fileID);
fprintf("\nPaths of Zero-Padding related files added in a txt file.\n")

end