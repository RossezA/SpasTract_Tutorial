function [conventionalCSImageReconstructed, coilSensitivityMapsCorrected] = conventionalCSReconstruction_V5 (undersampledKSpace, lambda1Array, lambda2Array, maxItrArray, correctImage, fftShiftCorrection, offsetCorrection, B0SamplingStrategy, nB0s, centerNormX, centerNormY, centerNormZ, cubeSize, savePathAndFileName, pathTxtDoc, showResults, saveResults, whichSliceX, whichSliceY, whichSliceZ, simplerSavePath, imgstoragePath)

    close all
 
    if (B0SamplingStrategy=="UndersampledB0")    
        fprintf("\nUndersampled B0 strategy chosen...\n")
        % --------------------------------------------------------------------------------------------
        % Coil sensitivity maps
        % --------------------------------------------------------------------------------------------
        for lambda1=lambda1Array
            for lambda2=lambda2Array
                for maxItr=maxItrArray             

                    sensitivityMapsAllVolumes = zeros(size(undersampledKSpace,1),size(undersampledKSpace,2),size(undersampledKSpace,3),size(undersampledKSpace,4),size(undersampledKSpace,5));


                    for volume=1:size(undersampledKSpace,5)
                        disp(['Calculating coil sensitivity maps... - Volume: ', num2str(volume)])
                        tic
                        KOrigBart_maps = bart('ecalib -c0. -m1', squeeze(undersampledKSpace(:,:,:,:,volume)));
                        coilSensMapTime = toc
                        sensitivityMapsAllVolumes(:,:,:,:,volume) = KOrigBart_maps;

                        %KOrigBart_maps = fftshift(sensitivityMapsAllVolumes,1);


                        if (showResults && saveResults && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized');
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';
                            saveas(gcf,strcat(savePathAndFileName, '_0_Coil-sensitivity-maps_Volume-',num2str(volume),  terminusName, imageExtension));
                        
                        elseif (showResults && not(saveResults) && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized');
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';
                        
                        elseif (not(showResults) && saveResults && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized',"Visible","off");
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';
                            saveas(gcf,strcat(savePathAndFileName, '_0_Coil-sensitivity-maps_Volume-',num2str(volume),  terminusName, imageExtension));
                        end

                    end


                    coilSensitivityMapsCorrected = sensitivityMapsAllVolumes;
                    if (correctImage)
                        fprintf("\nCorrecting coil sensitivity maps using fftshift...\n")
                        if (fftShiftCorrection(1))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,1);
                        end
                        if (fftShiftCorrection(2))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,2);
                        end
                        if (fftShiftCorrection(3))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,3);                                                                
                        end


                        coilSensitivityMapsCorrected = circshift(coilSensitivityMapsCorrected, offsetCorrection);

                    end
                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------
                    niftiwrite((single(abs(coilSensitivityMapsCorrected))),strcat(savePathAndFileName,'_coilSensitivityMaps_lamb1-',num2str(lambda1),'_lamb2-',num2str(lambda2),'_maxItr-',num2str(maxItr),'.nii'))
                    % -------------------------------------------------------------------------------------------------------------------------


                    % --------------------------------------------------------------------------------------------
                    % Apply 3D Compressed Sensing reconstruction
                    % --------------------------------------------------------------------------------------------
                    conventionalCSImageReconstructed = zeros(size(undersampledKSpace,1),size(undersampledKSpace,2),size(undersampledKSpace,3),size(undersampledKSpace,5));


                    for volume=1:size(undersampledKSpace,5)
                        disp(['CS image reconstruction... - Volume: ', num2str(volume)])
                        %-------------- reco avec GPU -g ----------------
                        % commandString = strcat('pics -g -l1 -r',num2str(lambda1),' -R T:7:0:',num2str(lambda2),' -i',num2str(maxItr),' -S')
                        
                        %-------------- reco sans GPU  ----------------
                        commandString = strcat('pics -l1 -r',num2str(lambda1),' -R T:7:0:',num2str(lambda2),' -i',num2str(maxItr),' -S')
                        tic
                        csReconstructedImage = bart(commandString, squeeze(undersampledKSpace(:,:,:,:,volume)), squeeze(sensitivityMapsAllVolumes(:,:,:,:,volume)));
                        csRecoTime = toc;
                        conventionalCSImageReconstructed(:,:,:,volume) = csReconstructedImage;


                        if (showResults && saveResults)    
                            figure;
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});

                            saveas(gcf,strcat(imgstoragePath, '/UnderSampledB0_3DCompressedSensingreconstruction-YZ', num2str(whichSliceX), '-XZ', num2str(whichSliceY),'-XY',num2str(whichSliceZ),'.jpg'));
                        
                        elseif (showResults && not(saveResults))    
                            figure;
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});

                        elseif (not(showResults) && saveResults)    
                            figure("Visible","off");
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});

                            saveas(gcf,strcat(imgstoragePath, '/UnderSampledB0_3DCompressedSensingreconstruction-YZ', num2str(whichSliceX), '-XZ', num2str(whichSliceY),'-XY',num2str(whichSliceZ),'.jpg'));
                        end

                        pause(1)
                    end

                    %[csReconstructedkSpace, csReconstructedImage] =  MultiChannelCSReconstruction (undersampledKSpace, '3D-CS', 1, 1, 1, showResults, whichSliceX, whichSliceY, whichSliceZ);
                    if (correctImage)
                        fprintf("\nCorrecting offset of reconstructed image...\n")
                        % Correcting image with FFTshift -------------------------------------------
                        if (fftShiftCorrection(1))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,1);
                        end
                        if (fftShiftCorrection(2))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,2);
                        end
                        if (fftShiftCorrection(3))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,3);
                        end
                        % --------------------------------------------------------------------------

                        % Correct offset --------------------------------------------------------------------------
                        conventionalCSImageReconstructed = circshift(conventionalCSImageReconstructed, offsetCorrection);
                        % -----------------------------------------------------------------------------------------

                    end


                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------                
                    % +Producing a copy that can go directly to Reco-CS-Only folder
                    niftiwrite((single(abs(conventionalCSImageReconstructed))),simplerSavePath);
 
                    pathAndfilename = strcat(savePathAndFileName,'_conventionalCS_lamb1-',num2str(lambda1),'_lamb2-',num2str(lambda2),'_maxItr-',num2str(maxItr),'.nii');
                    niftiwrite((single(abs(conventionalCSImageReconstructed))),pathAndfilename);
                    fprintf("\nReconstructed image saved as a NIfTI file.\n")
                    %%%%%%%%%%%%%%
                    % Changing NIfTI header [DONE FOR 50µm SETUP]
                    % mrconvert step
                    % Constructing the shell command
                    pathAndfilenameREOR = erase(pathAndfilename,".nii") + "_REOR.nii";
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
                    simplerSavePathREOR = erase(simplerSavePath,".nii") + "_REOR.nii";
                    cmdc = sprintf('cp %s %s', pathAndfilenameREOR, simplerSavePathREOR);
                    [status, result] = system(cmdc);
                    if status ~= 0
                        error('cp failed: %s', result);
                    else
                        disp('cp completed successfully.');
                    end
              
                    % -------------------------------------------------------------------------------------------------------------------------


                    % Pause step 3 ---
                    pause(3)
                    % ----------------

                    fileID = fopen(pathTxtDoc,'a');
                    fprintf(fileID,'%s',"#Reconstructed CS image RAW");
                    fprintf(fileID,'%s\r\n',pathAndfilename);
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
                    fprintf("\nPaths of CS reconstruction related files added in a txt file.\n")
                end
            end
        end
        
        
        
    elseif (B0SamplingStrategy=="FullySampledB0")
        fprintf("\nFully sampled B0 strategy chosen...\n")
        % --------------------------------------------------------------------------------------------
        % Coil sensitivity maps
        % --------------------------------------------------------------------------------------------
        for lambda1=lambda1Array
            for lambda2=lambda2Array
                for maxItr=maxItrArray             

                    sensitivityMapsAllVolumes = zeros(size(undersampledKSpace,1),size(undersampledKSpace,2),size(undersampledKSpace,3),size(undersampledKSpace,4),size(undersampledKSpace,5));


                    for volume=1:size(undersampledKSpace,5)
                        disp(['Calculating coil sensitivity maps... - Volume: ', num2str(volume)])
                        tic
                        KOrigBart_maps = bart('ecalib -c0. -m1', squeeze(undersampledKSpace(:,:,:,:,volume)));
                        coilSensMapTime = toc
                        sensitivityMapsAllVolumes(:,:,:,:,volume) = KOrigBart_maps;

                        %KOrigBart_maps = fftshift(sensitivityMapsAllVolumes,1);

                        % Show only the first and the last coil sensitivity maps, because they are almost the same ---------------------------------------
                        if (showResults && saveResults && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized');
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';                            
                            
                            saveas(gcf,strcat(savePathAndFileName, '_0_Coil-sensitivity-maps_Volume-',num2str(volume),  terminusName, imageExtension));
                        
                        elseif (showResults && not(saveResults) && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized');
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)                          
                            
                        elseif (not(showResults) && saveResults && ((volume==1 || volume==(size(undersampledKSpace,5)))))

                            figure('WindowState','maximized',"Visible","off");
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,1))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,2))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,3))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(KOrigBart_maps(whichSliceX,:,:,4))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,1))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,2))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,3))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(KOrigBart_maps(:,whichSliceY,:,4))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,1))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,2))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,3))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(KOrigBart_maps(:,:,whichSliceZ,4))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';                            
                            
                            saveas(gcf,strcat(savePathAndFileName, '_0_Coil-sensitivity-maps_Volume-',num2str(volume),  terminusName, imageExtension));    
                        end
                        %------------------------------------------------------------------------------------------------------------------------------
                    end


                    coilSensitivityMapsCorrected = sensitivityMapsAllVolumes;
                    if (correctImage)
                        fprintf("\nCorrectig coil sensitivity maps using fftshift...\n")
                        if (fftShiftCorrection(1))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,1);
                        end
                        if (fftShiftCorrection(2))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,2);
                        end
                        if (fftShiftCorrection(3))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,3);                                                                
                        end


                        coilSensitivityMapsCorrected = circshift(coilSensitivityMapsCorrected, offsetCorrection);

                    end
                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------
                    niftiwrite((single(abs(coilSensitivityMapsCorrected))),strcat(savePathAndFileName,'_coilSensitivityMaps_lamb1-',num2str(lambda1),'_lamb2-',num2str(lambda2),'_maxItr-',num2str(maxItr),'.nii'))
                    % -------------------------------------------------------------------------------------------------------------------------


                    % --------------------------------------------------------------------------------------------
                    % Apply 3D Compressed Sensing reconstruction
                    % --------------------------------------------------------------------------------------------
                    conventionalCSImageReconstructed = zeros(size(undersampledKSpace,1),size(undersampledKSpace,2),size(undersampledKSpace,3),size(undersampledKSpace,5));


                    for volume=1:size(undersampledKSpace,5)
                        disp(['CS image reconstruction... - Volume: ', num2str(volume)])
                        commandString = strcat('pics -g -l1 -r',num2str(lambda1),' -R T:7:0:',num2str(lambda2),' -i',num2str(maxItr),' -S');
                        tic
                        csReconstructedImage = bart(commandString, squeeze(undersampledKSpace(:,:,:,:,volume)), squeeze(sensitivityMapsAllVolumes(:,:,:,:,volume)));
                        csRecoTime = toc;
                        conventionalCSImageReconstructed(:,:,:,volume) = csReconstructedImage;


                        if (showResults && saveResults)    
                            figure;
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                            saveas(gcf,strcat(imgstoragePath, '/FullSampledB0_3DCompressedSensingreconstruction-YZ', num2str(whichSliceX), '-XZ', num2str(whichSliceY),'-XY',num2str(whichSliceZ),'.jpg'));
                        
                        elseif (showResults && not(saveResults))    
                            figure;
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});

                        elseif (not(showResults) && saveResults)    
                            figure("Visible","off");
                            subplot(1,3,1)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(conventionalCSImageReconstructed(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                            saveas(gcf,strcat(imgstoragePath, '/FullSampledB0_3DCompressedSensingreconstruction-YZ', num2str(whichSliceX), '-XZ', num2str(whichSliceY),'-XY',num2str(whichSliceZ),'.jpg'));    
                        end

                        pause(1)
                    end

                    
                    
                    
                    % Replacing B0 volumes reconstructed by only FFT transformed ones (there is no undersampling on B) volumes) ---
                    % Creating image ----------------------------------------
                    B0FullySampledImages = bart('fft -i 7', squeeze(undersampledKSpace(:,:,:,:,1:nB0s)));
                    B0FullySampledImageView = bart('rss 8', B0FullySampledImages);
                    B0FullySampledImageView = squeeze(B0FullySampledImageView);
                    % -------------------------------------------------------
                    
                    
                    %[csReconstructedkSpace, csReconstructedImage] =  MultiChannelCSReconstruction (undersampledKSpace, '3D-CS', 1, 1, 1, showResults, whichSliceX, whichSliceY, whichSliceZ);
                    if (correctImage)
                        fprintf("\nCorrecting offset of reconstructed image...\n")
                        % Correcting image with FFTshift -------------------------------------------
                        if (fftShiftCorrection(1))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,1);
                            B0FullySampledImageView = fftshift(B0FullySampledImageView,1);
                        end
                        if (fftShiftCorrection(2))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,2);
                            B0FullySampledImageView = fftshift(B0FullySampledImageView,2);
                        end
                        if (fftShiftCorrection(3))
                            conventionalCSImageReconstructed = fftshift(conventionalCSImageReconstructed,3);
                            B0FullySampledImageView = fftshift(B0FullySampledImageView,3);
                        end
                        % --------------------------------------------------------------------------

                        % Correct offset --------------------------------------------------------------------------
                        conventionalCSImageReconstructed = circshift(conventionalCSImageReconstructed, offsetCorrection);
                        B0FullySampledImageView = circshift(B0FullySampledImageView, offsetCorrection);
                        % -----------------------------------------------------------------------------------------

                    end                                        
                                                            
                    [csReconstructedImageOnlyNormalizedB0, B0RegularizationVector] =  imageNormalizationWrt3DReference(conventionalCSImageReconstructed(:,:,:,1:nB0s), B0FullySampledImageView, centerNormX, centerNormY, centerNormZ, cubeSize);
                    conventionalCSImageReconstructed(:,:,:,1:nB0s) = csReconstructedImageOnlyNormalizedB0;
                    % -------------------------------------------------------------------------------------------------------------

                    % Save .nii image in a specific folder and then modify with proper properties ------------------------------------------------------------------------------------                
                    % +Producing a copy that can go directly to Reco-CS-Only folder
                    niftiwrite((single(abs(conventionalCSImageReconstructed))),simplerSavePath);

                    pathAndfilename = strcat(savePathAndFileName,'_conventionalCS_lamb1-',num2str(lambda1),'_lamb2-',num2str(lambda2),'_maxItr-',num2str(maxItr),'.nii');
                    niftiwrite((single(abs(conventionalCSImageReconstructed))),pathAndfilename);
                    fprintf("\nReconstructed image saved as a NIfTI file.\n")
                    %%%%%%%%%%%%%%
                    % Changing NIfTI header [DONE FOR 50µm SETUP]
                    % mrconvert step
                    % Constructing the shell command
                    pathAndfilenameREOR = erase(pathAndfilename,".nii") + "_REOR.nii";
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
                    simplerSavePathREOR = erase(simplerSavePath,".nii") + "_REOR.nii";
                    cmdc = sprintf('cp %s %s', pathAndfilenameREOR, simplerSavePathREOR);
                    [status, result] = system(cmdc);
                    if status ~= 0
                        error('cp failed: %s', result);
                    else
                        disp('cp completed successfully.');
                    end
                    
                    % -------------------------------------------------------------------------------------------------------------------------


                    % Pause step 3 ---
                    pause(3)
                    % ----------------

                    fileID = fopen(pathTxtDoc,'a');
                    fprintf(fileID,'%s',"#Reconstructed CS image RAW");
                    fprintf(fileID,'%s\r\n',pathAndfilename);
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
                    fprintf("\nPaths of CS reconstruction related files added in a txt file.\n")
                end
            end
        end   
        
    end
end