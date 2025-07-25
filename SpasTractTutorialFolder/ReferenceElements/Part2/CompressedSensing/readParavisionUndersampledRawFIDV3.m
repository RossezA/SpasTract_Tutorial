function [kSpaceReconstructed, imageZPReconstructed, csMaskRetrieved] =  readParavisionUndersampledRawFIDV3 (fidPath, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, saveResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, imgstoragePath)

    % Split full path into parts
    parts = strsplit(fidPath, filesep);
    if length(parts) >= 1
        % Get last folder
        FIDChannel = parts{end};
    
    switch acquisitionType
        case 'CS-3D-FLASH'
            disp('Compressed Sensing 3D Flash selected')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+1i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires
            
            if (showResults && saveResults)
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_BrukerFID_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_BrukerFID_', FIDChannel,'.jpg'));
            end

           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim);
            
            % Open the file
            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            % Initialize a cell array to hold values (in case some are strings)
            allCell = {};
            
            %Scan the undersampling pattern file and place the content into a column vector
            % Read lines one by one
            while ~feof(fileID)
                line = strtrim(fgetl(fileID));
                if isempty(line) || startsWith(strtrim(line), '#')
                    continue;
                end

                    % Check if it contains a colon (e.g. "AcquisitionTypeID : 1")
                if contains(line, ':')
                    % Split on the colon and take the part after it
                    parts = strsplit(line, ':');
                    valStr = strtrim(parts{2});
                else
                    % Otherwise, treat the whole line as the value
                    valStr = line;
                end

                % Convert to numeric or logical if possible
                valNum = str2double(valStr);
                if ~isnan(valNum)
                    value = valNum;
                elseif strcmpi(valStr, 'true')
                    value = true;
                elseif strcmpi(valStr, 'false')
                    value = false;
                else
                    value = valStr; % Keep as string
                end

                % Store in list
                allCell{end+1} = value;
            end
            
            % Close the file
            fclose(fileID);
            
            allData = cell2mat(allCell);

            whichCase = allData(1); % useless ID that is not checked here to see which case is followed thereafter i.e "CS-3D-FLASH" 
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nVolumes = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12; %?          
            
            
            %% Data organization
            disp('ImagData size')
            size(imagData)            
            origReadoutDim
            disp('ImagData size / ReadOut dim')
            size(imagData,1)/origReadoutDim
            totalNumberOfLines
            originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines);
            
            % If fewer lines were acquired -----------------------------------------------------------
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, []);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2)) = originalOrganizedDataTemp;
            % ----------------------------------------------------------------------------------------           
            
            %imagData = imagData(1:(origReadoutDim*totalNumberOfLines))
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            
            if (showResults && saveResults)
                figure;
                imagesc(abs(squeeze(originalOrganizedData)).^0.25);
                title('Original lines acquired');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_Originallinesacquired_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                imagesc(abs(squeeze(originalOrganizedData)).^0.25);
                title('Original lines acquired');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                imagesc(abs(squeeze(originalOrganizedData)).^0.25);
                title('Original lines acquired');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_Originallinesacquired_', FIDChannel,'.jpg'));
            end
            
            
            %organizedData = originalOrganizedData(1:finalReadoutDim,:);

            posData = 1;
            globalPos = posInAllDataOfFirstSlice;
            for volume=1:nVolumes
                
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);

                        
            if (showResults && saveResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});             
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ))),[]);
                title({'Zero-filling image','X-Y plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-XY_', FIDChannel,'.jpg'));

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-XZ_', FIDChannel,'.jpg'));

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-YZ_', FIDChannel,'.jpg'));

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:)));
                title('Compressed Sensing mask retrieved');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_Mask_', FIDChannel,'.jpg'));
             
            elseif (showResults && not(saveResults))
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});             
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ))),[]);
                title({'Zero-filling image','X-Y plan'});

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:))),[]);
                title({'Zero-filling image','X-Z plan'});

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:))),[]);
                title({'Zero-filling image','Y-Z plan'});

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:)));
                title('Compressed Sensing mask retrieved');

            elseif (not(showResults) && saveResults)
                
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});             
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ))),[]);
                title({'Zero-filling image','X-Y plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-XY_', FIDChannel,'.jpg'));

                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-XZ_', FIDChannel,'.jpg'));

                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_UnderKspace-YZ_', FIDChannel,'.jpg'));

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure("Visible","off");
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:)));
                title('Compressed Sensing mask retrieved');
                saveas(gcf,strcat(imgstoragePath, '/3DCSFlash_Mask_', FIDChannel,'.jpg'));
                
            end
        
        
        case 'CS-2D-FLASH'
            disp('Compressed Sensing 2D-MS Flash selected but Work in Progress...')
        
        
        case 'CS-3D-SE-DTI_DWLoopFirst'
            disp('Compressed Sensing 3D DTI Spin Echo - Same undersampling pattern for all 3D volumes (Mono3D masks) - DW loop as the first one')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+1i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            if (showResults && saveResults)
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BrukerFID_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BrukerFID_', FIDChannel,'.jpg'));
            end


           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            % Open the file
            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            % Initialize a cell array to hold values (in case some are strings)
            allCell = {};
            
            %Scan the undersampling pattern file and place the content into a column vector
            % Read lines one by one
            while ~feof(fileID)
                line = strtrim(fgetl(fileID));
                if isempty(line) || startsWith(strtrim(line), '#')
                    continue;
                end

                    % Check if it contains a colon (e.g. "AcquisitionTypeID : 1")
                if contains(line, ':')
                    % Split on the colon and take the part after it
                    parts = strsplit(line, ':');
                    valStr = strtrim(parts{2});
                else
                    % Otherwise, treat the whole line as the value
                    valStr = line;
                end

                % Convert to numeric or logical if possible
                valNum = str2double(valStr);
                if ~isnan(valNum)
                    value = valNum;
                elseif strcmpi(valStr, 'true')
                    value = true;
                elseif strcmpi(valStr, 'false')
                    value = false;
                else
                    value = valStr; % Keep as string
                end

                % Store in list
                allCell{end+1} = value;
            end
            
            % Close the file
            fclose(fileID);
            
            allData = cell2mat(allCell);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            numberOfSlices = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12;
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            originalOrganizedData = permute(originalOrganizedData,[1 3 2]);
                        
            if (showResults && saveResults)
                figure;
                imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
                title('Original lines acquired');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_Originallinesacquired_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
                title('Original lines acquired');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
                title('Original lines acquired');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_Originallinesacquired_', FIDChannel,'.jpg'));
            end


            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            for volume=1:nVolumes
                posData = 1;
                globalPos = posInAllDataOfFirstSlice;
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData,volume);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);
        
            
            if (showResults && saveResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-XY_', FIDChannel,'.jpg'));
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-XZ_', FIDChannel,'.jpg'));

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-YZ_', FIDChannel,'.jpg'));

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_MaskYZVol1_', FIDChannel,'.jpg'));
                
            elseif (showResults && not(saveResults))
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});

                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                
            elseif (not(showResults) && saveResults)
                
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-XY_', FIDChannel,'.jpg'));
                
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-XZ_', FIDChannel,'.jpg'));

                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_UnderKspace-YZ_', FIDChannel,'.jpg'));

                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure("Visible","off");
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_MaskYZVol1_', FIDChannel,'.jpg'));

            end
            
        
        case 'CS-3D-SE-DTI_DWLoopLast'
            disp('Compressed Sensing 3D DTI Spin Echo - One undersampling pattern per 3D volume (Multi3D masks) - DW loop as the last one')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+1i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires
            
            if (showResults && saveResults)
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_BrukerFID_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_BrukerFID_', FIDChannel,'.jpg'));
            end


           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            % Open the file
            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            % Initialize a cell array to hold values (in case some are strings)
            allCell = {};
            
            %Scan the undersampling pattern file and place the content into a column vector
            % Read lines one by one
            while ~feof(fileID)
                line = strtrim(fgetl(fileID));
                if isempty(line) || startsWith(strtrim(line), '#')
                    continue;
                end

                    % Check if it contains a colon (e.g. "AcquisitionTypeID : 1")
                if contains(line, ':')
                    % Split on the colon and take the part after it
                    parts = strsplit(line, ':');
                    valStr = strtrim(parts{2});
                else
                    % Otherwise, treat the whole line as the value
                    valStr = line;
                end

                % Convert to numeric or logical if possible
                valNum = str2double(valStr);
                if ~isnan(valNum)
                    value = valNum;
                elseif strcmpi(valStr, 'true')
                    value = true;
                elseif strcmpi(valStr, 'false')
                    value = false;
                else
                    value = valStr; % Keep as string
                end

                % Store in list
                allCell{end+1} = value;
            end
            
            % Close the file
            fclose(fileID);
            
            allData = cell2mat(allCell);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nVolumes = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12;
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            %originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            %originalOrganizedData = reshape(imagData(1:(end-origReadoutDim)), origReadoutDim, totalNumberOfLines);
           
            %originalOrganizedData = permute(originalOrganizedData,[1 3 2]);
            
            % if (showResults && saveResults)
            %     figure;
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            %     saveas(gcf,strcat(imgstoragePath, '/DCS-SE-DTI_DWLoopLast_Originallinesacquired_', FIDChannel,'.jpg'));
            % elseif (showResults && not(saveResults))
            %     figure;
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            % elseif (not(showResults) && saveResults)
            %     figure("Visible","off");
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            %     saveas(gcf,strcat(imgstoragePath, '/DCS-SE-DTI_DWLoopLast_Originallinesacquired_', FIDChannel,'.jpg'));
            % end

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            

            posData = 1;
            globalPos = posInAllDataOfFirstSlice;
            for volume=1:nVolumes
                
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);
            
            if (showResults && saveResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});                
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-XY_', FIDChannel,'.jpg'));
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-XZ_', FIDChannel,'.jpg'));
                              
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-YZ_', FIDChannel,'.jpg'));
        
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_MaskYZVol1_', FIDChannel,'.jpg'));
                
            elseif (showResults && not(saveResults))
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                              
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
        
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
            
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});                
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-XY_', FIDChannel,'.jpg'));
                
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-XZ_', FIDChannel,'.jpg'));
                              
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_UnderKspace-YZ_', FIDChannel,'.jpg'));
        
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure("Visible","off");
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopLast_MaskYZVol1_', FIDChannel,'.jpg'));    
            end
            
            
        case 'CS-3D-SE-DTI_DWLoopFirst_BinaryGuideArray'
            disp('Compressed Sensing 3D DTI Spin Echo - One undersampling pattern per 3D volume (Multi3D masks) - DW loop as first one with binary guide array')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+1i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            if (showResults && saveResults)
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_BrukerFID_', FIDChannel,'.jpg'));
            elseif (showResults && not(saveResults))
                figure;
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
            elseif (not(showResults) && saveResults)
                figure("Visible","off");
                plot(1:size(imagData,1),abs(imagData));
                title('Bruker FID');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_BrukerFID_', FIDChannel,'.jpg'));
            end


           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            % Open the file
            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            % Initialize a cell array to hold values (in case some are strings)
            allCell = {};
            
            %Scan the undersampling pattern file and place the content into a column vector
            % Read lines one by one
            while ~feof(fileID)
                line = strtrim(fgetl(fileID));
                if isempty(line) || startsWith(strtrim(line), '#')
                    continue;
                end

                    % Check if it contains a colon (e.g. "AcquisitionTypeID : 1")
                if contains(line, ':')
                    % Split on the colon and take the part after it
                    parts = strsplit(line, ':');
                    valStr = strtrim(parts{2});
                else
                    % Otherwise, treat the whole line as the value
                    valStr = line;
                end

                % Convert to numeric or logical if possible
                valNum = str2double(valStr);
                if ~isnan(valNum)
                    value = valNum;
                elseif strcmpi(valStr, 'true')
                    value = true;
                elseif strcmpi(valStr, 'false')
                    value = false;
                else
                    value = valStr; % Keep as string
                end

                % Store in list
                allCell{end+1} = value;
            end
            
            % Close the file
            fclose(fileID);
            
            allData = cell2mat(allCell);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nB0s = allData(5);
            nDiffusionDirections = allData(6);
            nVolumes = nB0s + nDiffusionDirections;
            accelerationFactor = allData(7);
            centerSquareSize = allData(8);

            variableDensityOption = allData(9);
            ellipseOption = allData(10);
            seedValue = allData(11);

            totalNumberOfLines = allData(12);

            posInAllDataOfFirstSlice = 13; % et là 13 au lieu de 12 ?
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            %originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            %originalOrganizedData = reshape(imagData(1:(end-origReadoutDim)), origReadoutDim, totalNumberOfLines);
           
            %originalOrganizedData = permute(originalOrganizedData,[1 3 2]);

            % if (showResults && saveResults)
            %     figure;
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            %     saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_Originallinesacquired_', FIDChannel,'.jpg'));
            % elseif (showResults && not(saveResults))
            %     figure;
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            % elseif (not(showResults) && saveResults)
            %     figure("Visible","off");
            %     imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %     title('Original lines acquired');
            %     saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_Originallinesacquired_', FIDChannel,'.jpg'));
            % end

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            
            globalPos = posInAllDataOfFirstSlice;
            posData = 1;
            for whichSlice=1:finalPhase2Dim
                
                for whichLine=1:finalPhase1Dim
                    
                    for volume=1:nVolumes
                        
                        acquisitionMade = allData(globalPos);
                        
                        if (acquisitionMade==1)
                            
                            originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                            %fprintf("\n Row %d of slice %d of diffusion volume %d acquired.\n", whichLine, whichSlice, volume);
                            posData = posData + 1;
                        end
                        
                        
                        globalPos = globalPos + 1;
                    end
                        
                end
                
            end                                              
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);


            
             if (showResults && saveResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});             
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-XY_', FIDChannel,'.jpg'));
     
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-XZ_', FIDChannel,'.jpg'));
                               
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-YZ_', FIDChannel,'.jpg'));
         
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) ~= 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_MaskYZVol1_', FIDChannel,'.jpg'));
                
             elseif (showResults && not(saveResults))
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});             
     
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                               
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
         
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) ~= 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
            
             elseif (not(showResults) && saveResults)
                
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});             
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-XY_', FIDChannel,'.jpg'));
     
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-XZ_', FIDChannel,'.jpg'));
                               
                figure("Visible","off");
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});
                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_UnderKspace-YZ_', FIDChannel,'.jpg'));
         
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) ~= 0);
                figure("Visible","off");
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                saveas(gcf,strcat(imgstoragePath, '/3DCS-SE-DTI_DWLoopFirst_BinaryGuideArray_MaskYZVol1_', FIDChannel,'.jpg'));
                    
            end
            
            
       
        
        otherwise
            disp('No one selected')
            kSpaceReconstructed = NaN; 
            imageZPReconstructed = NaN;
            csMaskRetrieved = NaN;
    end



end