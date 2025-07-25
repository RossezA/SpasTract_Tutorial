function [acquisitionTypeID, inputReadoutDim, inputPhase1Dim, inputPhase2Dim, nEchos, outputReadoutDim, outputPhase1Dim, outputPhase2Dim, nChannels, nB0s, nDWI, nVolumes] = readAcquisitionFileV2(txtDocPath)
    
    % Open the file
    fileID = fopen(txtDocPath,'r');
    
    % Initialize a cell array to hold values (in case some are strings)
    allCell = {};
    %Scan the acquisition parameter file and place the content into a column vector
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
    
    %Attributing each line of allData to a parameter
    acquisitionTypeID = allData(1); % ID that is checked in the main script to see if the acquisition type is i.e "3D-FLASH" if ==1 or "3D-SE-DTI_DWLoopFirst" with ==2 or ... 
    inputReadoutDim= allData(2);
    inputPhase1Dim= allData(3);
    inputPhase2Dim = allData(4); 
    nEchos = allData(5);
    outputReadoutDim = allData(6);
    outputPhase1Dim = allData(7);
    outputPhase2Dim = allData(8);
    nChannels = allData(9);
    nB0s = allData(10);              
    nDWI = allData(11);
    nVolumes = nB0s + nDWI;
    
end
