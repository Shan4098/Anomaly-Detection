% Specify the main directory
mainDirectory = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\IDMT-ISA-ELECTRIC-ENGINE\train\';

% Get a list of subfolders in the main directory
subfolders = dir(mainDirectory);

% Iterate through each subfolder
for folderIndex = 1:length(subfolders)
    currentFolder = subfolders(folderIndex).name;
    
    % Skip '.' and '..' folders
    if strcmp(currentFolder, '.') || strcmp(currentFolder, '..')
        continue;
    end
    
    % Full path to the current subfolder
    currentFolderPath = fullfile(mainDirectory, currentFolder);
    
    % Check if the current item is a directory
    if isfolder(currentFolderPath)
        
        % Get all .wav files in the current subfolder
        wavFiles = dir(fullfile(currentFolderPath, '*.wav'));
        
        % Iterate through .wav files in the current subfolder
        for fileIndex = 1:length(wavFiles)
            currentFile = wavFiles(fileIndex).name;
            currentFilePath = fullfile(currentFolderPath, currentFile);
            
            % Read the original audio
            [y, Fs] = audioread(currentFilePath);
            
            % Resample the audio
            y_resamp = resample(y, 16000,Fs);
            
            % Create the output directory structure
            outputFolder = fullfile(mainDirectory, [currentFolder '_new']);
            
            % Create the output directory if it doesn't exist
            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end
            
            % Specify the output file name
            [~, baseFileName, ~] = fileparts(currentFile);
            outputFilePath = fullfile(outputFolder, [baseFileName '_new.wav']);
            
            % Write the resampled audio to the new subfolder
            audiowrite(outputFilePath, y_resamp, 16000);
        end
    end
end


% currentFilePath="C:\Users\91949\Downloads\test_haoml.wav";
% [y, Fs] = audioread(currentFilePath);
% y_resamp = resample(y, 16000,Fs);
% audiowrite(currentFilePath, y_resamp, 16000);
