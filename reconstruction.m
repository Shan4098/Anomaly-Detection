
%%                                                  Augmenting the Dataset

% Set the path to the input folder
input_folder = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\data\test_new\engine1_good_new';

% Set the path to the output folder
output_folder = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\data\test_new\engine1_good_new_augmented';

% Create the output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get a list of all .wav files in the input folder
wav_files = dir(fullfile(input_folder, '*.wav'));

% Parameters for time warping
warp_factor = 0.1; % Adjust the warp factor as needed

% Loop through each .wav file
for i = 1:length(wav_files)
    % Read the original audio file
    file_path = fullfile(input_folder, wav_files(i).name);
    [y, fs] = audioread(file_path);

    % Perform time warping
    y_warped = warpTime(y, warp_factor, fs);

    % Save warped file
    [~, file_name, ext] = fileparts(wav_files(i).name);
    output_path = fullfile(output_folder, [file_name '_warped' ext]);
    audiowrite(output_path, y_warped, fs);
end

% Function for time warping
function y_warped = warpTime(y, factor, fs)
    len = length(y);
    t_original = (0:len-1) / fs;

    % Generate a random warping function
    warping_function = 1 + factor * randn(1, len);

    % Interpolate the warped signal
    t_warped = cumsum(warping_function) / fs;
    y_warped = interp1(t_warped, y, t_original, 'linear', 'extrap');
end

%%                                            Combining the dataset

% 
% Specify the folder containing the .wav files
folderPath = "C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\IDMT-ISA-ELECTRIC-ENGINE\test_cut_new\engine1_good_new";

% Get a list of all .wav files in the folder
fileList = dir(fullfile(folderPath, '*.wav'));
filePaths = fullfile(folderPath, {fileList.name});

% Read and concatenate the audio files randomly
combinedAudio = [];
for i = randperm(length(filePaths))
    [audio, sampleRate] = audioread(filePaths{i});
    combinedAudio = [combinedAudio; audio];
end

% Save the combined audio to a new file
outputPath = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\IDMT-ISA-ELECTRIC-ENGINE\test_cut_new\Normal_test_3.wav';
audiowrite(outputPath, combinedAudio, sampleRate);




%% %                    Cutting the dataset into half 


% Input and output directories
inputDir = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\DATA_SET\Normal_test';
outputDir = 'C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\7551261\DATA_SET\Normal_test_2';

% Check if the output directory exists, and create it if necessary
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Get a list of all .wav files in the input directory
wavFiles = dir(fullfile(inputDir, '*.wav'));

% Process each .wav file
for i = 1:length(wavFiles)
    % Read the original .wav file
    originalFile = fullfile(inputDir, wavFiles(i).name);
    [waveform, sampleRate] = audioread(originalFile);

    % Calculate the midpoint index
    midpointIndex = floor(length(waveform) / 2);

    % Split the waveform into two halves
    waveform1 = waveform(1:midpointIndex);
    waveform2 = waveform(midpointIndex+1:end);

    % Save the two halves with modified file names
    [~, fileName, fileExt] = fileparts(originalFile);
    outputFileName1 = fullfile(outputDir, [fileName '_1' fileExt]);
    outputFileName2 = fullfile(outputDir, [fileName '_2' fileExt]);

    audiowrite(outputFileName1, waveform1, sampleRate);
    audiowrite(outputFileName2, waveform2, sampleRate);
end

disp('Splitting complete.');
