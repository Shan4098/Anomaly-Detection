% Read the .wav file
[audio, sampleRate] = audioread("C:\Users\91949\OneDrive\Desktop\Semester\HAOML\Project\DATA\IDMT-ISA-ELECTRIC-ENGINE\16_KHZ\train_cut_new\engine1_good_new\pure_0_new.wav");

% Normalize the audio data
audio = audio / max(abs(audio));

% Create and open the .h file
hFile = fopen('audio_data.h', 'w');

% Write array declaration and sample rate
fprintf(hFile, 'const int audioLength = %d;\n', length(audio));
fprintf(hFile, 'const int sampleRate = %d;\n', sampleRate);
fprintf(hFile, 'float audioData[] = {\n'); % Remove const

% Write the audio samples
fprintf(hFile, '  %f', audio(1));
for i = 2:length(audio)
    fprintf(hFile, ', %f', audio(i));
end

fprintf(hFile, '};\n');

% Close the .h file
fclose(hFile);

disp('Conversion complete.');
