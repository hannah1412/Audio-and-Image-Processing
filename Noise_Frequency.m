clc; 
close all;
clear all;

% Processing audio 1
filename = 'audio_in_noise3.wav';

[s_audio, s_rate] = audioread(filename); 
time = (0:length(s_audio)-1) / s_rate;  %time vector

% % %STEP 1:  Visualising audio in Time domain 
figure;
subplot(2,1,1);
plot(time, s_audio);
title('Audio_2 in Time domain (unfiltered)');
xlabel('Time (s)'),ylabel('Amplitude');
grid on;

n_samples = length(s_audio);
nfft = 2^nextpow2(n_samples);       %doule the amouon tof sample for reducing alias

% STEP 2: FFT with 0-paddings
fft_audio = fft(s_audio, nfft);
freq = (0:nfft-1) * (s_rate/nfft); %frequency vector

amplitude = abs(fft_audio)/ n_samples;
amplitude = amplitude(1:nfft/2+1);      %single sided 
freq = freq(1:nfft/2+1);

% % % Visualising Frequency domain 
figure, plot(freq, amplitude);
title('Audio1 signal in Frequency domain');
xlabel('Frequency (Hz)'), ylabel('|fft|');

% STEP 3: Find dominant frequency-> the noise 
max_amplitude = max(amplitude);
threshold = max_amplitude * 0.9;
[peak_amp, peak_loc] = findpeaks(amplitude, freq, 'MinPeakHeight', threshold, 'SortStr', 'descend', 'NPeaks',2);

disp(['Result for audio 3']);
disp(['Dominant Frequency : ', num2str(peak_loc), 'Hz']);
disp(['Dominant Ampitude: ', num2str(peak_amp)]);

% % % Visualising Dominant frequency -> noise 
figure, plot(freq, amplitude, 'b');
hold on;
plot(peak_loc, peak_amp, 'ro', 'MarkerSize',10, 'DisplayName', 'Dominant Freqency');
hold off; 
title('Dominant Frequency in Frequency domain(audio1)');
xlabel('Frequency (Hz)'), ylabel('|fft(x)|');


% STEP 4 - design Notch filter
filtered_audio = s_audio;
quality = 3;

for i = 1:length(peak_loc)
    notch_filter = peak_loc(i) /(s_rate/2);
    [b, a ] = iirnotch(notch_filter, notch_filter/quality);
    filtered_audio = filter(b,a,filtered_audio);
end

% % % Visualising filtered audio 
figure;
subplot(2, 1, 1), plot(time, s_audio), title('Audio in Time domain (unfiltered)');
xlabel('Time (s)'),ylabel('Amplitude');

subplot(2, 1, 2), plot(time, filtered_audio), title('Audio in Time domain (FILTERED)');
xlabel('Time (s)'),ylabel('Amplitude');

sound(filtered_audio, s_rate);
% Write to new files
audiowrite('filtered_audio3.wav', filtered_audio, s_rate);

% % % Visualising the Filtered audio in Frequency domain for further
% comparison
fft_filtered = fft(filtered_audio, nfft);
amplitude_filtered = abs(fft_filtered) / n_samples;
figure, subplot(2, 1, 1);
plot(freq, amplitude), title('OG audio in Frequency domain')
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
plot(freq(1:nfft/2+1), amplitude_filtered(1:nfft/2+1));
title('Filtered Audio Signal in Frequency domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

