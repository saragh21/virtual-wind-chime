% Step 1: Read the Audio File
[y, Fs] = audioread('wood.wav');
% Consider only the first 2 seconds for analysis
y = y(1:2*Fs);

% Generate a high-resolution spectrogram
nfft = 4096; % Increase or decrease depending on the resolution required
window = 4096; % Window size for the spectrogram
overlap = round(nfft * 0.9); % Overlap size
[S, F, T, P] = spectrogram(y, window, overlap, nfft, Fs);

% Convert power to decibels
PdB = 10*log10(abs(S)+eps); 

% Plot the spectrogram
figure;
spectrogram(y, window, overlap, nfft, Fs, 'yaxis');
title('Spectrogram of the First Two Seconds');
colormap('gray'); % This is a common colormap for spectrograms
colorbar; % Show the color scale

% Define the segment times for each tube hit
segment_times = [0.0, 0.2; 0.21, 0.3; 0.31, 0.4; 0.41, 0.59; 0.81, 0.89]; 
num_tubes = size(segment_times, 1);

% Initialize a table to store mode frequencies
mode_freqs = zeros(num_tubes, 5); 



% Analyze each segment to find the mode frequencies
for i = 1:num_tubes
    % Get the time indices for the current tube segment
    time_indices = T >= segment_times(i, 1) & T <= segment_times(i, 2);
    
    % Sum the power over time to get the frequency profile for this segment
    segment_power = sum(PdB(:, time_indices), 2);
    
    % DEBUG: Print the max power in the segment to check if there's signal
    fprintf('Max power in segment %d: %f dB\n', i, max(segment_power));
    
    % Find the peaks in the segment power spectrum
    [peaks, locs] = findpeaks(segment_power, 'MinPeakDistance', Fs/500); 
    
    % DEBUG: Print the number of peaks found
    fprintf('Number of peaks found in segment %d: %d\n', i, length(locs));
    
    % Ensure locs are within the bounds of F
    valid_locs = locs(locs <= length(F));
    
    % Get the frequency values of the peaks (up to the number of desired modes)
    num_modes_to_detect = min(length(valid_locs), size(mode_freqs, 2));
    if ~isempty(valid_locs) && num_modes_to_detect > 0
        mode_freqs(i, 1:num_modes_to_detect) = F(valid_locs(1:num_modes_to_detect));
    end
end

% Convert mode frequencies to a table for display
mode_table = array2table(mode_freqs, 'VariableNames', {'Mode1', 'Mode2', 'Mode3', 'Mode4', 'Mode5'}, ...
    'RowNames', {'Tube1', 'Tube2', 'Tube3', 'Tube4', 'Tube5'});

% Display the table
disp(mode_table);

% Number of points for FFT, use a power of 2 for speed
N = 2^nextpow2(length(y));
fft_result = fft(y, N);
fft_magnitude = abs(fft_result);
fft_freq = linspace(0, Fs, N);

% Preallocate the gains matrix
mode_gains = zeros(size(mode_freqs));

for tube = 1:size(mode_freqs, 1)
    for mode = 1:size(mode_freqs, 2)
        target_freq = mode_freqs(tube, mode);

        % Find the closest FFT bin to the target frequency
        [~, closest_bin] = min(abs(fft_freq - target_freq));

        % Take the magnitude of the FFT at this bin and its neighbors
        mag = fft_magnitude(closest_bin-1:closest_bin+1);

        % Frequencies for the bins
        freqs = fft_freq(closest_bin-1:closest_bin+1);

        % Fit a parabola to the three points (mag)
        p = polyfit(freqs - target_freq, mag, 2);

        % Vertex of the parabola gives the gain at the peak
        mode_gains(tube, mode) = -p(2) / (2*p(1));
    end
end

T60 = [2, 1, 0.4, 0.6, 0.16]; % Example T60 times in seconds for each mode

% Preallocate array for pole radii
R = zeros(size(T60));

% Calculate R for each T60
for i = 1:length(T60)
    R(i) = 10^(-6/(T60(i) * Fs));
end

% Display the pole radii
disp('Pole Radii:');
disp(R);


% Display the gains matrix
disp(mode_gains);
%soundsc(y,Fs);

duration = 25; % Total duration of the simulation in seconds
wind_force = simulate_wind_force(Fs, duration);

Rdecay = 0.9999; % Decay rate

% Update system energy
[E, dn] = update_system_energy(Fs, duration, Rdecay);

c = 19;

states = simulate_clapper_movement(E,c);


tube_modes = [...
    64.6,    1055.1, 2659.4, 3617.6, 4974.2;  % Tube1
    592.16,  1571.9, 2960.8, 4199,   5437.1;  % Tube2
    592.16,  1571.9, 3456.1, 5351,   6632.2;  % Tube3
    592.16,  1571.9, 2896.2, 3886.7, 5254.1;  % Tube4
    796.73,  1755,   2734.7, 4252.8, 5221.8    % Tube5
];

gains = [0.0550, 0.0877, 0.0719, 0.0586, 0.0331];  % Corrected values

pole_radii = [0.9998, 0.9997, 0.9992, 0.9995, 0.9980];

Rinput = 0.97;

a = sqrt(E + 0.1);

%final_sound = synthesize_tube_sound(Fs, duration, states, tube_modes, gains, pole_radii, Rinput, a);

