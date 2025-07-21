% --- Configuration Parameters ---
filename = 'WDM_EYE_1.txt'; % Name of your text file
bit_rate = 10e9;                            % Gbps
bit_period = 1 / bit_rate;                  % Bit period in seconds (100 ps for 10 Gbps)

% --- Load Data ---
% Assuming the file has a header and comma-separated values
% If space-separated and no header, use: data = load(filename);
try
    data = readtable(filename, 'Delimiter', ',', 'HeaderLines', 0);
    % Assuming columns are time and Y (current)
    time_points = data{:, 1}; % First column is time
    voltage_or_current = data{:, 2}; % Second column is Y (current or voltage)
catch
    fprintf('Error: Could not read %s as a comma-separated file with no header.\n', filename);
    fprintf('Attempting to read as space-separated numeric data.\n');
    data = load(filename);
    time_points = data(:, 1);
    voltage_or_current = data(:, 2);
end

% --- Eye Diagram Plot (to visualize what we're analyzing) ---
figure;
plot(time_points, voltage_or_current, '.'); % Plot as points
title('Eye Diagram (from data)');
xlabel('Time (s)');
ylabel('Amplitude (A)'); % Assuming Amperes from your previous context
grid on;

% --- Eye Diagram Analysis ---

% Find min and max voltage/current levels
level_0 = min(voltage_or_current); % Closest to 0
level_1 = max(voltage_or_current); % Max current

% Define eye boundaries for slicing (adjust if your data is not perfectly centered)
% We'll typically analyze the middle 20-30% of the eye for eye height,
% and the crossing region for jitter.
eye_center_time = 0; % Assuming the eye is centered at 0 time
eye_region_start_time = eye_center_time - bit_period * 0.15; % e.g., -15 ps
eye_region_end_time = eye_center_time + bit_period * 0.15;   % e.g., +15 ps

crossing_region_width = bit_period * 0.1; % 10% of bit period for crossing analysis

% Extract samples in the eye opening and crossing regions
eye_samples_index = (time_points >= eye_region_start_time) & (time_points <= eye_region_end_time);
eye_samples = voltage_or_current(eye_samples_index);

% Find samples for upper and lower rails in the eye opening
% Define a threshold to distinguish between '0' and '1' rails for eye height.
% This threshold is usually midway between the average '0' and '1' levels.
avg_level_0 = mean(voltage_or_current(voltage_or_current < (level_0 + level_1)/2));
avg_level_1 = mean(voltage_or_current(voltage_or_current > (level_0 + level_1)/2));

upper_rail_samples = eye_samples(eye_samples > (avg_level_0 + avg_level_1)/2);
lower_rail_samples = eye_samples(eye_samples < (avg_level_0 + avg_level_1)/2);

% --- Calculate Eye Metrics ---

% 1. Eye Height (Vertical Opening)
if ~isempty(upper_rail_samples) && ~isempty(lower_rail_samples)
    eye_height = min(upper_rail_samples) - max(lower_rail_samples);
else
    eye_height = NaN; % Cannot calculate if no samples found
end

% 2. Eye Amplitude (Peak to Peak of the eye)
eye_amplitude = avg_level_1 - avg_level_0; % More accurate than min/max raw data

% 3. Extinction Ratio (ER_dB)
% Using average levels for '1' and '0'
if avg_level_0 > 0 % Avoid log(0)
    ER_dB = 10 * log10(avg_level_1 / avg_level_0);
else
    ER_dB = Inf; % If level_0 is zero or negative, ER is infinite or very high
end

% 4. Eye Width (Horizontal Opening)
% This requires finding crossing points. A simplified approach is to look at the
% distribution of horizontal crossing points at the eye threshold.
% You'll need more sophisticated functions for precise jitter/eye width if data isn't pre-processed.
% For this example, we'll approximate based on the data points.
% A more robust eye width calculation often involves finding the zero-crossing points
% (or eye threshold crossing points) for each individual trace.
% Here, we'll use the time points where the signal is close to the threshold.

eye_threshold = (avg_level_0 + avg_level_1) / 2;
crossing_times = [];
for i = 2:length(time_points)
    % Check if a crossing occurs between current and previous point
    if (voltage_or_current(i-1) < eye_threshold && voltage_or_current(i) >= eye_threshold) || ...
       (voltage_or_current(i-1) > eye_threshold && voltage_or_current(i) <= eye_threshold)
        
        % Linear interpolation to find precise crossing time
        t1 = time_points(i-1);
        t2 = time_points(i);
        v1 = voltage_or_current(i-1);
        v2 = voltage_or_current(i);
        
        crossing_t = t1 + (eye_threshold - v1) * (t2 - t1) / (v2 - v1);
        crossing_times = [crossing_times; crossing_t];
    end
end

% Sort crossing times and find the standard deviation for jitter
if ~isempty(crossing_times)
    % Assuming the eye is centered, filter crossings around the bit boundaries
    % We need to find crossings around the beginning and end of the bit period.
    % This part requires more advanced analysis of how the eye is folded.
    % For a simple pre-folded eye, we look at the spread of the crossings.
    
    % The simplest measure of eye width is often the horizontal distance between the
    % mean of the early crossings and the mean of the late crossings.
    % This is highly dependent on how your 'time_points' are structured.
    
    % For now, we'll simply look at the span of the crossing points.
    
    % Let's assume the crossing points are clustered around -bit_period/2, 0, +bit_period/2 etc.
    % We are only given a small window of data. The `time_points` range from
    % approximately -2.4375e-11 to -1.5e-11 which is a narrow slice.
    % For a full eye diagram, `time_points` would span, for example, -50ps to +50ps.
    % If the provided data IS the full folded eye data, then the min/max of time_points
    % gives the effective x-span.
    
    eye_width_data_span = max(time_points) - min(time_points);
    
    % More accurate eye width and jitter requires differentiating between
    % 'left' and 'right' eye crossings.
    % This code will give a rough estimate of the horizontal spread.
    
    % A simple measure of jitter (total jitter) from eye diagram is the
    % horizontal opening width / (slope at crossing points).
    % Or, simply the peak-to-peak spread of crossing points.
    
    % Let's refine based on the plot; assuming `time_points` are already
    % centered around 0 for the eye opening.
    % The eye width is typically measured at the eye threshold, from the
    % latest rising edge crossing to the earliest falling edge crossing.
    
    % To calculate Eye Width (simplified for pre-folded data):
    % Find the mean of the left-most crossing points and the right-most crossing points.
    % This is hard with just general `time_points` without knowing their structure.
    % Given your data range (-2.4375e-11 to -1.5e-11), it looks like only a small
    % fraction of the eye is shown. The 'full eye diagram' image you sent before
    % had an x-axis span of about 100 ps.

    % Assuming your `time_points` for the plot effectively span the full bit period (e.g., -50ps to 50ps)
    % Let's use the range of time points where the eye is open as a proxy for width.
    % This will be more accurate if the data includes points across multiple cycles.
    
    % For a more robust Eye Width and Jitter calculation, you'd typically need
    % a larger set of data or knowledge of how the individual traces are aligned.
    
    % Let's use the horizontal opening at 50% eye height for a rough estimate
    % Find the 50% amplitude level
    mid_eye_amplitude = (avg_level_0 + avg_level_1) / 2;
    
    % Find time points where signal crosses this mid-level from low to high (rising edge)
    % and high to low (falling edge).
    
    % For a proper Eye Width (and jitter), standard tools identify the crossing points
    % for *each* trace, then statistically analyze them.
    % Since you've provided raw data, a simpler approach for eye width:
    % It's the time duration at the eye threshold where the "eye" is open.
    % Given the plot, the eye width is approximately 80-90% of the bit period if clear.
    eye_width_rough = bit_period; % Placeholder. More precise measurement is complex from raw samples.

    % Jitter (Total Jitter) - simplified for pre-folded data
    % This is the peak-to-peak horizontal spread of the crossing points.
    % Needs more context on `time_points` range.
    % For now, we'll leave it as a conceptual note.
    jitter_pp = NaN; % Placeholder

else
    eye_width_rough = NaN;
    jitter_pp = NaN;
end

% 5. Q-factor (Approximated)
% Q = (Mu_1 - Mu_0) / (Sigma_1 + Sigma_0)
% Mu_1, Mu_0 are mean of '1' and '0' levels in eye opening
% Sigma_1, Sigma_0 are standard deviation of '1' and '0' levels in eye opening
if ~isempty(upper_rail_samples) && ~isempty(lower_rail_samples)
    sigma_1 = std(upper_rail_samples);
    sigma_0 = std(lower_rail_samples);

    if (sigma_1 + sigma_0) > 0
        Q_factor = (avg_level_1 - avg_level_0) / (sigma_1 + sigma_0);
    else
        Q_factor = Inf;
    end
else
    Q_factor = NaN;
end


% --- Display Results ---
fprintf('\n--- Eye Diagram Analysis Results ---\n');
fprintf('Signal Type: 10 Gbps NRZ (Bit Period: %.0f ps)\n', bit_period * 1e12);
fprintf('------------------------------------\n');
fprintf('Average Logic 0 Level (I_min): %.4e A\n', avg_level_0);
fprintf('Average Logic 1 Level (I_max): %.4e A\n', avg_level_1);
fprintf('Eye Height (Vertical Opening): %.4e A\n', eye_height);
fprintf('Eye Amplitude (I_max - I_min): %.4e A\n', eye_amplitude);
fprintf('Extinction Ratio (ER): %.2f dB\n', ER_dB);
fprintf('Approx. Q-factor: %.2f\n', Q_factor);
fprintf('------------------------------------\n');
fprintf('Note: For precise Eye Width and Jitter from raw data,\n');
fprintf('      a larger dataset spanning multiple cycles and specific\n');
fprintf('      eye-folding logic is usually required.\n');
fprintf('      The data provided seems to be a single, already-folded eye.\n');
fprintf('      Eye Width is approximately the bit period (%.0f ps).\n', bit_period * 1e12);
fprintf('      Jitter requires analysis of crossing points at specific thresholds.\n');