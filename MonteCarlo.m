function [gain_samples, power_responses, time_outputs, power_envelope_max, power_envelope_min, time_envelope_max, time_envelope_min] = MonteCarlo(num_monte_carlo, coupling_coeffs, system_order, time_vector, num_points, nominal_reflectivity, reflectivity_tolerance, nominal_loss, loss_tolerance, effective_index, detector_tolerance, freq_difference, speed_of_light, ring_lengths, input_signal, reference_transfer_func)
% MonteCarlo - Perform Monte Carlo analysis for perturbed system parameters
%
% Inputs:
% num_monte_carlo - Number of Monte Carlo runs
% coupling_coeffs - Array of coupling coefficients
% system_order - System order
% time_vector - Time vector
% num_points - Number of points
% nominal_reflectivity - Nominal reflectivity parameters
% reflectivity_tolerance - Tolerance for reflectivity
% nominal_loss - Nominal loss parameters
% loss_tolerance - Tolerance for loss
% effective_index - Effective refractive index
% detector_tolerance - Tolerance for detector
% freq_difference - Frequency difference
% speed_of_light - Speed of light
% ring_lengths - Ring lengths array
% input_signal - Input signal
% reference_transfer_func - Reference transfer function
%
% Outputs:
% gain_samples - gain for each Monte Carlo run
% power_responses - All normalized power responses
% time_outputs - All time-domain outputs
% power_envelope_max - Maximum envelope of power responses
% power_envelope_min - Minimum envelope of power responses
% time_envelope_max - Maximum envelope of time outputs
% time_envelope_min - Minimum envelope of time outputs

% Initialize gain storage
gain_samples = zeros(num_monte_carlo, 1);

% Define gain measurement time window. We will average the gain over this window.
% Here we assume that the gain settles between 60 ns and 100 ns.
time_indices = find(time_vector >= 50e-9 & time_vector <= 60e-9);

% Preallocate storage for full distribution
time_outputs = zeros(num_monte_carlo, num_points);
power_responses = zeros(num_monte_carlo, num_points);

% Initialize envelope tracking
power_envelope_max = zeros(1, num_points);
power_envelope_min = Inf(1, num_points);
time_envelope_max = -Inf(1, num_points);
time_envelope_min = Inf(1, num_points);

% Main Monte Carlo loop
for run_idx = 1:num_monte_carlo
    % Generate perturbed parameters
    perturbed_reflectivity = nominal_reflectivity(1:system_order) .* ...
        (1 + reflectivity_tolerance * randn(1, system_order));

    perturbed_loss = nominal_loss(1:system_order) .* ...
        (1 - loss_tolerance * abs(randn(1, system_order)));

    perturbed_effective_index = effective_index * ...
        (1 + detector_tolerance * randn);

    
    % Build perturbed transfer function
    transfer_function = ones(1, num_points);
    propagation_constant = 2*pi*freq_difference / (speed_of_light/perturbed_effective_index);
    
    for ring_idx = 1:system_order
        ring_component = (1/coupling_coeffs(ring_idx)) .* ...
            ((1 - perturbed_reflectivity(ring_idx)^2) .* perturbed_loss(ring_idx) ./ ...
            (1 - perturbed_reflectivity(ring_idx)^2 .* perturbed_loss(ring_idx) .* ...
            exp(-1j*propagation_constant*ring_lengths(ring_idx))));
        
        % Compose the transfer functions
        transfer_function = transfer_function .* ring_component;
    end
    
    % Store normalized power response
    normalized_power = abs(transfer_function).^2 ./ max(abs(reference_transfer_func).^2);
    power_responses(run_idx, :) = normalized_power;
    
    % Compute time-domain output and store
    time_response = real(ifft(fftshift(input_signal .* transfer_function)));
    time_outputs(run_idx, :) = time_response;
    
    % Calculate gain for this run
    % We average the gain over the defined time window
    gain_samples(run_idx) = mean(time_response(time_indices));
    
    % Update running envelope tracking
    power_envelope_max = max(power_envelope_max, normalized_power);
    power_envelope_min = min(power_envelope_min, normalized_power);
    time_envelope_max = max(time_envelope_max, time_response);
    time_envelope_min = min(time_envelope_min, time_response);
end 
end