function [spatial_modes, temporal_modes, frequencies, periodicities] = dmd_eeg(eeg_data, rank, fs)
% Assume 'eeg_data' is a matrix of size (num_channels x num_time_points)
% with each row representing a channel pair and each column a time point.
    sampling_frequency=fs;

% Step 1: Create data matrices X and Y (lagged by 1 time step)
    X = eeg_data(:, 1:end-1);
    Y = eeg_data(:, 2:end);

    % Step 2: Perform SVD on X
    [U, S, V] = svd(X, 'econ');

    figure(1); plot(diag(S));
    U = U(:, 1:rank);
    S = S(1:rank, 1:rank);
    V = V(:, 1:rank);

    % Step 3: Calculate A_tilde matrix
    A_tilde = U' * Y * V / S;

    % Step 4: Perform eigendecomposition on A_tilde
    [W, D] = eig(A_tilde);

    % Step 5: Calculate spatial modes (Phi) and temporal modes (Psi)
    spatial_modes = Y * V / S * W;
%     spatial_modes = Y * V * inv(S) * W;
    temporal_modes = X' * spatial_modes;

    % Step 6: Extract frequencies and periodicities
    dt = 1 / sampling_frequency;
    eigenvalues = diag(D);
%     frequencies = abs(imag(log(eigenvalues))) / (2 * pi * dt); %in Hz
    frequencies = abs(log(eigenvalues)) / (2 * pi * dt); %in Hz
    periodicities = 1 ./ abs(frequencies);

%     eigenvalues = diag(D);
%     frequencies = abs(imag(log(eigenvalues))) / (2 * pi); % In Hz
%     periodicities = 1 ./ frequencies; % In seconds

for i = 1:size(temporal_modes, 2)
    [~, max_idx] = max(abs(temporal_modes(:, i)));
    sign_factor = sign(temporal_modes(max_idx, i));
    temporal_modes(:, i) = temporal_modes(:, i) * sign_factor;
end

for i = 1:size(spatial_modes, 2)
    [~, max_idx] = max(abs(spatial_modes(:, i)));
    sign_factor = sign(spatial_modes(max_idx, i));
    spatial_modes(:, i) = spatial_modes(:, i) * sign_factor;
end

end