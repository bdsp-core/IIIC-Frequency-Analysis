function [spatial_modes, temporal_modes, frequencies, periodicities] = logistic_dmd_eeg(eeg_data, rank, fs)
% Assume 'eeg_data' is a matrix of size (num_channels x num_time_points)
% with each row representing a channel pair and each column a time point.
    sampling_frequency=fs;

 % Define the logit and sigmoid functions
    logit = @(x) log(x ./ (1 - x));
    sigmoid = @(x) 1 ./ (1 + exp(-x));

    % Step 1: Transform the binary data using the logit function
    transformed_data = logit(eeg_data);
 
% Step 1: Create data matrices X and Y (lagged by 1 time step)
    X = eeg_data(:, 1:end-1);
    Y = eeg_data(:, 2:end);

    % Step 2: Perform SVD on X
    [U, S, V] = svd(X, 'econ');

%     figure; plot(diag(S));
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
    frequencies = abs(imag(log(eigenvalues))) / (2 * pi * dt); %in Hz
    periodicities = 1 ./ abs(frequencies);

%     eigenvalues = diag(D);
%     frequencies = abs(imag(log(eigenvalues))) / (2 * pi); % In Hz
%     periodicities = 1 ./ frequencies; % In 

 % Step 3: Transform the spatial_modes back to binary space using the sigmoid function
    spatial_modes = sigmoid(spatial_modes);


end