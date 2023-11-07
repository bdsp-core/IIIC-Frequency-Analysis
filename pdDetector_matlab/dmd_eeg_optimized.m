function [spatial_modes, temporal_modes, frequencies, periodicities] = dmd_eeg_optimized(eeg_data, rank, fs)
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

 % Step 4: Compute the low-rank approximation of X
    X_low_rank = U * S * V';

    % Step 5: Estimate the noise covariance matrix W
    noise = X - X_low_rank;
    W = cov(noise');

    % Step 6: Calculate A_tilde
%     A_tilde = U' * Y * V * inv(S) * inv(W);
    S_pinv_weighted = V * (S' * S + W)^-1 * S';
    A_tilde = U' * Y * S_pinv_weighted;
%     A_tilde = U' * Y * V * inv(S);

    % Step 7: Compute the eigendecomposition of A_tilde
    [W_eig, Lambda] = eig(A_tilde);
    
    % Step 8: Calculate spatial modes (Phi) and temporal modes (Psi)
    spatial_modes = Y * V * inv(S) * W_eig;
    temporal_modes = X' * spatial_modes;

    % Step 9: Calculate the frequencies and periodicities
%     frequencies = angle(Lambda) * sampling_frequency / (2 * pi);
%     periodicities = 1 ./ frequencies;

    dt = 1 / sampling_frequency;
    eigenvalues = diag(Lambda);
    frequencies = angle(eigenvalues) / (2 * pi * dt); %in Hz
    periodicities = 1 ./ abs(frequencies);

% 
% 
%     % Step 3: Calculate A_tilde matrix
%     A_tilde = U' * Y * V / S;
% 
%     % Step 4: Perform eigendecomposition on A_tilde
%     [W, D] = eig(A_tilde);
% 
%     % Step 5: Calculate spatial modes (Phi) and temporal modes (Psi)
%     spatial_modes = Y * V / S * W;
% %     spatial_modes = Y * V * inv(S) * W;
%     temporal_modes = X' * spatial_modes;
% 
%     % Step 6: Extract frequencies and periodicities
%     dt = 1 / sampling_frequency;
%     eigenvalues = diag(D);
% %     frequencies = abs(imag(log(eigenvalues))) / (2 * pi * dt); %in Hz
%     frequencies = abs(log(eigenvalues)) / (2 * pi * dt); %in Hz
%     periodicities = 1 ./ abs(frequencies);

%     eigenvalues = diag(D);
%     frequencies = abs(imag(log(eigenvalues))) / (2 * pi); % In Hz
%     periodicities = 1 ./ frequencies; % In seconds
end