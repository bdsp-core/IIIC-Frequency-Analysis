function [spatial_modes, spatial_modes_1D, temporal_modes, frequencies, periodicities] = dmd_eeg_2D(data, rank, sampling_frequency)
    % data: 25x25x2800 2D EEG data matrix
    % rank: the rank for truncating the SVD
    % sampling_frequency: the sampling frequency in Hz

    % Step 1: Reshape the data to unroll the spatial dimensions
    [n, m, l] = size(data);
    reshaped_data = reshape(data, [n * m, l]);

    % Step 2: Split the data into X and Y matrices
    X = reshaped_data(:, 1:end-1);
    Y = reshaped_data(:, 2:end);

    % Step 3: Compute the SVD of X
    [U, S, V] = svd(X, 'econ');

    % Step 4: Truncate the SVD matrices based on the rank
    U = U(:, 1:rank);
    S = S(1:rank, 1:rank);
    V = V(:, 1:rank);

    % Step 5: Calculate A_tilde
    A_tilde = U' * Y * V * inv(S);
    % Step 8: Calculate eigenvalues and eigenvectors of A_tilde
    [W, Lambda] = eig(A_tilde);

    % Step 6: Calculate spatial modes (Phi) and temporal modes (Psi)
    spatial_modes_1D = Y * V * inv(S) * W;
    temporal_modes = X' * spatial_modes_1D;

    % Step 7: Reshape spatial_modes_1D back to the original 2D spatial dimensions
    spatial_modes = reshape(spatial_modes_1D, [n, m, rank]);



    % Step 9: Calculate frequencies and periodicities
    dt = 1 / sampling_frequency;
    eigenvalues = diag(Lambda);
    frequencies = abs(imag(log(eigenvalues))) / (2 * pi * dt);
    periodicities = 1 ./ abs(frequencies);
end
