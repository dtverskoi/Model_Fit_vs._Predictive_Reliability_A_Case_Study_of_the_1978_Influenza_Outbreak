function hpd_region = empirical_hpd_4d(samples, alpha)
    % EMPIRICAL_HPD_4D computes the HPD region for 4D posterior samples.
    %
    % Inputs:
    %   samples - A matrix of posterior samples (n_samples x 4).
    %   alpha   - The credible interval level (e.g., 0.05 for 95% CI).
    %
    % Outputs:
    %   hpd_region - A structure containing:
    %                - region_samples: Samples within the HPD region.
    %                - density_threshold: The minimum density value in the HPD region.

    % Validate that the samples are 4D
    [n_samples, n_dimensions] = size(samples);
    if n_dimensions ~= 4
        error('The input samples must be a matrix of size n_samples x 4 for 4D data.');
    end
    
    % Parameters
    credible_mass = 1 - alpha; % Mass to include in the HPD region

    % Estimate posterior density using kernel density estimation (KDE)
    bins=std(samples)*(4/(n_samples*(2+n_dimensions)))^(1/(4+n_dimensions));
    [density, ~] = mvksdensity(samples, samples, 'Bandwidth', bins);

    % Combine densities with corresponding samples
    posterior_data = [samples, density];
    
    % Sort by density in descending order
    posterior_data = sortrows(posterior_data, -(n_dimensions + 1));
    
    % Calculate cumulative probability
    sorted_density = posterior_data(:, end);
    cumulative_prob = cumsum(sorted_density) / sum(sorted_density);

    % Find threshold where cumulative probability exceeds credible mass
    density_threshold = sorted_density(find(cumulative_prob >= credible_mass, 1));
    
    % Extract samples within the HPD region
    hpd_samples = posterior_data(posterior_data(:, end) >= density_threshold, 1:n_dimensions);
    
    % Output the HPD region and threshold
    hpd_region.region_samples = hpd_samples;
    hpd_region.density_threshold = density_threshold;
    
    % Display results
    fprintf('Density Threshold: %.4f\n', density_threshold);
    fprintf('Number of Samples in HPD Region: %d\n', size(hpd_samples, 1));
end
