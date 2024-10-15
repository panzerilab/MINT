% Step 1: Generate synthetic data
data = randn(50, 1000);

% Step 2: Call the PCA pipeline to reduce to 2 dimensions
ndims = 2;
[reduced_data, coeff] = pca_pipeline(data, ndims,23);

% Step 3: Plot the reduced data
figure;
scatter(reduced_data(1, :), reduced_data(2, :), 25, 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('Data Reduced to 2 Principal Components');
grid on;
