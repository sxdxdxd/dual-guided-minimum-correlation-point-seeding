function ExpressionData = calculateExpressivenessForStreamLines(X, V, ScalarField, PointPassedByStreamLines, Options)
if ndims(X) == 4
    x = X(:, :, :, 1); y = X(:, :, :, 2); z = X(:, :, :, 3);
    u = V(:, :, :, 1); v = V(:, :, :, 2); w = V(:, :, :, 3);
    N_x = size(x, 2); N_y = size(y, 1); N_z = size(z, 3);
    Xmin = x(1, 1, 1); Xmax = x(1, end, 1);
    Ymin = y(1, 1, 1); Ymax = y(end, 1, 1);
    Zmin = z(1, 1, 1); Zmax = z(1, 1, end);
    u_g = gpuArray(single(u));
    v_g = gpuArray(single(v));
    w_g = gpuArray(single(w));
    % ScalarField_g = gpuArray(single(ScalarField));
    OmegaField = calculateOmega3D(x, y, z, u, v, w);
else
    x = X(:, :, 1); y = X(:, :, 2);
    u = V(:, :, 1); v = V(:, :, 2);
    N_x = size(x, 2); N_y = size(y, 1);
    Xmin = x(1, 1); Xmax = x(1, end);
    Ymin = y(1, 1); Ymax = y(end, 1);
    u_g = gpuArray(single(u));
    v_g = gpuArray(single(v));
    % ScalarField_g = gpuArray(single(ScalarField));
    OmegaField = calculateOmega2D(x, y, u, v);
end
ExpressionData = struct;
if ndims(X) == 4
    Mask = generateMask3D(N_x, N_y, N_z, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, PointPassedByStreamLines);
    [U, V, W] = reconstructVectorField3D(Mask, u_g, v_g, w_g, ...
        Options.OvercorrectionFactor, Options.Mu, Options.DiffusionThreshold, Options.OptimizeMaxIterations);
    ThetaDifferenceField = calculateThetaDifference(cat(4, u, v, w), cat(4, U, V, W));
    ScalarField_R = reconstructScalarField3D(X, Mask, ScalarField);
    Velocity = cat(4, u, v, w);
    Velocity_R = cat(4, U, V, W);
    ExpressionData.PSNR_Phi = psnr(double(ScalarField_R), double(ScalarField), range(double(ScalarField), 'all'), 'DataFormat', 'SSS');
    ExpressionData.PSNR_V = psnr(double(Velocity_R), Velocity, range(Velocity, 'all'), 'DataFormat', 'SSSC');
    load("BinMap.mat", "BinMap");
    N_Bins = 360;
    BinField = calculateBinField3D(u, v, w, BinMap);
    BinField_R = calculateBinField3D(U, V, W, BinMap);
    TotalEntropy = calculateTotalEntropy3D(BinField, N_Bins);
    ExpressionData.ConditionalEntropy = calculateGlobalConditionalEntropy3D(BinField, BinField_R, N_Bins, TotalEntropy, 'H(X|Y)');
else
    Mask = generateMask2D(N_x, N_y, Xmin, Xmax, Ymin, Ymax, PointPassedByStreamLines);
    [U, V] = reconstructVectorField2D(Mask, u_g, v_g, ...
        Options.OvercorrectionFactor, Options.Mu, Options.DiffusionThreshold, Options.OptimizeMaxIterations);
    ThetaDifferenceField = calculateThetaDifference(cat(3, u, v),cat(3, U, V));
    ScalarField_R = reconstructScalarField2D(X, Mask, ScalarField);
    Velocity = cat(3, u, v);
    Velocity_R = cat(3, U, V);
    ExpressionData.PSNR_Phi = psnr(double(ScalarField_R), double(ScalarField), range(double(ScalarField), 'all'), 'DataFormat', 'SS');
    ExpressionData.PSNR_V = psnr(double(Velocity_R), Velocity, range(Velocity, 'all'), 'DataFormat', 'SSC');
    N_Bins = 60;
    BinField = calculateBinField2D(u, v, N_Bins);
    BinField_R = calculateBinField2D(U, V, N_Bins);
    TotalEntropy = calculateTotalEntropy2D(BinField, N_Bins);
    ExpressionData.ConditionalEntropy = calculateGlobalConditionalEntropy2D(BinField, BinField_R, N_Bins, TotalEntropy, 'H(X|Y)');
end

ExpressionData.AAD = mean(ThetaDifferenceField, 'all');
Error = abs((ScalarField_R - ScalarField) ./ ScalarField);
Error(ScalarField == 0) = 0;
Error = rmoutliers(Error(:), "gesd","ThresholdFactor", 0.001, "MaxNumOutliers", 100);
ExpressionData.ARE = mean(Error, 'all');%Relative Error
ExpressionData.Expressiveness = calculateExpressiveness(Options.alpha, ThetaDifferenceField, OmegaField, ScalarField_R, ScalarField, 'AddConst');
end
