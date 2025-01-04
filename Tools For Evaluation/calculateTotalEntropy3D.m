function TotalEntropy = calculateTotalEntropy3D(BinField, N_Bins)
N_x = size(BinField, 2);
N_y = size(BinField, 1);
N_z = size(BinField, 3);
Total_NUM = N_x * N_y * N_z;
assert(N_Bins > 0 && N_Bins <= intmax('int16'));
%% 计算全局直方图
%% Calculate the global histogram
Histogram = zeros(N_Bins, 1);
for I = 1:N_y
    for J = 1:N_x
        for K = 1:N_z
            Index = BinField(I, J, K);
            Histogram(Index) = Histogram(Index) + 1;
        end
    end
end
% [Histogram, ~] = histcounts(Bin_Field(:), 1:(N_Bins + 1));
ProbabilityDistribution = Histogram ./ Total_NUM;
TotalEntropy = calculateEntropy(ProbabilityDistribution);
end