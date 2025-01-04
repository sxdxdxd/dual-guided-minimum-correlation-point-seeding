function TotalEntropy = calculateTotalEntropy2D(BinField, N_Bins)
N_x = size(BinField, 2);
N_y = size(BinField, 1);
Total_NUM = N_x * N_y;
assert(N_Bins > 0 && N_Bins <= intmax('int8'));
%% 计算全局直方图
%% Calculate the global histogram
Histogram = zeros(N_Bins, 1);
for i = 1:N_y
    for j = 1:N_x
        Index = BinField(i, j);
        Histogram(Index) = Histogram(Index) + 1;
    end
end
ProbabilityDistribution = Histogram ./ Total_NUM;
TotalEntropy = calculateEntropy(ProbabilityDistribution);
end