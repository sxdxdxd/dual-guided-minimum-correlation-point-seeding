function Entropy = calculateEntropy(ProbabilityDistribution)
%% 计算熵场
%% Calculate the entropy field
ProbabilityDistributionIn = nonzeros(ProbabilityDistribution);
Log_Probability = -log2(ProbabilityDistributionIn);
Entropy = sum((ProbabilityDistributionIn .* Log_Probability), 'all');
end