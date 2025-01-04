function SeedPointIds = diffuseScalarField2D(x, y, Mask, ScalarField, ...
    TotalSeedNumToPut, N_threshold, Proportion, MaximumAspectRatio, UseTemplateIter, UseMask_MN)
N_x = size(y, 2);
N_y = size(y, 1);
%% 生成源场
%% Generate source field
if UseMask_MN
    Mask = log2(1 + Mask);
end
Index = (Mask > 0 );
IndexNegative = ~Index;
SourcePositive = ScalarField(Index) .* Mask(Index);
SourceNegative = ScalarField(IndexNegative);
TotalSource0 = sum([SourcePositive; -SourceNegative], 'all');
SourceNegative = SourceNegative.*( 1 + TotalSource0 ./ sum(SourceNegative, 'all'));
SourceField = gpuArray.zeros(N_y, N_x, 'single');
SourceField(Index) = SourcePositive;
SourceField(IndexNegative) = -SourceNegative;
%% 不鼓励在角点处放置流线种子点
%% Placing streamline seeds at corners is discouraged
CornerSource = max(SourcePositive, [], 'all');
SourceField([1, end], [1, end]) = CornerSource;
SourceField(Index) = SourceField(Index).*(1 - sum(SourceField, 'all') ./ sum(SourceField(Index), 'all'));
%% 释放用过的变量减小GPU内存占用
%% Release used variables to reduce GPU memory usage
clear Index IndexNegative SourcePositive DeltaPositive SourceNegative DeltaNegative CornerSource EdgeSource;
%% 求解poisson方程获得温度场
%% Solve the Poisson equation to obtain the temperature field
TemperatureField = solvePoissonEquation2D(x, y, -SourceField);
%% 全neumann边界条件poisson方程的解可以相差任意常数，这里加一个常数使得场中不存在负数
%% The solutions of the Poisson equation with full Neumann boundary conditions can differ by any constant
% Here, a constant is added so that there are no negative numbers in the field.
IntegrationConstant = abs(min(TemperatureField, [], 'all')) + max(min(abs(TemperatureField), [], 'all'), eps);
TemperatureField = TemperatureField + IntegrationConstant;
%% 放置种子点
%% Placing seeds
if TotalSeedNumToPut > N_threshold
    SeedNumForThisIter = TotalSeedNumToPut * Proportion;
    if SeedNumForThisIter > 9 || ~UseTemplateIter
        if UseTemplateIter
            Depth = max(ceil(log2(SeedNumForThisIter / 9)), 1);
        else
            Depth = max((floor(log2(SeedNumForThisIter))), 1);
        end
        Fields = splitFiled(TemperatureField, Depth, MaximumAspectRatio);
        MinIndex = zeros(size(Fields, 1), 1);
        for i = 1:size(Fields, 1)
            [~, I] = min(SplitFiledField(Fields{i}), [], 'all');
            MinIndex(i) = Fields{i}(I);
        end
    else
        [~, MinIndex] = min(TemperatureField, [], 'all');
    end
else
    [~, MinIndex] = min(TemperatureField, [], 'all'); 
end
[Iy, Ix] = ind2sub(size(x) - 2, MinIndex);
SeedPointIds = sub2ind(size(x), Iy + 1, Ix + 1);
end
