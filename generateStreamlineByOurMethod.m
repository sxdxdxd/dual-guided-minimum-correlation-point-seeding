function StreamlineSet = generateStreamlineByOurMethod(X, V, ScalarField, OmegaField, Options)
assert(Options.N_Seeds > 0 &&Options.Stepsize > 0 && Options.Max_Vertices > 0);
assert(Options.InitialRatio <= 1 && Options.Proportion <= 1 && Options.N_threshold < Options.N_Seeds);
assert((ndims(X) == 3 &&  ndims(V) == 3 && ndims(ScalarField) == 2) || (ndims(X) == 4 &&  ndims(V) == 4 && ndims(ScalarField) == 3));%只支持2维或3维的数据
assert(sum(isinf(X), "all") == 0 && sum(isinf(V), "all") == 0 && sum(isinf(ScalarField), "all") == 0);
assert(sum(isnan(X), "all") == 0 && sum(isnan(V), "all") == 0 && sum(isnan(ScalarField), "all") == 0);
assert(Options.alpha >= 0 && Options.alpha <= 1);
assert(any([Options.UseMask_01, Options.UseMask_0N, Options.UseMask_MN, Options.UseMaskGaussian]) ...
    && nnz([Options.UseMask_01, Options.UseMask_0N, Options.UseMask_MN, Options.UseMaskGaussian]) == 1);
%% 获取需要的参数
%% Get the required parameters
N_Iter = min(ceil((1 - Options.InitialRatio) * Options.N_Seeds), Options.N_Seeds - 1);
N_IterForV = floor((1 - Options.alpha) * N_Iter);
N_IterForPhi = N_Iter - N_IterForV;
N_Initial = Options.N_Seeds - N_Iter;
N_InitialForV = floor((1 - Options.alpha) * N_Initial);
N_InitialForPhi = N_Initial - N_InitialForV;
IsThreeDimensionalData = ndims(ScalarField) == 3;
%% 数据解包
%% data unpacking
if IsThreeDimensionalData
    x = X(:, :, :, 1); y = X(:, :, :, 2); z = X(:, :, :, 3);
    u = V(:, :, :, 1); v = V(:, :, :, 2); w = V(:, :, :, 3);
    N_x = size(x, 2); N_y = size(y, 1); N_z = size(z, 3);
    Xmin = x(1, 1, 1); Xmax = x(1, end, 1);
    Ymin = y(1, 1, 1); Ymax = y(end, 1, 1);
    Zmin = z(1, 1, 1); Zmax = z(1, 1, end);
else
    x = X(:, :, 1); y = X(:, :, 2);
    u = V(:, :, 1); v = V(:, :, 2);
    N_x = size(x, 2); N_y = size(y, 1);
    Xmin = x(1, 1); Xmax = x(1, end);
    Ymin = y(1, 1); Ymax = y(end, 1);
end
OmegaFlag = Options.alpha == 0;
%% Step 1: 处理标量场
%% Step 1: Process the scalar field
if min(ScalarField, [], 'all') < 0 && Options.alpha ~= 0
    if strcmp(Options.handlingNegative,'AddConst')
        ScalarField = ScalarField + abs(min(ScalarField, [], 'all'));
    else %strcmp(Options.handlingNegative,'Abs')
        ScalarField = abs(ScalarField);
    end
end
%% Step 2: 寻找标量场中的极大值点
%% Step 2: Find the maximum value point in the scalar field
if OmegaFlag || Options.alpha == 0
    ExtremePointList = findExtremePoints(N_Initial, OmegaField);
elseif Options.alpha == 1
    ExtremePointList = findExtremePoints(N_Initial, ScalarField);
else
    ExtremePointListV = findExtremePoints(N_InitialForV, OmegaField);
    ExtremePointListPhi = findExtremePoints(N_InitialForPhi, ScalarField);
end
%% Step 3:生成初始流线
%% Step 3: Generate initial streamline
%% 根据标量场选择初始种子点，生成初始流线，标记流线经过的点
%% Select the initial seed point based on the scalar field
if OmegaFlag || Options.alpha == 0
    StreamLines = generateInitialStreamlines(ExtremePointList, OmegaField, N_Initial);
elseif Options.alpha == 1
    StreamLines = generateInitialStreamlines(ExtremePointList, ScalarField, N_Initial);
else
    StreamLines_Phi = generateInitialStreamlines(ExtremePointListPhi, ScalarField, N_InitialForPhi);
    StreamLines_V = generateInitialStreamlines(ExtremePointListV, OmegaField, N_InitialForV);
    N_StreamlineV = size(StreamLines_V, 1);
    N_StreamlinePhi = size(StreamLines_Phi, 1);
    StreamLines = cat(1, StreamLines_Phi, StreamLines_V);
    clear StreamLines_V StreamLines_Phi;
end
%% 标记现有流线经过的点
%% marks the points passed by existing streamlines
MaskField = calculateStreamlineDensity(StreamLines);
%% 提前把变量传到GPU上
%% Transfer variables to GPU in advance
if OmegaFlag || Options.alpha == 0
    DiffuseField_g = gpuArray(single(OmegaField));
elseif Options.alpha == 1
    DiffuseField_g = gpuArray(single(ScalarField));
else%% Options.alpha ~= 0 && Options.alpha ~= 1
    ScalarField_g = gpuArray(single(ScalarField));
    OmegaField_g = gpuArray(single(OmegaField));
end
MaskField_g = gpuArray(single(MaskField));
x_g = gpuArray(single(x));
y_g = gpuArray(single(y));
if IsThreeDimensionalData
    z_g = gpuArray(single(z));
end
%% Step 4: 热扩散引导的最无关点种子重采样(基于物理的种子点采样)
%% Step 4: Physics-Based Seeding
%% 生成流线直到达到指定数量
%% Generate streamlines until the specified number is reached
while size(StreamLines, 1) < Options.N_Seeds
    %% 来回切换
    if ~OmegaFlag && Options.alpha ~= 0 && Options.alpha ~= 1
        if N_StreamlinePhi / N_StreamlineV <= Options.alpha / (1 - Options.alpha)
            DiffuseField_g = ScalarField_g;
            PhiFlag = true;
            N_SeedToPut = N_IterForPhi - N_StreamlinePhi;
        else
            DiffuseField_g = OmegaField_g;
            PhiFlag = false;
            N_SeedToPut = N_IterForV - N_StreamlineV;
        end
    else
         N_SeedToPut = Options.N_Seeds - size(StreamLines, 1);
         PhiFlag = false;
    end
    %% 选出与现有流线上的点最无关的一个（或N个）点
    %% Select the point (or N) that is least related to the points on the existing streamline
    if IsThreeDimensionalData
        NewSeedPointIds = diffuseScalarField3D(x_g, y_g, z_g, MaskField_g, DiffuseField_g, ...
            N_SeedToPut, Options.N_threshold, ...
            Options.Proportion, Options.MaximumAspectRatio, Options.UseTemplateIter, Options.UseMask_MN || Options.UseMaskGaussian);
        NewSeedPoints = [x(NewSeedPointIds), y(NewSeedPointIds), z(NewSeedPointIds)];
    else
        NewSeedPointIds = diffuseScalarField2D(x_g, y_g, MaskField_g, DiffuseField_g, ...
            N_SeedToPut, Options.N_threshold, ...
            Options.Proportion, Options.MaximumAspectRatio, Options.UseTemplateIter, Options.UseMask_MN || Options.UseMaskGaussian);
        NewSeedPoints = [x(NewSeedPointIds), y(NewSeedPointIds)];
    end
    if Options.UseTemplateIter
        if IsThreeDimensionalData
            NewSeedPoints = setSeedPointTemplate3D(x, y, z, NewSeedPointIds, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
        else
            NewSeedPoints = setSeedPointTemplate2D(x, y, NewSeedPointIds, Xmin, Xmax, Ymin, Ymax);
        end
    end
    %% 若此次迭代放置的所有种子点都和之前重复，就会陷入死循环，故直接退出
    %% If all the seed points placed in this iteration are repeated with the previous ones, it will fall into an infinite loop, so exit directly
    NewSeedPoints = setdiff(NewSeedPoints, cell2mat(StreamLines(:, 2)), 'rows');
    if isempty(NewSeedPoints)        
        warning("Streamline generation failed, please adjust the parameters and try again " + ...
            "(try to increase the ratio of the number of seed points placed to the number of remaining seed points each time the diffusion equation is solved).");
        break;
    end
    %% 生成新流线
    %% Generate new streamlines
    if IsThreeDimensionalData
        ExtraStreamLine = generate3DStreamlinesBySeeds(x, y, z, u, v, w, double(NewSeedPoints), Options.Stepsize, Options.Max_Vertices);
    else
        ExtraStreamLine = generate2DStreamlinesBySeeds(x, y, u, v, double(NewSeedPoints), Options.Stepsize, Options.Max_Vertices);
    end
    %% 筛选新流线
    %% Select new streamlines
    N_Select = max(floor((Options.N_Seeds - size(StreamLines, 1)) * Options.Proportion), 1);
    if PhiFlag
        SelectedExtraStreamLine = selectStreamlinesByScalarValue(ExtraStreamLine, N_Select, X, ScalarField);
    else
        SelectedExtraStreamLine = selectStreamlinesByScalarValue(ExtraStreamLine, N_Select, X, OmegaField);
    end
    %% 记录新增流线
    %% Record new streamlines
    Mask = calculateStreamlineDensity(SelectedExtraStreamLine);
    if Options.UseMask_01
        MaskField_g(Mask) = 1;
    else
        MaskField_g = MaskField_g + Mask;
    end
    StreamLines = cat(1, StreamLines, SelectedExtraStreamLine);
    if ~OmegaFlag && Options.alpha ~= 0 && Options.alpha ~= 1
        if PhiFlag
            N_StreamlinePhi = N_StreamlinePhi + size(SelectedExtraStreamLine, 1);
        else
            N_StreamlineV = N_StreamlineV + size(SelectedExtraStreamLine, 1);
        end
    end
end
StreamlineSet = StreamLines;
%% function1
    function ExtremePoints = findExtremePoints(N_Seeds, Field)
        if (IsThreeDimensionalData && N_Seeds > 27) || (~IsThreeDimensionalData && N_Seeds > 9)
            if Options.UseTemplateInit
                if IsThreeDimensionalData
                    TargetDepth = max(1, ceil(log2(N_Seeds / 27)));
                else
                    TargetDepth = max(1, ceil(log2(N_Seeds / 9)));
                end
            else
                TargetDepth = max(1, ceil(log2(N_Seeds)));
            end
            Fields = splitFiled(Field, TargetDepth, Options.MaximumAspectRatio);
            Fields(cellfun(@isempty,Fields(:, 1)), :) = [];
            ExtremePoints = zeros(size(Fields, 1), 1);
            for n = 1:size(Fields, 1)
                [~, I] = max(Field(Fields{n}), [], 'all');
                ExtremePoints(n, 1) = Fields{n}(I);
            end
        else
            [~, ExtremePoints] = max(Field, [], 'all');
        end
    end
%% function2
    function InitialStreamlins = generateInitialStreamlines(ExtremePoints, Field, N_Select)
        if N_Select == 0
            InitialStreamlins = {};
            return;
        end
        if Options.UseTemplateInit
            if IsThreeDimensionalData
                Seeds = setSeedPointTemplate3D(x, y, z, ExtremePoints, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
            else
                Seeds = setSeedPointTemplate2D(x, y, ExtremePoints, Xmin, Xmax, Ymin, Ymax);
            end
        else
            if IsThreeDimensionalData
                Seeds = [x(ExtremePoints), y(ExtremePoints), z(ExtremePoints)];
            else
                Seeds = [x(ExtremePoints), y(ExtremePoints)];
            end
        end
        if IsThreeDimensionalData
            TempStreamlins = generate3DStreamlinesBySeeds(x, y, z, u, v, w, Seeds, Options.Stepsize, Options.Max_Vertices);
        else
            TempStreamlins = generate2DStreamlinesBySeeds(x, y, u, v, Seeds, Options.Stepsize, Options.Max_Vertices);
        end
        InitialStreamlins = selectStreamlinesByScalarValue(TempStreamlins, N_Select, X, Field);
    end
%% function3
    function Mask = calculateStreamlineDensity(TargetStreamLines)
        if Options.UseMask_01
            if IsThreeDimensionalData
                Mask = generateMask3D(N_x, N_y, N_z, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, cell2mat(TargetStreamLines(:, 1)));
            else
                Mask = generateMask2D(N_x, N_y, Xmin, Xmax, Ymin, Ymax, cell2mat(TargetStreamLines(:, 1)));
            end
        elseif Options.UseMask_0N
            if IsThreeDimensionalData
                Mask = zeros(N_y, N_x, N_z);
                for n = 1:size(TargetStreamLines, 1)
                    Mask = Mask + generateMask3D(N_x, N_y, N_z, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, cell2mat(TargetStreamLines(n, 1)));
                end
            else
                Mask = zeros(N_y, N_x);
                for n = 1:size(TargetStreamLines, 1)
                    Mask = Mask + generateMask2D(N_x, N_y, Xmin, Xmax, Ymin, Ymax, cell2mat(TargetStreamLines(n, 1)));
                end
            end
        elseif Options.UseMask_MN
            Mask = generateMNMask(cell2mat(TargetStreamLines(:, 1)));
        elseif Options.UseMaskGaussian
            Mask = generateMNMask(cell2mat(TargetStreamLines(:, 1)));
            if IsThreeDimensionalData
                Mask = imgaussfilt3(double(Mask), Options.sigma, 'FilterSize', Options.KernelSize, 'Padding', 'symmetric');
            else
                Mask = imgaussfilt(double(Mask), Options.sigma, 'FilterSize', Options.KernelSize, 'Padding', 'symmetric');
            end
        end
    end
%% function4
    function Mask = generateMNMask(PointsPassedPyStreamline)
        if IsThreeDimensionalData
            Mask = zeros(N_y, N_x, N_z);
            [Index_X, Index_Y, Index_Z] = convertAbsoluteCoordinatesToGrid3D(PointsPassedPyStreamline, ...
                N_x, Xmax, Xmin, N_y, Ymax, Ymin, N_z, Zmax, Zmin);
            Index = sub2ind([N_y, N_x, N_z], Index_Y, Index_X, Index_Z);
            for n = 1:size(Index, 1)
                Mask(Index(n)) = Mask(Index(n)) + 1;
            end
        else
            Mask = zeros(N_y, N_x);
            [Index_X, Index_Y] = convertAbsoluteCoordinatesToGrid2D(PointsPassedPyStreamline, ...
                N_x, Xmax, Xmin, N_y, Ymax, Ymin);
            Index = sub2ind([N_y, N_x], Index_Y, Index_X);
            for n = 1:size(Index, 1)
                Mask(Index(n)) = Mask(Index(n)) + 1;
            end
        end
    end
end
