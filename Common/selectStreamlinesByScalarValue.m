function SelectedStreamlines = selectStreamlinesByScalarValue(StreamLine, N_Select, X, ScalarField)
N_Seeds = size(StreamLine, 1);
if ndims(ScalarField) == 3
    N_x = size(ScalarField, 2); N_y = size(ScalarField, 1); N_z = size(ScalarField, 3);
    Xmin = X(1, 1, 1, 1); Xmax = X(1, end, 1, 1);
    Ymin = X(1, 1, 1, 2); Ymax = X(end, 1, 1, 2);
    Zmin = X(1, 1, 1, 3); Zmax = X(1, 1, end, 3);
else
    N_x = size(ScalarField, 2); N_y = size(ScalarField, 1);
    Xmin = X(1, 1, 1); Xmax = X(1, end, 1);
    Ymin = X(1, 1, 2); Ymax = X(end, 1, 2);
end
if N_Seeds > N_Select && N_Select ~= 0
    %% 分配储存空间
    %% Allocate storage
    StreamlineWithScalarSum = cell(N_Seeds, 3);
    StreamlineWithScalarSum(:, 1:2) = StreamLine;
    %% 将流线按其经过点的标量之和排序
    %% Sort the streamlines by the sum of the scalars of the points they pass through
    if ndims(ScalarField) == 3
        for i = 1:N_Seeds
            [Index_X, Index_Y, Index_Z] = convertAbsoluteCoordinatesToGrid3D(StreamLine{i, 1}, ...
                N_x, Xmax, Xmin, N_y, Ymax, Ymin, N_z, Zmax, Zmin);
            Index = sub2ind([N_y, N_x, N_z], Index_Y, Index_X, Index_Z);
            StreamlineWithScalarSum{i, 3} = sum(ScalarField(Index));
        end
    else
        for i = 1:N_Seeds
            [Index_X, Index_Y] = convertAbsoluteCoordinatesToGrid2D(StreamLine{i, 1}, N_x, Xmax, Xmin, N_y, Ymax, Ymin);
            Index = sub2ind([N_y, N_x], Index_Y, Index_X);
            StreamlineWithScalarSum{i, 3} = sum(ScalarField(Index));
        end
    end
    StreamlineWithScalarSum = sortrows(StreamlineWithScalarSum, 3, 'descend');
    SelectedStreamlines = StreamlineWithScalarSum(1:N_Select, 1:2);
elseif N_Select ~= 0
    %% 若生成的流线数量不大于用户指定的种子点数量则直接返回
    %% If the number of generated streamlines is not greater than the number of seed points specified by the user, it will be returned directly
    SelectedStreamlines = StreamLine;
else
    SelectedStreamlines = {};
end
end