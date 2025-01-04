function Fileds = splitFiled(Field, TargetDepth, MaximumAspectRatio)
assert(ndims(Field) == 3 || ndims(Field) == 2);
%% 将场拆分为尽量均匀的2^TargetDepth份
%% Split the field into 2^TargetDepth parts as evenly as possible
if ndims(Field) == 3
    KDTree = buildKDTree_3D(Field, TargetDepth, MaximumAspectRatio);
else
    KDTree = buildKDTree_2D(Field, TargetDepth, MaximumAspectRatio);
end
LeafNodeIndex = find([KDTree.Depth] == TargetDepth);
Fileds1 = {KDTree(LeafNodeIndex).LeftRange};
Fileds2 = {KDTree(LeafNodeIndex).RightRange};
Fileds = [Fileds1, Fileds2]';
end
function kd_tree = buildKDTree_3D(Data, target_depth, MaximumAspectRatio)
assert(target_depth < 20 && target_depth > 0 && target_depth == fix(target_depth));
XLength = size(Data, 2);
YLength = size(Data, 1);
ZLength = size(Data, 3);
[XX, YY, ZZ] = meshgrid(1:XLength, 1:YLength, 1:ZLength);
IndexMat = sub2ind([YLength, XLength, ZLength], YY, XX, ZZ);
Index = 0;
NumOfNodes = 2^target_depth - 1;
kd_tree = struct('Father', [], 'LeftSon', [], 'RightSon', [], 'SplitDim', [], 'SplitPoint', [], 'Depth', [], ...
    'LeftRange', cell(NumOfNodes, 1), 'RightRange', cell(NumOfNodes, 1), 'LeftData', cell(NumOfNodes, 1), 'RightData', cell(NumOfNodes, 1));
buildKDTree_recurse(Data, IndexMat, []);
% 递归子函数
% Recursive subfunction
    function Currentindex = buildKDTree_recurse(CurrentData, CurrentIndexMat, PreviousIndex)
        Currentindex = Index + 1;
        Index = Index + 1;
        kd_tree(Currentindex).Father = PreviousIndex;
        if ~isempty(PreviousIndex)
            kd_tree(Currentindex).Depth = kd_tree(PreviousIndex).Depth + 1;
        else
            kd_tree(Currentindex).Depth = 1;
        end
        % 数据降维
        % Data dimensionality reduction
        XData = sum(CurrentData, [1, 3]);
        YData = sum(CurrentData, [2, 3]);
        ZData = sum(CurrentData, [1, 2]);
        SizeVec = [size(CurrentData, 1), size(CurrentData, 2), size(CurrentData, 3)];
        % 先划分长的那一边，如果差不多一样长则先按照最小误差原则选择划分维度
        % Divide the longer side first. If they are almost the same length, choose the dimension to divide according to the principle of minimum error.
        [Max, IMax] = max(SizeVec);
        [Min, ~] = min(SizeVec);
        if Max > Min * MaximumAspectRatio
            choose_dim = IMax;
            if IMax == 1
                SplitPoint = divideArrayIntoHalf(YData);
            elseif IMax == 2
                SplitPoint = divideArrayIntoHalf(XData);
            else
                SplitPoint = divideArrayIntoHalf(ZData);
            end
        else
            SplitY = divideArrayIntoHalf(YData);
            DeltaY = abs(sum(YData(1:SplitY - 1)) - sum(YData(SplitY:end)));
            SplitX = divideArrayIntoHalf(XData);
            DeltaX = abs(sum(XData(1:SplitX - 1)) - sum(XData(SplitX:end)));
            SplitZ = divideArrayIntoHalf(ZData);
            DeltaZ = abs(sum(ZData(1:SplitZ - 1)) - sum(ZData(SplitZ:end)));
            Split = [SplitY, SplitX, SplitZ];
            [~, choose_dim] = min([DeltaY, DeltaX, DeltaZ]);
            SplitPoint = Split(choose_dim);
        end
        % 求分割点
        % Find the split point
        if choose_dim == 1            
            left_data = CurrentData(1:SplitPoint - 1, :, :);
            right_data = CurrentData(SplitPoint:end, :, :);
            left_index = CurrentIndexMat(1:SplitPoint - 1, :, :);
            right_index = CurrentIndexMat(SplitPoint:end, :, :);
        elseif choose_dim == 2
            left_data = CurrentData(:, 1:SplitPoint - 1, :);
            right_data = CurrentData(:, SplitPoint:end, :);
            left_index = CurrentIndexMat(:, 1:SplitPoint - 1, :);
            right_index = CurrentIndexMat(:, SplitPoint:end, :);
        else
            left_data = CurrentData(:, :, 1:SplitPoint - 1);
            right_data = CurrentData(:, :, SplitPoint: end);
            left_index = CurrentIndexMat(:, :, 1:SplitPoint - 1);
            right_index = CurrentIndexMat(:, :, SplitPoint:end);
        end
        kd_tree(Currentindex).LeftData = left_data;
        kd_tree(Currentindex).RightData = right_data;
        kd_tree(Currentindex).LeftRange = left_index;
        kd_tree(Currentindex).RightRange = right_index;
        kd_tree(Currentindex).SplitDim = choose_dim;
        kd_tree(Currentindex).SplitPoint = SplitPoint;
        if kd_tree(Currentindex).Depth == target_depth
            return;
        end
        % 递归进行左右子树的构建
        % Recursively construct the left and right subtrees
        kd_tree(Currentindex).LeftSon = buildKDTree_recurse(left_data, left_index, Currentindex);
        kd_tree(Currentindex).RightSon = buildKDTree_recurse(right_data, right_index, Currentindex);
    end
end
%%
function kd_tree = buildKDTree_2D(Data, target_depth, MaximumAspectRatio)
assert(target_depth < 20 && target_depth > 0 && target_depth == fix(target_depth));
XLength = size(Data, 2);
YLength = size(Data, 1);
[XX, YY] = meshgrid(1:XLength, 1:YLength);
IndexMat = sub2ind([YLength, XLength], YY, XX);
Index = 0;
NumOfNodes = 2^target_depth - 1;
kd_tree = struct('Father', [], 'LeftSon', [], 'RightSon', [], 'SplitDim', [], 'SplitPoint', [], 'Depth', [], ...
    'LeftRange', cell(NumOfNodes, 1), 'RightRange', cell(NumOfNodes, 1), 'LeftData', cell(NumOfNodes, 1), 'RightData', cell(NumOfNodes, 1));
buildKDTree_recurse(Data, IndexMat, []);
% 递归子函数
% Recursive subfunction
    function Currentindex = buildKDTree_recurse(CurrentData, CurrentIndexMat, PreviousIndex)
        Currentindex = Index + 1;
        Index = Index + 1;
        kd_tree(Currentindex).Father = PreviousIndex;
        if ~isempty(PreviousIndex)
            kd_tree(Currentindex).Depth = kd_tree(PreviousIndex).Depth + 1;
        else
            kd_tree(Currentindex).Depth = 1;
        end
        % 数据降维
        % Data dimensionality reduction
        XData = sum(CurrentData, 1);
        YData = sum(CurrentData, 2);
        % 先划分长的那一边，如果差不多一样长则先按照最小误差原则选择划分维度
        % Divide the longer side first. If they are almost the same length, choose the dimension to divide according to the principle of minimum error.
        if size(CurrentData, 1) > size(CurrentData, 2) * MaximumAspectRatio
            choose_dim = 1;
            SplitPoint = divideArrayIntoHalf(YData);
        elseif size(CurrentData, 1) * MaximumAspectRatio < size(CurrentData, 2)
            choose_dim = 2;
            SplitPoint = divideArrayIntoHalf(XData);
        else
            SplitY = divideArrayIntoHalf(YData);
            DeltaY = abs(sum(YData(1:SplitY - 1)) - sum(YData(SplitY:end)));
            SplitX = divideArrayIntoHalf(XData);
            DeltaX = abs(sum(XData(1:SplitX - 1)) - sum(XData(SplitX:end)));
            Split = [SplitY, SplitX];
            [~, choose_dim] = min([DeltaY, DeltaX]);
            SplitPoint = Split(choose_dim);
        end
        % 求分割点
        % Find the split point
        if choose_dim == 1
            left_data = CurrentData(1:SplitPoint - 1, :, :);
            right_data = CurrentData(SplitPoint: end, :, :);
            left_index = CurrentIndexMat(1:SplitPoint - 1, :, :);
            right_index = CurrentIndexMat(SplitPoint: end, :, :);
        else
            left_data = CurrentData(:, 1:SplitPoint - 1, :);
            right_data = CurrentData(:, SplitPoint: end, :);
            left_index = CurrentIndexMat(:, 1:SplitPoint - 1, :);
            right_index = CurrentIndexMat(:, SplitPoint: end, :);
        end
        kd_tree(Currentindex).LeftData = left_data;
        kd_tree(Currentindex).RightData = right_data;
        kd_tree(Currentindex).LeftRange = left_index;
        kd_tree(Currentindex).RightRange = right_index;
        kd_tree(Currentindex).SplitDim = choose_dim;
        kd_tree(Currentindex).SplitPoint = SplitPoint;
        if kd_tree(Currentindex).Depth == target_depth
            return;
        end  
        % 递归进行左右子树的构建
        % Recursively construct the left and right subtrees
        kd_tree(Currentindex).LeftSon = buildKDTree_recurse(left_data, left_index, Currentindex);
        kd_tree(Currentindex).RightSon = buildKDTree_recurse(right_data, right_index, Currentindex);
    end
end
%% 将一维数组分成尽量均匀的两半（两个左闭右开区间，Mid处元素位于右区间内）
%% Divide the one-dimensional array into two halves as evenly as possible (two left closed and right open intervals, the element at Mid is in the right interval)
function Result = divideArrayIntoHalf(array)
L = length(array);
First = 1;
Last = L;
while First < Last
    Mid = First + floor((Last - First) / 2);
    if sum(array(1:Mid)) < sum(array(Mid:L))
        First = Mid + 1 ;
    else
        Last = Mid ;
    end
end
Result = First;
end