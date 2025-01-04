function S = generate2DStreamlinesBySeeds(x, y, u, v, SeedPointList, Step, N_VerticesMax)
%% 计算流线数据
%% Calculate streamline
StreamLine_positive = stream2(x, y, u, v, SeedPointList(:, 1), SeedPointList(:, 2), [Step, N_VerticesMax]);%正向流线计算
StreamLine_negative = stream2(x, y, -u, -v, SeedPointList(:, 1), SeedPointList(:, 2), [Step, N_VerticesMax]);%反向流线计算
%% 分配储存空间
%% Allocate storage
StreamLineNum = size(SeedPointList, 1);
StreamLine = cell(StreamLineNum, 2);
for i = 1:StreamLineNum
%% 将正向和反向数据合二为一
%% Combine forward and reverse streamline data
    if ~isempty(StreamLine_positive{i}) && ~isempty(StreamLine_negative{i})
    StreamLine{i, 1} = [[flipud(StreamLine_negative{i}(:, 1));StreamLine_positive{i}(:, 1)],...
        [flipud(StreamLine_negative{i}(:, 2));StreamLine_positive{i}(:, 2)]];%拼在一起
    elseif ~isempty(StreamLine_positive{i})  
            StreamLine{i, 1} = [StreamLine_positive{i}(:, 1),StreamLine_positive{i}(:, 2)];
    elseif ~isempty(StreamLine_negative{i})
            StreamLine{i, 1} = [StreamLine_negative{i}(:, 1),StreamLine_negative{i}(:, 2)];
    end
    %% 流线与其种子点一起储存
    %% Streamlines are stored with their seed points
     StreamLine{i, 2} = SeedPointList(i,:);
end
%% 去除空流线
%% Remove empty streamlines
StreamLine(cellfun(@isempty,StreamLine(:, 1)), :) = [];
S = StreamLine;
end