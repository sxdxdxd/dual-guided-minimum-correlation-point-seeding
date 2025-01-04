function S = generate3DStreamlinesBySeeds(x, y, z, u, v, w, SeedPointList, Step, N_VerticesMax)
StreamLineNum = size(SeedPointList, 1);
SeedPointListX = SeedPointList(:, 1);
SeedPointListY = SeedPointList(:, 2);
SeedPointListZ = SeedPointList(:, 3);
ThresholdNum = 1e4;
%% 计算流线数据
%% Calculate streamline
if StreamLineNum < ThresholdNum
    StreamLine_positive = stream3(x, y, z, u, v, w, SeedPointListX, SeedPointListY, SeedPointListZ, [Step, N_VerticesMax]);%正向流线计算
    StreamLine_negative = stream3(x, y, z, -u, -v, -w, SeedPointListX, SeedPointListY, SeedPointListZ, [Step, N_VerticesMax]);%反向流线计算
else
    ThreadNum = maxNumCompThreads;
    StreamLineNumPerThread = floor(StreamLineNum / ThreadNum);
    ParLength = StreamLineNumPerThread * maxNumCompThreads;
    SeedPointMatX = reshape(SeedPointListX(1:ParLength), ThreadNum, []);
    SeedPointMatY = reshape(SeedPointListY(1:ParLength), ThreadNum, []);
    SeedPointMatZ = reshape(SeedPointListZ(1:ParLength), ThreadNum, []);
    StreamLine_positive = cell(ThreadNum, StreamLineNumPerThread);
    StreamLine_negative = cell(ThreadNum, StreamLineNumPerThread);
    parfor i = 1:ThreadNum
        TempStreamLine_positive = stream3(x, y, z, u, v, w, SeedPointMatX(i, :), SeedPointMatY(i, :), SeedPointMatZ(i, :), [Step, N_VerticesMax]);%正向流线计算
        TempStreamLine_negative = stream3(x, y, z, -u, -v, -w, SeedPointMatX(i, :), SeedPointMatY(i, :), SeedPointMatZ(i, :), [Step, N_VerticesMax]);%反向流线计算
        [StreamLine_positive{i, :}] = deal(TempStreamLine_positive{:});
        [StreamLine_negative{i, :}] = deal(TempStreamLine_negative{:});
    end
    StreamLine_positive = reshape(StreamLine_positive, 1, []);
    StreamLine_negative = reshape(StreamLine_negative, 1, []);
    if ParLength < StreamLineNum
        t = ParLength + 1;
        TempStreamLine_positive0 = stream3(x, y, z, u, v, w, SeedPointListX(t:end), SeedPointListY(t:end), SeedPointListZ(t:end), [Step, N_VerticesMax]);%正向流线计算
        TempStreamLine_negative0 = stream3(x, y, z, -u, -v, -w, SeedPointListX(t:end), SeedPointListY(t:end), SeedPointListZ(t:end), [Step, N_VerticesMax]);%反向流线计算
        StreamLine_positive = [StreamLine_positive, TempStreamLine_positive0];
        StreamLine_negative = [StreamLine_negative, TempStreamLine_negative0];
    end
end
%% 将正向和反向数据合二为一,流线与其种子点一起储存
%% Combine forward and reverse streamline data (Streamlines are stored with their seed points)
StreamLine = cell(StreamLineNum, 2);
if StreamLineNum < ThresholdNum
    for i = 1:StreamLineNum
        if ~isempty(StreamLine_positive{i}) && ~isempty(StreamLine_negative{i})
            StreamLine{i, 1} = [[flipud(StreamLine_negative{i}(:, 1));StreamLine_positive{i}(:, 1)],...
                [flipud(StreamLine_negative{i}(:, 2));StreamLine_positive{i}(:, 2)]...
                [flipud(StreamLine_negative{i}(:, 3));StreamLine_positive{i}(:, 3)]];%拼在一起
        elseif ~isempty(StreamLine_positive{i})
            StreamLine{i, 1} = [StreamLine_positive{i}(:, 1), StreamLine_positive{i}(:, 2), StreamLine_positive{i}(:, 3)];
        elseif ~isempty(StreamLine_negative{i})
            StreamLine{i, 1} = [StreamLine_negative{i}(:, 1), StreamLine_negative{i}(:, 2), StreamLine_negative{i}(:, 3)];
        end
        StreamLine{i, 2} = SeedPointList(i,:);
    end
else
    parfor i = 1:StreamLineNum
        if ~isempty(StreamLine_positive{i}) && ~isempty(StreamLine_negative{i})
            StreamLine{i, 1} = [[flipud(StreamLine_negative{i}(:, 1));StreamLine_positive{i}(:, 1)],...
                [flipud(StreamLine_negative{i}(:, 2));StreamLine_positive{i}(:, 2)]...
                [flipud(StreamLine_negative{i}(:, 3));StreamLine_positive{i}(:, 3)]];%拼在一起
        elseif ~isempty(StreamLine_positive{i})
            StreamLine{i, 1} = [StreamLine_positive{i}(:, 1), StreamLine_positive{i}(:, 2), StreamLine_positive{i}(:, 3)];
        elseif ~isempty(StreamLine_negative{i})
            StreamLine{i, 1} = [StreamLine_negative{i}(:, 1), StreamLine_negative{i}(:, 2), StreamLine_negative{i}(:, 3)];
        end
    end
    parfor i = 1:StreamLineNum
        StreamLine{i, 2} = SeedPointList(i, :);
    end
end
%% 去除空流线
%% Remove empty streamlines
StreamLine(cellfun(@isempty,StreamLine(:, 1)), :) = [];
S = StreamLine;
end