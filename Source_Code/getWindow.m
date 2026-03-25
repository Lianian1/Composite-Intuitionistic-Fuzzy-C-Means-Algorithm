function window = getWindow(image, row, col, windowSize)
    % 确保传入的坐标和窗口尺寸是正整数
    if any([row, col, windowSize] < 1)
        error('Coordinates and window size must be positive integers.');
    end
    
    % 计算窗口大小的一半
    halfWindowSize = floor(windowSize/2);
    
    % 获取图像的尺寸
    [numRows, numCols] = size(image);

    % 计算窗口的起始和结束索引
    startX = max(1, row - halfWindowSize(1));
    endX = min(numRows, row + halfWindowSize(1));
    startY = max(1, col - halfWindowSize(2));
    endY = min(numCols, col + halfWindowSize(2));

    % 提取窗口
    window = image(startX:endX, startY:endY);
end
