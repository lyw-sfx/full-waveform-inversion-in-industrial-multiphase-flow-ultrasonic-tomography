function [D, A] = createTVMatrix(M, N, varargin)
    % CREATETVMATRIX 生成任意大小的二维Total Variation (TV)正则化矩阵。
    %
    % 输入:
    %   M, N - 图像的行数和列数。
    %   varargin - 可选参数。可以是 'isotropic' (默认) 或 'anisotropic'。
    %              - 'isotropic'：各向同性TV，组合x和y方向梯度 (D = [Dx; Dy])。
    %              - 'anisotropic'：各向异性TV，分别返回x和y方向梯度矩阵 (D = [Dx, Dy])。
    %
    % 输出:
    %   D - 梯度算子矩阵。
    %       如果是各向同性，D是一个 (M*N*2 - M - N) x (M*N) 的稀疏矩阵。
    %       如果是各向异性，D是一个细胞数组，D{1}=Dx, D{2}=Dy。
    %   A - TV正则化矩阵，即 D' * D。它是一个 (M*N) x (M*N) 的稀疏矩阵。
    %       当使用各向异性TV时，A也是一个细胞数组。
    %
    % 示例:
    %   % 生成一个3x3图像的各向同性TV矩阵
    %   [D, A] = createTVMatrix(3, 3);
    %
    %   % 生成一个5x5图像的各向异性TV矩阵
    %   [D_aniso, A_aniso] = createTVMatrix(5, 5, 'anisotropic');

    % 检查输入参数
    if nargin < 2
        error('至少需要输入图像的行数M和列数N。');
    end
    
    if ~isscalar(M) || ~isscalar(N) || M <= 0 || N <= 0 || M ~= round(M) || N ~= round(N)
        error('M和N必须是正整数。');
    end

    % 设置默认模式为各向同性
    mode = 'isotropic';
    if nargin >= 3
        mode = varargin{1};
        if ~strcmp(mode, 'isotropic') && ~strcmp(mode, 'anisotropic')
            error('可选参数只能是 ''isotropic'' 或 ''anisotropic''。');
        end
    end

    % 计算一维索引
    [I, J] = meshgrid(1:N, 1:M);
    idx = sub2ind([M, N], J, I);

    % ---------------------- 1. 构建x方向梯度矩阵 Dx ----------------------
    % Dx作用于图像向量x时，Dx*x的结果是一个向量，其中每个元素是对应位置的水平差分。
    % 即 (Dx*x)(k) = x(i,j+1) - x(i,j)，其中k是(i,j)的一维索引。
    % 注意，最后一列没有水平差分。
    rows_x = M * (N - 1);
    cols = M * N;
    Dx = sparse(rows_x, cols);
    
    % 非零元素的行索引
    row_idx_x = kron(1:rows_x, ones(1, 2)); 
    
    % 非零元素的列索引
    col_idx_x_1 = idx(:, 1:end-1);
    col_idx_x_2 = idx(:, 2:end);
    col_idx_x = [col_idx_x_1(:); col_idx_x_2(:)];
    
    % 非零元素的值 (+1 和 -1)
    vals_x = [-ones(rows_x, 1); ones(rows_x, 1)];
    
    % 构建稀疏矩阵Dx
    Dx = sparse(row_idx_x, col_idx_x, vals_x, rows_x, cols);

    % ---------------------- 2. 构建y方向梯度矩阵 Dy ----------------------
    % Dy作用于图像向量x时，Dy*x的结果是一个向量，其中每个元素是对应位置的垂直差分。
    % 即 (Dy*x)(k) = x(i+1,j) - x(i,j)，其中k是(i,j)的一维索引。
    % 注意，最后一行没有垂直差分。
    rows_y = (M - 1) * N;
    Dy = sparse(rows_y, cols);
    
    % 非零元素的行索引
    row_idx_y = kron(1:rows_y, ones(1, 2));
    
    % 非零元素的列索引
    col_idx_y_1 = idx(1:end-1, :);
    col_idx_y_2 = idx(2:end, :);
    col_idx_y = [col_idx_y_1(:); col_idx_y_2(:)];
    
    % 非零元素的值 (+1 和 -1)
    vals_y = [-ones(rows_y, 1); ones(rows_y, 1)];
    
    % 构建稀疏矩阵Dy
    Dy = sparse(row_idx_y, col_idx_y, vals_y, rows_y, cols);

    % ---------------------- 3. 根据模式组合矩阵 ----------------------
    if strcmp(mode, 'isotropic')
        % 各向同性TV: D = [Dx; Dy]
        D = [Dx; Dy];
        % 正则化矩阵 A = D' * D
        A = D' * D;
    else % anisotropic
        % 各向异性TV: D是一个包含Dx和Dy的细胞数组
        D = {Dx, Dy};
        % 正则化矩阵A也是一个细胞数组
        A = {Dx' * Dx, Dy' * Dy};
    end

end