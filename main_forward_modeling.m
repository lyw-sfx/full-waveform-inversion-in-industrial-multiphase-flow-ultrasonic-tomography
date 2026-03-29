%% 初始化 FWI forward
close all; 
clear; 
clc;
% set(groot, 'DefaultFigureColormap', jet);
%% 正问题 网格结构体

c_min = 1200;           % 设置介质最小声速
points_per_length = 3;  % 每波长网格点数3 (非均匀介质)
f = 1e6;              % 超声最大频率
diameter = 0.15;        % 直径150mm (单位:mm)

FW_dx = c_min/(points_per_length*f);  %网格步长 一个波长内至少3个网格点
FW_dx = FW_dx ;
Fw_Nx = round(diameter/FW_dx);                % 网格数取整
FW_N = 410;                                       % 选择2或3的倍数加匹配层
FW_kgrid = kWaveGrid(FW_N, FW_dx, FW_N, FW_dx);     % 初始化k-wave网格对象

%% 正问题 介质属性结构体

% 初始化全区域先设为水
% f = 2e6;                  % 超声中心频率 (Hz)
T = 20;                        % 温度
% c_water = waterSoundSpeed(T);  % 水中声速 (m/s)
c_water = 1500;
rho_water = waterDensity(T) ;  % 水的密度 (kg/m3)  
alpha_water_ = waterAbsorption(f/1e6, T);      % 水中衰减系数
alpha_water = alpha_water_/(f/1e6^1.43) ;    % 单位转换(dB/(MHz^y)*cm)

medium_refw.sound_speed = c_water * ones(FW_N, FW_N); % 声速分布矩阵m
medium_refw.density = rho_water * ones(FW_N, FW_N);   % 密度分布矩阵kg 
medium_refw.alpha_coeff = alpha_water* ones(FW_N, FW_N);   
                                                      % 吸收系数
medium_refw.alpha_power = 1.43;                       % 吸收指数
                                    % 吸收指数在0-3之间且不等于1
                                    % 可以默认选1.42
                                    % 如果过大则计算错误
% figure;
% imagesc(medium.sound_speed);
% colormap;

%% 正问题 设置仿真时间步长

% 设置结束时间: 波传播对角线距离+n个最长周期
% 好像不用加 50 个周期 10
% n = 0;
% t_end = FW_N*sqrt(2)*FW_dx / c_min + n/f;
t_end = 1.5e-4; % 设置为150us 与实验一致
% 使用 maketime 创建时间数组并设置结束时间
% 可以先设置完整个 medium 再设置 t_array 
% FW_kgrid.t_array = makeTime(FW_kgrid,medium_refw.sound_speed,0.3,t_end);
% 因为正逆问题网格数不一致，所以这里不同maketimetime
FW_kgrid.t_array = 0:1/(20*f):t_end; % 50MHz采样频率
% FW_kgrid.t_array = 0:1/(50*e6):t_end; % 50MHz采样频率
%% 正问题 声源结构体

% 设置声源位置
sensor_radius = 0.075;   % 半径[m]
num_sensor = 32;  % 激励传感器数量
num_rec = 7;      % 接收传感器数量
source_cart = makeCartCircle(sensor_radius, ...
                             num_sensor, ...
                             [0,0], ...
                             2*pi, ...
                             false);
% 脉冲周期数
periods = 1;
% 设置信号峰值 手动放缩设置峰值
target_amplitude = 10;
% 设置采样频率，和kgrid保持一致，满足奈奎斯特采样定理
fs = 1 / (FW_kgrid.dt);
% 生成正弦脉冲信号（k-wave函数，基于采样率、频率、周期数）
% 信号长度大于会自动补0
p = target_amplitude*toneBurst(10*fs, f, periods,'Plot', false);
% 'SignalLength','auto'
% plot_single_sided_amplitude_spectrum(p, 10*fs);
%% 正问题 定义传感器结构体

sensor_radius = 0.075;   % [m]
sensor_positions = makeCartCircle(sensor_radius, ...
                             num_sensor, ...
                             [0,0], ...
                             2*pi, ...
                             false);

%% 正问题 ref_water

% 硬件加速暂时不要用 
fw_input_args = {'DataCast' , 'gpuArray-single',...  % 'gpuArray-single' 
                 'PlotSim'  , false,...   % 仿真过程
                 'PMLSize'  , 20,...      % 完美匹配层大小
                 'PMLInside', false...    % 设置完美匹配层位置
                 'PlotPML'  , false,...
                 'RecordMovie',false
                };

%%
% 计算纯水模型
measured_data_refw = run_simulation_forward( ...
                     FW_kgrid, medium_refw, ...
                     source_cart, p, ...
                     sensor_positions, num_rec, num_sensor,...
                     fw_input_args,1);

data_show(measured_data_refw,num_rec, ...
    'Sensor Position','Time Step','全水场域');

%% 正问题 油滴模型

% 构建网格掩码
x_vec = FW_kgrid.x_vec;              
y_vec = FW_kgrid.y_vec;              
[X, Y] = meshgrid(x_vec, y_vec);

c_oil = 1200;                    % 油中声速 (m/s)
rho_oil = 800;                   % 油的密度 (kg/m3)
% alpha_power_oil = 1.43;
% alpha_oil_ = 0.0001;                    % 油中衰减系数
% alpha_oil = alpha_oil_ /(f/1e6^1.43); % 单位转换(dB/(MHz^y)*cm)
alpha_oil = 5*alpha_water;
medium_oil = medium_refw;

oil_center = [0.0, 0.0];      % 油滴中心
oil_radius = 0.01;              % 油滴半径
% 生成包含物掩码 判断网格点是否在包含物内部
inclusion_dis = ((X - oil_center(1)).^2 + (Y - oil_center(2)).^2);
inclusion_mask = inclusion_dis <= oil_radius^2;
% 更新包含物区域的介质属性
medium_oil.sound_speed(inclusion_mask) = c_oil;      % 包含物声速
medium_oil.density(inclusion_mask) = rho_oil;        % 包含物密度
medium_oil.alpha_coeff(inclusion_mask) = alpha_oil;  % 包含物衰减系数

% oil_center = [-0.025, 0];      % 油滴中心
% oil_radius = 0.015;              % 油滴半径
% 生成包含物掩码 判断网格点是否在包含物内部
% inclusion_dis = ((X - oil_center(1)).^2 + (Y - oil_center(2)).^2);
% inclusion_mask = inclusion_dis <= oil_radius^2;
% 更新包含物区域的介质属性
% medium_oil.sound_speed(inclusion_mask) = c_oil;      % 包含物声速
% medium_oil.density(inclusion_mask) = rho_oil;        % 包含物密度
% medium_oil.alpha_coeff(inclusion_mask) = alpha_oil;  % 包含物衰减系数

% 计算油滴模型
measured_data_oil = run_simulation_forward(FW_kgrid, medium_oil, ...
                        source_cart, p, ...
                        sensor_positions, num_rec, num_sensor,...
                        fw_input_args,0);

[row_measured_data,col_measured_data] = size(measured_data_oil);
%%
data_show(measured_data_oil,num_rec, ...
    'Sensor Position','Time Step','油水混合场域');

%% fwi 逆问题初始化

%% fwi 逆问题网格结构体初始化
% 用于验证网格无关性-验证没通过所以网格数一样
INV_dx = FW_dx;
INV_N = FW_N;          
INV_kgrid = kWaveGrid(INV_N, INV_dx, INV_N, INV_dx);

%%
% 设置仿真初始场域
medium_rec = medium_oil;

%%
% 定义逆问题传感器
% sensor_allfield = zeros(INV_N,INV_N);
% 直接赋逻辑值
% sa_index = round(INV_N/6)+1;
% sensor_allfield(sa_index:INV_N-sa_index,sa_index:INV_N-sa_index) = true;
sensor_allfield = ones(INV_N,INV_N);
%%
% 设置结束时间: 波传播对角线距离+n个最长周期
% 好像不用加 50 个周期 10
% n = 5;
% t_end = FW_N*sqrt(2)*FW_dx / c_min + n/f;
% 使用 maketime 创建时间数组并设置结束时间
% 可以先设置完整个 medium 再设置 t_array 
% INV_kgrid.t_array = makeTime(INV_kgrid,medium_rec.sound_speed,0.3,t_end);
% INV_kgrid.t_array = 0:1/(10*f):t_end;  % 生成时间数组

% 设置结束时间
t_end = 1.2e-4;                       % 设置为150us 与实验一致
INV_kgrid.t_array = 0:1/(20*f):t_end; % 50MHz采样频率

%%
% 设置可选参数
iv_input_args = {'DataCast', 'single',...   % 'off''gpuArray-single'
                 'PlotSim', false,...    % 仿真过程
                 'PMLSize', 20,...       % 完美匹配层大小
                 'PMLInside', false,...   % 设置完美匹配层位置
                 'PlotPML', false,...     % 
                 'PlotLayout', false,...
                 'DisplayMask', 'off'
                };

%%

num_sensor_i = 1;
medium_rec_adj = medium_rec;
    
% 构建正向传播声源结构体
source_cart_temp = source_cart( : , num_sensor_i);
forward_source.p_mask = cart2grid(INV_kgrid, source_cart_temp);
forward_source.p = p;
% 正向传播用于梯度计算 需记录整个场声压变化
forward_sensor_1.mask = sensor_allfield;

% 正向传播k-wave仿真
forward_data_1 = kspaceFirstOrder2D(INV_kgrid, medium_rec_adj, ...
                     forward_source, forward_sensor_1, ...
                     iv_input_args{:});
%%
% save_matrix_500step_cols(forward_data_1, 'F:\net_practice\k_wave\data_0112\o_15', 'phan3', 'single')
%%
% for j=556800:6:900000
%     iner_show = forward_data_1(j,:);
%     plot(iner_show);
% 
% end
% filtered_img = mean_filter_2d(adjoint_grad, 3);
% max_adj = max(filtered_img);
% min_adj = min(filtered_img);
% nomal_adj = (filtered_img - min_adj)./(max_adj- min_adj);
% filtered_img = mean_filter_2d(adjoint_grad, 5);
%% k-wave 正问题仿真函数

% 串行激励并行接收 fwi正演函数
function measured_data = run_simulation_forward(kgrid, medium, ...
                         source_position, source_p, ...
                         sensor_positions, num_rec, num_sensor,...
                         input_args, z)
% 输入参数：
%   kgrid            - 场分割参数 结构体
%   medium           - 介质参数 结构体
%   source_position  - 声源位置物理坐标 非结构体 
%   source_p         - 声源性质 非结构体
%   sensor_positions - 传感器位置物理坐标 非结构体
%   num_rec          - 接收传感器数量
%   num_sensor       - 总传感器数量
%   input_args       - k-wave可选参数

% 输出参数：
%   measured_data    - 所有测量值

for i = 1 : 1 : num_sensor
    % 构建接-收传感器索引
    index_rec_center = mod(i + num_sensor/2, num_sensor);
    if (index_rec_center == 0)
        index_rec_center = num_sensor;
    end
    index_rec_num = (num_rec-1)/2;
    
    index_1 = index_rec_center-index_rec_num ;
    index_2 = index_rec_center+index_rec_num ;
    index_rec = mod(index_1:index_2,num_sensor) ;
    index_rec(index_rec == 0) = num_sensor;
    % 构建声源结构体
    source_cart_temp = source_position( : , i);
    source_temp.p_mask = cart2grid(kgrid, source_cart_temp);
    source_temp.p = source_p;
    % 构建传感器结构体
    sensor_positions_temp = sensor_positions(:,index_rec);
    sensor_temp.mask = sensor_positions_temp;
    % sensor_temp.record = {'p','p_min_all'};

    % 仿真
    sensor_data_gpu = kspaceFirstOrder2D(kgrid, medium, ...
                                     source_temp, sensor_temp, ...
                                     input_args{:});
    % 定位'DataCast'键的位置（忽略大小写）
    data_cast_idx = find(strcmpi(input_args, 'DataCast'));

    % 判断是否存在该键且值为'gpuArray-single'
    if ~isempty(data_cast_idx)&&strcmp(input_args{data_cast_idx+1},...
            'gpuArray-single')
        disp(" 有GPU加速 ");
        sensor_data = gather(sensor_data_gpu);
        measured_data(:, (i-1)*num_rec+1:num_rec*i) = sensor_data';
        clear sensor_data_gpu;
    else
        disp(" 无GPU加速 ");
        measured_data(:, (i-1)*num_rec+1:num_rec*i) = sensor_data_gpu';
    end
    %disp("无GPU加速");
    %measured_data(:, (i-1)*5+1:5*i) = sensor_data';
    if (z == 1)
        fprintf("仅计算一个角度");
        break;
    end
end

end

%% 接收波形展示函数

function data_show(data,num_rec,y_name,x_name,pic_name)
% data 每列为一个传感器接收
% num_rec 接收传感器个数
% y轴标签 默认'Sensor Position'
% x轴标签 默认'Time Step'
figure('Name',pic_name);
imagesc(data(:,1:num_rec)');
%colormap(getColorMap); %别用这个函数一用就闪退
ylabel(y_name);
xlabel(x_name);
colormap("jet");
colorbar;

end

%% 二阶差分函数

function forward_ddt = secondOrderDiff(forward_in, dt)
    % secondOrderDiff: 快速计算一维/二维序列的二阶差分
    % 输入：
    %   forward_in: 时间序列（2维数组[行×时间列] 每行对应一个时间序列）
    %   dt: 时间步长（采样间隔，标量）
    % 输出：
    %   forward_ddt: 二阶差分结果（与forward_in同维度、同数据类型）
    
    % 步骤1：预处理参数，统一数组维度，避免重复计算
    [M, N] = size(forward_in);  % M=行数（1维时M=1），N=时间列数
    dt_sq = dt^2;  % 提前计算dt的平方，减少冗余运算
    forward_ddt = zeros(size(forward_in), class(forward_in)); 
    
    % 步骤2：分情况处理（1维/2维，保留原计算逻辑）
    if N >= 3
        % 核心优化1：二维数组批量向量化切片（替代外层行循环+内层点循环）
        % 对所有行的内部点（列2~N-1）一次性完成中心差分计算
        forward_ddt(:, 2:N-1) = (forward_in(:, 3:end) - ...
            2*forward_in(:, 2:end-1) + ...
            forward_in(:, 1:end-2)) / dt_sq;
        
        % 核心优化2：批量处理所有行的首尾边界点（无循环）
        forward_ddt(:, 1) = (forward_in(:, 3) - ...
            2*forward_in(:, 2) + ...
            forward_in(:, 1)) / dt_sq;
        forward_ddt(:, end) = (forward_in(:, end) - ...
            2*forward_in(:, end-1) + forward_in(:, end-2)) / dt_sq;
    else
        % 序列过短（N<3），用diff的二次差分近似，支持二维数组
        if M == 1
            % 1维序列，直接按原逻辑处理
            ddt_temp = diff(diff(forward_in)) / dt_sq;
            forward_ddt(1:length(ddt_temp)) = ddt_temp;
        else
            % 2维序列，逐行进行二次差分（此时N<3，行数多也无太大开销）
            for i = 1:M
                ddt_temp = diff(diff(forward_in(i, :))) / dt_sq;
                forward_ddt(i, 1:length(ddt_temp)) = ddt_temp;
            end
        end
    end
end
%% 信号滤波函数

function filtered_signal = filterSignal(signal_in, fs, filter_type, ...
    cutoff_freq, order)
    % 对矩阵信号的每一列进行滤波处理
    % 输入参数:
    %   signal      - 输入信号 (M x N 矩阵, M 为信号长度, N 为通道数)
    %   fs          - 采样频率 (Hz)
    %   filter_type - 滤波类型: 'low'(低通), 'high'(高通), 
    %                          'bandpass'(带通), 'bandstop'(带阻)
    %   cutoff_freq - 截止频率 (Hz):
    %                 - 低通/高通: 标量 (如 100 表示 100Hz)
    %                 - 带通/带阻: 二元素向量 (如 [50, 200] 表示 50-200Hz)
    %   order       - 滤波器阶数 (正整数, 默认 4)
    %
    % 输出参数:
    %   filtered_signal - 滤波后的信号 (M x N 矩阵, 与输入信号维度相同)

    % 参数校验与默认值设置
    signal = double(signal_in);
    if nargin < 5
        order = 4; % 默认 4 阶滤波器
    end
    if ~ismatrix(signal) || ~isnumeric(signal)
        error('输入信号必须是数值矩阵');
    end
    if fs <= 0
        error('采样频率必须为正数');
    end
    valid_types = {'low', 'high', 'bandpass', 'bandstop'};
    if ~ismember(filter_type, valid_types)
        error(['滤波类型必须是 ''low'', ''high'',' ...
            ' ''bandpass'' 或 ''bandstop''']);
    end

    % 设计滤波器 (只需要设计一次)
    nyquist = fs / 2;
    Wn = cutoff_freq / nyquist;
    [b, a] = butter(order, Wn, filter_type);

    % 获取输入信号的维度
    [num_samples, num_channels] = size(signal);

    % 初始化输出矩阵
    filtered_signal = zeros(size(signal), class(signal));

    % 对每一列信号进行滤波
    for i = 1:num_channels
        % 提取第 i 列信号
        col_signal = signal(:, i);

        % 使用 filtfilt 进行零相位滤波
        filtered_col = filtfilt(b, a, col_signal);

        % 将滤波后的列存入输出矩阵
        filtered_signal(:, i) = filtered_col;
    end

end

%% 图像滤波/平滑函数

function filtered_img = mean_filter_2d(input_img, filter_size, varargin)
    % 二维均值滤波函数，使用 filter2 实现
    % 输入参数：
    %   input_img   - 输入图像（灰度图：M×N 矩阵；RGB图：M×N×3 矩阵）
    %   filter_size - 滤波器大小，指定为标量 n（生成 n×n 滤波器）或 
    %                                      [m, n]（生成 m×n 滤波器）
    %   varargin    - 可选参数，格式为 'shape', shape_type
    %                 shape_type: 'same'（默认，输出与输入大小相同）、
    %                             'valid'（仅输出完全重叠区域）、
    %                             'full'（输出完整卷积结果）
    % 输出参数：
    %   filtered_img - 滤波后的图像，与输入图像类型、维度一致
    
    % ---------------------- 1. 输入参数验证与处理 ----------------------
    % 验证输入图像合法性
    if ~ismatrix(input_img) && ~ndims(input_img) == 3
        error('输入图像必须是二维（灰度图）或三维（RGB图）矩阵');
    end
    if ndims(input_img) == 3 && size(input_img, 3) ~= 3
        error('RGB图像必须是 M×N×3 矩阵');
    end
    
    % 验证并处理滤波器大小
    if isscalar(filter_size)
        % 标量输入：生成 n×n 方阵
        n = filter_size;
        % 修改点：检查值是否为整数，而非类型
        if n <= 0 || mod(n, 1) ~= 0
            error('滤波器大小必须是正整数');
        end
        filter_rows = n;
        filter_cols = n;
    elseif isvector(filter_size) && length(filter_size) == 2
        % 向量输入：生成 m×n 矩阵
        [filter_rows,filter_cols]=deal(filter_size(1),filter_size(2));
        % 修改点：检查值是否为整数，而非类型
        if filter_rows<=0||filter_cols<=0||~all(mod(filter_size, 1) == 0)
            error('滤波器大小必须是正整数向量 [m, n]');
        end
    else
        error('滤波器大小必须是正整数标量或 [m, n] 向量');
    end
    
    % 处理可选参数（边缘处理方式）
    if nargin < 3
        shape = 'same';  % 默认输出与输入大小相同
    else
        if ~strcmp(varargin{1}, 'shape') || ~ismember(varargin{2}, ...
                {'same', 'valid', 'full'})
            error('可选参数格为 ''shape'', ''same''/''valid''/''full''');
        end
        shape = varargin{2};
    end
    
    % ---------------------- 2. 生成均值滤波器核 ----------------------
    % 均值滤波器核：每个元素为 1/(m×n)，总和为 1
    mean_kernel=ones(filter_rows,filter_cols)/(filter_rows*filter_cols);
    
    % ---------------------- 3. 对图像进行均值滤波 ---------------------
    if ismatrix(input_img)
        % 灰度图：直接滤波
        img_double = double(input_img);  % 转换为double避免精度损失
        filtered_double = filter2(mean_kernel, img_double, shape);
        % 转换回原图像类型（如uint8、uint16等）
        filtered_img = cast(filtered_double, class(input_img));
    else
        % RGB图：分通道滤波
        filtered_img = zeros(size(input_img), class(input_img));
        for c = 1:3
            img_double = double(input_img(:, :, c));
            filtered_double = filter2(mean_kernel, img_double, shape);
            filtered_img(:,:,c)=cast(filtered_double,class(input_img));
        end
    end
    
    % ---------------------- 4. 可选：裁剪到原图像大小（若shape为'full'）
    % 若用户选择'full'但希望输出与输入大小一致，可取消以下注释
    % if strcmp(shape, 'full')
    %     [m, n, ~] = size(input_img);
    %     filtered_img = imcrop(filtered_img, [1, 1, n, m]);
    % end
end

%% 箱式限制函数

function clipped_matrix = clip_matrix_values(original_matrix, ...
    lower_bound, upper_bound)
    % 输入参数：
    %   original_matrix - 待处理的原始矩阵（支持任意维度的数值矩阵）
    %   lower_bound     - 下限阈值（数值型）
    %   upper_bound     - 上限阈值（数值型）
    % 输出参数：
    %   clipped_matrix  - 裁剪后的矩阵（与原始矩阵维度、数据类型一致）
    
    % ---------------------- 1. 输入参数验证 ----------------------
    % 验证原始矩阵是否为数值型
    if ~isnumeric(original_matrix)
        error('输入矩阵必须是数值型矩阵');
    end
    % 验证上下限阈值是否为数值
    if ~isnumeric(lower_bound) || ~isnumeric(upper_bound)
        error('上下限阈值必须是数值');
    end
    % 验证下限是否小于等于上限
    if lower_bound > upper_bound
        error('下限阈值不能大于上限阈值');
    end
    
    % ---------------------- 2. 执行矩阵值裁剪 ----------------------
    % 方法1：使用逻辑索引（高效，适合大型矩阵）
    clipped_matrix = original_matrix;  % 复制原始矩阵
    clipped_matrix(clipped_matrix < lower_bound) = lower_bound;  
    % 小于下限的设为下限
    clipped_matrix(clipped_matrix > upper_bound) = upper_bound;  
    % 大于上限的设为上限
    
    % ---------------------- 3. 保持数据类型一致性 --------------------
    % 将裁剪后的矩阵转换回原始矩阵的数据类型（避免类型转换导致的精度问题）
    clipped_matrix = cast(clipped_matrix, class(original_matrix));
end

%% 二值化函数

function result_matrix = threshold_assign(original_matrix, ...
    threshold, value_above, value_below)
    % 根据阈值对矩阵元素进行二值化赋值：大于阈值的元素赋value_above，
    %                                小于阈值的元素赋value_below
    % 输入参数：
    %   original_matrix - 待处理的原始矩阵（支持任意维度的数值矩阵）
    %   threshold       - 阈值（数值型，用于判断元素大小的基准）
    %   value_above     - 大于阈值的元素所要赋予的值（数值型）
    %   value_below     - 小于阈值的元素所要赋予的值（数值型）
    % 输出参数：
    %   result_matrix   - 处理后的矩阵（与原始矩阵维度、数据类型一致）
    
    % ---------------------- 1. 输入参数校验 ----------------------
    % 验证原始矩阵是否为数值型
    if ~isnumeric(original_matrix)
        error('输入的 original_matrix 必须是数值型矩阵');
    end
    % 验证阈值、赋值是否为数值
    if (~isnumeric(threshold)|| ...
            ~isnumeric(value_above)|| ...
            ~isnumeric(value_below))
        error('threshold、value_above、value_below 必须是数值型');
    end
    % 验证阈值是否为标量（避免多值阈值导致逻辑混乱）
    if ~isscalar(threshold)
        error('threshold 必须是标量数值');
    end
    
    % ---------------------- 2. 初始化结果矩阵 ----------------------
    % 复制原始矩阵，确保结果矩阵与原始矩阵维度、数据类型完全一致
    result_matrix = original_matrix;
    
    % ---------------------- 3. 根据阈值进行元素赋值 
    % 逻辑索引：找到大于阈值的元素，赋值为 value_above
    result_matrix(original_matrix > threshold) = value_above;
    % 逻辑索引：找到小于阈值的元素，赋值为 value_below
    result_matrix(original_matrix < threshold) = value_below;
    
    
end

%% 灰度化成像结果

% 全局阈值分割二值化（支持任意数值矩阵）
% 功能：将任意数值范围的矩阵按阈值分割为二值矩阵（高值/低值）

function binary_matrix = global_thresholding(any_matrix, ...
    varargin)
    % 输入参数：
    %   any_matrix - 任意数值范围的矩阵（2D/3D，支持double、single等类型
    %   varargin   - 可选参数：
    %              - 第1个可选参数：手动阈值，不提供则自动计算
    %              - 第2个可选参数：高值（大于阈值的元素赋值，默认1）
    %              - 第3个可选参数：低值（小于等于阈值的元素赋值，默认0）
    
    % ---------------------- 1. 输入参数校验与处理 ---------------------
    % 校验输入矩阵合法性
    if ~isnumeric(any_matrix) || isempty(any_matrix)
        error('输入必须是non-empty的数值矩阵');
    end
    % 处理2D/3D矩阵（3D矩阵按通道独立处理，默认最后一维为通道）
    if ndims(any_matrix) > 3
        error('暂不支持3D以上矩阵');
    end
    [m, n, c] = size(any_matrix);
    if c == 1
        c = [];  % 2D矩阵标记为单通道
    end
    
    % 解析可选参数（手动阈值、高低值）
    manual_threshold = [];
    high_val = 1;
    low_val = 0;
    if nargin >= 2 && ~isempty(varargin{1})
        manual_threshold = varargin{1};
        if ~isscalar(manual_threshold) || ~isnumeric(manual_threshold)
            error('手动阈值必须是标量数值');
        end
    end
    if nargin >= 3
        high_val = varargin{2};
        if ~isscalar(high_val) || ~isnumeric(high_val)
            error('高值必须是标量数值');
        end
    end
    if nargin >= 4
        low_val = varargin{3};
        if ~isscalar(low_val) || ~isnumeric(low_val)
            error('低值必须是标量数值');
        end
    end
    
    % ---------------------- 2. 计算阈值（自动/手动） 
    % 处理2D矩阵（单通道）
    if isempty(c)
        % 提取矩阵的数值范围（用于Otsu归一化和结果显示）
        min_val = min(any_matrix(:));
        max_val = max(any_matrix(:));
        if min_val == max_val
           warning('矩阵所有元素值相同，直接返回全高值矩阵');
           binary_matrix=ones(size(any_matrix), ...
               class(any_matrix))*high_val;
           return;
        end
        
        % 计算阈值
        if isempty(manual_threshold)
            % 自动阈值：Otsu方法（适配任意数值分布，先归一化到0-1）
            normalized_matrix=(any_matrix-min_val)/(max_val-min_val);
            level = graythresh(normalized_matrix);  % 0-1范围的阈值
            threshold = min_val+level*(max_val-min_val);
            fprintf('自动阈值：%.0f(原数据：[%.0f, %.0f])\n', ...
                threshold, min_val, max_val);
        else
            % 手动阈值
            threshold = manual_threshold;
            fprintf('手动阈值：%.4f(原数据：[%.4f, %.4f])\n', ...
                threshold, min_val, max_val);
        end
        
        % 执行二值化：大于阈值→高值，小于等于→低值
        binary_matrix = zeros(size(any_matrix), class(any_matrix));
        binary_matrix(any_matrix > threshold) = high_val;
        binary_matrix(any_matrix <= threshold) = low_val;
    
    % 处理3D矩阵（多通道，按通道独立二值化）
    else
        binary_matrix = zeros(size(any_matrix), class(any_matrix));
        for ch = 1:c
            channel_matrix = any_matrix(:, :, ch);
            min_val = min(channel_matrix(:));
            max_val = max(channel_matrix(:));
            if min_val == max_val
            warning('第%d通道所有元素值相同，返回全高值', ch);
            binary_matrix(:,:,ch)=ones(m,n,class(any_matrix))*high_val;
            continue;
            end
            
            % 计算阈值
            if isempty(manual_threshold)
               normalized=(channel_matrix-min_val)/(max_val-min_val);
               level = graythresh(normalized);
               threshold = min_val + level * (max_val - min_val);
               fprintf('第%d通道Otsu阈值：%.4f（范围：[%.4f, %.4f]）\n',...
                   ch, threshold, min_val, max_val);
            else
               threshold = manual_threshold;
               fprintf('第%d通道手动阈值：%.4f（范围：[%.4f, %.4f]）\n',...
                   ch, threshold, min_val, max_val);
            end
            
            % 二值化
            binary_matrix(:, :, ch) = low_val;
            binary_matrix(channel_matrix > threshold, ch) = high_val;
        end
    end
    
    % ---------------------- 3. 结果可视化（可选，仅2D矩阵或3D单通道）
    if isempty(c) || c == 1
       figure('Color', 'white', 'Position', [100, 100, 800, 400]);
       
       % 显示原始矩阵（用jet colormap适配任意数值）
       subplot(1, 2, 1);
       % org_matrix = any_matrix;
       imagesc(any_matrix);
       colormap(jet);
       colorbar;
       title('原始矩阵（任意数值范围）', 'FontSize', ...
           12, 'FontWeight', 'bold');
       xlabel('列'); ylabel('行');
       
       % 显示二值化结果（黑白colormap）
       subplot(1, 2, 2);
       imagesc(binary_matrix);
       colormap(gray);
       colorbar;
       title(sprintf('二值化结果（阈值：%.4f，高值：%.1f，低值：%.1f）',...
           threshold, high_val, low_val), ...
           'FontSize', 12, ...
           'FontWeight', 'bold');
       xlabel('列'); ylabel('行');
    else
       fprintf('3D矩阵二值化完成，共%d个通道\n', c);
    end
    
end

%%

function plot_single_sided_amplitude_spectrum(signal, fs)
    % 仅计算并可视化信号的单边幅频谱（无其他冗余可视化）
    % 输入参数：
    %   signal - 输入信号（数值型，自动转为列向量）
    %   fs     - 采样频率（Hz，正数值）
    
    % ---------------------- 1. 输入参数校验 ----------------------
    % 校验信号合法性
    if ~isnumeric(signal) || isempty(signal)
        error('输入信号必须为非空的数值矩阵');
    end
    if ~iscolumn(signal)
        signal = signal(:);  % 强制转为列向量
        % warning('输入信号非列向量，已自动转置为列向量');
    end
    
    % 校验采样频率
    if fs <= 0 || ~isnumeric(fs)
        error('采样频率fs必须为正数值');
    end
    
    % ---------------------- 2. 计算单边幅频谱 ----------------------
    N = length(signal);  % 信号长度
    N_fft = 2^nextpow2(N);  % 补零优化FFT效率
    fft_result = fft(signal, N_fft);  % FFT计算
    
    % 频率轴（单边：0 ~ fs/2 Hz）
    freq_axis = (0:N_fft/2) * (fs / N_fft);  % 包含Nyquist频率点
    
    % 幅值谱归一化（保证与实际信号幅值一致）
    amp_spectrum = abs(fft_result(1:N_fft/2 + 1)) / N;
    amp_spectrum(2:end-1) = 2 * amp_spectrum(2:end-1);  % 修正非DC分量幅值
    
    % ---------------------- 3. 可视化单边幅频谱（仅这一个图） ------------
    figure('Color', 'white', 'Position', [100, 100, 800, 500]);
    
    % 绘制幅频谱曲线
    plot(freq_axis,amp_spectrum,'LineWidth',1.5,'Color',[0.2,0.6,0.8]);
    title('信号单边幅频谱', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('频率 (Hz)', 'FontSize', 12);
    ylabel('幅值', 'FontSize', 12);
    grid on;
    set(gca, 'GridAlpha', 0.3);  % 网格半透明，不遮挡曲线
    xlim([0, fs/2]);  % 固定显示到Nyquist频率
    ylim([0, max(amp_spectrum) * 1.1]);  % 预留10%顶部空间，避免峰值贴顶
    
    % 标注主要峰值频率（便于识别关键成分）
    [peak_vals, peak_idxs] = findpeaks(amp_spectrum, freq_axis, ...
        'MinPeakHeight', max(amp_spectrum)*0.1);  %仅标注幅值≥10%峰值的频
    if ~isempty(peak_idxs)
        text(peak_idxs, peak_vals, ...
        sprintf('%.1f Hz', peak_idxs), ...
         'FontSize',10,'VerticalAlignment', ...
         'bottom','Color',[0.8,0.2,0.2]);
    end
    
    % ---------------------- 4. 输出关键信息（命令窗口） ----------
    fprintf('采样频率：%.1f Hz\n', fs);
    fprintf('Nyquist频率：%.1f Hz\n', fs/2);
    fprintf('频谱分辨率：%.3f Hz\n', fs/N_fft);
    if ~isempty(peak_idxs)
        fprintf('主要频率成分（幅值≥10%%峰值）：');
        for i = 1:length(peak_idxs)
            fprintf('%.1f Hz (幅值：%.2f) ', peak_idxs(i), peak_vals(i));
        end
        fprintf('\n');
    end
end

function save_matrix_500step_cols(matrix, save_path, ...
    base_filename, save_mode)
% SAVE_MATRIX_500STEP_COLS 提取矩阵中500、1000、1500…列并保存为mat文件
% 输入参数：
%   matrix       - 待提取的原始矩阵（每行是样本，每列是维度/传感器）
%   save_path    - 保存路径（如'F:\k_wave\data\'，自动创建不存在的路径）
%   base_filename - 保存文件的基础名称（如'forward_data'）
%   save_mode    - 保存模式（可选，'merge'=合并所有目标列保存；
%   'single'=每个列单独保存，默认'merge'）

% 1. 输入参数校验与默认值设置
% 检查矩阵有效性
if ~ismatrix(matrix)
    error('输入的不是有效矩阵！');
end
n_rows = size(matrix, 1);   % 矩阵行数
n_cols = size(matrix, 2);   % 矩阵总列数

% 检查保存路径
if nargin < 2 || isempty(save_path)
    save_path = pwd;  % 默认当前工作目录
else
    if save_path(end) ~= '\' && save_path(end) ~= '/'
        save_path = [save_path, '\'];
    end
    if ~exist(save_path, 'dir')
        mkdir(save_path);
        fprintf('创建保存路径：%s\n', save_path);
    end
end

% 基础文件名默认值
if nargin < 3 || isempty(base_filename)
    base_filename = 'matrix_500step_cols';
end

% 保存模式默认值（merge=合并保存，single=单列单独保存）
if (nargin < 4 || isempty(save_mode) ...
    || ~ismember(save_mode, {'merge','single'}))
    save_mode = 'merge';
    fprintf('未指定保存模式，默认按「合并所有目标列」保存\n');
end

% 2. 计算需要提取的列索引（500、1000、1500…）
start_col = 501;    % 起始列
step_col = 500;     % 步长500
col_indices = start_col:step_col:n_cols;  % 生成目标列索引

% 边界校验：若无符合条件的列（如矩阵不足500列）
if isempty(col_indices)
    warning('矩阵总列数（%d列）不足500列，无需要保存的列！', n_cols);
    return;
end
fprintf('矩阵共%d列，将提取第%s列（共%d列）\n', ...
    n_cols, join(string(col_indices), '、'), length(col_indices));

% 3. 提取目标列并保存
% 提取所有目标列
target_cols = matrix(:, col_indices);

if strcmp(save_mode, 'merge')
    % 方式1：合并所有目标列保存为一个mat文件
    save_filename = [save_path, base_filename, '_500step_cols.mat'];
    save(save_filename, 'target_cols', 'col_indices');  % 同时保存列索引
    fprintf('合并保存完成：%s\n', save_filename);
    fprintf('保存的列索引：%s\n', join(string(col_indices), '、'));

elseif strcmp(save_mode, 'single')
    % 方式2：每个目标列单独保存为mat文件
    for idx = 1:length(col_indices)
        curr_col = col_indices(idx);  % 当前列号
        curr_data = matrix(:, curr_col);  % 提取当前列
        % 生成文件名（包含列号，如forward_data_col500.mat）
        save_filename = [save_path, base_filename, '_col', num2str(curr_col), '.mat'];
        save(save_filename, 'curr_data', 'curr_col');  % 保存列数据+列号
        fprintf('第%d列保存完成：%s\n', curr_col, save_filename);
    end
    fprintf('所有单列保存完成！共保存%d个文件\n', length(col_indices));
end

end





