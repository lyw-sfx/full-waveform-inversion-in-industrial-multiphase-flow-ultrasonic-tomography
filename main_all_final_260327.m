%% 初始化 FWI - 声速和衰减系数联合重构（归一化约束版本）
% 基于伴随法的全波形反演进行工业油水两相流超声层析成像
%
% 本程序实现了基于伴随法的全波形反演(FWI)，用于重建油水两相流的
% 声速分布和衰减系数分布。采用交替更新策略，分别优化声速和衰减系数。
%
% 主要特点：
% 1. 联合重构：同时重建声速(c)和衰减系数(alpha)
% 2. TV正则化：使用全变分(Total Variation)正则化保持边缘
% 3. 联合正则化：引入声速-衰减结构相似性约束
% 4. 自适应线搜索：动态调整学习率保证收敛
% 5. 归一化约束：使用Min-Max归一化将场域映射到物理合理范围（非硬截断）
%
% 与上一版本的主要区别：
% - 使用归一化约束代替箱式限制（clip）
% - 保持场域的相对分布，整体缩放到目标范围
% - 避免硬截断导致的梯度不连续问题
%
% 作者：[您的名字]
% 日期：2026-03-27
% 版本：3.0

close all; 
clear; 
clc;

%% ==================== 正问题设置 ====================
% 正问题：根据已知的介质分布，计算传感器接收到的声压信号

%% 正问题 网格结构体设置
% 计算网格参数，确保满足Nyquist采样定理
% 每波长至少3个网格点，避免数值频散

c_min = 1200;           % 介质最小声速 [m/s]，用于确定最大网格尺寸
points_per_length = 3;  % 每波长网格点数（非均匀介质建议≥3）
f = 1e6;                % 超声中心频率 [Hz]
diameter = 0.15;        % 成像区域直径 [m] = 150mm

% 计算网格步长：dx = c_min / (points_per_length * f)
% 确保一个波长内至少有points_per_length个网格点
FW_dx = c_min / (points_per_length * f);
Fw_Nx = round(diameter / FW_dx);  % 理论网格数
FW_N = 410;                        % 实际使用网格数（选择2或3的倍数便于计算）

% 初始化k-wave网格对象
% kWaveGrid自动计算时间步长dt，满足CFL稳定性条件
FW_kgrid = kWaveGrid(FW_N, FW_dx, FW_N, FW_dx);

%% 正问题 介质属性结构体设置
% 定义水和油的基本物理参数
% 注意：k-wave使用幂律吸收模型 alpha = alpha_0 * f^y

T = 20;                        % 环境温度 [°C]
c_water = 1500;                % 水中声速 [m/s]（理论值约1482m/s@20°C）
rho_water = waterDensity(T);   % 水的密度 [kg/m³]

% 计算水的衰减系数
% waterAbsorption返回单位 [dB/(MHz^y·cm)]
alpha_water_ = waterAbsorption(f/1e6, T);
% 转换为k-wave使用的单位 [dB/(MHz^y·cm)]，y=1.43
alpha_water = alpha_water_ / (f/1e6)^1.43;

% 初始化参考介质（纯水）
medium_refw.sound_speed = c_water * ones(FW_N, FW_N);     % 声速场 [m/s]
medium_refw.density = rho_water * ones(FW_N, FW_N);       % 密度场 [kg/m³]
medium_refw.alpha_coeff = alpha_water * ones(FW_N, FW_N); % 吸收系数 [dB/(MHz^y·cm)]
medium_refw.alpha_power = 1.43;                           % 吸收指数 y（0<y<3, y≠1）

%% 正问题 仿真时间设置
% 设置仿真总时长，确保声波有足够时间传播到最远距离
% 对于直径150mm的圆形区域，最大传播距离为对角线长度

t_end = 1.5e-4;                       % 仿真结束时间 [s] = 150μs
FW_kgrid.t_array = 0:1/(20*f):t_end;  % 时间数组，采样率20*f = 20MHz

%% 正问题 声源设置
% 使用环形阵列，32个传感器均匀分布在圆周上
% 采用对向发射-接收模式（一个发射，对侧7个接收）

sensor_radius = 0.075;    % 传感器阵列半径 [m] = 75mm
num_sensor = 32;          % 传感器总数
num_rec = 7;              % 每个发射对应的接收传感器数（奇数，包含对侧）

% 生成传感器位置（笛卡尔坐标）
source_cart = makeCartCircle(sensor_radius, num_sensor, [0,0], 2*pi, false);

% 声源信号参数
periods = 5;              % 正弦脉冲周期数
target_amplitude = 10;    % 信号峰值幅度 [Pa]
fs = 1 / (FW_kgrid.dt);   % 实际采样频率 [Hz]

% 生成tone burst信号（加窗正弦波）
% 10*fs是为了在k-wave内部插值时保持精度
p = target_amplitude * toneBurst(10*fs, f, periods, 'Plot', false);

%% 正问题 传感器结构体
sensor_positions = makeCartCircle(sensor_radius, num_sensor, [0,0], 2*pi, false);

%% 正问题 参考场仿真（纯水模型）
% 计算纯水情况下的接收信号，用于后续差分成像

fw_input_args = {'DataCast', 'gpuArray-single',...  % 使用GPU加速（单精度）
                 'PlotSim', false,...                % 不显示仿真过程
                 'PMLSize', 20,...                   % PML厚度（网格点数）
                 'PMLInside', false,...              % PML在网格外部
                 'PlotPML', false};                  % 不显示PML

% 运行正演仿真
measured_data_refw = run_simulation_forward(FW_kgrid, medium_refw, ...
                     source_cart, p, sensor_positions, num_rec, num_sensor, ...
                     fw_input_args, 0);

% 显示参考场接收波形
data_show(measured_data_refw, num_rec, 'Sensor Position', 'Time Step', '全水场域');

%% 正问题 油滴模型设置
% 在纯水中添加油滴 inclusion，模拟油水两相流

x_vec = FW_kgrid.x_vec;              
y_vec = FW_kgrid.y_vec;              
[X, Y] = meshgrid(x_vec, y_vec);     % 生成网格坐标

% 油的物理参数
c_oil = 1200;              % 油中声速 [m/s]（比水低）
rho_oil = 800;             % 油的密度 [kg/m³]（比水低）
alpha_power_oil = 1.43;    % 油的吸收指数（与水相同）
alpha_oil_ = 0.0001;       % 油的衰减系数 [dB/(MHz^y·cm)]（比水低）
alpha_oil = alpha_oil_ / (f/1e6)^1.43;  % 单位转换

% 从参考介质复制，然后修改油滴区域
medium_oil = medium_refw;

% 定义油滴位置和大小
oil_center = [-0.02, -0.00];   % 油滴中心坐标 [m]
oil_radius = 0.0075;           % 油滴半径 [m] = 7.5mm

% 生成油滴掩码（判断网格点是否在油滴内部）
inclusion_dis = ((X - oil_center(1)).^2 + (Y - oil_center(2)).^2);
inclusion_mask = inclusion_dis <= oil_radius^2;

% 更新油滴区域的介质参数
medium_oil.sound_speed(inclusion_mask) = c_oil;
medium_oil.density(inclusion_mask) = rho_oil;
medium_oil.alpha_coeff(inclusion_mask) = alpha_oil;

% 定义油滴位置和大小
oil_center = [0.02, -0.00];   % 油滴中心坐标 [m]
oil_radius = 0.0075;          % 油滴半径 [m] = 7.5mm

% 生成油滴掩码（判断网格点是否在油滴内部）
inclusion_dis = ((X - oil_center(1)).^2 + (Y - oil_center(2)).^2);
inclusion_mask = inclusion_dis <= oil_radius^2;

% 更新油滴区域的介质参数
medium_oil.sound_speed(inclusion_mask) = c_oil;
medium_oil.density(inclusion_mask) = rho_oil;
medium_oil.alpha_coeff(inclusion_mask) = alpha_oil;

%% 正问题 油滴模型仿真
measured_data_oil = run_simulation_forward(FW_kgrid, medium_oil, ...
                        source_cart, p, sensor_positions, num_rec, num_sensor, ...
                        fw_input_args, 0);

[row_measured_data, col_measured_data] = size(measured_data_oil);

data_show(measured_data_oil, num_rec, 'Sensor Position', 'Time Step', '油水混合场域');

%% ==================== FWI 逆问题设置 ====================
% 逆问题：根据接收到的声压信号，反演介质参数分布

%% FWI 逆问题网格设置
INV_dx = FW_dx;           % 逆问题网格步长（与正问题相同）
INV_N = FW_N;             % 逆问题网格数
INV_kgrid = kWaveGrid(INV_N, INV_dx, INV_N, INV_dx);

%% FWI 初始场设置
% 初始猜测：纯水介质（与真实值不同，验证算法收敛性）

medium_rec = struct();
medium_rec.sound_speed = c_water * ones(INV_N, INV_N);      % 声速初值
medium_rec.density = rho_water * ones(INV_N, INV_N);        % 密度初值
medium_rec.alpha_coeff = alpha_water * ones(INV_N, INV_N);  % 衰减系数初值
medium_rec.alpha_power = 1.43;

%% FWI 传感器掩码设置
% 定义需要记录全场声压的区域（用于伴随法梯度计算）

sensor_allfield = zeros(INV_N, INV_N);
sa_index = round(INV_N/4) + 1;  % 重建区域边界索引
% 重建区域为内部区域（排除PML边界）
sensor_allfield(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index) = true;

%% FWI 时间设置
INV_kgrid.t_array = 0:1/(20*f):t_end;

%% FWI 仿真参数
iv_input_args = {'DataCast', 'gpuArray-single',...
                 'PlotSim', false,...
                 'PMLSize', 20,...
                 'PMLInside', false,...
                 'PlotPML', false,...
                 'PlotLayout', false,...
                 'DisplayMask', 'off'};

%% FWI 重建区域参数
rec_mesh_num = INV_N - 2*sa_index + 1;  % 重建区域网格数
num_k = rec_mesh_num;

% 初始化梯度向量
grad_speed = zeros(rec_mesh_num^2, 1);  % 声速梯度
grad_alpha = zeros(rec_mesh_num^2, 1);  % 衰减系数梯度

%% TV正则矩阵构建
% 使用独立的createTVMatrix函数构建TV正则化矩阵
% D: 差分算子矩阵，A = D'*D: TV梯度矩阵

[D, A] = createTVMatrix(rec_mesh_num, rec_mesh_num);

%% ==================== 迭代参数设置 ====================

num_iter = 0;              % 当前迭代计数
max_iter = 10;             % 最大迭代次数

% 学习率（步长）设置
lr_speed = 1e-3;           % 声速学习率
lr_alpha = 1e-3;           % 衰减系数学习率

% TV正则权重
beta_tv_speed = 1e-4;      % 声速TV正则权重
beta_tv_alpha = 1e-1;      % 衰减TV正则权重

% 联合正则化参数
beta_joint = 1e-2;         % 联合正则化权重（初始值，将自适应调整）
joint_reg_adaptive = true; % 是否使用自适应权重
joint_reg_window = 5;      % 结构相似性计算窗口大小

% 参数物理限制（归一化目标范围）
% 注意：使用Min-Max归一化将场域映射到以下范围
% 归一化公式：x_norm = min + (x - x_min) / (x_max - x_min) * (max - min)
c_min_limit = 1200;                    % 声速归一化目标最小值 [m/s]
c_max_limit = 1500;                    % 声速归一化目标最大值 [m/s]
alpha_min_limit = min(alpha_oil,alpha_water) ;   % 衰减系数归一化目标最小值
alpha_max_limit = max(alpha_oil,alpha_water) ;   % 衰减系数归一化目标最大值

% 损失函数历史记录
loss_history = zeros(max_iter, 1);
data_fidelity_history = zeros(max_iter, 1);
tv_reg_speed_history = zeros(max_iter, 1);
tv_reg_alpha_history = zeros(max_iter, 1);
joint_reg_history = zeros(max_iter, 1);

% 交替更新标志：1=更新声速，2=更新衰减系数
update_flag = 1;

fprintf('========================================\n');
fprintf('    FWI 声速-衰减联合重构\n');
fprintf('========================================\n');
fprintf('网格尺寸: %d x %d\n', FW_N, FW_N);
fprintf('重建区域: %d x %d\n', rec_mesh_num, rec_mesh_num);
fprintf('传感器数: %d（发射）x %d（接收）\n', num_sensor, num_rec);
fprintf('最大迭代: %d\n', max_iter);
fprintf('声速学习率: %.2e, 衰减学习率: %.2e\n', lr_speed, lr_alpha);
fprintf('TV权重(声速): %.2e, TV权重(衰减): %.2e\n', beta_tv_speed, beta_tv_alpha);
fprintf('联合正则权重: %.2e (自适应: %s)\n', beta_joint, iif(joint_reg_adaptive, '是', '否'));
fprintf('========================================\n\n');

%% ==================== 主迭代循环 ====================
% 使用交替更新策略：
% 奇数次迭代更新声速，偶数次迭代更新衰减系数
% 每次更新包含：前向仿真 → 残差计算 → 伴随仿真 → 梯度计算 → 线搜索更新

while num_iter < max_iter
    
    fprintf('\n========== 第 %d 次迭代 ==========\n', num_iter + 1);
    fprintf('当前更新参数: %s\n', iif(update_flag == 1, '声速 (c)', '衰减系数 (α)'));
    
    % 重置当前参数的梯度
    if update_flag == 1
        grad_speed = zeros(rec_mesh_num^2, 1);
    else
        grad_alpha = zeros(rec_mesh_num^2, 1);
    end
    
    % 初始化前向仿真数据存储
    forward_data_2_all = zeros(size(measured_data_oil));
    
    %% 内循环：遍历所有传感器对
    % 对于每个发射传感器，计算对侧接收传感器的响应
    
    for num_sensor_i = 1:num_sensor
        
        medium_rec_adj = medium_rec;
        
        % 计算接收传感器索引（对侧±3个传感器）
        index_rec_center = mod(num_sensor_i + num_sensor/2, num_sensor);
        if index_rec_center == 0
            index_rec_center = num_sensor;
        end
        index_rec_num = (num_rec-1)/2;    
        index_1 = index_rec_center - index_rec_num;
        index_2 = index_rec_center + index_rec_num;
        index_rec = mod(index_1:index_2, num_sensor);
        index_rec(index_rec == 0) = num_sensor;
        
        % 数据索引范围
        num_index_1 = (num_sensor_i-1)*num_rec + 1;
        num_index_2 = num_sensor_i*num_rec;
        
        %% 前向传播仿真（全场记录）
        % 用于后续梯度计算（需要全场声压）
        
        source_cart_temp = source_cart(:, num_sensor_i);
        forward_source.p_mask = cart2grid(INV_kgrid, source_cart_temp);
        forward_source.p = p;
        forward_sensor_1.mask = sensor_allfield;
        
        forward_data_1_gpu = kspaceFirstOrder2D(INV_kgrid, medium_rec_adj, ...
                             forward_source, forward_sensor_1, iv_input_args{:});
        
        % GPU数据传输
        data_cast_idx = find(strcmpi(iv_input_args, 'DataCast'));
        if ~isempty(data_cast_idx) && strcmp(iv_input_args{data_cast_idx+1}, 'gpuArray-single')
            forward_data_1 = gather(forward_data_1_gpu);
            clear forward_data_1_gpu;
        else
            forward_data_1 = forward_data_1_gpu;
        end
        
        % 低通滤波（保留主要频率成分）
        forward_data_1_ = forward_data_1';
        forward_data_1_ = filterSignal(forward_data_1_, fs, 'low', 1.1*f);
        forward_data_1 = forward_data_1_';
        
        %% 前向传播仿真（传感器记录）
        % 用于计算残差
        
        forward_sensor_2.mask = sensor_positions(:, index_rec);
        forward_data_2_gpu = kspaceFirstOrder2D(INV_kgrid, medium_rec_adj, ...
                             forward_source, forward_sensor_2, iv_input_args{:});
        
        if ~isempty(data_cast_idx) && strcmp(iv_input_args{data_cast_idx+1}, 'gpuArray-single')
            forward_data_2 = gather(forward_data_2_gpu);
            clear forward_data_2_gpu;
        else
            forward_data_2 = forward_data_2_gpu;
        end
        
        % 低通滤波
        forward_data_2_ = forward_data_2';
        forward_data_2_ = filterSignal(forward_data_2_, fs, 'low', 1.1*f);
        forward_data_2 = forward_data_2_';
        
        % 存储仿真数据
        forward_data_2_all(:, num_index_1:num_index_2) = forward_data_2_;
        
        %% 残差计算
        % 残差 = 实测数据 - 仿真数据
        
        measured_data = measured_data_oil(:, num_index_1:num_index_2);
        adjoint_compute = measured_data' - forward_data_2;
        
        % 残差滤波
        adjoint_compute_ = adjoint_compute';
        adjoint_compute_ = filterSignal(adjoint_compute_, fs, 'low', 1.1*f);
        adjoint_compute = adjoint_compute_';
        
        % 时间反转（伴随法要求）
        adjoint_compute_flip = fliplr(adjoint_compute);
        adjoint_compute_flip_ = adjoint_compute_flip';
        
        %% 伴随源设置
        % 将反转残差作为伴随问题的源
        
        adjoint_source = struct();
        adjoint_source_position = sensor_positions(:, index_rec);
        adjoint_source.p_mask = cart2grid(INV_kgrid, adjoint_source_position);
        adjoint_source.p = adjoint_compute_flip;
        
        adjoint_sensor = struct();
        adjoint_sensor.mask = sensor_allfield;
        
        %% 伴随仿真
        % 计算伴随场，用于梯度计算
        
        adjoint_data_gpu = kspaceFirstOrder2D(INV_kgrid, medium_rec_adj, ...
                                      adjoint_source, adjoint_sensor, iv_input_args{:});
        
        if ~isempty(data_cast_idx) && strcmp(iv_input_args{data_cast_idx+1}, 'gpuArray-single')
            adjoint_data = gather(adjoint_data_gpu);
            clear adjoint_data_gpu;
        else
            adjoint_data = adjoint_data_gpu;
        end
        
        % 伴随场滤波
        adjoint_data_ = adjoint_data';
        adjoint_data_ = filterSignal(adjoint_data_, fs, 'low', 1.1*f);
        adjoint_data = adjoint_data_';
        
        % 二次反转（k-wave时间步进的要求）
        adjoint_data = fliplr(adjoint_data);
        
        %% 梯度计算
        [row_forward_data_1, ~] = size(forward_data_1);
        dt = INV_kgrid.dt;
        
        % 声速梯度计算
        if update_flag == 1
            % 声速梯度公式：∂J/∂c = -2/c³ ∫(∂²p/∂t² · q) dt
            % 其中p是正向场，q是伴随场
            
            forward_data_1_ddt = secondOrderDiff(forward_data_1, dt);
            speed_field_matrix = medium_rec_adj.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
            speed_field_vector = speed_field_matrix(:);
            
            for i = 1:row_forward_data_1
                grad_1 = sum(forward_data_1_ddt(i,:) .* adjoint_data(i,:));
                grad_2 = speed_field_vector(i)^3;
                grad_speed(i) = grad_1 / grad_2 + grad_speed(i);
            end
        end
        
        % 衰减系数梯度计算
        if update_flag == 2
            % 衰减梯度公式：∂J/∂α = ∫(∂p/∂t · q) dt
            % 使用累积积分近似时间导数
            
            alpha_field_matrix = medium_rec_adj.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
            alpha_field_vector = alpha_field_matrix(:);
            
            for i = 1:row_forward_data_1
                % 时间积分近似
                grad_alpha_temp = sum(cumtrapz(forward_data_1(i,:)) .* adjoint_data(i,:));
                grad_alpha(i) = grad_alpha_temp + grad_alpha(i);
            end
        end
    end
    
    %% 损失函数计算
    % 总损失 = 数据保真度 + TV正则(声速) + TV正则(衰减) + 联合正则
    
    field_speed = medium_rec.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
    field_alpha = medium_rec.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
    
    % 数据保真度项（L2范数）
    data_fidelity = norm(forward_data_2_all - measured_data_oil, 'fro')^2;
    
    % TV正则项
    tv_reg_speed = norm(D * field_speed(:), 2);
    tv_reg_alpha = norm(D * field_alpha(:), 2);
    
    % 联合正则项：声速和衰减的结构相似性
    [joint_reg, ssim_val] = compute_joint_regularization(field_speed, field_alpha, joint_reg_window);
    
    % 自适应调整联合正则权重
    if joint_reg_adaptive && num_iter > 0
        % 根据SSIM值调整权重：SSIM越低，权重越大（鼓励相似性）
        beta_joint = 1e-2 * (1 - ssim_val) * 10;
        beta_joint = max(min(beta_joint, 1e-1), 1e-3);  % 限制范围
    end
    
    % 总损失
    current_loss = data_fidelity + ...
                   beta_tv_speed * tv_reg_speed + ...
                   beta_tv_alpha * tv_reg_alpha + ...
                   beta_joint * joint_reg;
    
    % 记录历史
    loss_history(num_iter + 1) = current_loss;
    data_fidelity_history(num_iter + 1) = data_fidelity;
    tv_reg_speed_history(num_iter + 1) = tv_reg_speed;
    tv_reg_alpha_history(num_iter + 1) = tv_reg_alpha;
    joint_reg_history(num_iter + 1) = joint_reg;
    
    % 输出损失信息
    fprintf('数据保真项:     %.4e\n', data_fidelity);
    fprintf('TV正则(声速):   %.4e (权重: %.2e)\n', tv_reg_speed, beta_tv_speed);
    fprintf('TV正则(衰减):   %.4e (权重: %.2e)\n', tv_reg_alpha, beta_tv_alpha);
    fprintf('联合正则:       %.4e (权重: %.2e, SSIM: %.3f)\n', joint_reg, beta_joint, ssim_val);
    fprintf('总损失:         %.4e\n', current_loss);
    
    %% 线搜索和参数更新
    if update_flag == 1
        % 更新声速（包含联合正则梯度）
        [medium_rec, lr_speed, success] = line_search_update_speed(...
            medium_rec, grad_speed, lr_speed, field_speed, ...
            beta_tv_speed, beta_joint, A, D, c_min_limit, c_max_limit, ...
            sa_index, INV_N, field_alpha, joint_reg_window, ...
            forward_data_2_all, measured_data_oil);
        
        if success
            fprintf('声速更新成功，新学习率: %.2e\n', lr_speed);
        else
            fprintf('声速更新失败，学习率过小\n');
        end
        
        % 切换到衰减更新
        update_flag = 2;
        
    else
        % 更新衰减系数（包含联合正则梯度）
        [medium_rec, lr_alpha, success] = line_search_update_alpha(...
            medium_rec, grad_alpha, lr_alpha, field_alpha, ...
            beta_tv_alpha, beta_joint, A, D, alpha_min_limit, alpha_max_limit, ...
            sa_index, INV_N, field_speed, joint_reg_window, ...
            forward_data_2_all, measured_data_oil);
        
        if success
            fprintf('衰减系数更新成功，新学习率: %.2e\n', lr_alpha);
        else
            fprintf('衰减系数更新失败，学习率过小\n');
        end
        
        % 切换回声速更新，并增加迭代计数
        update_flag = 1;
        num_iter = num_iter + 1;
    end
    
    %% 显示当前重建结果
    if mod(num_iter, 1) == 0
        figure(100);
        set(gcf, 'Name', sprintf('FWI重建过程 - 迭代 %d', num_iter));
        
        subplot(2, 2, 1);
        imagesc(field_speed);
        title(sprintf('声速重建 (迭代 %d)', num_iter));
        colorbar;
        colormap(jet);
        axis image;
        
        subplot(2, 2, 2);
        imagesc(field_alpha);
        title(sprintf('衰减系数重建 (迭代 %d)', num_iter));
        colorbar;
        colormap(jet);
        axis image;
        
        subplot(2, 2, 3);
        imagesc(medium_oil.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
        title('真实声速');
        colorbar;
        colormap(jet);
        axis image;
        
        subplot(2, 2, 4);
        imagesc(medium_oil.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
        title('真实衰减系数');
        colorbar;
        colormap(jet);
        axis image;
        
         drawnow;
    end
end

%% ==================== 结果显示 ====================

%% 损失函数收敛曲线
figure('Name', '损失函数收敛曲线', 'Position', [100, 100, 1200, 400]);

subplot(1, 4, 1);
plot(1:num_iter, loss_history(1:num_iter), 'b-o', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('总损失');
title('总损失函数');
grid on;

subplot(1, 4, 2);
plot(1:num_iter, data_fidelity_history(1:num_iter), 'r-s', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('数据保真度');
title('数据保真度项');
grid on;

subplot(1, 4, 3);
plot(1:num_iter, tv_reg_speed_history(1:num_iter), 'g-^', 'LineWidth', 2);
hold on;
plot(1:num_iter, tv_reg_alpha_history(1:num_iter), 'm-v', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('TV正则');
title('TV正则项');
legend('声速', '衰减');
grid on;

subplot(1, 4, 4);
plot(1:num_iter, joint_reg_history(1:num_iter), 'c-d', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('联合正则');
title('联合正则项');
grid on;

%% 最终重建结果
figure('Name', '最终重建结果', 'Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
imagesc(medium_rec.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('声速重建结果');
colorbar;
colormap(jet);
axis image;
clim([c_min_limit, c_max_limit]);

subplot(1, 2, 2);
imagesc(medium_rec.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('衰减系数重建结果');
colorbar;
colormap(jet);
axis image;

%% 与真实值对比
figure('Name', '声速对比', 'Position', [100, 100, 1200, 300]);

subplot(1, 3, 1);
imagesc(medium_oil.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('真实声速');
colorbar;
colormap(jet);
axis image;
clim([c_min_limit, c_max_limit]);

subplot(1, 3, 2);
imagesc(medium_rec.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('重建声速');
colorbar;
colormap(jet);
axis image;
clim([c_min_limit, c_max_limit]);

subplot(1, 3, 3);
speed_error = medium_rec.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index) - ...
              medium_oil.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
imagesc(speed_error);
title('声速误差');
colorbar;
colormap(jet);
axis image;

% 计算声速重建质量指标
speed_rmse = sqrt(mean(speed_error(:).^2));
speed_mae = mean(abs(speed_error(:)));
speed_maxe = max(abs(speed_error(:)));
fprintf('\n声速重建误差:\n');
fprintf('  RMSE: %.2f m/s\n', speed_rmse);
fprintf('  MAE:  %.2f m/s\n', speed_mae);
fprintf('  MaxE: %.2f m/s\n', speed_maxe);

figure('Name', '衰减系数对比', 'Position', [100, 100, 1200, 300]);

subplot(1, 3, 1);
imagesc(medium_oil.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('真实衰减系数');
colorbar;
colormap(jet);
axis image;

subplot(1, 3, 2);
imagesc(medium_rec.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index));
title('重建衰减系数');
colorbar;
colormap(jet);
axis image;

subplot(1, 3, 3);
alpha_error = medium_rec.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index) - ...
              medium_oil.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index);
imagesc(alpha_error);
title('衰减系数误差');
colorbar;
colormap(jet);
axis image;

% 计算衰减系数重建质量指标
alpha_rmse = sqrt(mean(alpha_error(:).^2));
alpha_mae = mean(abs(alpha_error(:)));
alpha_maxe = max(abs(alpha_error(:)));
fprintf('\n衰减系数重建误差:\n');
fprintf('  RMSE: %.4e\n', alpha_rmse);
fprintf('  MAE:  %.4e\n', alpha_mae);
fprintf('  MaxE: %.4e\n', alpha_maxe);

fprintf('\n========== FWI联合重构完成 ==========\n');

%% ==================== 辅助函数 ====================

%% 联合正则化计算函数
% 计算声速和衰减系数之间的结构相似性
% 使用Min-Max归一化到[0,1]后，计算局部结构相似性

function [joint_reg, ssim_val] = compute_joint_regularization(field_speed, field_alpha, window_size)
    % 输入:
    %   field_speed - 声速场矩阵
    %   field_alpha - 衰减系数场矩阵
    %   window_size - 局部窗口大小（奇数）
    %
    % 输出:
    %   joint_reg - 联合正则项值（越小表示越相似）
    %   ssim_val  - 结构相似性指数（0-1，越大越相似）
    
    % Min-Max归一化到[0,1]
    c_norm = (field_speed - min(field_speed(:))) / (max(field_speed(:)) - min(field_speed(:)) + eps);
    alpha_norm = (field_alpha - min(field_alpha(:))) / (max(field_alpha(:)) - min(field_alpha(:)) + eps);
    
    % 计算局部均值
    kernel = ones(window_size) / (window_size^2);
    mu_c = filter2(kernel, c_norm, 'same');
    mu_alpha = filter2(kernel, alpha_norm, 'same');
    
    % 计算局部方差和协方差
    sigma_c_sq = filter2(kernel, c_norm.^2, 'same') - mu_c.^2;
    sigma_alpha_sq = filter2(kernel, alpha_norm.^2, 'same') - mu_alpha.^2;
    sigma_c_alpha = filter2(kernel, c_norm .* alpha_norm, 'same') - mu_c .* mu_alpha;
    
    % SSIM参数
    C1 = 0.01^2;  % 亮度比较常数
    C2 = 0.03^2;  % 对比度比较常数
    
    % 计算SSIM图
    ssim_map = ((2 * mu_c .* mu_alpha + C1) .* (2 * sigma_c_alpha + C2)) ./ ...
               ((mu_c.^2 + mu_alpha.^2 + C1) .* (sigma_c_sq + sigma_alpha_sq + C2));
    
    % 平均SSIM
    ssim_val = mean(ssim_map(:));
    
    % 联合正则项：1 - SSIM（使优化目标为最小化）
    % 当SSIM=1（完全相同）时，joint_reg=0
    % 当SSIM=0（完全不同）时，joint_reg=1
    joint_reg = 1 - ssim_val;
end

%% 联合正则化梯度计算函数
% 计算联合正则项对声速或衰减的梯度

function grad_joint = compute_joint_gradient(field_primary, field_secondary, window_size, is_speed)
    % 输入:
    %   field_primary   - 当前更新的场（声速或衰减）
    %   field_secondary - 另一个场（衰减或声速）
    %   window_size     - 局部窗口大小
    %   is_speed        - true表示计算声速梯度，false表示计算衰减梯度
    %
    % 输出:
    %   grad_joint - 联合正则项的梯度
    
    % Min-Max归一化
    c_norm = (field_primary - min(field_primary(:))) / (max(field_primary(:)) - min(field_primary(:)) + eps);
    alpha_norm = (field_secondary - min(field_secondary(:))) / (max(field_secondary(:)) - min(field_secondary(:)) + eps);
    
    % 归一化导数（用于链式法则）
    range_primary = max(field_primary(:)) - min(field_primary(:)) + eps;
    
    % 简化计算：直接计算归一化场的差异梯度
    % 目标是使两个归一化场尽可能相似
    diff_norm = c_norm - alpha_norm;
    
    % 平滑差异（局部平均）
    kernel = ones(window_size) / (window_size^2);
    diff_smooth = filter2(kernel, diff_norm, 'same');
    
    % 梯度：差异越大，梯度越大（推动场向相似方向变化）
    grad_joint = 2 * diff_smooth / range_primary;
end

%% 更新声速（包含联合正则，使用归一化约束）
function [medium_rec, new_lr, success] = line_search_update_speed(...
    medium_rec, grad_speed, lr, field_org, beta_tv, beta_joint, A, D, ...
    c_min, c_max, sa_index, INV_N, field_alpha, joint_window, ...
    forward_data_2_all, measured_data_oil)
    % 使用Min-Max归一化将声速场约束到指定范围
    % 归一化公式：x_norm = c_min + (x - x_min) / (x_max - x_min) * (c_max - c_min)
    
    num_k = size(field_org, 1);
    grad_matrix = reshape(grad_speed, [num_k, num_k]);
    
    % 平滑梯度
    grad_matrix = mean_filter_2d(grad_matrix, 5);
    
    % 计算联合正则梯度
    grad_joint_matrix = compute_joint_gradient(field_org, field_alpha, joint_window, true);
    
    max_attempts = 5;
    success = false;
    new_lr = lr;
    
    for attempt = 1:max_attempts
        % 计算TV梯度
        tv_grad = reshape(A * field_org(:), [num_k, num_k]);
        
        % 总梯度 = 数据梯度 + TV正则梯度 + 联合正则梯度
        total_grad = grad_matrix + beta_tv * tv_grad + beta_joint * grad_joint_matrix;
        
        % 梯度下降更新
        field_update = field_org - new_lr * total_grad;
        
        % 使用Min-Max归一化将场域映射到[c_min, c_max]范围
        % 保持场域的相对分布，整体缩放到目标范围内
        field_update = normalize_to_range(field_update, c_min, c_max);
        
        % 检查更新量
        update_norm = norm(field_update - field_org, 'fro');
        relative_update = update_norm / (norm(field_org, 'fro') + eps);
        
        if relative_update < 0.5
            % 接受更新
            medium_rec.sound_speed(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index) = field_update;
            success = true;
            new_lr = min(new_lr * 1.1, 1e0);
            fprintf('  声速更新: 相对变化 %.4f, 学习率 %.2e\n', relative_update, new_lr);
            break;
        else
            new_lr = new_lr / 2;
            if new_lr < 1e-4
                fprintf('  声速更新: 学习率过小，跳过更新\n');
                break;
            end
        end
    end
end

%% 更新衰减系数（包含联合正则，使用归一化约束）
function [medium_rec, new_lr, success] = line_search_update_alpha(...
    medium_rec, grad_alpha, lr, field_org, beta_tv, beta_joint, A, D, ...
    alpha_min, alpha_max, sa_index, INV_N, field_speed, joint_window, ...
    forward_data_2_all, measured_data_oil)
    % 使用Min-Max归一化将衰减系数场约束到指定范围
    % 归一化公式：x_norm = alpha_min + (x - x_min) / (x_max - x_min) * (alpha_max - alpha_min)
    
    num_k = size(field_org, 1);
    grad_matrix = reshape(grad_alpha, [num_k, num_k]);
    
    % 平滑梯度
    grad_matrix = mean_filter_2d(grad_matrix, 5);
    
    % 计算联合正则梯度
    grad_joint_matrix = compute_joint_gradient(field_org, field_speed, joint_window, false);
    
    max_attempts = 5;
    success = false;
    new_lr = lr;
    
    for attempt = 1:max_attempts
        % 计算TV梯度
        tv_grad = reshape(A * field_org(:), [num_k, num_k]);
        
        % 总梯度 = 数据梯度 + TV正则梯度 + 联合正则梯度
        total_grad = grad_matrix + beta_tv * tv_grad + beta_joint * grad_joint_matrix;
        
        % 梯度下降更新
        field_update = field_org - new_lr * total_grad;
        
        % 使用Min-Max归一化将场域映射到[alpha_min, alpha_max]范围
        % 保持场域的相对分布，整体缩放到目标范围内
        field_update = normalize_to_range(field_update, alpha_min, alpha_max);
        
        % 检查更新量
        update_norm = norm(field_update - field_org, 'fro');
        relative_update = update_norm / (norm(field_org, 'fro') + eps);
        
        if relative_update < 0.5
            % 接受更新
            medium_rec.alpha_coeff(sa_index:INV_N-sa_index, sa_index:INV_N-sa_index) = field_update;
            success = true;
            new_lr = min(new_lr * 1.1, 1e-1);
            fprintf('  衰减更新: 相对变化 %.4f, 学习率 %.2e\n', relative_update, new_lr);
            break;
        else
            new_lr = new_lr / 2;
            if new_lr < 1e-6
                fprintf('  衰减更新: 学习率过小，跳过更新\n');
                break;
            end
        end
    end
end

%% k-wave 正问题仿真函数
function measured_data = run_simulation_forward(kgrid, medium, ...
                         source_position, source_p, ...
                         sensor_positions, num_rec, num_sensor, ...
                         input_args, z)
    % 串行激励、并行接收的FWI正演仿真
    %
    % 输入:
    %   kgrid           - k-wave网格对象
    %   medium          - 介质参数结构体
    %   source_position - 声源位置坐标
    %   source_p        - 声源信号
    %   sensor_positions- 传感器位置坐标
    %   num_rec         - 接收传感器数量
    %   num_sensor      - 总传感器数量
    %   input_args      - k-wave可选参数
    %   z               - 0=全角度，1=单角度（测试用）
    %
    % 输出:
    %   measured_data   - 所有传感器的接收信号 [时间样本 x 传感器数]

    for i = 1:num_sensor
        % 计算接收传感器索引（对侧±num_rec/2）
        index_rec_center = mod(i + num_sensor/2, num_sensor);
        if index_rec_center == 0
            index_rec_center = num_sensor;
        end
        index_rec_num = (num_rec-1)/2;
        
        index_1 = index_rec_center - index_rec_num;
        index_2 = index_rec_center + index_rec_num;
        index_rec = mod(index_1:index_2, num_sensor);
        index_rec(index_rec == 0) = num_sensor;
        
        % 构建声源
        source_cart_temp = source_position(:, i);
        source_temp.p_mask = cart2grid(kgrid, source_cart_temp);
        source_temp.p = source_p;
        
        % 构建传感器
        sensor_positions_temp = sensor_positions(:, index_rec);
        sensor_temp.mask = sensor_positions_temp;

        % 运行k-wave仿真
        sensor_data_gpu = kspaceFirstOrder2D(kgrid, medium, ...
                                         source_temp, sensor_temp, ...
                                         input_args{:});
        
        % 处理GPU数据
        data_cast_idx = find(strcmpi(input_args, 'DataCast'));
        if ~isempty(data_cast_idx) && strcmp(input_args{data_cast_idx+1}, 'gpuArray-single')
            sensor_data = gather(sensor_data_gpu);
            measured_data(:, (i-1)*num_rec+1:num_rec*i) = sensor_data';
            clear sensor_data_gpu;
        else
            measured_data(:, (i-1)*num_rec+1:num_rec*i) = sensor_data_gpu';
        end
        
        % 单角度模式（用于快速测试）
        if z == 1
            fprintf("仅计算一个角度\n");
            break;
        end
    end
end

%% 接收波形展示函数
function data_show(data, num_rec, y_name, x_name, pic_name)
    % 显示传感器接收到的波形
    %
    % 输入:
    %   data     - 接收数据 [时间样本 x 总通道数]
    %   num_rec  - 每个发射对应的接收通道数
    %   y_name   - Y轴标签
    %   x_name   - X轴标签
    %   pic_name - 图窗标题

    figure('Name', pic_name);
    imagesc(data(:, 1:num_rec)');
    ylabel(y_name);
    xlabel(x_name);
    colormap("jet");
    colorbar;
    title(pic_name);
end

%% 二阶差分函数
function forward_ddt = secondOrderDiff(forward_in, dt)
    % 计算时间序列的二阶差分（用于声速梯度计算）
    %
    % 输入:
    %   forward_in - 输入信号 [空间点数 x 时间样本]
    %   dt         - 时间步长
    %
    % 输出:
    %   forward_ddt - 二阶时间导数

    [M, N] = size(forward_in);
    dt_sq = dt^2;
    forward_ddt = zeros(size(forward_in), class(forward_in));
    
    if N >= 3
        % 内部点：中心差分
        forward_ddt(:, 2:N-1) = (forward_in(:, 3:end) - ...
            2*forward_in(:, 2:end-1) + ...
            forward_in(:, 1:end-2)) / dt_sq;
        
        % 边界点：向前/向后差分
        forward_ddt(:, 1) = (forward_in(:, 3) - ...
            2*forward_in(:, 2) + ...
            forward_in(:, 1)) / dt_sq;
        forward_ddt(:, end) = (forward_in(:, end) - ...
            2*forward_in(:, end-1) + forward_in(:, end-2)) / dt_sq;
    else
        % 序列太短，使用简单差分
        if M == 1
            ddt_temp = diff(diff(forward_in)) / dt_sq;
            forward_ddt(1:length(ddt_temp)) = ddt_temp;
        else
            for i = 1:M
                ddt_temp = diff(diff(forward_in(i, :))) / dt_sq;
                forward_ddt(i, 1:length(ddt_temp)) = ddt_temp;
            end
        end
    end
end

%% 信号滤波函数
function filtered_signal = filterSignal(signal_in, fs, filter_type, cutoff_freq, order)
    % 对信号进行Butterworth滤波
    %
    % 输入:
    %   signal_in    - 输入信号 [样本 x 通道]
    %   fs           - 采样频率 [Hz]
    %   filter_type  - 滤波类型: 'low', 'high', 'bandpass', 'bandstop'
    %   cutoff_freq  - 截止频率 [Hz]
    %   order        - 滤波器阶数（默认4）
    %
    % 输出:
    %   filtered_signal - 滤波后信号

    signal = double(signal_in);
    if nargin < 5
        order = 4;
    end
    if ~ismatrix(signal) || ~isnumeric(signal)
        error('输入信号必须是数值矩阵');
    end
    if fs <= 0
        error('采样频率必须为正数');
    end
    valid_types = {'low', 'high', 'bandpass', 'bandstop'};
    if ~ismember(filter_type, valid_types)
        error('滤波类型必须是 low, high, bandpass 或 bandstop');
    end

    % 设计Butterworth滤波器
    nyquist = fs / 2;
    Wn = cutoff_freq / nyquist;
    [b, a] = butter(order, Wn, filter_type);

    % 对每个通道进行零相位滤波
    [num_samples, num_channels] = size(signal);
    filtered_signal = zeros(size(signal), class(signal));

    for i = 1:num_channels
        col_signal = signal(:, i);
        filtered_col = filtfilt(b, a, col_signal);
        filtered_signal(:, i) = filtered_col;
    end
end

%% 图像滤波/平滑函数
function filtered_img = mean_filter_2d(input_img, filter_size, varargin)
    % 二维均值滤波
    %
    % 输入:
    %   input_img    - 输入图像
    %   filter_size  - 滤波器大小（标量n表示nxn，或[m,n]）
    %   varargin     - 可选参数 'shape', 'same'/'valid'/'full'
    %
    % 输出:
    %   filtered_img - 滤波后图像

    if ~ismatrix(input_img) && ~ndims(input_img) == 3
        error('输入图像必须是二维或三维矩阵');
    end
    if ndims(input_img) == 3 && size(input_img, 3) ~= 3
        error('RGB图像必须是 M×N×3 矩阵');
    end
    
    % 解析滤波器大小
    if isscalar(filter_size)
        n = filter_size;
        if n <= 0 || mod(n, 1) ~= 0
            error('滤波器大小必须是正整数');
        end
        filter_rows = n;
        filter_cols = n;
    elseif isvector(filter_size) && length(filter_size) == 2
        [filter_rows, filter_cols] = deal(filter_size(1), filter_size(2));
        if filter_rows <= 0 || filter_cols <= 0 || ~all(mod(filter_size, 1) == 0)
            error('滤波器大小必须是正整数向量 [m, n]');
        end
    else
        error('滤波器大小必须是正整数标量或 [m, n] 向量');
    end
    
    % 解析形状参数
    if nargin < 3
        shape = 'same';
    else
        if ~strcmp(varargin{1}, 'shape') || ~ismember(varargin{2}, {'same', 'valid', 'full'})
            error('可选参数格式为 shape, same/valid/full');
        end
        shape = varargin{2};
    end
    
    % 构建均值滤波核
    mean_kernel = ones(filter_rows, filter_cols) / (filter_rows * filter_cols);
    
    % 应用滤波
    if ismatrix(input_img)
        img_double = double(input_img);
        filtered_double = filter2(mean_kernel, img_double, shape);
        filtered_img = cast(filtered_double, class(input_img));
    else
        filtered_img = zeros(size(input_img), class(input_img));
        for c = 1:3
            img_double = double(input_img(:, :, c));
            filtered_double = filter2(mean_kernel, img_double, shape);
            filtered_img(:, :, c) = cast(filtered_double, class(input_img));
        end
    end
end

%% Min-Max归一化函数
function normalized_matrix = normalize_to_range(original_matrix, lower_bound, upper_bound)
    % 使用Min-Max线性归一化将矩阵值映射到指定范围
    % 归一化公式：x_norm = lower + (x - x_min) / (x_max - x_min) * (upper - lower)
    % 特点：保持场域的相对分布，整体缩放到目标范围内
    %
    % 输入:
    %   original_matrix - 输入矩阵
    %   lower_bound     - 目标范围下限值
    %   upper_bound     - 目标范围上限值
    %
    % 输出:
    %   normalized_matrix - 归一化后的矩阵

    if ~isnumeric(original_matrix)
        error('输入矩阵必须是数值型矩阵');
    end
    if ~isnumeric(lower_bound) || ~isnumeric(upper_bound)
        error('上下限阈值必须是数值');
    end
    if lower_bound > upper_bound
        error('下限阈值不能大于上限阈值');
    end
    
    % 计算当前场域的最小值和最大值
    x_min = min(original_matrix(:));
    x_max = max(original_matrix(:));
    
    % 如果场域已经是常数，直接返回范围中值
    if abs(x_max - x_min) < eps
        normalized_matrix = ones(size(original_matrix)) * (lower_bound + upper_bound) / 2;
        normalized_matrix = cast(normalized_matrix, class(original_matrix));
        return;
    end
    
    % Min-Max归一化：线性映射到[lower_bound, upper_bound]
    % 公式：x_norm = lower + (x - x_min) / (x_max - x_min) * (upper - lower)
    normalized_matrix = lower_bound + ...
        (original_matrix - x_min) / (x_max - x_min) * (upper_bound - lower_bound);
    
    % 保持原始数据类型
    normalized_matrix = cast(normalized_matrix, class(original_matrix));
end

%% 箱式限制函数（保留用于兼容性，但不再在更新函数中使用）
function clipped_matrix = clip_matrix_values(original_matrix, lower_bound, upper_bound)
    % 将矩阵值截断到指定范围内（硬约束）
    % 注意：此函数会改变场域的相对分布，仅用于特殊场景
    %
    % 输入:
    %   original_matrix - 输入矩阵
    %   lower_bound     - 下限值
    %   upper_bound     - 上限值
    %
    % 输出:
    %   clipped_matrix  - 裁剪后的矩阵

    if ~isnumeric(original_matrix)
        error('输入矩阵必须是数值型矩阵');
    end
    if ~isnumeric(lower_bound) || ~isnumeric(upper_bound)
        error('上下限阈值必须是数值');
    end
    if lower_bound > upper_bound
        error('下限阈值不能大于上限阈值');
    end
    
    clipped_matrix = original_matrix;
    clipped_matrix(clipped_matrix < lower_bound) = lower_bound;
    clipped_matrix(clipped_matrix > upper_bound) = upper_bound;
    clipped_matrix = cast(clipped_matrix, class(original_matrix));
end

%% 条件表达式辅助函数
function result = iif(condition, true_val, false_val)
    % 内联if函数（类似C语言的?:运算符）
    %
    % 输入:
    %   condition - 条件表达式
    %   true_val  - 条件为真时的返回值
    %   false_val - 条件为假时的返回值
    %
    % 输出:
    %   result    - 根据条件返回的值

    if condition
        result = true_val;
    else
        result = false_val;
    end
end
