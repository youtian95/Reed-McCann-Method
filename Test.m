%% 定义
type = 2;
pga = 0.1:0.1:2;
N_LHS = 20;
Seed = 1;
N_MC = 10^6;

if type == 1
    % 三个相同相邻结构同时失效，example 1
    cutset = [0,0,0];
    MedianCapacity = [0.7,0.7,0.7];
    EpistemicLogSTD = ...
        [0.3,0.15,0.15; ...
        0.15,0.3,0.15; ...
        0.15,0.15,0.3];
    AleatoryLogSTD = ...
        [0.25,0.15,0.15; ...
        0.15,0.25,0.15; ...
        0.15,0.15,0.25];
elseif type == 2
    % 三个相同结构，两个相邻，一个在低楼层，example 2
    cutset = [0,0,0];
    MedianCapacity = [0.7,0.7,1];
    EpistemicLogSTD = ...
        [0.3,0.15,0.12; ...
        0.15,0.3,0.12; ...
        0.12,0.12,0.4];
    AleatoryLogSTD = ...
        [0.25,0.15,0.1; ...
        0.15,0.25,0.1; ...
        0.1,0.1,0.3];
end

obj = ReedMcCann(MedianCapacity,EpistemicLogSTD,AleatoryLogSTD);

%% 计算
% Reed-McCann方法
P = obj.FailureCutset(pga,cutset,N_LHS,Seed); %(i_pga, i_LHS)
P = mean(P,2)';
% MC模拟概率
% [P_MC,~] = obj.MonteCarlo(pga,cutset,N_MC,Seed);

%% 绘图
plot(pga,P_MC);
