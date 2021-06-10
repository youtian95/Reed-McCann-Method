classdef ReedMcCann < handle
    % Reed-McCann方法计算一个割集失效的概率
    %
    % 2021年6月2日19:43:42
    % 
    % 参考文献
    % [1] Budnitz RJ. Correlation of seismic performance in similar SSCs 
    %   (structures, systems, and components). United States Nuclear Regulatory 
    %   Commission, Office of Nuclear Regulatory …; 2017.
    
    properties (SetAccess = private)
        N % 结构数量
        theta_m % 能力中值，向量
        beta_U % 认知不确定性，(i,j)，矩阵
        beta_R % 随机不确定性，(i,j)，矩阵
        theta_U % 考虑认知不确定性随机生成的能力值，(k,i)，k表示第k次LHS抽样
    end
    
    properties (Access = private)
        x_range = 0.1:0.1:10; %数值积分自变量范围
    end
    
    methods (Access = public)
        function obj = ReedMcCann(MedianCapacity,EpistemicLogSTD,AleatoryLogSTD)
            % 输入：  
            % MedianCapacity - 一系列结构的抗震能力中值，向量，theta_m(i)
            % EpistemicLogSTD - 认知不确定性的对数标准差，beta_U(i,j)，矩阵。
            %       i≠j时，每个元素为i和j共有的随机变量的对数标准差
            %       i=j时，元素为i结构的认知不确定性的完整的对数标准差
            % AleatoryLogSTD - 随机不确定性的对数标准差，beta_R(i,j)，矩阵。
            %       i≠j时，每个元素为i和j共有的随机变量的对数标准差
            %       i=j时，元素为i结构的随机不确定性的完整的对数标准差
            
            obj.N = size(MedianCapacity,2);
            
            assert(size(MedianCapacity,1) == 1) ;
            assert(size(EpistemicLogSTD,1) == obj.N);
            assert(size(EpistemicLogSTD,2) == obj.N);
            assert(size(AleatoryLogSTD,1) == obj.N);
            assert(size(AleatoryLogSTD,2) == obj.N);
            assert(issymmetric(EpistemicLogSTD));
            assert(issymmetric(AleatoryLogSTD));
            
            
            obj.theta_m = MedianCapacity;
            
            % obj.beta_U, 认知不确定性，(i,j)，矩阵
            obj.beta_U = EpistemicLogSTD;
            obj.beta_U(logical(eye(obj.N))) = 0;
            % 对角元素列向量
            temp = vecnorm(EpistemicLogSTD(logical(eye(obj.N))),2,2).^2 ...
                - vecnorm(EpistemicLogSTD.*(ones(obj.N)-eye(obj.N)),2,2).^2;
            assert(all(temp>=0));
            temp = sqrt(temp);
            obj.beta_U = obj.beta_U + diag(temp);
            
            % obj.beta_R, 随机不确定性，(i,j)，矩阵
            obj.beta_R = AleatoryLogSTD;
            obj.beta_R(logical(eye(obj.N))) = 0;
            % 对角元素列向量
            temp = vecnorm(AleatoryLogSTD(logical(eye(obj.N))),2,2).^2 ...
                - vecnorm(AleatoryLogSTD.*(ones(obj.N)-eye(obj.N)),2,2).^2;
            assert(all(temp>=0));
            temp = sqrt(temp);
            obj.beta_R = obj.beta_R + diag(temp);
        end
        
        function P = FailureCutset(obj,pga,cutset,N_LHS,Seed)
            % 计算一个割集在某PGA下的失效的概率
            %
            % 输入：
            % pga - PGA向量
            % cutset - 逻辑向量，表示各个结构是否失效，0为失效，1为不失效
            % N_LHS - LHS抽样次数
            % seed - 随机数种子
            %
            % 输出：
            % P - 割集失效概率, (i_pga, i_LHS)
            
            rng(Seed);
            
            % Stage 1 - step 1
            p = zeros(N_LHS,obj.N);
            for j=1:obj.N
                p_int = randperm(N_LHS)';
                p_int = p_int-1;
                p_decimal = rand(N_LHS,1);
                p(:,j) = p_int + p_decimal;
            end
            p = p./N_LHS;
            theta_U_reduced = zeros(N_LHS,obj.N); %每行为一次抽样
            for j=1:obj.N
                % 每列一个结构
                theta_U_reduced(:,j) = icdf('Lognormal',p(:,j), ...
                    log(obj.theta_m(j)),obj.beta_U(j,j));
            end
            
            % Stage 1 - step 2
            p = zeros(obj.N,obj.N,N_LHS);
            for i=1:obj.N
                for j=(i+1):obj.N
                    p_int = randperm(N_LHS)';
                    p_int = p_int-1;
                    p_decimal = rand(N_LHS,1);
                    p(i,j,:) = p_int + p_decimal;
                    p(j,i,:) = p(i,j,:);
                end
            end
            p = p./N_LHS;
            CF = zeros(obj.N,obj.N,N_LHS);
            for i=1:obj.N
                for j=(i+1):obj.N
                    CF(i,j,:) = icdf('Lognormal',p(i,j,:), ...
                        0,obj.beta_U(i,j));
                    CF(j,i,:) = CF(i,j,:);
                end
            end
            
            % Stage 1 - step 3
            obj.theta_U = theta_U_reduced; %每行为一次抽样 (k,i)
            for k=1:N_LHS
                for i=1:obj.N
                    for j=1:obj.N
                        if i~=j
                            obj.theta_U(k,i) = obj.theta_U(k,i) ...
                                * CF(i,j,k);
                        end
                    end
                end
            end
            
            % Stage 2
            f = waitbar(0,'Please wait...');
            P = zeros(numel(pga),N_LHS); % (i_pga, i_LHS)
            for i_LHS=1:N_LHS
                waitbar(i_LHS./N_LHS,f,'数值积分...');
                for i_pga=1:numel(pga)
                    PGA = pga(i_pga);
                    if PGA<=0
                        P(i_pga) = 0;
                        continue;
                    end
                    %自变量, 数量为 (N^2-N)/2, 顺序为从上到下，从左到右
                    X_N = {}; 
                    expression = '[';
                    for i=1:((obj.N^2-obj.N)/2-1)
                        expression = [expression,'X_N{'];
                        expression = [expression,num2str(i)];
                        expression = [expression,'},'];
                    end
                    expression = [expression,'X_N{(obj.N^2-obj.N)/2}'];
                    expression = [expression,']=ndgrid(obj.x_range);'];
                    eval(expression); % [X1,X2,...,Xn] = ndgrid(x1,x2,...,xn)
                    % 积分函数 F
                    % 第1部分
                    F = ones(size(X_N{1}));
                    for i=1:obj.N
                        % 中间变量
                        SumLgX = 0; 
                        for j=1:obj.N
                            if i~=j
                                SumLgX = SumLgX + ...
                                    log(X_N{ReedMcCann.ij2vec(i,j,obj.N)});
                            end
                        end
                        temp = log(PGA) - log(obj.theta_U(i_LHS,i)) ...
                            -SumLgX;
                        temp = temp ./ obj.beta_R(i,i);
                        if ~cutset(i) %失效
                            F = F .* normcdf(temp);
                        else
                            F = F .* (1-normcdf(temp));
                        end
                    end
                    % 第2部分
                    for i=1:obj.N
                        for j=(i+1):obj.N
                            x = X_N{obj.ij2vec(i,j,obj.N)};
                            F = F ./ obj.beta_R(i,j) .* normpdf(log(x) ...
                                ./ obj.beta_R(i,j)) ./ x;
                        end
                    end
                    % 数值积分
                    for i=1:((obj.N^2-obj.N)/2)
                        F = trapz(obj.x_range,F,i);
                    end
                    P(i_pga, i_LHS) = F;
                end
            end
            close(f);
        end
        
        function [P,Sim_MC] = MonteCarlo(obj,pga,cutset,N_MC,Seed)
            % 采用蒙特卡罗模拟直接计算一个割集在某PGA下的失效的概率
            %
            % 输入：
            % pga - PGA向量
            % cutset - 逻辑向量，表示各个结构是否失效，0为失效，1为不失效
            % N_MC - MC抽样次数
            % seed - 随机数种子
            %
            % 输出：
            % P - 割集失效概率, (i_pga)
            % Sim_MC - MC模拟的结果, (i_structure, i_pga, i_MC), 
            %       每一个元素为各个结构MC模拟失效的结果, 0为失效，1为未失效
            
            rng(Seed);
            
            % 能力 Cap(i_structure, i_pga, i_MC)
            Cap = zeros(obj.N,numel(pga),N_MC);
            
            Sim_MC = zeros(obj.N,numel(pga),N_MC);
            
            for i_pga = 1:numel(pga)
                PGA = pga(i_pga);
                % 认知不确定性
                normU = zeros(obj.N,obj.N,N_MC);
                % 随机不确定性
                normR = zeros(obj.N,obj.N,N_MC);
                for i=1:obj.N
                    for j=i:obj.N
                        normU(i,j,:) = randn(1,1,N_MC);
                        normU(j,i,:) = normU(i,j,:);
                        normR(i,j,:) = randn(1,1,N_MC);
                        normR(j,i,:) = normR(i,j,:);
                    end
                end
                normU = normU .* obj.beta_U;
                normR = normR .* obj.beta_R;
                % 能力
                Cap(:,i_pga,:) = log(obj.theta_m') + ...
                    sum(normU,2) + sum(normR,2);
                Sim_MC(:,i_pga,:) = (exp(Cap(:,i_pga,:)) > PGA);
            end
            
            % 概率
            P = ~xor(Sim_MC,cutset');
            P = all(P,1);
            P = sum(P,3)./N_MC;
        end
    end
    
    methods (Static)
        function [i,j] = vec2ij(i_vec,N)
            % 将对称矩阵的线性排列（从上到下，从左到右，除了对角线的元素）索引
            %       对应到(i,j)索引
            %
            % 输入：
            % i_vec - 线性排列索引
            % N - 矩阵行列的数量
            %
            % 输出：
            % i,j - (i,j)索引
            
            assert(i_vec <= (N^2-N)/2);
            
            N_temp = i_vec;
            for row = 1:N
                N_temp = N_temp - (N - row);
                if N_temp <= 0
                    i = row;
                    j = N_temp + (N - row) + 1;
                    break;
                end
            end
        end
        function i_vec = ij2vec(i,j,N)
            % 将对称矩阵的(i,j)索引对应到线性排列（从上到下，从左到右，
            %   除了对角线的元素）索引
            %
            % 输入：
            % i,j - (i,j)索引
            % N - 矩阵行列的数量
            %
            % 输出：
            % i_vec - 线性排列索引
            
            i_vec = 0;
            if j<i % 左下角换到右上对称位置
                temp = j;
                j = i;
                i = temp;
            end
            for k=1:(i-1)
                i_vec = i_vec + (N - k);
            end
            i_vec = i_vec + j - i;
        end
    end
end

