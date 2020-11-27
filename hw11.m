clear all;

% グローバル変数：obj_two_period.m, obj_projection.mと変数を共有
global w grid_w beta gamma rent

%% *** カリブレーション ***
beta  = 0.985.^30;     % 割引因子
gamma = 2.0;           % 相対的危険回避度
rent  = 1.025.^30-1.0; % 純利子率
wH = 1.0;
wL = 0.1;
p = 0.5;

% *** パラメータ ***
na1    =  10;   % 所得グリッドの数
a1_max = 1.0;   % 所得グリッドの最大値
a1_min = 0.1;   % 所得グリッドの最小値
na2    = 401;   % 貯蓄グリッドの数
a2_max = 1.0;   % 貯蓄グリッドの最大値
a2_min = -0.1; % 貯蓄グリッドの最小値
%a2_min = 0.0; % 貯蓄グリッドの最小値

%% グリッドポイントを計算
grid_a1 = linspace(a1_min, a1_max, na1)';
grid_a2 = linspace(a2_min, a2_max, na2)';

%% あらゆる(w,a)の組み合わせについて生涯効用を計算
obj = zeros(na1, na2);

for i = 1:na1
    for j = 1:na2
        c1 = grid_a1(i) - grid_a2(j);
        if c1 > 0.0
            c2H = (1.0+rent)*grid_a2(j)+wH;
            c2L = (1.0+rent)*grid_a2(j)+wL;
%             if c2L > 0.0
                obj(j, i) = CRRA(c1, gamma) ...
                    + beta*p*CRRA(c2H, gamma) ...
                    + beta*(1-p)*CRRA(c2L, gamma);
%             else
%             % 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
%             obj(j, i) = -10000.0;
%             end                
        else
            % 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
            obj(j, i) = -10000.0;
        end
    end
end

%% 効用を最大にする操作変数を探し出す：政策関数
pol = zeros(na1,1);

% 各wについて生涯効用を最大にするaを探す
for i = 1:na1
    [maxv, maxl] = max(obj(:,i));
    pol(i) = grid_a2(maxl);
end

figure;
plot(grid_a1, pol, 'bo-', 'color', 'blue', 'MarkerEdgeColor', 'b', 'MarkerSize', 12, 'linewidth', 3);
xlabel('若年期の所得：a1', 'Fontsize', 16);
ylabel('若年期の貯蓄：a2', 'Fontsize', 16);
xlim([0, 1]);
%ylim([0, 0.5]);
set(gca, 'Fontsize', 16);
grid on;