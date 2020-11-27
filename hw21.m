clear all;

% グローバル変数：obj_two_period.m, obj_projection.mと変数を共有
global w grid_w beta gamma rent

%% *** カリブレーション ***
beta  = 0.985; %0.985.^30;     % 割引因子
gamma = 2.0;           % 相対的危険回避度
rent  = 1.01; %1.025.^30-1.0; % 純利子率
wH = 1.0;
wL = 0.1;
p = 0.5;

% *** パラメータ ***
na1    =  10;   % 所得グリッドの数
a1_max = 1.0;   % 所得グリッドの最大値
a1_min = 0.1;   % 所得グリッドの最小値
na2    = 401;   % 貯蓄グリッドの数
a2_max = 1.0;   % 貯蓄グリッドの最大値
a2_min = 0.0; % 貯蓄グリッドの最小値
%a2_min = 0.0; % 貯蓄グリッドの最小値

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
maxit = 1000;    % 繰り返し計算の最大値
tol  = 1.0e-005; % 許容誤差(STEP 2)
dif1 = 1;        % 価値関数の繰り返し誤差
dif2 = 1.0;      % 政策関数の繰り返し誤差
count = 1;

%% グリッドポイントを計算
grid_a1 = linspace(a1_min, a1_max, na1)';
grid_a2 = linspace(a2_min, a2_max, na2)';

%% STEP 3: 効用関数の組み合わせ
util = -10000.0*ones(na1, na2, 2);

for i = 1:na1
    for j = 1:na2
        cH = (1+rent)*grid_a1(i) - grid_a2(j) + wH;
        cL = (1+rent)*grid_a1(i) - grid_a2(j) + wL;
        if cH > 0.0
            util(j,i,1) = CRRA(cH, gamma);
        end
        if cL > 0.0
            util(j,i,2) = CRRA(cL, gamma);
        end
    end
end

%% STEP 4: 価値関数を繰り返し計算
while it < maxit && dif1 > tol

    % ベルマン方程式: V(a,w;a')
    for i = 1:nk
        vkp(:,i) = util(:,i) + beta.*vfcn;
    end
    
    % 最適化: 各kについてV(k;k')を最大にするk'を探す
    [Tvfcn, ploc] = max(vkp);
    Tvfcn = Tvfcn';
    Tpfcn = kgrid(ploc);
    
    % 繰り返し計算誤差を確認
    dif1 = max(abs((Tvfcn-vfcn)./vfcn));

    % 価値関数・政策関数をアップデート
    vfcn = Tvfcn;
    pfcn = Tpfcn;
    fprintf('iteration index: %i, iteration diff of value: %d, iteration diff of policy: %d \n', it, dif1, dif2);
    
    it = it + 1;

end

% %% 効用を最大にする操作変数を探し出す：政策関数
% pol = zeros(na1,1);
% 
% % 各wについて生涯効用を最大にするaを探す
% for i = 1:na1
%     [maxv, maxl] = max(obj(:,i));
%     pol(i) = grid_a2(maxl);
% end
% 
% figure;
% plot(grid_a1, pol, 'bo-', 'color', 'blue', 'MarkerEdgeColor', 'b', 'MarkerSize', 12, 'linewidth', 3);
% xlabel('若年期の所得：a1', 'Fontsize', 16);
% ylabel('若年期の貯蓄：a2', 'Fontsize', 16);
% xlim([0, 1]);
% %ylim([0, 0.5]);
% set(gca, 'Fontsize', 16);
% grid on;