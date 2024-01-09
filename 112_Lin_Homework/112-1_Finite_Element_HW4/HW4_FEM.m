clear;close all;clc

%% Define Element

nele = 24; % 元素個數 Can be modified
nodes = nele+1; % 節點個數
L = 8; % 總長度 Can be modified
C = 1; % 本構關係 Can be modified
h = L/nele; % 元素長
J = h/1; % Jacobin (原長/轉換後長)

% 建立空矩陣
K = zeros(nodes,nodes);
F = zeros(nodes,1);
U = zeros(nodes,1);
Ud = zeros(nodes,1);

%% Define Body Force 邊界條件

b = @(x) 0.1; % Body Force Can be modified 無直接槓掉
BC_0 = 1; % 起點邊界條件 Can be modified 
BC_L = -1; % 終點邊界條件 Can be modified 若無輸入555
BC_tl = 0; % 終點邊界條件 Can be modified 

%% Define Mesh

% Mesh(1,1)第一個元素的第一個節點
% Mesh(1,2)第一個元素的第二個節點

Mesh = zeros(nele:2);

for i = 1:nele

    Mesh(i,1) = i;
    Mesh(i,2) = i+1;

end

% NodeLocation 節點的座標

NodeLocation = zeros(nodes,1);

for i = 1:nele+1

    NodeLocation(i) = (i-1)*h;

end    
%% Define Shape Function

N1 = @(k) (1-k);
N2 = @(k) (k);
Ni = @(k) [N1(k),N2(k)];
B = [-1,1]*1/J; % 對Ni微分

n = 0.001;
csi = 0:n:1;


%% Compute K F Matrix 

for i = 1:nele

    Ke = C*transpose(B)*B*J; % C*B^T*B*J
    K(Mesh(i,1),Mesh(i,1)) = K(Mesh(i,1),Mesh(i,1))+Ke(1,1);
    K(Mesh(i,1),Mesh(i,2)) = K(Mesh(i,1),Mesh(i,2))+Ke(1,2);
    K(Mesh(i,2),Mesh(i,1)) = K(Mesh(i,2),Mesh(i,1))+Ke(2,1);
    K(Mesh(i,2),Mesh(i,2)) = K(Mesh(i,2),Mesh(i,2))+Ke(2,2);

end


for i = 1:nele
    b_argument = NodeLocation(Mesh(i,1))+csi*h; % 將元素長(x)離散化轉成(csi)
    b_trapz = b(b_argument); % 將外力b(x)轉換成csi(x)

    % 定義積分函數
    Intergrand1 = csi.*b_trapz*h;
    Intergrand2 = (1-csi).*b_trapz*h;

    % Trapezoidal Rule 積分計算
    fe1 = trapz(csi,Intergrand1);
    fe2 = trapz(csi,Intergrand2);

    % 組合F矩陣
    F(Mesh(i,1)) = F(Mesh(i,1))+fe1;
    F(Mesh(i,2)) = F(Mesh(i,2))+fe2;
end

%% Compute Dispacement and Stress 加上邊界條件 

if BC_L == 555 % 當終點無邊界條件時

    U(1) = BC_0;
    F(nodes) = F(nodes)+BC_tl; % 加上桿端外力
    u = inv(K(2:nodes,2:nodes))*F(2:nodes);
    U(2:nodes) = U(2:nodes)+u;
    
    for i = 1:nele
        ud = (U(i+1)-U(i))/h;
        Ud(i) = Ud(i)+ud;
    end
    
    sigma = C*Ud;
    
    figure
    plot(NodeLocation,U);grid on;
    xlabel('Length')
    ylabel('Displacement')
    legend('Displacement')
    
    figure
    stairs(NodeLocation,sigma)
    xlabel('Length')
    ylabel('Stress')
    legend('Stress')

else % 當終點有邊界條件時

    % 定義固定邊界條件
    fixed_nodes = [1, nodes]; % 第一個節點和最後一個節點
    fixed_values = [BC_0, BC_L]; % 對應的位移值
    
    % 條件數值積分法處理固定邊界條件
    penalty = 10000; % 懲罰因子（可以調整）
    for i = 1:length(fixed_nodes)
        index = fixed_nodes(i);
        K(index, :) = 0; % 將這些節點對應的整行清零
        K(index, index) = penalty; % 對角線加上懲罰因子
        F(index) = penalty * fixed_values(i); % 將對應的力修改為懲罰值乘上固定的位移值
    end
    
    % 計算變位
    U = K \ F;
    
    % 計算節點處的位移梯度

    for i = 1:nele
  
        ud = (U(i+1)-U(i))/h;
        Ud(i) = Ud(i)+ud;
    
    end
    
    sigma = C*Ud; % 求應力分布
    
    % 繪製位移圖,應力圖
    plot(NodeLocation, U);
    xlabel('Length');
    ylabel('Displacement');
    legend('Displacement');
    grid on;
    
    figure
    stairs(NodeLocation,sigma);
    xlabel('Length')
    ylabel('Stress')
    legend('Stress')

end