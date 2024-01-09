clear;close all;clc;

% 定義問題參數
L = 2;          % 區間的長度
nelements = 2;  % 元素的數量
nnodes = nelements + 1; % 節點的數量
h = L / nelements;  % 元素的大小

% 定義線性形狀函數
N1 = @(x) (1 - x) / h;
N2 = @(x) x / h;

% 建立全局勁度矩陣和荷載向量
K = zeros(nnodes, nnodes);
F = zeros(nnodes, 1);
U = zeros(nnodes, 1);
% 組裝全局勁度矩陣和荷載向量
for e = 1:nelements
    % 元素的端點
    x1 = (e - 1) * h;
    x2 = e * h;
    
    % 元素的勁度矩陣
    Ke = [1/h+1/2, -1/h-1/2; -1/h+1/2, 1/h-1/2];
    
    % 元素的荷載向量
    fb = h/2*[8/12*(e-1)+4/12*e;4/12*(e-1)+8/12*e];
    
    % 組裝到全局矩陣和向量中
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + Ke;
    F(e:e+1) = F(e:e+1) + fb;
end

K

% 加上桿端外力fp
fp = [0;2];
F(2:3)=F(2:3)+fp

%%

% 加入初始條件u(0)=0,u'(2)=2, 算 U = inv(K)*F

u = inv(K(2:nnodes,2:nnodes))*F(2:nnodes);
U(2:nnodes) = U(2:nnodes)+u;

disp(['Stiffness matrix = ',newline]);disp(K)
disp(['Force vector = ',newline]);disp(F)
disp(['Displacement of each node = ',newline]);disp(U)

