clear;close all;clc
%%
% 設參數
n_ele = 2;              % 元素個數
nodes = n_ele+1;        % 節點個數
L = 2;                  % 總長
h = L/n_ele;            % 個元素長
r_xk = h/2;             % dx/dξ 值

% 建立空矩陣 KU=F
K = zeros(nodes,nodes); % K矩陣
F = zeros(nodes,1);     % F向量
U = zeros(nodes,1);     % U向量

% 定義 Shape function
N1 = @(k) (1-k)/2;
N2 = @(k) (1+k)/2;
Ni = @(k) [N1(k); N2(k)];
Ni_pram = [-1/2; 1/2]*1/r_xk;

for n = 1:n_ele
        
        x_h = [(n-1)*h;n*h];                     % x_hat 每elem.的節點xi及xi+1值
        
        % K矩陣運算
        A = @(k) (Ni_pram-Ni(k))*Ni_pram'*r_xk;  % (N'-N)N'*dx/dξ
        Ke = integral(A,-1,1,'ArrayValue',true); % 將上式積分
        K(n:n+1,n:n+1)=K(n:n+1,n:n+1)+Ke;        % 將各elem.之Ke存入K矩陣

        % F向量運算
        B = @(k) Ni(k).*Ni(k)'*x_h*r_xk;         % N*N'*x_h*dx/dξ
        fb = integral(B,-1,1,'ArrayValue',true); % 將上式積分
        F(n:n+1) =F(n:n+1)+fb ;                  % 將各elem.之fb存入F矩陣
end

% 加入邊界條件u(0)=0,u'(2)=2, 算 U = inv(K)*F
F(nodes) = F(nodes)+2;
u = inv(K(2:nodes,2:nodes))*F(2:nodes);
U(2:nodes) = U(2:nodes)+u;

disp(['Stiffness matrix = ',newline]);          disp(K)
disp(['Force vector = ',newline]);              disp(F)
disp(['Displacement value of each node = ',newline]); disp(U)

plot(0:h:L,U)
xlabel('Value of Xi')
ylabel('Value of Ui')
