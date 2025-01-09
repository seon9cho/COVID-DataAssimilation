%{

SIR model

%}
clc; clear all; close all; beep off; 

% -------------------------------------------------------------------------
% actual solution
% -------------------------------------------------------------------------

plot_it = 'on';

p.k = 1;
p.q = 0.3;
p.I_0 = 0.1;
p.R_0 = 0;

% SIR_ODE 
% S = y(1), I = y(2), R = y(3)
SIR_ODE = @(x,y)[-p.k*y(1)*y(2);p.k*y(1)*y(2)-p.q*y(2);p.q*y(2)];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol_actual = ode15s(SIR_ODE,[0,20],[1-p.I_0-p.R_0,p.I_0,p.R_0],options);

if strcmp(plot_it,'on')
    figure; hold on;
    plot(sol_actual.x,sol_actual.y,'-b','LineWidth',2);
end

% -------------------------------------------------------------------------
% get data
% -------------------------------------------------------------------------

t = linspace(0,5,300);

temp = deval(sol_actual,t);

I_data = temp(2,:)+normrnd(zeros(1,length(t)),0.05);

R_data = temp(3,:)+normrnd(zeros(1,length(t)),0.05);

if strcmp(plot_it,'on')
    plot(t,I_data,'.c','MarkerSize',8)
    plot(t,R_data,'.b','MarkerSize',8)
end


% -------------------------------------------------------------------------
% interpolation
% -------------------------------------------------------------------------

poly_I = polyfit(t,I_data,3);
poly_R = polyfit(t,R_data,3);
Idfun = @(t) poly_I(1)*t.^3+poly_I(2)*t.^2+poly_I(3)*t+poly_I(4);
Rdfun = @(t) poly_R(1)*t.^3+poly_R(2)*t.^2+poly_R(3)*t+poly_R(4);

% poly_I = polyfit(t,I_data,2);
% poly_R = polyfit(t,R_data,2);
% Idfun = @(t) poly_I(1)*t.^2+poly_I(2)*t+poly_I(3);
% Rdfun = @(t) poly_R(1)*t.^2+poly_R(2)*t+poly_R(3);

% poly_I = polyfit(t,I_data,1);
% poly_R = polyfit(t,R_data,1);
% Idfun = @(t) poly_I(1)*t+poly_I(2);
% Rdfun = @(t) poly_R(1)*t+poly_R(2);

if strcmp(plot_it,'on')
    plot(t,Idfun(t),'-m','LineWidth',2);
    plot(t,Rdfun(t),'-m','LineWidth',2);
end

% -------------------------------------------------------------------------
% gradient descent
% -------------------------------------------------------------------------

% initial guess
p.k = 1.1;
p.q = 0.25;
p.I_0 = 0.03;
p.R_0 = 0.01;

% Require S_0+I_0+R_0 = 1
p.S_0 = 1-p.I_0-p.R_0;
if p.S_0 < 0
    error('S_0 < 0')
end

for cnt = 1:30

    % -------------------------------------------------------------------------
    % Solve SIR for guess
    % -------------------------------------------------------------------------
    
    [sol_SIR,sol_adjoint,Sfun,Ifun,Rfun] = solve_ODEs(p,Idfun,Rdfun,'off',20);
    
    % -------------------------------------------------------------------------
    % gradients
    % -------------------------------------------------------------------------
    
    P0 = deval(sol_adjoint,0);
    
    PS = @(x)[1,0,0]*deval(sol_adjoint,x);
    PI = @(x)[0,1,0]*deval(sol_adjoint,x);
    PR = @(x)[0,0,1]*deval(sol_adjoint,x);
    
    fun1 = @(x) (PS(x)-PI(x)).*Sfun(x).*Ifun(x);
    J_k = integral(fun1,0,5);
    
    fun2 = @(x) (PI(x)-PR(x)).*Ifun(x);
    J_q = integral(fun2,0,5);
    
    J_S0 = -P0(1);
    J_I0 = -P0(2);
    J_R0 = -P0(3);

    % find the direction where the derivative is greatest.
    vec = randn(1,4);
    best_vec = vec/norm(vec);
    best_der = abs(best_vec*[J_k;J_q;J_I0-J_S0;J_R0-J_S0]);
    for q = 1:1000
        vec = randn(1,4);
        vec = vec/norm(vec);
        der = abs(vec*[J_k;J_q;J_I0-J_S0;J_R0-J_S0]);
        if der > best_der
            best_vec = vec;
            best_der = der;
        end
    end

    alpha = [0.01*(rand(10,1)-0.5);0]; % values for line search
    
    J = zeros(size(alpha));
    for j = 1:length(alpha)
    
        new_params = [p.k,p.q,p.I_0,p.R_0]+alpha(j)*best_vec;
    
        p1.k = new_params(1);
        p1.q = new_params(2);
        p1.I_0 = new_params(3);
        p1.R_0 = new_params(4);
    
        p1.S_0 = 1-p1.I_0-p1.R_0;
  
        [new_SIR,new_adjoint,new_Sfun,new_Ifun,new_Rfun] = solve_ODEs(p1,Idfun,Rdfun,'off',5);
    
        f3 = @(x) (new_Ifun(x)-Idfun(x)).^2+(new_Rfun(x)-Rdfun(x)).^2;
    
        J(j) = integral(f3,0,5);
    
    end
    
    ind = find(J==min(J));

    min_J = J(ind);
    disp(min_J);

    new_params = [p.k,p.q,p.I_0,p.R_0]+alpha(ind)*best_vec;
    
    p.k = new_params(1);
    p.q = new_params(2);
    p.I_0 = new_params(3);
    p.R_0 = new_params(4);
    
    % Require S_0+I_0+R_0 = 1
    p.S_0 = 1-p.I_0-p.R_0;

end

disp(p);

[SIR,adjoint,new_Sfun,new_Ifun,new_Rfun] = solve_ODEs(p,Idfun,Rdfun,plot_it,20);

if strcmp('plot_it','on')
        plot(SIR.x,SIR.y,'-g','LineWidth',2);
end






