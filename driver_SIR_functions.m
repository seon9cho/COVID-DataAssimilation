%{

SIR model

%}
clc; clear all; close all; beep off; 

% -------------------------------------------------------------------------
% actual solution
% -------------------------------------------------------------------------

plot_it = 'on';

p.k = @(t) 1-t.^2/100;
p.q = @(t) 1-0.1*tanh(t);
p.I_0 = 0.05;
p.R_0 = 0;

p_true = p;

% SIR_ODE 
% S = y(1), I = y(2), R = y(3)
SIR_ODE = @(t,y)[-p.k(t)*y(1)*y(2);p.k(t)*y(1)*y(2)-p.q(t)*y(2);p.q(t)*y(2)];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol_actual = ode15s(SIR_ODE,[0,20],[1-p.I_0-p.R_0,p.I_0,p.R_0],options);

if strcmp(plot_it,'on')
    figure; hold on;
    plot(sol_actual.x,sol_actual.y,'-b','LineWidth',2);
end

% -------------------------------------------------------------------------
% get data
% -------------------------------------------------------------------------

t = linspace(0,3,300);

temp = deval(sol_actual,t);

I_data = temp(2,:)+normrnd(zeros(1,length(t)),0.001);

if strcmp(plot_it,'on')
    plot(t,I_data,'.c','MarkerSize',8)
end

% -------------------------------------------------------------------------
% interpolation
% -------------------------------------------------------------------------

poly_I = polyfit(t,I_data,3);
Idfun = @(t) poly_I(1)*t.^3+poly_I(2)*t.^2+poly_I(3)*t+poly_I(4);

if strcmp(plot_it,'on')
    plot(t,Idfun(t),'-m','LineWidth',2);
    drawnow;
end

% -------------------------------------------------------------------------
% gradient descent
% -------------------------------------------------------------------------

% initial guess
p.k = @(t) ones(1,length(t));
p.q = @(t) ones(1,length(t));
p.I_0 = Idfun(0);
p.R_0 = 0;

% Require S_0+I_0+R_0 = 1
p.S_0 = 1-p.I_0-p.R_0;
if p.S_0 < 0
    error('S_0 < 0')
end


for cnt = 1:50

    cnt

    % -------------------------------------------------------------------------
    % Solve SIR for guess
    % -------------------------------------------------------------------------
    
    [sol_SIR,sol_adjoint,Sfun,Ifun,Rfun] = solve_ODEs(p,Idfun,0,'off',5);
    
    % -------------------------------------------------------------------------
    % gradients
    % -------------------------------------------------------------------------
    
    P0 = deval(sol_adjoint,0);
    
    PS = @(x)[1,0,0]*deval(sol_adjoint,x);
    PI = @(x)[0,1,0]*deval(sol_adjoint,x);
    PR = @(x)[0,0,1]*deval(sol_adjoint,x);
    
    J_S0 = -P0(1);
    J_I0 = -P0(2);
    J_R0 = -P0(3);

    % find the direction where the derivative is greatest.
    best_vec = randn(1,10);
    best_der = 1e20;

    % polynomials of the form $p(x) = a_0+a_1*x+a_2*x^2+a_3*x^3$
    for q = 1:15

        vec = randn(1,10);
        
        delta_k = @(x) (vec(1:4)*[ones(1,length(x));x;x.^2;x.^3]);
        delta_q = @(x) (vec(5:8)*[ones(1,length(x));x;x.^2;x.^3]);

        fun1 = @(x) ((PS(x)-PI(x)).*Sfun(x).*Ifun(x)).*delta_k(x);
        J_k = integral(fun1,0,5);
    

        fun2 = @(x) ((PI(x)-PR(x)).*Ifun(x).*delta_q(x));
        J_q = integral(fun2,0,5);

        der = abs(J_k+J_q+(J_I0-J_S0)*vec(9)+(J_R0-J_S0)*vec(10));

        fk = @(x) (delta_k(x).^2);
        fq = @(x) (delta_q(x).^2);

        vec_norm_squared = integral(fk,0,5)+integral(fq,0,5)+vec(9)^2+vec(10)^2;
        der = der/sqrt(vec_norm_squared);

        if der > best_der
            best_vec = vec;
            best_der = der;
        end
    end

    alpha = [0.01*(rand(10,1)-0.5);0]; % values for line search
    
    delta_k = @(x) (best_vec(1:4)*[ones(1,length(x));x;x.^2;x.^3]);
    delta_q = @(x) (best_vec(5:8)*[ones(1,length(x));x;x.^2;x.^3]);

    J = zeros(size(alpha));
    for j = 1:length(alpha)
    
        p1.k = @(t) p.k(t)+alpha(j)*delta_k(t);
        p1.q = @(t) p.q(t)+alpha(j)*delta_q(t);
        p1.I_0 = p.I_0+alpha(j)*best_vec(9);
        p1.R_0 = p.R_0+alpha(j)*best_vec(10);


    
        p1.S_0 = 1-p1.I_0-p1.R_0;
  
        [new_SIR,new_adjoint,new_Sfun,new_Ifun,new_Rfun] = solve_ODEs(p1,Idfun,0,'off',5);
    
        f3 = @(x) (new_Ifun(x)-Idfun(x)).^2;
    
        J(j) = integral(f3,0,5);
    
    end
    
    ind = find(J==min(J));

    min_J = J(ind);
    disp(min_J);

    % delta_k = @(x) (vec(1:4)*[ones(1,length(x));x;x.^2;x.^3]);
    % delta_q = @(x) (vec(5:8)*[ones(1,length(x));x;x.^2;x.^3]);

    p.k = @(t) (p.k(t)+alpha(ind)*delta_k(t));
    p.q = @(t) (p.q(t)+alpha(ind)*delta_q(t));
    p.I_0 = p.I_0+alpha(ind)*best_vec(9);
    p.R_0 = p.R_0+alpha(ind)*best_vec(10);
    
    % Require S_0+I_0+R_0 = 1
    p.S_0 = 1-p.I_0-p.R_0;


    SIR_ODE = @(t,y)[-p.k(t)*y(1)*y(2);p.k(t)*y(1)*y(2)-p.q(t)*y(2);p.q(t)*y(2)];

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    sol = ode15s(SIR_ODE,[0,20],[1-p.I_0-p.R_0,p.I_0,p.R_0],options);

    % if strcmp(plot_it,'on')
    %     figure; hold on;
    %     plot(sol_actual.x,sol_actual.y,'-b','LineWidth',2);
    %     plot(sol.x,sol.y,'-g','LineWidth',2);
    %     plot(t,I_data,'.c','MarkerSize',8);
    %     plot(t,Idfun(t),'-m','LineWidth',2);
    %     drawnow;
    % 
    % end

    

end

disp(p);

[SIR,adjoint,new_Sfun,new_Ifun,new_Rfun] = solve_ODEs(p,Idfun,0,plot_it,20);

if strcmp(plot_it,'on')
    figure; hold on;
    plot(sol_actual.x,sol_actual.y,'-b','LineWidth',2);
    plot(sol.x,sol.y,'-g','LineWidth',2);
    plot(t,I_data,'.c','MarkerSize',8);
    plot(t,Idfun(t),'-m','LineWidth',2);
    drawnow;
end






