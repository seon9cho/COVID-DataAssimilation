function [sol_SIR,sol_adjoint,Sfun,Ifun,Rfun] = solve_ODEs(p,Idfun,Rdfun,plot_it,T)



% SIR_ODE 
% S = y(1), I = y(2), R = y(3)
SIR_ODE = @(x,y)[-p.k*y(1)*y(2);p.k*y(1)*y(2)-p.q*y(2);p.q*y(2)];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol_SIR = ode15s(SIR_ODE,[0,T],[1-p.I_0-p.R_0,p.I_0,p.R_0],options);

Sfun = @(x)[1,0,0]*deval(sol_SIR,x);
Ifun = @(x)[0,1,0]*deval(sol_SIR,x);
Rfun = @(x)[0,0,1]*deval(sol_SIR,x);

if strcmp(plot_it,'on')
    plot(sol_SIR.x,sol_SIR.y,'-k','LineWidth',2);
end

% temp = deval(sol_SIR,t);
% Idat = I_data-temp(2,:);
% Rdat = R_data-temp(3,:);
% 
% figure;
% plot(t,Idat,'.b')


% -------------------------------------------------------------------------
% Solve the adjoint equations
% -------------------------------------------------------------------------

adjoint_ODE = @(x,P)[p.k*Ifun(x)*(P(1)-P(2)); ...
    2*(Ifun(x)-Idfun(x))+p.k*Sfun(x)*(P(1)-P(2))+p.q*(P(2)-P(3));...
    2*(Rfun(x)-Rdfun(x))];


sol_adjoint = ode15s(adjoint_ODE,[T,0],[0,0,0],options);





