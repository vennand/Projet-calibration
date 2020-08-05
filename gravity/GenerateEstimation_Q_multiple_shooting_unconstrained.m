function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting_unconstrained(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

N_cardinal_coor = data.nCardinalCoor;

tau_base = SX.zeros(3,1);
x = SX.sym('x', model.nx);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));

G = SX.sym('G',N_cardinal_coor);
forDyn = @(x,G)[  x(model.idx_v)
    FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), tau_base, G )];
f = Function('f', {x, G}, {forDyn(x,G)});

fJx = Function('fJx', {x}, {L(x)});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
Jx = {};
g={};
lbg = [];
ubg = [];

G = MX.sym('G',N_cardinal_coor);

kalman_q = data.kalman_q;
kalman_v = data.kalman_v;

X_kalman = vertcat(kalman_q, kalman_v);

Xk = MX.sym(['X_' '0'], model.nx);
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

Jx = {Jx{:}, fJx(X_kalman(:,1) - Xk)};

M = 4;
DT = dN/M;
for k=0:Nint-1
%     Xk = RK4('x0',Xk,'p',Uk);
%     Xk = Xk.xf;
    for j=1:M
        k1 = f(Xk, G);
        k2 = f(Xk + DT/2 * k1, G);
        k3 = f(Xk + DT/2 * k2, G);
        k4 = f(Xk + DT * k3, G);

        Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xkend = Xk;
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Jx = {Jx{:}, fJx(X_kalman(:,k+2) - Xk)};
    
    g = {g{:}, Xkend - Xk};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
end

bounds = rotx(data.gravityRotationBound) * vertcat([0;0;0],data.gravity);

w = {w{:}, G};
lbw = [lbw; [-abs(bounds(5)); -abs(bounds(5)); data.gravity(3)]];
ubw = [ubw; [ abs(bounds(5));  abs(bounds(5)); bounds(6)]];

% g = {g{:}, G'*G - data.gravity(:)' * data.gravity(:)};
% lbg = [lbg; -1e-12];
% ubg = [ubg;  1e-12];

Jx = vertcat(Jx{:});
w = vertcat(w{:});
g = vertcat(g{:});
prob = struct('f', sum(Jx), 'x', w, 'g', g);

if nargout > 5
    objFunc = Function('J',  {w}, {Jx});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end