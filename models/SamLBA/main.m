function main()

% Irreversible stoichiometric matrix
S = [-2 0 -1 0 0;
    -1 0 0 0 0;
    0 -1 0 -2 0;
    1 -1 0 0 -2;
    0 1 -2 0 0;
    0 0 2 -1 0;
    0 0 1 -1 0;
    0 0 0 2 0;
    0 0 0 1 -1;
    0 0 0 0 1];

% Reversible stoichiometric matrix plus uptake
S = [S -S [eye(3) ; zeros(7,3)]];

% Number of metabolites
n = size(S,1);

% Number of reactions
m = size(S,2);

% Doubling time (hours)
T = 1;

% Growth rate
gr = log(2)/T;

% Biomass stoichiometry vector
p = 10;
b_mean = [0 0 0 0 3 1 0 2 1 0]';
b = [];
for i = 1 : p
    b = [b b_mean+(~(b_mean==zeros(n,1))).*(.5*randn(n,1))];
end
b = b./kron(ones(n,1),sum(b));
b_av = zeros(n,1);

% Stoichiometrix matrix plus biomass reaction
S_b_av = [S -gr*b_av -eye(n)];
for i = 1 : p
    S_b{i} = [S -gr*b(:,i) -eye(n)];
end

% Flux bounds
v_max = [100*ones(1,10) .1 .1 .1 100*ones(1,11)]';

% Robust Optimization
r = 1.5;
e_b = zeros(n+m+1,1); e_b(m+1) = 1;
e_vec = [0:.025:1];%[.4 : .025 : .85];

count = 0;
for e = e_vec
    count = count + 1;
    
    cvx_begin
    variable v(n+m+1,1);
    variable l(1,1);
    variable t(1,1);
    minimize (e*l-(1-e)*e_b'*v);
    subject to;
    t >= 0;
    r^2*t <= l;
    v >= 0;
    v <= v_max;
    [l-r^2*t zeros(1,p) (S_b_av*v)';
        zeros(p,1) t*eye(p) [S_b{1}*v S_b{2}*v S_b{3}*v S_b{4}*v S_b{5}*v S_b{6}*v S_b{7}*v S_b{8}*v S_b{9}*v S_b{10}*v]';
        S_b_av*v [S_b{1}*v S_b{2}*v S_b{3}*v S_b{4}*v S_b{5}*v S_b{6}*v S_b{7}*v S_b{8}*v S_b{9}*v S_b{10}*v] eye(n)] -1e-2*eye(n+p+1) == semidefinite(n+p+1);
    cvx_end
    
    for i = 1 : p
        constraint(count,i) = norm(S_b{i}*v);
    end
    objective(count) = e_b'*v;
    lambda(count) = l;
    v_opt{count} = v;
    
end

% Display the net (irreversible) fluxes
[v(1:5)-v(6:10) ; v(11:24)]

% Plots
hold on;
for i = 1 : p
    h3 = semilogy(e_vec,constraint(:,i),'-b');
end
h1 = semilogy(e_vec,objective,'-r','LineWidth',1);
h2 = semilogy(e_vec,lambda,'--k','LineWidth',1);
axis([min(e_vec) max(e_vec) 0 1.01*max(objective)]);
h = legend([h1 h2 h3],'$e_b^T v_b$','$\lambda$','$\|S_{b_i}v_b\|_2$');
set(h,'Interpreter','Latex','FontSize',24,'Box','off','Location','SouthWest')
grid on;
box on;
hold off;

end