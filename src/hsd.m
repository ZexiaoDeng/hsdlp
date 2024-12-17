function [x, y, phi, z, w, psi, fval, status] = hsd(A, b, c)
% function : Implementation of the Homogeneous Self-Dual Long Step Method
%
% primal problem:
%   max   c'x
%   s.t.  Ax <= b
%         x >= 0
% 
% dual problem:
%   min   b'y
%   s.t.  A'y >= c
%           y >= 0
%
% dx, dy, dz : step direction
% fx, fy, gx, gy :
% phi, psi, dphi, dpsi :
% D, E : diagonal matrices
% gamma, beta, delta, mu, theta : parameters
%

EPS = 1.0e-12;
MAX_ITER = 600;

status = 5;

% Allocate memory for arrays.
n = size(c(:), 1);
m = size(b(:), 1);

% Initialization
x = ones(n, 1);
z = ones(n, 1);
y = ones(m, 1);
w = ones(m, 1);

phi = 1.0;
psi = 1.0;

fprintf('HSD method\n');
fprintf('-----------------------------------------------------------------------------\n');
fprintf('%4s   | %22s     | %22s     | %8s |\n', '', 'PRIMAL', 'DUAL', '');
fprintf('%4s   | %14s    %8s | %14s    %8s | %8s |\n', 'iter', 'Obj Value',...
    'Infeas', 'Obj Value', 'Infeas', 'mu');
fprintf('-----------------------------------------------------------------------------\n');

% Iteration.
% beta = 0.80;
% delta = 2 * (1 - beta);

for iter = 1: MAX_ITER,
    % STEP 1: Compute mu.
    mu = (z'*x + w'*y + psi*phi) / (n + m + 1);
    if (mod(iter, 2) ~= 0)
        delta = 0;
    else
        delta = 1.0;
    end
    % STEP 1: Compute primal and dual objective function values.
    primal_obj = c'*x;
    dual_obj = b'*y;
    % STEP 2: Check stopping rule.
    if (mu < EPS),
        if (phi > EPS), status = 0; break;            % optimal
        elseif (dual_obj < 0.0),  status = 2; break;  % primal infeasible
        elseif (primal_obj > 0.0), status = 4; break; % dual infeasible
        else status = 7; break;                       % numerical problem
        end
    end
    % STEP 3: Compute infeasibilities.
    rho = A*x - b*phi + w;
    normr = norm(rho, 2)/phi;
    rho_hat = -(1 - delta)*rho + w - delta*mu./y;
    sigma = -A'*y + c*phi + z;
    norms = norm(sigma, 2)/phi;
    sigma_hat = -(1 - delta)*sigma + z - delta*mu./x;
    gamma_hat = -(1 - delta)*(dual_obj - primal_obj + psi) + psi - delta*mu/phi;

    fprintf('%4d   | %14.7e    %8.1e | %14.7e    %8.1e | %8.1e |\r',...
        iter, primal_obj/phi, normr, dual_obj/phi, norms, mu);

    % STEP 4: Compute step directions.
    X = diag(x); Y = diag(y); Z = diag(z); W = diag(w);
    D = X^(-1)*Z;
    E = Y^(-1)*W;
    quasiD = [-E, A; A', D];
    fyfx = quasiD\[rho_hat; -sigma_hat;];
    gygx = quasiD\[-b; -c];
    
    fy = fyfx(1: m, 1);
    fx = fyfx(m + 1: end, 1);
    gy = gygx(1: m, 1);
    gx = gygx(m + 1: end, 1);
    
    dphi = (c'*fx - b'*fy + gamma_hat)/(c'*gx - b'*gy - psi/phi);
    
    dx = fx - gx*dphi;
    dy = fy - gy*dphi;
    dz = delta*mu./x - z - D*dx;
    dw = delta*mu./y - w - E*dy;
    dpsi = delta*mu/phi - psi - (psi/phi)*dphi;
    
    % STEP 5: Compute step length.
    if (mod(iter,2) == 0),
        theta = 1.0;
        continue;
    else
        theta = 0.0;
        % Attention!!!!!
        theta = max([theta; -dx./x; -dy./y; -dz./z; -dw./w; -dphi/phi; -dpsi/psi]);
        theta = min(0.95/theta, 1.0);
    end
    
%     if (theta < 4*beta/(n+m+1)),
%         printf('ratio = %10.3e \n', theta*(n+m+1)/(4*beta));
%         status = 7;
%         break;
%     end
%     if (theta < 1.0), theta = theta * 0.9999; end
    
    % STEP 6: Update the iteration.
    x = x + theta*dx;
    y = y + theta*dy;
    z = z + theta*dz;
    w = w + theta*dw;
    phi = phi + theta*dphi;
    psi = psi + theta*dpsi;
    
end
x = x/phi;
z = z/phi;
y = y/phi;
w = w/phi;
fval = primal_obj/phi;
end
% End of solver

