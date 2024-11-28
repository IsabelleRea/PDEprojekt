L = pi;
T = 1; % tau>0
M = 50; % antal x-steg
N = 100; % antal tidssteg

dx = L / M-1;
dt = T / N;
s = dt / dx^2; % stabilitetskonstant

x = linspace(0, L, M);  % Rumssteg
tau = linspace(0, T, N);  % Tidssteg

% Initialisera u (lösningen)
u = zeros(M, N);  % u[n, m] representerar u(x_n, tau_m)

gx_1 = @(x) x .* (L - x);
gx_2 = @(x) x .* (x < L/2); % (x < L/2) returnerar ett om det uppfylls, annars noll

% g = gx_1(x);
g = gx_2(x);

% initial- och randvillkor
u(:, 1) = g;         % u(x, 0) = g(x)
u(1, :) = 0;         % u(0, tau) = 0
u(end, :) = 0;       % u(L, tau) = 0

for m=2:M-1
    for n=1:N-1 % ej randpunkter
        u(m, n + 1) = s*(u(m+1, n) + u(m-1, n)) + (1-2*s)*u(m, n);
    end
end


figure;
surf(tau, x, u)
grid on;

