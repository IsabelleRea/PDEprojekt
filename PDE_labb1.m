format long

M = 1000;
k = 2*pi;
N = 6;
K = 10;

x1 = rand(M, 1);
x2 = rand(M, 1);

xp = [x1 x2];

[n1, n2] = ndgrid(-N:N, -N:N);

f_hatt = zeros(size(n1));
v1_hatt = zeros(size(n1));
v2_hatt = zeros(size(n1));

f_square_sum = zeros(size(M));
v1_square_sum = zeros(size(M));
v2_square_sum = zeros(size(M));

f_reconstructed = zeros(size(M));
v1_reconstructed = zeros(size(M));
v2_reconstructed = zeros(size(M));

% Loop over samples
for p = 1:M
    % Compute f(x1(p), x2(p))
    f = f_function(x1(p), x2(p));
    [v1, v2] = v(x1(p), x2(p));
    
    % Compute the exponential term
    exp_term = exp(-1i * (n1 * x1(p) + n2 * x2(p)));
    
    % Accumulate f_hatt
    f_hatt = f_hatt + (1 / M) * f * exp_term;
    v1_hatt = v1_hatt + (1 / M) * v1 * exp_term;
    v2_hatt = v2_hatt + (1 / M) * v2 * exp_term;
end

% Compute f_square_sum
for p = 1:M
    % Compute f(x1(p), x2(p))
    f = f_function(x1(p), x2(p));
    [v1, v2] = v(x1(p), x2(p));
    
    % Compute the exponential term
    exp_term = exp(-1i * (n1 * x1(p) + n2 * x2(p)));
    
    % Reconstruct f from f_hatt
    f_reconstructed_p = sum(f_hatt .* exp_term, 'all');
    v1_reconstructed_p = sum(v1_hatt .* exp_term, 'all');
    v2_reconstructed_p = sum(v2_hatt .* exp_term, 'all');

    f_reconstructed(p) = f_reconstructed_p;
    v1_reconstructed(p) = v1_reconstructed_p;
    v2_reconstructed(p) = v2_reconstructed_p;
    
    % Compute quadratic difference
    f_square_sum(p) = abs(f - f_reconstructed_p).^2;
    v1_square_sum(p) = abs(v1 - v1_reconstructed_p).^2;
    v2_square_sum(p) = abs(v2 - v2_reconstructed_p).^2;
end

min_square_sum_f = min(f_square_sum);
index_of_min_f = find(f_square_sum == min_square_sum_f);

f_streck = f_reconstructed(index_of_min_f)

min_square_sum_v1 = min(v1_square_sum);
index_of_min_v1 = find(v1_square_sum == min_square_sum_v1);

v1_streck = v1_reconstructed(index_of_min_v1)

min_square_sum_v2 = min(v2_square_sum);
index_of_min_v2 = find(v2_square_sum == min_square_sum_v2);

v2_streck = v2_reconstructed(index_of_min_v2)

v1_array = zeros(K);
v2_array = zeros(K);



for i=1:K
    v1_array(i) = euler(v1_array(i-1), prim1, 1/K);
    v2_array(i) = euler(v2_array(i-1), prim2, 1/K);
end



%% 
function f = f_function(x, y)
    if 0.4 <= x && x <= 0.6 && 0.4 <= y && y <= 0.6
        f = 1;
    else
        f = 0;
    end
end


%%
function [v1, v2] = v(x, y)
    v1= y;
    v2 = 1-x;
end


%%
function u = euler(h, u0, u_prim)
    u = u0 + h*u_prim;
end






