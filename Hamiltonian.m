% Parameters
n = 7; % Chain length (number of sites)
J = 1; % Coupling constant
h = 0.5; % Transverse field

% Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_z = [1, 0; 0, -1];

% Initialize an empty 2^n x 2^n Hamiltonian matrix
H = zeros(2^n);

% Generate Hamiltonian matrix
for i = 1:n
    % Identify the interacting spins
    j = mod(i, n) + 1; % For periodic boundary conditions
    
    % Tensor product of identity matrices before the i-th site
    pre_site = kron(eye(2^(i-1)), sigma_z);
    
    % Tensor product of identity matrices between the i-th and j-th sites
    if j - i > 1
        mid_site = kron(eye(2^(j-i-1)), sigma_z);
        % Tensor product of identity matrices after the j-th site
        post_site = eye(2^(n-j));
        % Add the contribution of the σ_i^z * σ_j^z term
        H = H - J * kron(kron(pre_site, mid_site), post_site);
    else
        % Tensor product of identity matrices after the j-th site
        post_site = eye(2^(n-j));
        % Add the contribution of the σ_i^z * σ_j^z term
        H = H - J * kron(pre_site, post_site);
    end
    
    % Add the contribution of the σ_i^x term
    pre_site_x = kron(eye(2^(i-1)), sigma_x);
    H = H - h * kron(pre_site_x, eye(2^(n-i)));
end
