n = 4; % Number of sites in the 1D chain
h = 1.0; % Transverse field strength

H = tfimHamiltonian(n,h);
imagesc(H);

function x = pauliX(n,i)
    x = kron(kron(eye(2^(n-i-1)),[0,1;1,0]),eye(2^i));
end

function z = pauliZ(n,i)
    z = kron(kron(eye(2^(n-i-1)),[1,0;0,-1]),eye(2^i));
end

function H = tfimHamiltonian(n,h)
    H = zeros(2^n,2^n);
    for i = 1:n
        H = H - pauliZ(n,i-1)*pauliZ(n, mod(i,n));
        H = H - h*pauliX(n,i-1);
    end
end