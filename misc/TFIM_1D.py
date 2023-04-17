import matplotlib.pyplot as plt
import numpy as np

def pauli_x(n, i):
    """Return the Pauli X operator at site i for an n-site system."""
    return np.kron(np.kron(np.eye(2 ** (n - i - 1)), np.array([[0, 1], [1, 0]])), np.eye(2 ** i))

def pauli_z(n, i):
    """Return the Pauli Z operator at site i for an n-site system."""
    return np.kron(np.kron(np.eye(2 ** (n - i - 1)), np.array([[1, 0], [0, -1]])), np.eye(2 ** i))

def tfim_hamiltonian(n, h):
    """Return the Hamiltonian matrix for a transverse field Ising model with n sites and transverse field strength h."""
    H = np.zeros((2 ** n, 2 ** n))
    
    for i in range(n):
        H -= np.matmul(pauli_z(n, i), pauli_z(n, (i + 1) % n))
        H -= h * pauli_x(n, i)
    
    return H

n = 18  # Number of sites in the 1D chain
h = 10  # Transverse field strength

H = tfim_hamiltonian(n, h)
H = H + ((1+h)*n)*np.eye()
print("Hamiltonian matrix:")
plt.imshow(H, cmap='coolwarm', origin='upper')
plt.colorbar(label='Matrix element value')
plt.title('Transverse field Ising model Hamiltonian')
plt.xlabel('Column index')
plt.ylabel('Row index')
plt.show()
