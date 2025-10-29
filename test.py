import sympy as sp
from sympy import Matrix, eye

def analyze_matrix(A):
    A = Matrix(A)
    n = A.shape[0]

    print("Input Matrix:")
    sp.pprint(A)

    eig_data = A.eigenvects()

    eigenvalues = []
    eigenvectors = []

    for val, mult, vecs in eig_data:
        eigenvalues.extend([val] * mult)
        eigenvectors.extend(vecs)

    print("\nEigenvalues:")
    print(", ".join(str(ev) for ev in eigenvalues))

    total_eigvecs = len(eigenvectors)
    if total_eigvecs < n:
        for val, mult, vecs in eig_data:
            I = eye(n)
            M = A - val * I
            gen_vecs = []
            while len(vecs) + len(gen_vecs) < mult:
                for k in range(2, n + 2):
                    Nk = (M ** k).nullspace()
                    if len(Nk) > len(vecs) + len(gen_vecs):
                        new_vecs = [v for v in Nk if v not in vecs + gen_vecs]
                        gen_vecs.extend(new_vecs)
                        break
                if not gen_vecs:
                    break
            eigenvectors.extend(gen_vecs)

    if len(eigenvectors) < n:
        I = eye(n)
        for i in range(n):
            if len(eigenvectors) >= n:
                break
            e = I[:, i]
            if e not in eigenvectors:
                eigenvectors.append(e)

    print("\nEigenvectors (columns):")
    V = Matrix.hstack(*eigenvectors[:n])
    sp.pprint(V)

    print("\nJordan Block:")
    J = sp.zeros(n)
    for i in range(n):
        J[i, i] = eigenvalues[i]
        if i < n - 1 and eigenvalues[i] == eigenvalues[i + 1]:
            J[i, i + 1] = 1
    sp.pprint(J)


A1 = [
    [5, 3, 2],
    [-1, 8, 1],
    [-1, 2, 8]
]

# A2 = [
#     [3,1,-7],
#     [0,3,14],
#     [0,0,5]
# ]

analyze_matrix(A1)
