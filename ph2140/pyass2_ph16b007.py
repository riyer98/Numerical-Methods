#Q1
def comm(A,B):
    """this function takes 2 matrices or arrays A and B as input.
    then it prints the commutator AB-BA and anticommutator AB+BA.
    it uses the dot function from numpy for the matrix multipliations.

    the matrices must be multipliable both ways with AB and BA having the
    same dimensionality. thus we can conclude that the input must only be
    square matrices of the same dimension or else the output will be invalid."""
    import numpy as np
    print("commutator\n",np.dot(A,B)-np.dot(B,A))
    print("anticommutator\n",np.dot(A,B)+np.dot(B,A))


#Q2
def mat_op():
    """this function takes no input. it randomly generates a 3x3 matrix using
    the np.random.random command. then it prints the conjugate, transpose,
    hermitian adjoint, determinant, inverse and 100th power of the matrix.

    it then prints the product A.A inverse, which comes out to be the identity
    matrix of order 3. due to limited prescion in python there maybe non-zero
    terms of the order of 10^-17 to 10^-16."""
    import numpy as np
    A=np.random.random((3,3))+1j*np.random.random((3,3))
    print("random matrix\n",A)
    print("conjugate\n",A.conj())
    print("transpose\n",A.T)
    print("hermitian adjoint\n",A.conj().T)
    print("determinant\n",np.linalg.det(A))
    B=np.linalg.inv(A)
    print("inverse\n",B)
    print("100th power\n",np.linalg.matrix_power(A,100))
    print("product of A and A inverse\n",np.dot(A,B))


#Q3
def real_p(n,k):
    """this program takes the dimensionality of matrix n and number of
    realisations k as input. it generates a matrix with components being
    standard normal numbers, k times.

    in the array of eigenvalues, if the program encounters even a single
    complex number, it will increment an index j by 1. then the probability
    that at least 1 eigenvalue is complex is given by j/k (value of j after
    k iterations).
    thus the probablility that no eigenvalues are complex, i.e., all are real,
    is given by 1-j/k.

    to determine if an eigenvalue is real, the magnitude of its imaginary
    component must be less than some r>0, which is the presicion limit of
    python. the number here is ideally 2^-52 (as 52 bits are for fraction) which
    is about 2.2*10^-16. to be on the safe side i have taken 3*10^-16 as the
    error can get magnified sometimes, while performing operations with complex
    numbers."""
    import numpy as np
    j=0
    for i in range(k):
        A=np.random.standard_normal(size=(n,n))
        eig=np.linalg.eigvals(A)
        for l in eig:
            if abs(l.imag)>=3*10**-16:
                j+=1
                break
            
    return (1-(j/k))


#Q4
def eplot(h):
    """this program plots the energy as a function of theta (for a given input
    h). the eigenvalues of the hamiltonian matrix are appended to lists, as is
    the value of theta.
    the energy values are independent of theta for any given h and come out to
    be a straight line parallel to the x-axis. the magnitude of the energies
    depends on the value of h."""
    import numpy as np
    import pylab as plt
    sx=np.array([[0,1],[1,0]])
    sy=np.array([[0,-1j],[1j,0]])
    sz=np.array([[1,0],[0,-1]])
    th=[]
    E1=[]
    E2=[]
    for t in np.linspace(-np.pi,np.pi):
        H=sz + h*np.cos(t)*sx + h*np.sin(t)*sy
        eval1=np.linalg.eigvalsh(H)
        th.append(t)
        E1.append(eval1[0])
        E2.append(eval1[1])

    plt.plot(th,E1)
    plt.plot(th,E2)
    plt.xlabel('theta')
    plt.ylabel('energy')
    plt.title('energy as function of theta')
    plt.show()

