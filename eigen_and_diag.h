struct vector;

//B-matrix(B, L) takes in an array of vector pointers B, and a matrix L. It
//   returns the B-matrix of L through a heap allocated matrix pointer if 
//   possible(client must free the pointer using matrix_destroy). Otherwise
//   it outputs an error message and returns NULL. 
//requires: B and L are not NULL;
//          number of pointers in B matches the dimension of L
//          pointers in B are not NULL.
//effects: may allocate heap memory
//         may print message
struct matrix *B_matrix(const struct vector * const B[], 
                        const struct matrix * const L);

//eigenvalue_2x2(A, lambda1, lambda2) takes in a pointer to a 2 x 2 matrix A.
//   It modifies *lambda1 and *lambda2 to the real eigenvalues of *A if 
//   possible, or outputs an error message and modifies *lambda1 and *lambda2
//   to INT_MIN.
//requires: A, lambda1, lambda2 are not NULL.
//effects: prints message
//         modifies *lambda1 and *lambda2
void eigenvalue_2x2(const struct matrix * const A, long double * const lambda1,
                    long double * const lambda2);

//eigenvalue_3x3(A, lambda1, lambda2, lambda3) takes in a pointer to a 3 x 3
//   matrix A. It modifies *lambda1, *lambda2 and *lambda3 to the real 
//   eigenvalues of *A if possible, or outputs an error message and modifies
//   *lambda1, *lambda2 and *lambda3 INT_MIN.
//requires: A, lambda1, lambda2, lambda3 are not NULL.
//effects: prints message
//         modifies *lambda1, *lambda2 and *lambda3
void eigenvalue_3x3(const struct matrix * const A, long double * const lambda1,
                    long double * const lambda2, long double * const lambda3);

//eigenvectors_2x2(A, v1, v2) takes in a pointer to a 2 x 2 matrix A, and two
//   pointers to struct vector pointers. It changes as many vector pointers as 
//   possible, (starting from *v1, then *v2) to reflect basis of eigenspaces of
//   A. Then it returns the sum of g(lambda) of A. If no eigenvectors can be
//   found, it outputs an error message.
//requires: A, v1, v2 are not NULL;
//effects: may allocate heap memory
//         may modify *v1 and *v2
//         may print output.
int eigenvectors_2x2(const struct matrix * const A, struct vector **v1,
                     struct vector ** const v2);

//eigenvectors_3x3(A, v1, v2, v3) takes in a pointer to a 3 x 3 matrix A, and
//   three struct vector pointers. It changes as many vector pointers as 
//   possible, (starting from v1, then v2, then v3) to reflect basis of 
//   eigenspaces of A. Then it returns the sum of g(lambda) of A. If no 
//   eigenvectors can be found, it outputs an error message.
//requires: A, v1, v2, v3 are not NULL;
//effects: may allocate heap memory
//         may modify *v1, *v2 and *v3
//         may print output
int eigenvectors_3x3(const struct matrix * const A, struct vector ** const v1,
                     struct vector ** const v2, struct vector ** const v3);

//diagonalize_2x2(A, P, D, P_inv) takes in four struct matrix pointers where A
//   is a 2 x 2 matrix. If possible, it modifies P, D and P_inv such that P 
//   diagonalizes A, or (P_inv)(D)(P) = A. Otherwise it outputs an error
//   message and leave P, D and P_inv unchanged.
//requires: A is not NULL;
//          *A is a 2 x 2 matrix.
//effects: may modify *P, *D and *P_inv
//         may allocate heap memory
void diagonalize_2x2(const struct matrix * const A, struct matrix ** const P,
                     struct matrix ** const D, struct matrix ** const P_inv);


//diagonalize_3x3(A, P, D, P_inv) takes in four struct matrix pointers where A
//   is a 3 x 3 matrix. If possible, it modifies P, D and P_inv such that P 
//   diagonalizes A, or (P_inv)(D)(P) = A. Otherwise it outputs an error
//   message and leave P, D and P_inv unchanged.
//requires: A is not NULL;
//          *A is a 3 x 3 matrix.
//effects: may modify *P, *D and *P_inv
//         may allocate heap memory
void diagonalize_3x3(const struct matrix * const A, struct matrix ** const P,
                     struct matrix ** const D, struct matrix ** const P_inv);

