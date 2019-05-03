//matrix_cof(A, i ,j) returns the cofactor A[ij] if possible. Otherwise it 
//   outputs an error message and returns INT_MIN.
//requires: A is not NULL;
//effects: may print message
long double matrix_cof(const struct matrix * const A, const int i, 
                       const int j);

//cof_matrix(A) returns the cofactor matrix of A through a heap allocated 
//   matrix pointer if possible (the client must free the pointer using
//   matrix_destroy).Otherwise it outputs an error message and returns NULL.
//requires: A is not NULL;
//effects: may print message
struct matrix *cof_matrix(const struct matrix * const A);


//adj_matrix(A) returns the adjugate matrix of A through a heap allocated 
//   matrix pointer if possible (the client must free the pointer using
//   matrix_destroy).Otherwise it outputs an error message and returns NULL.
//requires: A is not NULL;
//effects: may print message
struct matrix *adj_matrix(const struct matrix * const A);

//matrix_inverse(A) returns the inverse of A through a heap allocated 
//   matrix pointer if possible (the client must free the pointer using
//   matrix_destroy).Otherwise it outputs an error message and returns NULL.
//requires: A is not NULL;
//effects: may print message
struct matrix *matrix_inverse(const struct matrix * const A);

//matrix_det(A) returns the determinant of A if possible. Otherwise it 
//   outputs an error message and returns INT_MIN.
//requires: A is not NULL;
//effects: may print message
long double matrix_det(const struct matrix * const A);



