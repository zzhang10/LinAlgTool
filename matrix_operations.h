#include <stdbool.h>

//matrix_add(A, B) takes in two struct matrix pointers and outputs the 
//   result of A+B through a matrix pointer if possible (client must free the
//   pointer). Otherwise it prints an error message and returns NULL.
//requires: A, B are not NULL.
//effects: may print output
//         may allocate heap memory
struct matrix *matrix_add(const struct matrix * const A, 
                          const struct matrix * const B);

//matrix_mult_scalar(A, c) takes in a struct matrix pointer and a long double.
//   It returns c(A) through a matrix pointer if possible (client must free 
//   the pointer). Otherwise it prints an error message and returns NULL.
//requires: A, v1 are not NULL.
//effects: may print output
//         may allocate heap memory
struct matrix *matrix_mult_scalar(const struct matrix * const A,
                                  const long double c);


//matrix_mult_vector(A, v1) takes in a struct matrix pointer and a struct 
//   vector pointer. It returns A(v1) through a vector pointer if possible 
//   (client must free the pointer). Otherwise it prints an error message and 
//   returns NULL.
//requires: A, v1 are not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *matrix_mult_vector(const struct matrix * const A,
                                  const struct vector * const v1);

//matrix_mult_matrix(A, B) takes in two struct matrix pointers and outputs the 
//   result of AB through a matrix pointer if possible (client must free the
//   pointer). Otherwise it prints an error message and return NULL.
//requires: A, B are not NULL.
//effects: may print output
//         may allocate heap memory
struct matrix *matrix_mult_matrix(const struct matrix * const A, 
                                  const struct matrix * const B);

//rotation_matrix(theta) takes in a long double as the angle in RADIANS. Then 
//  it outputs the rotation matrix R[theta] through a matrix pointer if 
//   possible (client must free the pointer). Otherwise it prints an error 
//   message and returns NULL.
//effects: may print output
//         may allocate heap memory
struct matrix *rotation_matrix(const long double theta);

//is_RREF(A) returns true if A represents a matrix in RREF, false otherwise.
//requires: A is not NULL, *A is not empty.
//effects: may print output
bool is_RREF(const struct matrix * const A);

//RREF(A)takes in a struct matrix pointer A, and returns the RREF of A through
//   a matrix pointer if possible (client must free the pointer). Otherwise
//   it prints an error message and returns NULL. 
//requires: A is not NULL;
//effects: may print output
//         may allocate heap memory
struct matrix *RREF(const struct matrix * const A);


//matrix_transpose(A) takes in a struct matrix pointer A, and returns the 
//   transposeof A through a matrix pointer if possible (client must free 
//   the pointer). Otherwise it prints an error message and returns NULL.
//requires: A is not NULL;
//effects: may print output
//         may allocate heap memory
struct matrix *matrix_transpose(const struct matrix * const A);

//matrix_rank(A) takes in a struct matrix pointer A, and returns the 
//   rank of A if possible, Otherwise it prints an error message and returns 
//   INT_MIN.
//requires: A is not NULL;
//effects: may print output
int matrix_rank(const struct matrix * const A);

//matrix_power(A, n) takes in a struct matrix pointer A and an integer n, and
//   returns A^n through a heap_allocated struct matrix pointer if possible (
//   the client must free the pointer with matrix_destroy). Otherwise it 
//   outputs and error message and returns NULL.
//requires: A is not NULL.
//effects: may print output
//         may allocate heap memory
struct matrix *matrix_power(const struct matrix * const A, const int n);

