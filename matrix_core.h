//You have all seen a vector before, but now...
struct vector;
//Whoa, a matrix.
struct matrix;

//matrix_create() returns a heap-allocated struct matrix pointer that caller 
//   must free using matrix_destroy().
//effects: allocates heap memory
struct matrix *matrix_create();

//quick_matrix_input(values, m, n) takes an array of long doubles and two 
//   integers m and n. If possible, it returns a heap-allocated pointer to an
//   m by n matrix formed with the numbers; the caller must free this pointer
//   with matrix_destroy(). If matrix creation is not possible, it output an
//   error message, and returns NULL.
//requires: values is not NULL
//          for positive m and n, values contains exactly m * n long doubles.
//effects: allocates heap memory 
struct matrix *quick_matrix_input(const long double values[], const int m,
                                   const int n);

//matrix_size(A, m, n) takes in a struct matrix pointer and two int pointers.
//  it then modifies *m and *n to reflect the number of rows and columns of 
//  A.
//requires: A, m, n are not NULL
//effects: may modify *m and *n
void matrix_size(const struct matrix *const A, int * const m, int * const n);

//matrix_add_row(A, v1) takes in a struct vector pointer and a matrix pointer. 
//   It then attempts to add the corresponding vector as a row in the matrix.
//   It outputs an error if the operation cannot be done.
//requires: A and v1 are not NULL;
//effects: may modify *A
//         may print message
void matrix_add_row(struct matrix * const A, const struct vector * const v1);

//matrix_del_row(A, m) takes in a matrix pointer, and removes row m of 
//   *A if possible. Otherwise it outputs an error message.
//requires: A is not NULL;
//effects: may modify *A
//         may print message
void matrix_del_row(struct matrix * const A, const int m);

//matrix_replace_row(A, index, v1) takes in a matrix pointer A, an int index 
//   and a vector v1. It replaces the index-th row of *A with *v1 if possible.
//   otherwise it outputs an error message.
//requires: A and v1 are not NULL;
//effects: may modify *A
//         may print message
void matrix_replace_row(struct matrix * const A, const int index,
                        const struct vector * const v1);

//matrix_dupe_row(A, index) returns the index-th row of the matrix through
//   a vector pointer if possible. Otherwise it outputs and error message and
//   returns NULL.
//effects: may allocate heap memory
//         may print message
struct vector *matrix_dupe_row(const struct matrix * const A, const int index);


//matrix_swap_row(A, r1, r2) swaps row r1 and r2 in matrix A if possible. 
//   Otherwise it outputs an error message.
//effects: may print message
void matrix_swap_row(struct matrix * const A, const int r1, const int r2);

//matrix_sum_row(A, r1, r2) adds row r2 to r1 of matrix A if possible. 
//   Otherwise it outputs an error message.
//effects: may print message
void matrix_sum_row(struct matrix * const A, const int r1, const int r2);

//matrix_mult_row(A, r1, c) multiplies row r1 in A by c if possible. 
//   Otherwise it outputs and error message.
//effects: may print message
void matrix_mult_row(struct matrix * const A, const int r1, 
                     const long double c);

//matrix_add_mult_row(A, r1, r2, c) adds c times row r2 to row r1 in A
//   if possible. Otherwise it outputs and error message.
//effects: may print message
void matrix_add_mult_row(struct matrix * const A, const int r1, const int r2, 
                         const long double c);



//matrix_add_col(A, v1) takes in a struct vector pointer and a matrix pointer. 
//   It then attempts to add the corresponding vector as a column in the 
//   matrix. It outputs an error if the operation cannot be done.
//requires: A and v1 are not NULL;
//effects: may modify *A
//         may print message
void matrix_add_col(struct matrix * const A, const struct vector * const v1);

//matrix_del_col(A, n) takes in a matrix pointer, and removes column n of 
//   *A if possible. Otherwise it outputs an error message.
//requires: A is not NULL;
//effects: may modify *A
//         may print message
void matrix_del_col(struct matrix * const A, const int n);

//matrix_replace_col(A, index, v1) takes in a matrix pointer A, an int index 
//   and a vector v1. It replaces the index-th column of *A with *v1 if
//   possible. Otherwise it outputs an error message.
//requires: A and v1 are not NULL;
//effects: may modify *A
//         may print message
void matrix_replace_col(struct matrix * const A, const int index, 
                        const struct vector * const v1);

//matrix_dupe_col(A, index) returns the index-th column of the matrix through
//   a vector pointer if possible. Otherwise it outputs and error message and
//   returns NULL.
//requires: A is not NULL;
//effects: may allocate heap memory
//         may print message
struct vector *matrix_dupe_col(const struct matrix * const A, const int index);

//matrix_dupe(A) returns a new heap-allocated pointer to a duplicate copy of
//   the matrix passed in if possible. Otherwise it outputs an error message.
//requires: A is not NULL;
//effects: may allocate heap memory
//         may print message
struct matrix *matrix_dupe(const struct matrix * const A);

//matrix_elem(A, m, n) returns A[mn] if possible, and an error message if not.
//requires: A is not NULL;
//effects: may allocate heap memory
//         may print message
long double matrix_elem(const struct matrix * const A, const int m,
                        const int n);


//matrix_print(A) takes in a pointer to a matrix and prints the matrix.
//effects: prints output
void matrix_print(const struct matrix * const A);

//matrix_destroy(A) frees heap memory allocated to A if it is not NULL
//requires: A is not NULL
void matrix_destroy(struct matrix * const A);
 



