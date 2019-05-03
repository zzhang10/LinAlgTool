#include <stdbool.h>

//vector_list_print(vector_list, n) takes in an array of struct vector pointers
//   and prints the first n vectors in the array if possible. Otherwise it 
//   outputs an error message.
//requires: vector_list is not NULL;
//          there are at least n vector pointers in vector_list
//          first n pointers in vector_list are not NULL;
//effects: prints output
void vector_list_print(const struct vector * const vector_list[], const int n);


//linearly_independent(vector_list, n) takes in an array of struct vector 
//   pointers and an integer n. It returns true if the first n vectors in
//   the list of vectors are linearly independent, and false otherwise.
//requires: vector_list is not NULL;
//          there are at least n vector pointers in vector_list
//          first n pointers in vector_list are not NULL;   
//          n >= 0
bool linearly_independent(const struct vector * const vector_list[], 
                          const int n);

//is_basis(vector_list, dim, n) takes in an array of struct vector pointers 
//   and an integer n. It returns true if the first n vectors in vector_list
//   form a basis for R[dim], and false otherwise.
//requires: vector_list is not NULL;
//          there are at least n vector pointers in vector_list
//          first n pointers in vector_list are not NULL; 
//          n >= 0
//          dim > 0
bool is_basis(const struct vector * const vector_list[], const int dim,
              const int n);

//in_span(vector_list, n, v1) takes in an array of struct vector pointers,
//   an integer n and a vector pointer v1. It returns true if v1 is in the
//   span of the first n vectors, or false otherwise.
//requires: vector_list and v1 are not NULL;
//          there are at least n vector pointers in vector_list
//          first n pointers in vector_list are not NULL; 
//          n >= 0
bool in_span(const struct vector * const vector_list[], const int n, 
             const struct vector * const v1);

//find_basis(vector_list, n) takes in an array of struct vector pointers 
//   and an integer n. If possible, it returns a heap-allocated pointer to a
//   matrix whose column vectors form a basis for the span of first n vectors.
//   The client must free this pointer with matrix_destroy(). If a basis cannot
//   be found, it outputs an error message and return NULL.
//requires: vector_list is not NULL;
//          there are at least n vector pointers in vector_list
//          first n pointers in vector_list are not NULL; 
//          n >= 0
//effects: may print message
//         may allocate heap memory
struct matrix *find_basis(const struct vector * const vector_list[], 
                          const int n);

//B_coord(basis, n, v1) takes in an array of struct vector pointers basis, an
//   integer n, and a struct vector pointer v1. If possible, it returns the 
//   B-coordinate of v1 where the basis is the first n vectors in basis,
//   through a heap allocated struct vector pointer that the client must free 
//   using vector_destroy(). If B-coordinate cannot be found, it prints an
//   error message and returns NULL.
//requires: basis and v1 are not NULL
//          there are at least n pointers in basis
//          first n pointers in basis are not NULL
//effects: may print message
//         may allocate heap memory
struct vector *B_coord(const struct vector * const basis[], const int n,
                       const struct vector * const v1);

//chance_of_coord_matrix(B1, B2, n) takes in two arrays of struct vector 
//   pointers B1 and B2 and an integer n. It returns a heap-allocated pointer 
//   to the change of coordinates matrix from the first n vectors of B2 to the 
//   first n vectors of B1 if possible. The client must free the memory with 
//   matrix_destroy(). If the change of coordinate matrix cannot be found, it 
//   prints and error message and returns NULL.
//requires: B1 and B2 are not NULL
//          B1 and B2 contain at least n pointers
//          first n pointers in B1 and B2 are not NULL
//effects: may print message
//         may allocate heap memory
struct matrix *change_of_coord_matrix(const struct vector * const B1[], 
                                      const struct vector * const B2[],
                                      const int n);




