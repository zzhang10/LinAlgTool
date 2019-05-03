#include "vector_core.h"
#include "matrix_core.h"
#include "matrix_operations.h"
#include "vector_space.h"
#include <assert.h>
#include <stdio.h>

//See header file for documentation

void vector_list_print(const struct vector * const vector_list[], 
                       const int n) {
  assert(vector_list);
  for (int i = 0; i < n; i++) {
    assert(vector_list[i]);
  }
  if (n <= 0) {
    printf("Nothing to be printed when n is non-positive.\n");
  } else {
    for (int i = 0; i < n; i++) {
      printf("Vector %d of %d: ", i + 1, n);
      vector_print(vector_list[i]);
    }
  }
}

//vector_list_valid(vector_list, n) returns true if vector_list contains
//   pointers to n vectors with same dimesion. Otherwise it outputs an error
//   message and returns false.
//requires: vector_list is not NULL
//          first n pointers in vector_list is not NULL
//effects: prints message
static bool vector_list_valid(const struct vector * const vector_list[], 
                              const int n) {
  assert(vector_list);
  for (int i = 0; i < n; i++) {
    assert(vector_list[i]);
  }
  if (n <= 0) {
    printf("Invalid input. n must be greater than 0.\n");
    return false;
  } else {
    const int dim = vector_dim(vector_list[0]);
    if (dim < 1) {
      printf("Invalid input. At least one vector has a dimension of less ");
      printf("than 1.\n");
      return false;
    }
    for (int i = 1; i < n; i++) {
      if (vector_dim(vector_list[i]) != dim) {
        printf("Invalid input. All vectors must have same dimension.\n");
        return false;
      }
    }
    return true;
  }
}


bool linearly_independent(const struct vector * const vector_list[], 
                          const int n) {
  if (vector_list_valid(vector_list, n)) {
    struct matrix *vectors = matrix_create();
    for (int i = 0; i < n; i++) {
      matrix_add_col(vectors, vector_list[i]);
    }
    const int rank = matrix_rank(vectors);
    matrix_destroy(vectors);
    if (rank == n) {
      return true;
    }
  }
  return false;
}


bool is_basis(const struct vector * const vector_list[], const int dim,
              const int n) {
  if ((linearly_independent(vector_list, n)) && (n == dim)) {
    return true;
  }
  return false;
}


bool in_span(const struct vector * const vector_list[], const int n, 
             const struct vector * const v1) {
  assert(v1);
  if (vector_list_valid(vector_list, n)) {
    if (vector_dim(v1) != vector_dim(vector_list[0])) {
      printf("Invalid input. v1 must have same dimension has vectors in");
      printf("the list of vectors.\n");
      return false;
    }
    struct matrix *vectors = matrix_create();
    for (int i = 0; i < n; i++) {
      matrix_add_col(vectors, vector_list[i]);
    }
    const int coef_rank = matrix_rank(vectors);
    matrix_add_col(vectors, v1);
    const int augm_rank = matrix_rank(vectors);
    matrix_destroy(vectors);
    if (coef_rank == augm_rank) {
      return true;
    } 
  }
  return false;
}

struct matrix *find_basis(const struct vector * const vector_list[], 
                          const int n) {
  if (vector_list_valid(vector_list, n)) {
    struct matrix *basis = matrix_create();
    int basis_dim = 0;
    for (int i = 0; i < n; i++) {
      matrix_add_col(basis, vector_list[i]);
      if (matrix_rank(basis) != basis_dim + 1) {
        matrix_del_col(basis, basis_dim + 1);
      } else {
        basis_dim++;
      }
    }
    return basis;
  }
  return NULL;
}

struct vector *B_coord(const struct vector * const basis[], const int n,
                       const struct vector * const v1) {
  assert(v1);
  if (!linearly_independent(basis, n)) {
    printf("Invalid input. the first n vectors in basis are not linearly");
    printf("independent.\n");
  } else if (!in_span(basis, n, v1)) {
    printf("Invalid input. Vector is not in span.\n");
  } else {
    struct matrix *vectors = matrix_create();
    for (int i = 0; i < n; i++) {
      matrix_add_col(vectors, basis[i]);
    }
    matrix_add_col(vectors, v1);
    struct matrix *result = RREF(vectors);
    int height = vector_dim(basis[0]);
    while (height > n) {
      matrix_del_row(result, height);
      height --;
    }
    struct vector *Bcoords = matrix_dupe_col(result, n + 1);
    matrix_destroy(vectors);
    matrix_destroy(result);
    return Bcoords;
  }
  return NULL;
}


struct matrix *change_of_coord_matrix(const struct vector * const B1[], 
                                      const struct vector * const B2[],
                                      const int n) {
  if (vector_list_valid(B1,n) && vector_list_valid(B2, n)) {
    if (!(linearly_independent(B1, n) && linearly_independent(B2, n))) {
      printf("Invalid input. At least one of the vector sets is not a");
      printf("basis.\n");
    } else if (vector_dim(B1[0]) != vector_dim(B2[0])) {
      printf("Invalid input. The vector sets must be of same dimension.\n");
    } else {
      for (int i = 0; i < n; i++) {
        if (!in_span(B1, n, B2[i])) {
          printf("Invalid input. The vector sets are not basis of the same");
          printf("vector space.\n");
          return NULL;
        } 
      }
      struct matrix *result = matrix_create();
      for (int i = 0; i < n; i++) {
        struct vector *temp = B_coord(B1, n, B2[i]);
        matrix_add_col(result, temp);
        vector_destroy(temp);
      }
      return result;
    }
  }
  return NULL;
}







