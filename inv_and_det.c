#include "vector_core.h"
#include "matrix_core.h"
#include "matrix_operations.h"
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "inv_and_det.h"
#include "settings.h"


//See header file for documentation

long double matrix_det(const struct matrix * const A) {
  assert(A);
  int m, n = 0;
  matrix_size(A, &m, &n);
  if ((m != n) || (m < 1)) {
    printf("Invalid input. Matrix must be n x n where n is positive.\n");
  } else if (m == 1) {
    return matrix_elem(A, 1, 1);
  } else {
    int det_total = 0;
    for (int i = 1; i <= m; i++) {
      struct matrix *submatrix = matrix_dupe(A);
      matrix_del_col(submatrix, 1);
      matrix_del_row(submatrix, i);
      if (i % 2 == 1) {
        det_total += matrix_elem(A, i, 1) * matrix_det(submatrix);
      } else {
        det_total -= matrix_elem(A, i, 1) * matrix_det(submatrix);
      }
      matrix_destroy(submatrix);
    }
    return det_total;
  }
  return INT_MIN;
}


long double matrix_cof(const struct matrix * const A, const int i, 
                       const int j) {
  assert(A);
  int m, n = 0;
  matrix_size(A, &m, &n);
  if ((m != n) || (m < 2)) {
    printf("Invalid input. Matrix must be n x n where n >= 2.\n");
  } else if (matrix_elem(A, i, j) != INT_MIN) {
    struct matrix *submatrix = matrix_dupe(A);
    matrix_del_row(submatrix, i);
    matrix_del_col(submatrix, j);
    const long double det = matrix_det(submatrix);
    matrix_destroy(submatrix);
    if ((i + j) % 2 == 0) {
      return det;
    } else {
      return -det;
    }
  }
  return INT_MIN;
}

struct matrix *cof_matrix(const struct matrix * const A) {
  assert(A);
  int m, n = -1;
  matrix_size(A, &m, &n);
  if ((m < 2) || (m != n)) {
    printf("Invalid input. Matrix must be n by n where n >= 2. \n");
  } else {
    struct matrix *cof = matrix_create();
    for (int i = 1; i <= n; i++) {
      struct vector *temp = vector_create();
      for (int j = 1; j <= n; j++) {
        vector_add_elem(temp, matrix_cof(A, i, j));
      }
      matrix_add_row(cof, temp);
      vector_destroy(temp);
    }
    return cof;
  }
  return NULL;
}

struct matrix *adj_matrix(const struct matrix * const A) {
  assert(A);
  struct matrix *temp = cof_matrix(A);
  if (temp) {
    struct matrix *adj = matrix_transpose(temp);
    matrix_destroy(temp);
    return adj;
  } else {
    return NULL;
  }
}

struct matrix *matrix_inverse(const struct matrix * const A) {
  assert(A);
  if (matrix_det(A) != INT_MIN) {
    if ((-PRECISION < matrix_det(A)) && (matrix_det(A) < PRECISION)) {
      printf("The matrix is not invertible.\n");
    } else {
      struct matrix *adj = adj_matrix(A);
      struct matrix *inv = matrix_mult_scalar(adj, 1/matrix_det(A));
      matrix_destroy(adj);
      return inv;
    }
  }
  return NULL;
}

