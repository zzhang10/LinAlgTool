#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_operations.h"
#include "vector_core.h"
#include "matrix_operations.h"
#include "matrix_core.h"
#include "settings.h"


//See header file for documentation

//valid_matrix(A) returns true if A is a valid, non-empty matrix. It returns 
//   false and prints error pessage if not.
//requires: A is not NULL;
//effects: may print output.
static bool valid_matrix(const struct matrix * const A) {
  assert(A);
  int m = 0;
  int n = 0;
  matrix_size(A, &m, &n);
  if ((m == 0) || (n == 0)) {
    printf("The input matrix must not be empty.\n");
    return false;
  }
  return true;
}



struct matrix *matrix_add(const struct matrix * const A,
                          const struct matrix * const B) {
  assert(A);
  assert(B);
  int m1, n1, m2, n2 = 0;
  matrix_size(A, &m1, &n1);
  matrix_size(B, &m2, &n2);
  if ((m1 == 0) || (n1 == 0) || (m2 == 0) || (n2 == 0)) {
    printf("The input matrix must not be empty.\n");
    return NULL;
  } else if ((m1 != m2) || (n1 != n2)) {
    printf("The input matrices must be of the same size.\n");
    return NULL;
  } else {
    struct matrix *dupe = matrix_create();
    for (int i = 1; i <= m1; i++) {
      struct vector *dupe_row_A = matrix_dupe_row(A, i);
      struct vector *dupe_row_B = matrix_dupe_row(B, i);
      struct vector *dupe_row_sum = vector_add(dupe_row_A, dupe_row_B);
      matrix_add_row(dupe, dupe_row_sum);
      vector_destroy(dupe_row_A);
      vector_destroy(dupe_row_B);
      vector_destroy(dupe_row_sum);
    }
    return dupe;
  }
}

struct matrix *matrix_mult_scalar(const struct matrix * const A,
                                  const long double c) {
  if (valid_matrix(A)) {
    struct matrix *dupe = matrix_create();
    int m, n = 0;
    matrix_size(A, &m, &n);
    for (int i = 1; i <= m; i++) {
      struct vector *dupe_row = matrix_dupe_row(A, i);
      struct vector *dupe_row_mult = vector_mult(dupe_row, c);
      matrix_add_row(dupe, dupe_row_mult);
      vector_destroy(dupe_row);
      vector_destroy(dupe_row_mult);
    }
    return dupe;
  }
  return NULL;
}


struct vector *matrix_mult_vector(const struct matrix * const A, 
                                  const struct vector * const v1) {
  if (valid_matrix(A)) {
    int m, n = 0;
    matrix_size(A, &m, &n);
    if (vector_dim(v1) != n) {
      printf("The height of the vector must match the width of the matrix.\n");
    } else {
      struct vector *result = vector_create();
      for (int i = 1; i <= m; i++) {
        struct vector *dupe_row = matrix_dupe_row(A, i);
        long double entry = vector_dot(v1, dupe_row);
        vector_add_elem(result, entry);
        vector_destroy(dupe_row);
      }
      return result;
    }
  }
  return NULL;
}


struct matrix *matrix_mult_matrix(const struct matrix * const A, 
                                  const struct matrix * const B) {
  if (valid_matrix(A) && valid_matrix(B)) {
    int m1, n1, m2, n2 = 0;
    matrix_size(A, &m1, &n1);
    matrix_size(B, &m2, &n2);
    if (m2 != n1) {
      printf("The height of the second matrix must match the width of the");
      printf(" first matrix.\n");
    } else {
      struct matrix *result = matrix_create();
     
      for (int i = 1; i <= n2; i++) {
        
        struct vector *duped_col = matrix_dupe_col(B, i);
        struct vector *duped_col_mult = matrix_mult_vector(A, duped_col);

        matrix_add_col(result, duped_col_mult);
        vector_destroy(duped_col);
        vector_destroy(duped_col_mult);
      }
      return result;
    }
  }
  return NULL;
}


struct matrix *rotation_matrix(const long double theta) {
  struct matrix *result = matrix_create();
  struct vector *row1 = vector_create();
  struct vector *row2 = vector_create();
  long double cosine = cos(theta);
  long double sine = sin(theta);
  vector_add_elem(row1, cosine);
  vector_add_elem(row1, -sine);
  vector_add_elem(row2, sine);
  vector_add_elem(row2, cosine);
  matrix_add_row(result, row1);
  matrix_add_row(result, row2);
  vector_destroy(row1);
  vector_destroy(row2);
  return result;
}

struct matrix *matrix_transpose(const struct matrix * const A) {
  if (valid_matrix(A)) {
    int m, n = 0;
    matrix_size(A, &m, &n);
    struct matrix *result = matrix_create();
    for (int i = 1; i <= m; i++) {
      struct vector *row = matrix_dupe_row(A, i);
      matrix_add_col(result, row);
      vector_destroy(row);
    }
    return result;
  }
  return NULL;
}




//is_leading(A, m, n, precision) returns true if A[mn] is a leading number, 
//   fale otherwise.
//requires: m and n are in bound of A
static bool is_leading (struct matrix * const A, const int m, const int n) {
  assert(A);
  if ((-PRECISION < matrix_elem(A, m, n)) && 
      (PRECISION > matrix_elem(A, m, n))) {
    return false;
  }
  bool leading = true;
  for (int i = 1; i < n; i ++) {
    if ((-PRECISION > matrix_elem(A, m, i)) || 
        (PRECISION < matrix_elem(A, m, i))) {
      leading = false;
    }
  }
  return leading;
}


struct matrix *RREF(const struct matrix * const A) {
  if (valid_matrix(A)) {
    struct matrix *result = matrix_dupe(A);
    int rows, cols = 0;
    matrix_size(result, &rows, &cols);
    if (rows <= 1) {
      return result;
    }
    int leading_row = 1;
    for (int i = 1; i <= cols; i++) {
      for (int j = leading_row; j <= rows; j++) {
        if (is_leading(result, j, i)) {
          matrix_mult_row(result, j, 1 / matrix_elem(result, j, i));
          for (int k = 1; k <= rows; k++) {
            if ((-PRECISION > matrix_elem(result, k, i) ||
                 PRECISION < matrix_elem(result, k, i)) && k != j) {
              matrix_add_mult_row(result, k, j, -matrix_elem(result, k, i) /
                                  matrix_elem(result, j, i));
            }
          }
          matrix_swap_row(result, leading_row, j);
          leading_row++;
          break;
        }
      }
    }
    return result;
  }
  return NULL;
}

bool is_RREF(const struct matrix * const A) {
  if (valid_matrix(A)) {
    int m, n = 0;
    matrix_size(A, &m, &n);
    struct matrix *RREF_A = RREF(A);
    for (int i = 1; i <= m; i++) {
      for (int j = 1; j <= n; j++) {
        if ((matrix_elem(A, i, j) + PRECISION < matrix_elem(RREF_A, i, j)) ||
            (matrix_elem(A, i, j) - PRECISION > matrix_elem(RREF_A, i, j))) {
          matrix_destroy(RREF_A);
          return false;
        }
      }
    }
    matrix_destroy(RREF_A);
    return true;
  }
  return false;
}

int matrix_rank(const struct matrix * const A) {
  if (valid_matrix(A)) {
    int m, n = 0;
    matrix_size(A, &m, &n);
    struct matrix *RREF_A = RREF(A);
    int rank = 0;
    for (int i = 1; i <= m; i++) {
      for (int j = 1; j <= n; j++) {
        if (is_leading(RREF_A, i, j)) {
          rank++;
          break;
        }
      }
    }
    matrix_destroy(RREF_A);
    return rank;
  }
  return INT_MIN;
}

struct matrix *matrix_power(const struct matrix * const A, const int n) {
  struct matrix *result = matrix_dupe(A);
  for(int i = 1; i < n; i++) {
    struct matrix *temp = result;
    result = matrix_mult_matrix(result, A);
    matrix_destroy(temp);
  }
  return result;
}
    









