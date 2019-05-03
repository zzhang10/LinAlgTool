#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector_core.h"
#include "vector_operations.h"
#include "settings.h"

//see header file for documentation



struct matrix {
  int width;
  int height;
  int maxheight;
  struct vector **rows;
};


struct matrix *matrix_create() {
  struct matrix *current = malloc(sizeof(struct matrix));
  current->width = 0;
  //no maxwidth parameter since vector rows dynamically realloc themselves
  current->height = 0;
  current->maxheight = 1;
  current->rows = malloc(sizeof(struct vector *));
  return current;
}


void matrix_size(const struct matrix * const A, int * const m, int * const n) {
  assert(A);
  assert(n);
  assert(m);
  *m = A->height;
  *n = A->width;
}


void matrix_add_row(struct matrix * const A, const struct vector * const v1) {
  assert(A);
  assert(v1);
  if(A->height == 0) {
    struct vector *temp = vector_dupe(v1);
    A->rows[0] = temp;
    A->height ++;
    A->width = vector_dim(v1);
    return;
  } else if (A->width != vector_dim(v1)) {
    printf("A vector with %d elements cannot be added as a ", vector_dim(v1));
    printf("row of a matrix with %d columns.\n", A->width);
    return;
  } else if (A->height == A->maxheight) {
    A->maxheight *= 2;
    A->rows = realloc(A->rows, A->maxheight * sizeof(struct vector *));
  }
  struct vector *temp = vector_dupe(v1);
  A->rows[A->height] = temp;
  A->height ++;
}



void matrix_replace_row(struct matrix * const A, const int index, 
                        const struct vector * const v1) {
  assert(A);
  assert(v1);
  if (index <= 0 || index > A->height) {
    printf("Row %d does not exist in a matrix with %d rows.\n", index,
           A->height);
  } else if (A->width != vector_dim(v1)) {
    printf("A vector with %d elements cannot be a ", vector_dim(v1));
    printf("replacement of a row in a matrix with %d columns.\n", A->width);
    return;
  } else {
    struct vector *temp = A->rows[index - 1];
    struct vector *temp2 = vector_dupe(v1);
    A->rows[index - 1] = temp2;
    vector_destroy(temp);
  }
}



struct vector *matrix_dupe_row(const struct matrix * const A, 
                               const int index) {
  assert(A);
  if (index <= 0 || index > A->height) {
    printf("Row %d does not exist in a matrix with %d rows.\n", index,
           A->height);
    return NULL;
  } else {
    struct vector *dupe = vector_dupe(A->rows[index - 1]);
    return dupe;
  }
}

void matrix_del_row(struct matrix * const A, const int m) {
  assert(A);
  if (A->height == 0) {
    printf("The matrix has no rows to remove.\n");
  } else if ((m > A->height) || (m <= 0)) {
    printf("Row %d cannot be found in a matrix with %d rows.\n", m,
           A->height);
  } else {
    struct vector *temp = A->rows[m - 1];
    for (int i = m - 1; i < A->height - 1; i++) {
      A->rows[i] = A->rows[i + 1];
    }
    A->height --;
    vector_destroy(temp);
  }
}

void matrix_swap_row(struct matrix * const A, const int r1, const int r2) {
  assert(A);
  if (r1 <= 0 || r1 > A->height || r2 <= 0 || r2 > A->height) {
    printf("Rows %d and %d cannot both be found in a matrix with %d rows.\n",
           r1, r2, A->height);
    return;
  } else {
    struct vector *temp = A->rows[r1 - 1];
    A->rows[r1 - 1] = A->rows[r2 - 1];
    A->rows[r2 - 1] = temp;
  }
}


void matrix_sum_row(struct matrix * const A, const int r1, const int r2) {
  assert(A);
  if (r1 <= 0 || r1 > A->height || r2 <= 0 || r2 > A->height) {
    printf("Rows %d and %d cannot both be found in a matrix with %d rows.\n",
           r1, r2, A->height);
    return;
  } else {
    struct vector *sum = vector_add(A->rows[r1 - 1], A->rows[r2 - 1]);
    struct vector *temp = A->rows[r1 - 1];
    A->rows[r1 - 1] = sum;
    vector_destroy(temp);
  }
}


void matrix_mult_row(struct matrix * const A, const int r1,
                     const long double c) {
  assert(A);
  if (r1 <= 0 || r1 > A->height) {
    printf("Rows %d cannot be found in a matrix with %d rows.\n",
           r1, A->height);
    return;
  } else {
    struct vector *product = vector_mult(A->rows[r1 - 1], c);
    struct vector *temp = A->rows[r1 - 1];
    A->rows[r1 - 1] = product;
    vector_destroy(temp);
  }
}

void matrix_add_mult_row(struct matrix * const A, const int r1, const int r2, 
                         const long double c) {
  assert(A);
  if (r1 <= 0 || r1 > A->height || r2 <= 0 || r2 > A->height) {
    printf("Rows %d and %d cannot be found in a matrix with %d rows.\n",
           r1, r2, A->height);
    return;
  } else {
    struct vector *product = vector_mult(A->rows[r2 - 1], c);
    struct vector *sum_product = vector_add(A->rows[r1 - 1], product);
    struct vector *temp = A->rows[r1 - 1];
    A->rows[r1 - 1] = sum_product;
    vector_destroy(temp);
    vector_destroy(product);
  }
}


struct matrix *quick_matrix_input(const long double values[], const int m,
                                  const int n) {
  assert(values);
  if ((m < 0) || (n < 0)) {
    printf("A matrix cannot have negative width or height.\n");
    return NULL;
  } else {
    int entry = 0;
    struct matrix *current = matrix_create();
    for (int i = 0; i < m; i++) {
      struct vector *temp = vector_create();
      for (int j = 0; j < n; j++) {
        vector_add_elem(temp, values[entry]);
        entry++;
      }
      matrix_add_row(current, temp);
      vector_destroy(temp);
    }
    return current;
  }
}


void matrix_add_col(struct matrix * const A, const struct vector * const v1) {

  assert(A);
  assert(v1);
  if(A->width == 0) {
    for (int i = 1; i <= vector_dim(v1); i++) {
      struct vector *temp = vector_create();
      vector_add_elem(temp, vector_elem(v1, i));
      matrix_add_row(A, temp);
      vector_destroy(temp);
    }
    A->height = vector_dim(v1);
  } else if (A->height != vector_dim(v1)) {
    printf("A vector with %d elements cannot be added as a ", vector_dim(v1));
    printf("column of a matrix with %d rows.\n", A->height);
  } else {
    for (int i = 0; i < A->height; i++) {
      vector_add_elem(A->rows[i], vector_elem(v1,i + 1));
    }
    A->width ++;
  }
}



void matrix_replace_col(struct matrix * const A, const int index,
                        const struct vector * const v1) { 
  assert(A);
  assert(v1);
  if (index <= 0 || index > A->width) {
    printf("Column %d does not exist in a matrix with %d columns.\n", index,
           A->width);
  } else if (A->height != vector_dim(v1)) {
    printf("A vector with %d elements cannot be a ", vector_dim(v1));
    printf("replacement of a column in a matrix with %d rows.\n", A->height);
    return;
  } else {
    for (int i = 0; i < A->height; i++) {
      vector_replace(A->rows[i], index, vector_elem(v1,i + 1));
    }
  }
}




struct vector *matrix_dupe_col(const struct matrix * const A, 
                               const int index) {
  assert(A);
  if (index <= 0 || index > A->width) {
    printf("Column %d does not exist in a matrix with %d rows.\n", index,
           A->height);
    return NULL;
  } else {
    struct vector *dupe = vector_create();
    for (int i = 0; i < A->height; i++) {
      vector_add_elem(dupe, vector_elem((A->rows)[i], index));
    }    
    return dupe;
  }
}

void matrix_del_col(struct matrix * const A, const int n) {
  assert(A);
  if (A->width == 0) {
    printf("The matrix has no columns to remove.\n");
  } else if ((n > A->width) || (n <= 0)) {
    printf("Columns %d cannot be found in a matrix with %d columns.\n",n,
           A->width);
  } else {
    for (int i = n; i < A->width; i++) {
      struct vector *temp = matrix_dupe_col(A, i + 1);
      matrix_replace_col(A, i, temp);
      vector_destroy(temp);
    }
    A->width --;
    for (int i = 0; i < A->height; i++) {
      vector_remove(A->rows[i]);
    }
  }
}

struct matrix *matrix_dupe(const struct matrix * const A) {
  assert(A);
  struct matrix *result = matrix_create();
  for (int i = 0; i < A->height; i++) {
    struct vector *row = vector_dupe(A->rows[i]);
    matrix_add_row(result, row);
    vector_destroy(row);
  }
  return result;
}

long double matrix_elem(const struct matrix * const A, const int m, 
                        const int n) {
  assert(A);
  int m1, n1 = 0;
  matrix_size(A, &m1, &n1);
  if (m <= 0 || n <= 0 || m > m1 || n > n1) {
    printf("Entry %d, %d does not exist in a %d by %d matrix.\n",
           m, n, m1, n1);
    return INT_MIN;
  } else {
    return vector_elem(A->rows[m - 1], n);
  }
}

void matrix_print(const struct matrix * const A) {
  if (!A) {
    return;
  }
  if (A->height == 0) {
    printf("[Empty]\n");
  } else {
    for (int i = 0; i < A->height; i++) {
      vector_print_bracket(A->rows[i], MATRIX_BRACKET_LEFT, 
                           MATRIX_BRACKET_RIGHT);
    }
  }
  printf("\n");
}

void matrix_destroy(struct matrix * const A) {
  if (!A) {
    return;
  } else {
    for (int i = A->height - 1; i >= 0; i--) {
      vector_destroy(A->rows[i]);
    }
    free(A->rows);
    free(A);
  }
}

