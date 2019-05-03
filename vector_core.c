#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "vector_core.h"
#include "settings.h"


//see header file for documentation

struct vector {
  int dim;
  int maxdim;
  long double *value;
};

struct vector *vector_create() {
  struct vector *current = malloc(sizeof(struct vector));
  current->dim = 0;
  current->maxdim = 1;
  current->value = malloc(sizeof(long double));
  return current;
}

struct vector *quick_vector_input(const long double values[], const int n) {
  assert(values);
  if (n < 0) {
    printf("A vector cannot have negative number of elements.\n");
    return NULL;
  } else {
    struct vector *current = vector_create();
    for (int i = 0; i < n; i++) {
      vector_add_elem(current, values[i]);
    }
    return current;
  }
}


int vector_dim(const struct vector * const v1) {
  assert(v1);
  return v1->dim;
}


void vector_add_elem(struct vector * const v1, const long double x) {
  assert(v1);
  if (v1->dim == v1->maxdim) {
    v1->maxdim *= 2;
    v1->value = realloc(v1->value, v1->maxdim * sizeof(long double));
  }
  v1->value[v1->dim] = x;
  v1->dim ++;
}


struct vector *vector_dupe(const struct vector * const v1) {
  struct vector *duped = vector_create();
  for (int i = 0; i < v1->dim; i++) {
    vector_add_elem(duped, v1->value[i]);
  }
  return duped;
}

void vector_remove(struct vector * const v1) {
  assert(v1);
  if (v1->dim == 0) {
    printf("The vector has no coordinates to remove.\n");
  } else {
    v1->dim --;
  }
}

long double vector_elem(const struct vector * const v1, const int index) {
  assert(v1);
  if ((index <= 0) || (index > v1->dim)) {
    printf("Element %d does not exist in a vector with %d elements\n",
           index, v1->dim);
    return INT_MIN;
  } else {
    return v1->value[index - 1];
  }
}

void vector_replace(const struct vector * const v1, const int index, 
                    const long double x) {
  assert(v1);
  if (index <= 0 || index > v1->dim) {
    printf("Element %d does not exist in a vector with %d elements\n",
           index, v1->dim);
  } else {
    v1->value[index - 1] = x;
  }
}

void vector_print_bracket(const struct vector * const v1, const char left, 
                          const char right) {
  if (!v1) {
    printf("The vector is currently null (uninitialized).\n");
    return;
  }
  printf("%c", left);
  for (int i = 0; i < v1->dim; i++) {
    if ((-PRECISION < v1->value[i]) && (PRECISION > v1->value[i])) {
      const long double x = 0;
      printf("%14.5Lf", x);
    } else {
      printf("%14.5Lf", v1->value[i]);
    }
    if (i != v1->dim - 1) {
      printf(" ");
    }
  }
  printf("%c\n", right);
}

void vector_print(const struct vector * const v1) {
  vector_print_bracket(v1, VECTOR_BRACKET_LEFT, VECTOR_BRACKET_RIGHT);
}


void vector_destroy(struct vector * const v1) {
  if (!v1) {
    return;
  } else {
    free(v1->value);
    free(v1);
  }
}


