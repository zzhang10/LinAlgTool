#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_operations.h"
#include "vector_core.h"

//see header file for documentation


//valid_vectors(v1, v2) returns true if v1 and v2 are non-empty and of the
//   same length. It returns false and prints error pessage if not.
//requires: v1, v2 are not NULL;
//effects: may print output.
static bool valid_vectors(const struct vector * const v1,
                          const struct vector * const v2) {
  assert(v1);
  assert(v2);
  if ((vector_dim(v1) == 0) || (vector_dim(v2) == 0)) {
    printf("Invalid input. A vector is empty.\n");
    return false;
  } else if (vector_dim(v1) != vector_dim(v2)) {
    printf("Invalid input. The vectors are not of the same dimension.\n");
    return false;
  }
  return true;
}

struct vector *vector_mult(const struct vector * const v1,
                           const long double c) {
  assert(v1);
  if (vector_dim(v1) == 0) {
    printf("Invalid input. The vector is empty.\n");
    return NULL;
  } else {
    struct vector *new = vector_create();
    for (int i = 1; i <= vector_dim(v1); i++) {
      vector_add_elem(new, c * vector_elem(v1, i));
    }
    return new;
  }
}

struct vector *vector_add(const struct vector * const v1,
                          const struct vector * const v2) {
  if (!valid_vectors(v1,v2)) {
    return NULL;
  } else {
    struct vector *new = vector_create();
    for (int i = 1; i <= vector_dim(v1); i++) {
      vector_add_elem(new, vector_elem(v1, i) + vector_elem(v2, i));
    }
    return new;
  }
}

struct vector *vector_cross(const struct vector * const v1,
                            const struct vector * const v2) {
  assert(v1);
  assert(v2);
  if ((vector_dim(v1) != 3) || (vector_dim(v2) != 3)) {
    printf("Invalid input. Cross product in this course is defined for");
    printf("vectors in R(3) only.");
    return NULL;
  } else {
    struct vector *new = vector_create();
    long double x1 = ((vector_elem(v1, 2) * vector_elem(v2, 3)) -
                      (vector_elem(v1, 3) * vector_elem(v2, 2)));
    long double x2 = ((vector_elem(v1, 3) * vector_elem(v2, 1)) -
                      (vector_elem(v1, 1) * vector_elem(v2, 3)));
    long double x3 = ((vector_elem(v1, 1) * vector_elem(v2, 2)) -
                      (vector_elem(v1, 2) * vector_elem(v2, 1)));
    vector_add_elem(new, x1);
    vector_add_elem(new, x2);
    vector_add_elem(new, x3);
    return new;
  }
}

long double vector_dot(const struct vector * const v1,
                       const struct vector * const v2) {
  if (!valid_vectors(v1,v2)) {
    return INT_MIN;
  } else {
    long double total = 0;
    for (int i = 1; i <= vector_dim(v1); i++) {
      total += vector_elem(v1, i) * vector_elem(v2, i);
    }
    return total;
  }
}

struct vector *vector_proj(const struct vector * const v1,
                           const struct vector * const v2) {
  if (!valid_vectors(v1,v2)) {
    return NULL;
  } else {
    const long double factor = vector_dot(v1,v2) / vector_dot(v1,v1);
    return vector_mult(v1, factor);
  }
}


struct vector *vector_perp(const struct vector * const v1, 
                           const struct vector * const v2) {
  if (!valid_vectors(v1,v2)) {
    return NULL;
  } else {
    struct vector *temp1 = vector_proj(v1, v2);
    struct vector *temp2 = vector_mult(temp1, -1);
    struct vector *new = vector_add(v2, temp2);
    vector_destroy(temp1);
    vector_destroy(temp2);
    return new;
  }
}


long double vector_norm(const struct vector * const v1) {
  assert(v1);
  if (vector_dim(v1) == 0) {
    printf("Invalid input. The vector is empty.\n");
    return INT_MIN;
  }
  return pow(vector_dot(v1,v1), 0.5);
}

long double vector_angle(const struct vector * const v1,
                         const struct vector * const v2) {
  if (!valid_vectors(v1,v2)) {
    return INT_MIN;
  } else if (vector_norm(v1) * vector_norm(v2) == 0) {
    printf("At least one of the vectors is the zero vector.");
    printf("The angle is not defined.\n");
    return INT_MIN;
  } else {
    return acos(vector_dot(v1, v2) / (vector_norm(v1) * vector_norm(v2)));
  }
}

void scalar_equation(const struct vector * const v1, 
                     const struct vector * const v2, const long double x1,
                     const long double x2, const long double x3) {
  assert(v1);
  assert(v2);
  if ((vector_dim(v1) != 3) || (vector_dim(v2) != 3)) {
    printf("The current version only supports scalar equations in R(3)n\n");
    return;
  } else {
    struct vector *n = vector_cross(v1,v2);
    long double c = (vector_elem(n, 1) * x1) + (vector_elem(n, 2) * x2) +
      (vector_elem(n, 3) * x3);
    printf("The scalar equation is:\n%Lf X1 + ", vector_elem(n, 1));
    printf("%Lf X2 + %Lf X3 = %Lf\n", vector_elem(n, 2), vector_elem(n, 3), c);
    vector_destroy(n);
  }
}






