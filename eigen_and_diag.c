#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include "eigen_and_diag.h"
#include "matrix_operations.h"
#include "matrix_core.h"
#include "vector_core.h"
#include "inv_and_det.h"
#include "settings.h"


//See header file for documentation

struct matrix *B_matrix(const struct vector * const B[], 
                        const struct matrix * const L) {
  assert(L);
  int m, n = -1;
  matrix_size(L, &m, &n);
  for (int i = 0; i < n; i++) {
    assert(B[i]);
    if (vector_dim(B[i]) != n) {
      printf("Invalid input. Height of vectors in B must be the same as the");
      printf("matrix L.\n");
      return NULL;
    }
  }
  if ((m != n) || (m < 1)) {
    printf("Invalid input. [L] must be an n x n matrix where n > 0.\n");
    return NULL;
  } else {
    struct matrix *temp1 = matrix_create();
    for (int i = 0; i < n; i++) {
      matrix_add_col(temp1, B[i]);
    }
    struct matrix *temp1_inv = matrix_inverse(temp1);
    struct matrix *temp2 = matrix_mult_matrix(L, temp1);
    struct matrix *Bmatrix = matrix_mult_matrix(temp1_inv, temp2);
    matrix_destroy(temp1);
    matrix_destroy(temp1_inv);
    matrix_destroy(temp2);
    return Bmatrix;
  }
  return NULL;
}


void eigenvalue_2x2(const struct matrix * const A, long double * const lambda1,
                    long double * const lambda2) {
  assert(A);
  int m, n = 0;
  matrix_size(A, &m, &n);
  if (m != 2 || n != 2) {
    printf("Invalid input. The matrix is not 2 x 2.\n");
    *lambda1 = INT_MIN;
    *lambda2 = INT_MIN;
    return;
  } else {
    const long double a = matrix_elem(A, 1, 1);
    const long double b = matrix_elem(A, 1, 2);
    const long double c = matrix_elem(A, 2, 1);
    const long double d = matrix_elem(A, 2, 2);
    const long double power1 = -a - d;
    const long double constant = a * d - b * c;
    if ((power1 * power1 - 4 * constant) < 0) {
      printf("The matrix has non-real eigenvalues.\n");
      printf("The function currently supports only real eigenvalues.\n");
      *lambda1 = INT_MIN;
      *lambda2 = INT_MIN;
      return;
    } else {
      const long double root1 = 
        (-power1 + sqrt(power1 * power1 - 4 * constant)) / 2;
      const long double root2 = 
        (-power1 - sqrt(power1 * power1 - 4 * constant)) / 2;
      *lambda1 = root1;
      *lambda2 = root2;
    }
  }
}

//cubic_roots(a, b, c, x0, x1, x2) takes in three long doubles as the 2nd 
//   degree, 1st degree and constant coefficients of a degree 3 polynomial, 
//   where the 3rd degree coefficient is 1. Then it updates *x0, *x1 and *x2
//   to reflect the real roots of the polynomial if they exist, and return the
//   number of real roots.
//requires: x0, x1, x2 are not NULL.
//effects: may modify *x0, *x1, and *x2.

//IMPORTANT NOTE: this function is adapted from the gsl_poly_solve_cubic 
//   function from the GNU Scientific Library for C. The changes made were: 
//   altered all doubles into long doubles; altered a function call to fabs() 
//   into fabsl() to match long double type; changed code format to make it 
//   more compact. This function is copyrighted under the GNU General Public 
//   License v3.0.

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)
#define M_PI 3.14159265358979323846264338327950288

static int cubic_roots(long double a, long double b, long double c, 
                       long double *x0, long double *x1, long double *x2) {
  long double q = (a * a - 3 * b);
  long double r = (2 * a * a * a - 9 * a * b + 27 * c);
  long double Q = q / 9;
  long double R = r / 54;
  long double Q3 = Q * Q * Q;
  long double R2 = R * R;
  long double CR2 = 729 * r * r;
  long double CQ3 = 2916 * q * q * q;
  if (R == 0 && Q == 0) {
    *x0 = - a / 3 ;
    *x1 = - a / 3 ;
    *x2 = - a / 3 ;
    return 3;
  } else if (CR2 == CQ3) {
    long double sqrtQ = sqrt (Q);
    if (R > 0) {
      *x0 = -2 * sqrtQ  - a / 3;
      *x1 = sqrtQ - a / 3;
      *x2 = sqrtQ - a / 3;
    } else {
      *x0 = - sqrtQ  - a / 3;
      *x1 = - sqrtQ - a / 3;
      *x2 = 2 * sqrtQ - a / 3;
    }
    return 3 ;
  } else if (R2 < Q3) {
    long double sgnR = (R >= 0 ? 1 : -1);
    long double ratio = sgnR * sqrt (R2 / Q3);
    long double theta = acos (ratio);
    long double norm = -2 * sqrt (Q);
    *x0 = norm * cos (theta / 3) - a / 3;
    *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
    *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
    if (*x0 > *x1) SWAP(*x0, *x1) ;
    if (*x1 > *x2) {
      SWAP(*x1, *x2) ;
      if (*x0 > *x1) SWAP(*x0, *x1) ;
    }
    return 3;
  } else {
    long double sgnR = (R >= 0 ? 1 : -1);
    long double A = -sgnR * pow (fabsl (R) + sqrt (R2 - Q3), 1.0/3.0);
    long double B = Q / A ;
    *x0 = A + B - a / 3;
    return 1;
  }
}


void eigenvalue_3x3(const struct matrix * const A, long double * const lambda1,
                    long double * const lambda2, long double * const lambda3) {
  assert(A);
  int m, n = 0;
  matrix_size(A, &m, &n);
  if (m != 3 || n != 3) {
    printf("Invalid input. The matrix is not 3 x 3.\n");
    *lambda1 = INT_MIN;
    *lambda2 = INT_MIN;
    *lambda3 = INT_MIN;
    return;
  } else {
    const long double a = matrix_elem(A, 1, 1);
    const long double b = matrix_elem(A, 1, 2);
    const long double c = matrix_elem(A, 1, 3);
    const long double d = matrix_elem(A, 2, 1);
    const long double e = matrix_elem(A, 2, 2);
    const long double f = matrix_elem(A, 2, 3);
    const long double g = matrix_elem(A, 3, 1);
    const long double h = matrix_elem(A, 3, 2);
    const long double i = matrix_elem(A, 3, 3);
    const long double power2 = a + e + i;
    const long double power1 = -(a * e) - (a * i) + (c * g) + (b * d) - (e * i)
      + (h * f);
    const long double constant = (a * e * i) - (a * f * h) - (b * d * i) +
      (b * f * g ) - (c * e * g ) + (c * d * h);
    long double eigenvalue1, eigenvalue2, eigenvalue3 = INT_MIN;
    int root_number = cubic_roots(-power2, -power1, -constant, &eigenvalue1,
                                  &eigenvalue2, &eigenvalue3);
    if (root_number != 3) {
      printf("The matrix has non-real eigenvalues.\n");
      printf("The function currently supports only real eigenvalues.\n");
      *lambda1 = INT_MIN;
      *lambda2 = INT_MIN;
      *lambda3 = INT_MIN;
      return;
    } else {
      *lambda1 = eigenvalue1;
      *lambda2 = eigenvalue2;
      *lambda3 = eigenvalue3;
    }
  }
}


int eigenvectors_2x2(const struct matrix * const A, struct vector ** const v1,
                     struct vector ** const v2) {
  long double lambda1, lambda2 = 0;
  eigenvalue_2x2(A, &lambda1, &lambda2);
  if (lambda1 == INT_MIN && lambda2 == INT_MIN) {
    return 0;
  } else {
    long double entries[] = 
    {matrix_elem(A, 1, 1) - lambda1, matrix_elem(A, 1, 2),
     matrix_elem(A, 2, 1), matrix_elem(A, 2, 2) - lambda1};           
    struct matrix *temp = quick_matrix_input(entries, 2, 2);
    struct matrix *rref = RREF(temp);
    matrix_destroy(temp);
    if ((-PRECISION < matrix_elem(rref, 1, 1)) &&
        (matrix_elem(rref, 1, 1) < PRECISION) &&
        (-PRECISION < matrix_elem(rref, 1, 2)) &&
        (matrix_elem(rref, 1, 2) < PRECISION)) {
      matrix_destroy(rref);
      long double v_1[] = {1, 0};
      long double v_2[] = {0, 1};
      *v1 = quick_vector_input(v_1, 2);
      *v2 = quick_vector_input(v_2, 2);
      return 2;
    } else if ((-PRECISION < matrix_elem(rref, 1, 1)) &&
               (matrix_elem(rref, 1, 1) < PRECISION)) {
      //algebraic of 2, geometric of 1
      matrix_destroy(rref);
      long double v_1[] = {1, 0};
      *v1= quick_vector_input(v_1, 2);
      return 1;
    } else {
      long double v_1[] = 
      {-matrix_elem(rref, 1, 2) / matrix_elem(rref, 1, 1), 1};
      *v1 = quick_vector_input(v_1, 2);
      matrix_destroy(rref);
      long double entries[] = 
      {matrix_elem(A, 1, 1) - lambda2, matrix_elem(A, 1, 2),
       matrix_elem(A, 2, 1), matrix_elem(A, 2, 2) - lambda2};
      struct matrix *temp = quick_matrix_input(entries, 2, 2);
      struct matrix *rref = RREF(temp);
      matrix_destroy(temp);
      long double v_2[] = 
      {-matrix_elem(rref, 1, 2) / matrix_elem(rref, 1, 1), 1};
      *v2 = quick_vector_input(v_2, 2);
      matrix_destroy(rref);
      return 2;
    }
  }
}

int eigenvectors_3x3(const struct matrix * const A, struct vector ** const v1,
                     struct vector ** const v2, struct vector ** const v3) {
  long double lambda1, lambda2, lambda3 = 0;
  eigenvalue_3x3(A, &lambda1, &lambda2, &lambda3);
  if (!(lambda1 == INT_MIN && lambda2 == INT_MIN && lambda3 == INT_MIN)) {
    const long double a = matrix_elem(A, 1, 1);
    const long double b = matrix_elem(A, 1, 2);
    const long double c = matrix_elem(A, 1, 3);
    const long double d = matrix_elem(A, 2, 1);
    const long double e = matrix_elem(A, 2, 2);
    const long double f = matrix_elem(A, 2, 3);
    const long double g = matrix_elem(A, 3, 1);
    const long double h = matrix_elem(A, 3, 2);
    const long double i = matrix_elem(A, 3, 3);
    long double a1, b1, c1, e1, f1 = 0;
    //making sure any lambda with multiplicity of 2 stay in the front:
    if (lambda1 == lambda3) {
      long double temp = lambda2;
      lambda2 = lambda3;
      lambda3 = temp;
    }
    if (lambda2 == lambda3) {
      long double temp = lambda3;
      lambda3 = lambda1;
      lambda1 = temp;
    }
    printf("%Lf, %Lf, %Lf\n", lambda1,lambda2, lambda3);
    int first_part_total = 0;
    long double entries[] = {a - lambda1, b, c, d, e - lambda1, f, g, h,
                             i - lambda1};
    struct matrix *temp = quick_matrix_input(entries, 3, 3);
    struct matrix *rref = RREF(temp);
    int rref_rank = matrix_rank(rref);
    a1 = matrix_elem(rref, 1, 1);
    b1 = matrix_elem(rref, 1, 2);
    c1 = matrix_elem(rref, 1, 3);
    e1 = matrix_elem(rref, 2, 2);
    f1 = matrix_elem(rref, 2, 3);
    matrix_destroy(temp);
    matrix_destroy(rref);
    if (rref_rank == 0) {// 3 identical eigenvalues
      long double v_1[] = {1, 0, 0};
      long double v_2[] = {0, 1, 0};
      long double v_3[] = {0, 0, 1};
      *v1 = quick_vector_input(v_1, 3);
      *v2 = quick_vector_input(v_2, 3);
      *v3 = quick_vector_input(v_3, 3);
      //diagonalizable, it has g(lambda) = a(lambda) = 3
      return 3;
    } else if (rref_rank == 1) {//may be diagonalizable
      if ((a1 > PRECISION) || (a1 < -PRECISION)) {
        long double v_1[] = {-b1 / a1, 1, 0};
        long double v_2[] = {-c1 / a1, 0, 1};
        *v1 = quick_vector_input(v_1, 3);
        *v2 = quick_vector_input(v_2, 3);
      } else if ((b1 > PRECISION) || (b1 < -PRECISION)) {
        long double v_1[] = {1, -a1 / b1, 0};
        long double v_2[] = {0, -c1 / b1, 1};
        *v1 = quick_vector_input(v_1, 3);
        *v2 = quick_vector_input(v_2, 3);
      } else {
        long double v_1[] = {1, 0, 0};
        long double v_2[] = {0, 1, 0};
        *v1 = quick_vector_input(v_1, 3);
        *v2 = quick_vector_input(v_2, 3);
      }
      first_part_total = 2;
      if ((lambda1 == lambda2) && (lambda2 == lambda3)) {
        //in this subcase, not diagonalizable, also no more eigenvectors
        return 2;
      }
    } else if (rref_rank == 2) {//not diagonalizable
      if (a1 > PRECISION || a1 < -PRECISION) {
        if (e1 > PRECISION || e1 < -PRECISION) {
          long double v_1[] = {-c1 / a1, -f1 / e1, 1};
          *v1 = quick_vector_input(v_1, 3);
        } else {
          long double v_1[] = {-b1 / a1, 1, 0};
          *v1 = quick_vector_input(v_1, 3);
        }
      } else {
        long double v_1[] = {1, 0, 0};
        *v1 = quick_vector_input(v_1, 3);
      }
      first_part_total = 1;
      if (lambda2 == lambda3) {//in this subcase, no more eigenvectors
        return 1;
      }
    }
    //at this point we must have lambda1 == lambda2 != lambda3, or
    //   lambda1 != lambda2 != lambda3. Otherwise it would have returned
    //   previously.
    long double entries2[] = {a - lambda3, b, c, d, e - lambda3, f, g, h,
                              i - lambda3};
    temp = quick_matrix_input(entries2, 3, 3);
    rref = RREF(temp);
    a1 = matrix_elem(rref, 1, 1);
    b1 = matrix_elem(rref, 1, 2);
    c1 = matrix_elem(rref, 1, 3);
    e1 = matrix_elem(rref, 2, 2);
    f1 = matrix_elem(rref, 2, 3);
    matrix_destroy(rref);
    matrix_destroy(temp);
    //rref must have rank 2, since alg/geom multiplicity of lambda3 must be 1
    if (a1 > PRECISION || a1 < -PRECISION) {
      if (e1 > PRECISION || e1 < -PRECISION) {
        long double v_3[] = {-c1 / a1, -f1 / e1, 1};
        *v3 = quick_vector_input(v_3, 3);
      } else {
        long double v_3[] = {-b1 / a1, 1, 0};
        *v3 = quick_vector_input(v_3, 3);
      }
    } else {
      long double v_3[] = {1, 0, 0};
      *v3 = quick_vector_input(v_3, 3);
    }
    if (lambda1 == lambda2) {//in this subcase, no more eigenvectors
      return first_part_total + 1;
    }
    long double entries3[] = {a - lambda2, b, c, d, e - lambda2, f, g, h,
                              i - lambda2};
    temp = quick_matrix_input(entries3, 3, 3);
    rref = RREF(temp);
    a1 = matrix_elem(rref, 1, 1);
    b1 = matrix_elem(rref, 1, 2);
    c1 = matrix_elem(rref, 1, 3);
    e1 = matrix_elem(rref, 2, 2);
    f1 = matrix_elem(rref, 2, 3);
    matrix_destroy(rref);
    matrix_destroy(temp);
    if (a1 > PRECISION || a1 < -PRECISION) {
      if (e1 > PRECISION || e1 < -PRECISION) {
        long double v_2[] = {-c1 / a1, -f1 / e1, 1};
        *v2 = quick_vector_input(v_2, 3);
      } else {
        long double v_2[] = {-b1 / a1, 1, 0};
        *v2 = quick_vector_input(v_2, 3);
      }
    } else {
      long double v_2[] = {1, 0, 0};
      *v2 = quick_vector_input(v_2, 3);
    }
    return 3;
  }
  return 0;
}


void diagonalize_2x2(const struct matrix * const A, struct matrix ** const P,
                     struct matrix ** const D, struct matrix ** const P_inv) {
  struct vector *v1, *v2 = NULL;
  int eigenvector_num = eigenvectors_2x2(A, &v1, &v2);
  if (eigenvector_num == 1) {
    printf("The matrix is not diagonalizable.\n");
  } else if (eigenvector_num == 2) {
    long double lambda1, lambda2 = 0;
    eigenvalue_2x2(A, &lambda1, &lambda2);
    struct matrix *p = matrix_create();
    matrix_add_col(p, v1);
    matrix_add_col(p, v2);
    *P = p;
    struct matrix *p_inv = matrix_inverse(p);
    *P_inv = p_inv;
    vector_destroy(v1);
    vector_destroy(v2);
    long double entries[] = {lambda1, 0, 0, lambda2};
    struct matrix *d = quick_matrix_input(entries, 2, 2);
    *D = d;
  }
}



void diagonalize_3x3(const struct matrix * const A, struct matrix ** const P,
                     struct matrix ** const D, struct matrix ** const P_inv) {
  struct vector *v1, *v2, *v3 = NULL;
  int eigenvector_num = eigenvectors_3x3(A, &v1, &v2, &v3);
  if ((eigenvector_num == 1) || (eigenvector_num == 2)) {
    printf("The matrix is not diagonalizable.\n");
  } else if (eigenvector_num == 3) {
    long double lambda1, lambda2, lambda3 = 0;
    eigenvalue_3x3(A, &lambda1, &lambda2, &lambda3);
    struct matrix *p = matrix_create();
    matrix_add_col(p, v1);
    matrix_add_col(p, v2);
    matrix_add_col(p, v3);
    *P = p;
    struct matrix *p_inv = matrix_inverse(p);
    *P_inv = p_inv;
    vector_destroy(v1);
    vector_destroy(v2);
    vector_destroy(v3);
    long double entries[] = {lambda1, 0, 0, 0, lambda2, 0, 0, 0, lambda3};
    struct matrix *d = quick_matrix_input(entries, 3, 3);
    *D = d;
  }
}


