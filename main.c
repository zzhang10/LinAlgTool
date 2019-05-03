
/*/////////////////////////////////////////////////////////////////////////////

                       Linear Algebra I Toolbox 1.0

                               By Zack Zhang


/////////////////////////////////////////////////////////////////////////////*/
//Welcome to Linear Algebra Toolbox 1.0. This is a mini library written in 
//   C99 for linear algebra I. It covers most of the common functions taught
//   in class, including vector operations, matrix operations, RREF, inverses, 
//   determinants and diagonalization of up to 3 x 3 matrices. Please note
//   the program supports only real numbers.
//Before you start using the toolbox:

//Check that you have the required libraries:
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Read documentation of functions in their header files:
#include "matrix_core.h"
#include "vector_core.h"
#include "vector_space.h"
#include "vector_operations.h"
#include "matrix_operations.h"
#include "inv_and_det.h"
#include "eigen_and_diag.h"

//Modify calculation precision and printing configs in settings.c.
//Description of the parameters are in settings.h:
#include "settings.h"

///////////////////////////////////////////////////////////////////////////////

//Here's a framework to make calulations easier. You may enter numbers like
//   you write a matrix or a vector. Below are some examples that you may 
//   modify to fit your own calulations, but make sure the width, height and
//   dim match the matrices and vectors.
const long double A1_entries[] = {1, 0, 0, 0,
                                  0, 1, 0, 0,
                                  0, 0, 1, 0,
                                  0, 0, 0, 1};
const int A1_width = 4;
const int A1_height = 4;


const long double A2_entries[] = {1, -2, 0.33, 12120, 1,
                                  9, -0.229, 1, 10, -1.2};
const int A1_width = 5;
const int A1_height = 2;


const long double A3_entries[] = {2};
const int A3_width = 1;
const int A3_height = 1;


const long double v1_entries[] = {0, 1, 2, 3, 4};
const int v1_dim = 5;


const long double v2_entries[] = {1};
const int v2_dim = 1;



int main(void) {
  //Generating matrices and vectors according to input above:
  struct matrix *A1 = quick_matrix_input(A1_entries, A1_width, A1_height);
  struct matrix *A2 = quick_matrix_input(A2_entries, A2_width, A2_height);
  struct matrix *A3 = quick_matrix_input(A3_entries, A3_width, A3_height);
  struct vector *v1 = quick_vector_input(v1_entries, v1_dim);
  struct vector *v2 = quick_vector_input(v2_entries, v2_dim);
  /////////////////////////////////////////////////////////////////////////////
  //Enter your calculations below:











  /////////////////////////////////////////////////////////////////////////////
  //Matrices and vectors are dynamically allocated. Make sure to destroy them
  //   to prevent memory leaks.
  matrix_destroy(A1);
  matrix_destroy(A2);
  matrix_destroy(A3);
  vector_destroy(v1);
  vector_destroy(v2);
}
