# Linear Algebra I Toolbox

# Important: Currently optimized for Seashell environment. May not compile/link in other environments

### What tools does it contain?
You may find functions that cover most of the stuff you learn in Linear Algebra 1, such as vector and matrix operations, RREF, inverses, determinants and diagonalization (up to 3 x 3). Contact me if you would like to see anything added.

### I am just learning linear algebra and I don’t know a lot about these functions…
Linear algebra can be a little daunting at first. Try out some of the functions with different inputs and learn as you go! If you give a function invalid inputs, it will almost always tell you why the input is invalid before returning. The only exception is when you input NULL pointers that can cause segmentation faults.

### Doesn’t C accumulate errors when doing calculations with floating points?
Unfortunately, matrices are not limited to only integer entries, so we will have to deal with floating point inaccuracies. There is a “precision macro” defined in the program that helps with this problem when trying to determine whether a value is exact (like leading ones in a matrix).

#### Note: The program uses the following C libraries: assert.h, limits.h, stdbool.h, stdio.h, stdlib.h and math.h.
####       The value INT_MIN is a sentinel value. Matrices and vectors with INT_MIN as their entries may cause undefined behavior.

