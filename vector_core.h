//A struct vector represents a vector in Euclidean space.
struct vector;

//vector_create() returns a heap allocated struct vector pointer that caller 
//   must free using vector_destroy().
//effects: allocates heap memory
struct vector *vector_create();

//quick_vector_input(values, n) takes an array of long doubles and an
//   integer n. If possible, it returns a heap-allocated pointer to a
//   vector in R[n] formed with the array of numbers; the caller must free 
//   this pointer with vector_destroy(). If vector creation is not possible,
//   it outputs an error message, and returns NULL.
//requires: values is not NULL
//          for positive n, values contains at least n long doubles.
//effects: may allocate heap memory 
//         may print message
//note: values may contain more than n long doubles, but the extraneous
//      values are disregarded.
struct vector *quick_vector_input(const long double values[], const int n);


//vector_dim(v1) returns the number of dimensions (or elements) of *v1
//requires: v1 is not NULL
int vector_dim(const struct vector * const v1);


//vector_add_elem(v1, x) takes in a vector pointer v1 and a long double x. It 
//   adds x to *v1 as the last element.
//requires: v1 is not NULL
//effects: modifies *v1
void vector_add_elem(struct vector * const v1, const long double x);


//vector_remove(v1) removes the last element of the vector in v1. It displays
//   an error message if the vector is empty.
//requires: v1 is not NULL
//effects: may modify *v1
//         may print output
void vector_remove(struct vector * const v1);

//vector_elem(v1, index) returns the index-th element of the struct vector 
//   that v1 points to. It will return INT_MIN if the index is out of range.
long double vector_elem(const struct vector * const v1, const int index);

//vector_replace(v1, index, x) takes in a struct vector pointer, an integer
//   index, and a long double x. If possible, it replaces the index-th 
//   coordinate  in the vector with x. Otherwise it displays an error message.
//requires: v1 is not NULL
//effects: may modify *v1
//         may print output
void vector_replace(const struct vector * const v1, const int index, 
                    const long double x);

//vector_dupe(v1) returns a new struct vector pointer with identical vector as
//   the one passed in.
//effects: may allocate heap memory
struct vector *vector_dupe(const struct vector * const v1);

//vector_print(v1) prints the vector that v1 points to.
//effects: prints output
void vector_print(const struct vector * const v1);

//vector_print_bracket(v1, left, right) prints the vector that v1 points to
//   with the given characters left and right as brackets
//effects: prints output
void vector_print_bracket(const struct vector * const v1, const char left, 
                          const char right);

//vector_destroy(v1) frees heap memory allocated to v1 if it is not NULL.
//effects: may free heap memory
void vector_destroy(struct vector * const v1);


