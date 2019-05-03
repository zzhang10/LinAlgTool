//vector_mult(v1, c) takes in a struct vector pointer and a long double c, and 
//   returns the result of c(v1) through a vector pointer if possible (the 
//   caller must free the pointer). If not, it will output an error message 
//   and return NULL.
//requires: v1 is not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *vector_mult(const struct vector * const v1, 
                           const long double c);

//vector_add(v1, v2) takes in two struct vector pointers, and returns the 
//   result of v1 + v2 through a vector pointer if possible (the caller must
//   free the pointer).If not, it will output an error message and return NULL.
//requires: v1 and v2 are not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *vector_add(const struct vector * const v1, 
                          const struct vector * const v2);

//vector_cross(v1, v2) takes in two struct vector pointers, and returns the 
//   result of v1 x v2 through a vector pointer if possible (the caller must
//   free the pointer). If not, it will output an error message and return
//   NULL.
//requires: v1 and v2 are not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *vector_cross(const struct vector * const v1, 
                            const struct vector * const v2);

//vector_dot(v1, v2) takes in two struct vector pointers, and returns the 
//   result of v1 dot v2 if possible. If not, it will output an error message 
//   and return INT_MIN.
//requires: v1 and v2 are not NULL.
//effects: prints output.
long double vector_dot(const struct vector * const v1,
                       const struct vector * const v2);

//vector_proj(v1, v2) takes in two struct vector pointers, and returns the 
//   result of projecting v2 onto v1 through a vector pointer if possible (the
//   caller must free the pointer). If not, it will output an error message 
//   and return NULL.
//requires: v1 and v2 are not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *vector_proj(const struct vector * const v1, 
                           const struct vector * const v2);

//vector_perp(v1, v2) takes in two struct vector pointers, and returns the
//   result of taking the perpendicular of v1 onto v2 with a vector pointer
//   if possible (the caller must free the pointer). If not, it will output an
//   error message and return NULL.
//requires: v1 and v2 are not NULL.
//effects: may print output
//         may allocate heap memory
struct vector *vector_perp(const struct vector * const v1,
                           const struct vector * const v2);

//vector_norm(v1) takes in a struct vector pointer, and returns the norm (or
//   length) of the vector if possible. If the vector is empty, it will output
//   an error message and return INT_MIN.
//effects: may print output
long double vector_norm(const struct vector * const v1);

//vector_angle(v1, v2) takes in two struct vector pointers, and returns the
//   angle (in radians) between the two vectors if possible. If not, it will 
//   output an error message and return INT_MIN.
//requires: v1 and v2 are not NULL.
//effects: may print output
long double vector_angle(const struct vector * const v1, 
                         const struct vector * const v2);

//scalar_equation(v1, v2, x1, x2, x3) takes in two struct vector pointers and 
//   three long doubles representing a point on the plane. It will print the
//   scalar equation of the plane if possible, or output an error message.
//requires: v1 and v2 are not NULL.
//effects: may print output
void scalar_equation(const struct vector * const v1, 
                     const struct vector * const v2, const long double x1,
                     const long double x2, const long double x3);


