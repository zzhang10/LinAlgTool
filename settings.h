//Long doubles are very precise, but calculation errors still build up over
//   time. When checking critical values (such as leading ones for rank), 
//   functions use a given precision. For example, a number within 
//   (1-PRECISION, 1+PRECISION) will be treated as 1 when calculating rank.
extern const long double PRECISION;


//The following parameters control the brackets of vectors and matrices. For
//   example, you may want to print a vector with different side brackets as 
//   a matrix with one row. 
extern const char VECTOR_BRACKET_LEFT;
extern const char VECTOR_BRACKET_RIGHT;
extern const char MATRIX_BRACKET_LEFT;
extern const char MATRIX_BRACKET_RIGHT;


