#ifndef AMP_UTILITIES_H
#define AMP_UTILITIES_H

#define AMP_SUCCESS 0
#define AMP_ALLOC_MEMOEY_FAILURE 1

/* Return the AMP error string for a given error number.  */
char *amp_error_string(int error);
// 用于给AMP生成的低精度的变量申请内存空间
int ampMalloc(void **p, int size);

#endif