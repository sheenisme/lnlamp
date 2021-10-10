#include <stdio.h>
#include <stdlib.h>

#include "amp_utilities.h"

/* Return the AMP error string for a given error number.
 */
char *amp_error_string(int error)
{
    int errorCount;
    int index;

    char *errorString[] = {
        [AMP_SUCCESS] = "AMP_SUCCESS",                          //0代表成功
        [AMP_ALLOC_MEMOEY_FAILURE] = "AMP_ALLOC_MEMOEY_FAILURE" //1代表申请内存失败
    };

    errorCount = sizeof(errorString) / sizeof(errorString[0]);
    index = error;

#ifdef DEBUG
    printf("error ID is %d \n", error);
#endif

    return (index >= 0 && index < errorCount) ? errorString[index] : "Unspecified Error";
}

// 用于给AMP生成的低精度的变量申请内存空间
int ampMalloc(void **p, int size)
{
#ifdef DEBUG
    printf("will malloc the size of memory is : %d \n", size);
#endif

    if ((size > 0) && (p != NULL))
    {
        *p = malloc(size);
    }
    else
    {
#ifdef DEBUG
        printf("ampMalloc meets some Errors, may the size <= 0 or the pointer is null, please check.\n");
#endif
        return AMP_ALLOC_MEMOEY_FAILURE;
    }
    return AMP_SUCCESS;
}