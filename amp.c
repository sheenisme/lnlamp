#include <isl/aff.h>
#include <isl/ast.h>

static __isl_give isl_printer *print_amp_macros(__isl_take isl_printer *p)
{
    /* const char *macros =
        "#define ampCheckReturn(ret) \\\n"
        "  do { \\\n"
        "    cudaError_t cudaCheckReturn_e = (ret); \\\n"
        "    if (cudaCheckReturn_e != cudaSuccess) { \\\n"
        "      fprintf(stderr, \"CUDA error: %s\\n\", "
        "cudaGetErrorString(cudaCheckReturn_e)); \\\n"
        "      fflush(stderr); \\\n"
        "    } \\\n"
        "    assert(cudaCheckReturn_e == cudaSuccess); \\\n"
        "  } while(0)\n\n";
    */
    const char *macros = "test\\n";

    p = isl_printer_print_str(p, macros);
    return p;
}

/* Set the names of the macros that may appear in a printed isl AST.
 */
__isl_give isl_printer *amp_print_macros(__isl_take isl_printer *p)
{
    p = print_amp_macros(p);

    return p;
}