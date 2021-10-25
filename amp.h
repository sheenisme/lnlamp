#ifndef AMP_H
#define AMP_H

#include <isl/ast.h>

/* An access to an outer array element or an iterator.
 * Accesses to iterators have an access relation that maps to an unnamed space.
 * An access may be both read and write.
 * If the access relation is empty, then the output dimension may
 * not be equal to the dimension of the corresponding array.
 */
struct amp_stmt_access
{
    /* Access reads elements */
    int read;
    /* Access writes elements */
    int write;
    /* All writes are definite writes. */
    int exact_write;
    /* Is a single, fixed element being accessed? */
    isl_bool fixed_element;
    /* The number of index expressions specified in the access. */
    int n_index;

    /* May access relation */
    isl_map *access;
    /* May access relation with as domain a mapping from iteration domain
	 * to a reference identifier.
	 */
    isl_map *tagged_access;
    /* The reference id of the corresponding pet_expr. */
    isl_id *ref_id;

    struct amp_stmt_access *next;
};

/* A representation of a user statement.
 * "stmt" points to the corresponding pet statement.
 * "id" is the identifier of the instance set of the statement.
 * "accesses" is a linked list of accesses performed by the statement.
 * If the statement has been killed, i.e., if it will not be scheduled,
 * then this linked list may be empty even if the actual statement does
 * perform accesses.
 */
struct amp_stmt
{
    isl_id *id;
    struct pet_stmt *stmt;

    struct amp_stmt_access *accesses;
};
/* Represents an outer array possibly accessed by a amp_prog.
 */
struct amp_array_info
{
    /* The array data space. */
    isl_space *space;
    /* Element type. */
    char *type;
    /* Element size. */
    int size;
    /* Name of the array. */
    char *name;
    /* Declared extent of original array. */
    isl_set *declared_extent;
    /* AST expression for declared size of original array. */
    isl_ast_expr *declared_size;
    /* Extent of the array that needs to be copied. */
    isl_set *extent;
    /* Number of indices. */
    unsigned n_index;
    /* For each index, a bound on "extent" in that direction. */
    isl_multi_pw_aff *bound;
    /* The corresponding access AST expression, if the array needs
	 * to be allocated on the device.
	 */
    isl_ast_expr *bound_expr;

    /* All references to this array; point to elements of a linked list. */
    int n_ref;
    struct amp_stmt_access **refs;

    /* Is this array accessed at all by the program? */
    int accessed;

    /* Is this a scalar that is read-only within the entire program? */
    int read_only_scalar;

    /* Are the elements of the array structures? */
    int has_compound_element;

    /* Are the elements only accessed through constant index expressions? */
    int only_fixed_element;

    /* Is the array local to the scop? */
    int local;
    /* Is the array local and should it be declared on the host? */
    int declare_local;

    /* Is the corresponding global device memory accessed in any way? */
    int global;

    /* Should the array be linearized? */
    int linearize;

    /* Order dependences on this array.
	 * Only used if live_range_reordering option is set.
	 * It is set to NULL otherwise.
	 */
    isl_union_map *dep_order;
};

/* "read" and "write" contain the original access relations, possibly
 * involving member accesses.
 *
 * The elements of "array", as well as the ranges of "copy_in" and "copy_out"
 * only refer to the outer arrays of any possible member accesses.
 */
typedef struct amp_prog
{
    isl_ctx *ctx;

    struct ppcg_scop *scop;

    /* Set of parameter values */
    isl_set *context;

    /* All potential read accesses in the entire program */
    isl_union_map *read;

    /* All potential write accesses in the entire program */
    isl_union_map *may_write;
    /* All definite write accesses in the entire program */
    isl_union_map *must_write;
    /* All tagged definite kills in the entire program */
    isl_union_map *tagged_must_kill;

    /* The set of inner array elements that may be preserved. */
    isl_union_set *may_persist;

    /* A mapping from all innermost arrays to their outer arrays. */
    isl_union_map *to_outer;
    /* A mapping from the outer arrays to all corresponding inner arrays. */
    isl_union_map *to_inner;
    /* A mapping from all intermediate arrays to their outer arrays,
	 * including an identity mapping from the anonymous 1D space to itself.
	 */
    isl_union_map *any_to_outer;

    /* Order dependences on non-scalars. */
    isl_union_map *array_order;

    /* Array of statements */
    int n_stmts;
    struct amp_stmt *stmts;

    int n_array;
    struct amp_array_info *array;
} amp_prog;

amp_prog *amp_prog_alloc(__isl_take isl_ctx *ctx, struct ppcg_scop *scop);
void *amp_prog_free(amp_prog *prog);
__isl_give isl_ast_node *amp_build_array_bounds(__isl_take isl_ast_node *node, amp_prog *prog, __isl_keep isl_ast_build *build);

__isl_give isl_printer *amp_print_macros(__isl_take isl_printer *p);
__isl_give isl_printer *declare_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog);
__isl_give isl_printer *allocate_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog);

__isl_give isl_schedule *amp_schedule_again(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_keep isl_schedule *schd);
#endif