#ifndef AMP_H
#define AMP_H

#include <isl/ast.h>

/* A group of array references in a kernel that should be handled together.
 * If private_tile is not NULL, then it is mapped to registers.
 * Otherwise, if shared_tile is not NULL, it is mapped to shared memory.
 * Otherwise, it is accessed from global memory.
 * Note that if both private_tile and shared_tile are set, then shared_tile
 * is only used inside group_common_shared_memory_tile.
 */
struct amp_array_ref_group
{
    /* The references in this group access this local array. */
    struct amp_local_array_info *local_array;
    /* This is the corresponding array. */
    struct amp_array_info *array;
    /* Position of this group in the list of reference groups of array. */
    int nr;

    /* The following fields are use during the construction of the groups.
	 * access is the combined access relation relative to the private
	 * memory tiling.  In particular, the domain of the map corresponds
	 * to the first thread_depth dimensions of the kernel schedule.
	 * write is set if any access in the group is a write.
	 * exact_write is set if all writes are definite writes.
	 * slice is set if there is at least one access in the group
	 * that refers to more than one element
	 * "min_depth" is the minimum of the tile depths and thread_depth.
	 */
    isl_map *access;
    int write;
    int exact_write;
    int slice;
    int min_depth;

    /* The shared memory tile, NULL if none. */
    struct amp_array_tile *shared_tile;

    /* References in this group; point to elements of a linked list. */
    int n_ref;
    struct amp_stmt_access **refs;
};

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

/* Internal data structure for gpu_group_references.
 *
 * scop represents the input scop.
 * kernel_depth is the schedule depth where the kernel launch will
 * be introduced, i.e., it is the depth of the band that is mapped
 * to blocks.
 * shared_depth is the schedule depth at which the copying to/from
 * shared memory is computed.  The copy operation may then
 * later be hoisted to a higher level.
 * thread_depth is the schedule depth where the thread mark is located,
 * i.e., it is the depth of the band that is mapped to threads and also
 * the schedule depth at which the copying to/from private memory
 * is computed.  The copy operation may then later be hoisted to
 * a higher level.
 * n_thread is the number of schedule dimensions in the band that
 * is mapped to threads.
 * privatization lives in the range of thread_sched (i.e., it is
 * of dimension thread_depth + n_thread) and encodes the mapping
 * to thread identifiers (as parameters).
 * host_sched contains the kernel_depth dimensions of the host schedule.
 * shared_sched contains the first shared_depth dimensions of the
 * kernel schedule.
 * copy_sched contains the first thread_depth dimensions of the
 * kernel schedule.
 * thread_sched contains the first (thread_depth + n_thread) dimensions
 * of the kernel schedule.
 * full_sched is a union_map representation of the entire kernel schedule.
 * The schedules are all formulated in terms of the original statement
 * instances, i.e., those that appear in the domains of the access
 * relations.
 */
struct amp_group_data
{
    struct ppcg_scop *scop;
    int kernel_depth;
    int shared_depth;
    int thread_depth;
    int n_thread;
    isl_set *privatization;
    isl_union_map *host_sched;
    isl_union_map *shared_sched;
    isl_union_map *copy_sched;
    isl_union_map *thread_sched;
    isl_union_map *full_sched;
};

/* Represents an outer array accessed by a ppcg_kernel, localized
 * to the context of this kernel.
 *
 * "array" points to the corresponding array in the gpu_prog.
 * The "n_group" "groups" are the reference groups associated to the array.
 * If "force_private" is set, then the array (in practice a scalar)
 * must be mapped to a register.
 * "global" is set if the global device memory corresponding
 * to this array is accessed by the kernel.
 * "bound" is equal to array->bound specialized to the current kernel.
 * "bound_expr" is the corresponding access AST expression.
 */
struct amp_local_array_info
{
    struct amp_array_info *array;

    int n_group;
    struct amp_array_ref_group **groups;

    int global;

    unsigned n_index;
    isl_multi_pw_aff *bound;
    isl_ast_expr *bound_expr;
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

    /* Identifier of the next kernel. */
    int kernel_id;
} amp_prog;

enum ppcg_group_access_type
{
    ppcg_access_global,
    ppcg_access_shared
};

/* Representation of a local variable in a amp kernel.
 */
struct amp_ppcg_kernel_var
{
    struct amp_array_info *array;
    enum ppcg_group_access_type type;
    char *name;
    isl_vec *size;
};

/* Representation of a kernel.
 *
 * prog describes the original code from which the kernel is extracted.
 *
 * id is the sequence number of the kernel.
 *
 * block_ids contains the list of block identifiers for this kernel.
 * thread_ids contains the list of thread identifiers for this kernel.
 *
 * the first n_grid elements of grid_dim represent the specified size
 * of the grid.
 * the first n_block elements of block_dim represent the specified or
 * effective size of the block.
 * Note that in the input file, the sizes of the grid and the blocks
 * are specified in the order x, y, z, but internally, the sizes
 * are stored in reverse order, so that the last element always
 * refers to the x dimension.
 *
 * grid_size reflects the effective grid size.
 * grid_size_expr contains a corresponding access AST expression, built within
 * the context where the launch appears.
 *
 * context contains the values of the parameters and outer schedule dimensions
 * for which any statement instance in this kernel needs to be executed.
 *
 * n_sync is the number of synchronization operations that have
 * been introduced in the schedule tree corresponding to this kernel (so far).
 *
 * core contains the spaces of the statement domains that form
 * the core computation of the kernel.  It is used to navigate
 * the tree during the construction of the device part of the schedule
 * tree in gpu_create_kernel.
 *
 * expanded_domain contains the original statement instances,
 * i.e., those that appear in the domains of access relations,
 * that are involved in the kernel.
 * contraction maps those original statement instances to
 * the statement instances that are active at the point
 * in the schedule tree where the kernel is created.
 *
 * arrays is the set of possibly accessed outer array elements.
 *
 * space is the schedule space of the AST context.  That is, it represents
 * the loops of the generated host code containing the kernel launch.
 *
 * n_array is the total number of arrays in the input program and also
 * the number of element in the array array.
 * array contains information about each array that is local
 * to the current kernel.  If an array is not used in a kernel,
 * then the corresponding entry does not contain any information.
 *
 * any_force_private is set if any array in the kernel is marked force_private
 *
 * block_filter contains constraints on the domain elements in the kernel
 * that encode the mapping to block identifiers, where the block identifiers
 * are represented by "n_grid" parameters with as names the elements
 * of "block_ids".
 *
 * thread_filter contains constraints on the domain elements in the kernel
 * that encode the mapping to thread identifiers, where the thread identifiers
 * are represented by "n_block" parameters with as names the elements
 * of "thread_ids".
 *
 * copy_schedule corresponds to the schedule dimensions of
 * the (tiled) schedule for this kernel that have been taken into account
 * for computing private/shared memory tiles.
 * The domain corresponds to the original statement instances, i.e.,
 * those that appear in the leaves of the schedule tree.
 * copy_schedule_dim is the dimension of this schedule.
 *
 * sync_writes contains write references that require synchronization.
 * Each reference is represented by a universe set in a space [S[i,j] -> R[]]
 * with S[i,j] the statement instance space and R[] the array reference.
 */
struct amp_ppcg_kernel
{
    isl_ctx *ctx;
    struct ppcg_options *options;

    struct amp_prog *prog;

    int id;

    isl_ast_expr *size_expr;
    isl_set *context;

    isl_union_set *core;
    isl_union_set *arrays;

    isl_union_pw_multi_aff *contraction;
    isl_union_set *expanded_domain;

    isl_space *space;

    int n_array;
    struct amp_local_array_info *array;

    int n_var;
    struct amp_ppcg_kernel_var *var;

    isl_union_set *thread_filter;
    isl_union_pw_multi_aff *copy_schedule;
    int copy_schedule_dim;

    isl_ast_node *tree;
};

/* The current index is such that if you add "shift",
 * then the result is always a multiple of "stride",
 * where "stride" may be equal to 1.
 * Let D represent the initial tile->depth dimensions of the computed schedule.
 * The spaces of "lb" and "shift" are of the form
 *
 *	D -> [b]
 */
struct amp_array_bound
{
    isl_val *size;
    isl_aff *lb;

    isl_val *stride;
    isl_aff *shift;
};

/* A tile of an outer array.
 *
 * requires_unroll is set if the schedule dimensions that are mapped
 * to threads need to be unrolled for this (private) tile to be used.
 *
 * "depth" reflects the number of schedule dimensions that affect the tile.
 * The copying into and/or out of the tile is performed at that depth.
 *
 * n is the dimension of the array.
 * bound is an array of size "n" representing the lower bound
 *	and size for each index.
 *
 * tiling maps a tile in the global array to the corresponding
 * shared/private memory tile and is of the form
 *
 *	{ [D[i] -> A[a]] -> T[(a + shift(i))/stride - lb(i)] }
 *
 * where D represents the initial "depth" dimensions
 * of the computed schedule.
 */
struct amp_array_tile
{
    isl_ctx *ctx;
    int requires_unroll;
    int depth;
    int n;
    struct amp_array_bound *bound;
    isl_multi_aff *tiling;
};

amp_prog *amp_prog_alloc(__isl_take isl_ctx *ctx, struct ppcg_scop *scop);
void *amp_prog_free(amp_prog *prog);
__isl_give isl_ast_node *amp_build_array_bounds(__isl_take isl_ast_node *node, amp_prog *prog, __isl_keep isl_ast_build *build);

__isl_give isl_printer *amp_print_macros(__isl_take isl_printer *p);
__isl_give isl_printer *declare_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog);
__isl_give isl_printer *allocate_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog);
char *amp_get_lower_precision_type(char *type);

__isl_give isl_schedule *amp_schedule_again(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_keep isl_schedule *schd);

int amp_array_is_read_only_scalar(struct amp_array_info *array);
int amp_array_is_scalar(struct amp_array_info *array);

struct amp_array_tile *amp_array_ref_group_tile(struct amp_array_ref_group *group);
#endif