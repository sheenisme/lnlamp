#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include <isl/polynomial.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/map_to_basic_set.h>
#include <isl/set.h>
#include <isl/aff.h>
#include <isl/ilp.h>
#include <isl/flow.h>
#include <isl/schedule.h>
#include <isl/schedule_node.h>
#include <isl/options.h>
#include <isl/ast.h>
#include <isl/val.h>
#include <isl/ast_build.h>

#include "amp.h"
#include "hybrid.h"
#include "schedule.h"
#include "ppcg_options.h"
#include "print.h"
#include "util.h"

// dump amp_stmt_access
void amp_stmt_access_dump(struct amp_stmt_access *access)
{
    if (access)
    {
        fprintf(stderr, "       the   amp_stmt_access  is on the below :\n");
        fprintf(stderr, "          the       read      is %d \n", access->read);
        fprintf(stderr, "          the       write     is %d \n", access->write);
        fprintf(stderr, "          the    exact_write  is %d \n", access->exact_write);
        fprintf(stderr, "          the   fixed_element is %d \n", access->fixed_element);
        fprintf(stderr, "          the      n_index    is %d \n", access->n_index);
        fprintf(stderr, "          the    exact_write  is %d \n", access->read);
        fprintf(stderr, "          the     access      is:   \n");
        isl_map_dump(access->access);
        fprintf(stderr, "          the  tagged_access  is: \n");
        isl_map_dump(access->tagged_access);
        fprintf(stderr, "          the     ref_id      is: \n");
        isl_id_dump(access->ref_id);
        fprintf(stderr, "          the      next       is: \n");
        if (access->next)
        {
            amp_stmt_access_dump(access->next);
        }
        else
            fprintf(stderr, "NULL\n");
    }
    else
    {
        fprintf(stderr, "    @WARN:          sorry, the amp_stmt_access is NULL !!! \n");
    }
}

// dump the amp_group_data
void amp_group_data_dump(struct amp_group_data *data)
{
    if (data)
    {
        fprintf(stderr, "       the     amp_group_data    is on the below : \n");
        fprintf(stderr, "          the   kernel_depth     is %d \n", data->kernel_depth);
        fprintf(stderr, "          the   shared_depth     is %d \n", data->shared_depth);
        fprintf(stderr, "          the   thread_depth     is %d \n", data->thread_depth);
        fprintf(stderr, "          the     n_thread       is %d \n", data->n_thread);
        // fprintf(stderr, "          the   privatization    is:   \n");
        // isl_set_dump(data->privatization);
        fprintf(stderr, "          the    host_sched     is:   \n");
        isl_union_map_dump(data->host_sched);
        fprintf(stderr, "          the   shared_sched     is:   \n");
        isl_union_map_dump(data->shared_sched);
        fprintf(stderr, "          the    copy_sched      is:   \n");
        isl_union_map_dump(data->copy_sched);
        fprintf(stderr, "          the   thread_sched     is:   \n");
        isl_union_map_dump(data->thread_sched);
        fprintf(stderr, "          the    full_sched      is:   \n");
        isl_union_map_dump(data->full_sched);
    }
    else
    {
        fprintf(stderr, "    @WARN:\n       sorry, the amp_group_data is NULL !!! \n");
    }
}

// dump the amp_array_bound
void amp_array_bound_dump(struct amp_array_bound *bound)
{
    if (bound)
    {
        fprintf(stderr, "       the amp_array_bound is on the below :\n");
        fprintf(stderr, "          the      size    is: \n");
        if (bound->size)
            isl_val_dump(bound->size);
        fprintf(stderr, "          the      lb      is: \n");
        if (bound->lb)
            isl_aff_dump(bound->lb);
        fprintf(stderr, "          the    stride    is: \n");
        if (bound->stride)
            isl_val_dump(bound->stride);
        fprintf(stderr, "          the    shift     is: \n");
        if (bound->shift)
            isl_aff_dump(bound->shift);
    }
    else
    {
        fprintf(stderr, "   @WARN: \n       sorry, the amp_array_bound is NULL !!! \n");
    }
}

// dump the amp_array_tile
void amp_array_tile_dump(struct amp_array_tile *tile)
{
    if (tile)
    {
        fprintf(stderr, "       the amp_array_tile is on the below :\n");
        fprintf(stderr, "          the   depth     is %d \n", tile->depth);
        fprintf(stderr, "          the     n       is %d \n", tile->n);
        fprintf(stderr, "          the   tiling    is:   \n");
        if (tile->tiling)
            isl_multi_aff_dump(tile->tiling);
        fprintf(stderr, "          the   bound     is:   \n");
        amp_array_bound_dump(tile->bound);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_array_tile is NULL !!! \n");
    }
}

// dump the amp_array_info
void amp_array_info_dump(struct amp_array_info *info)
{
    if (info)
    {
        fprintf(stderr, "       the       amp_array_info     is on the below :\n");
        fprintf(stderr, "          the         type          is %s \n", info->type);
        fprintf(stderr, "          the         name          is %s \n", info->name);
        fprintf(stderr, "          the         size          is %d \n", info->size);
        fprintf(stderr, "          the        space          is :  \n");
        isl_space_dump(info->space);
        fprintf(stderr, "          the    declared_extent    is :  \n");
        isl_set_dump(info->declared_extent);
        fprintf(stderr, "          the     declared_size     is :  \n");
        isl_ast_expr_dump(info->declared_size);
        fprintf(stderr, "          the        extent         is :  \n");
        isl_set_dump(info->extent);
        fprintf(stderr, "          the        bound          is :  \n");
        isl_multi_pw_aff_dump(info->bound);
        fprintf(stderr, "          the      bound_expr       is :  \n");
        isl_ast_expr_dump(info->bound_expr);
        fprintf(stderr, "          the         n_ref         is %d \n", info->n_ref);
        for (int r = 0; r < info->n_ref; r++)
        {
            fprintf(stderr, "          the         refs[%d]      is:   \n", r);
            amp_stmt_access_dump(info->refs[r]);
        }
        fprintf(stderr, "          the        n_index        is %d \n", info->n_index);
        fprintf(stderr, "          the       accessed        is %d \n", info->accessed);
        fprintf(stderr, "          the   read_only_scalar    is %d \n", info->read_only_scalar);
        fprintf(stderr, "          the has_compound_element  is %d \n", info->has_compound_element);
        fprintf(stderr, "          the  only_fixed_element   is %d \n", info->only_fixed_element);
        fprintf(stderr, "          the         local         is %d \n", info->local);
        fprintf(stderr, "          the    declare_local      is %d \n", info->declare_local);
        fprintf(stderr, "          the         global        is %d \n", info->global);
        fprintf(stderr, "          the      linearize        is %d \n", info->linearize);
        fprintf(stderr, "          the      dep_order        is :  \n");
        isl_union_map_dump(info->dep_order);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_array_info is NULL !!! \n");
    }
}

// dump the amp_array_ref_group
void amp_array_ref_group_dump(struct amp_array_ref_group *group)
{
    if (group)
    {
        fprintf(stderr, "       the amp_array_ref_group is on the below :\n");
        fprintf(stderr, "          the      nr          is %d \n", group->nr);
        fprintf(stderr, "          the    write         is %d \n", group->write);
        fprintf(stderr, "          the  exact_write     is %d \n", group->exact_write);
        fprintf(stderr, "          the    slice         is %d \n", group->slice);
        fprintf(stderr, "          the   min_depth      is %d \n", group->min_depth);
        fprintf(stderr, "          the     n_ref        is %d \n", group->n_ref);
        fprintf(stderr, "          the    access        is: \n");
        isl_map_dump(group->access);
        fprintf(stderr, "          the  shared_tile     is: \n");
        amp_array_tile_dump(group->shared_tile);
        fprintf(stderr, "          the     array        is: \n");
        amp_array_info_dump(group->array);
        // fprintf(stderr, "          the    local_array   is: \n");
        // amp_local_array_info_dump(group->local_array);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_array_ref_group is NULL !!! \n");
    }
}

// dump the amp_local_array_info
void amp_local_array_info_dump(struct amp_local_array_info *info)
{
    if (info)
    {
        fprintf(stderr, "       the amp_local_array_info is on the below :\n");
        fprintf(stderr, "          the     global        is %d \n", info->global);
        fprintf(stderr, "          the     n_index       is %d \n", info->n_index);
        fprintf(stderr, "          the      bound        is :  \n");
        isl_multi_pw_aff_dump(info->bound);
        fprintf(stderr, "          the     bound_expr    is :  \n");
        isl_ast_expr_dump(info->bound_expr);
        fprintf(stderr, "          the     n_group       is %d \n", info->n_group);
        for (int i = 0; i < info->n_group; i++)
        {
            fprintf(stderr, "          the    groups[%d]     is :  \n", i);
            amp_array_ref_group_dump(info->groups[i]);
        }
        fprintf(stderr, "          the     array         is :  \n");
        amp_array_info_dump(info->array);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_local_array_info is NULL !!! \n");
    }
}

// dump amp_stmt
void amp_stmt_dump(struct amp_stmt *stmts)
{
    if (stmts)
    {
        fprintf(stderr, "       the      amp_stmt   is on the below :\n");
        fprintf(stderr, "          the      id      is: \n");
        isl_id_dump(stmts->id);
        fprintf(stderr, "          the   accesses   is: \n");
        amp_stmt_access_dump(stmts->accesses);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_stmt is NULL !!! \n");
    }
}

// dump amp_prog
void amp_prog_dump(struct amp_prog *prog)
{
    if (prog)
    {
        fprintf(stderr, "       the      amp_prog       is on the below :\n");
        fprintf(stderr, "          the      context     is: \n");
        isl_set_dump(prog->context);
        fprintf(stderr, "          the      read        is: \n");
        isl_union_map_dump(prog->read);
        fprintf(stderr, "          the    may_write     is: \n");
        isl_union_map_dump(prog->may_write);
        fprintf(stderr, "          the    must_write    is: \n");
        isl_union_map_dump(prog->must_write);
        fprintf(stderr, "          the tagged_must_kill is: \n");
        isl_union_map_dump(prog->tagged_must_kill);
        fprintf(stderr, "          the   may_persist   is: \n");
        isl_union_set_dump(prog->may_persist);
        fprintf(stderr, "          the      to_outer   is: \n");
        isl_union_map_dump(prog->to_outer);
        fprintf(stderr, "          the      to_inner   is: \n");
        isl_union_map_dump(prog->to_inner);
        fprintf(stderr, "          the  any_to_outer   is: \n");
        isl_union_map_dump(prog->any_to_outer);
        fprintf(stderr, "          the    array_order  is: \n");
        isl_union_map_dump(prog->array_order);
        fprintf(stderr, "          the   n_stmts   is: %d \n", prog->n_stmts);
        for (int s = 0; s < prog->n_stmts; s++)
        {
            fprintf(stderr, "          the  stmts[%d]  is: \n", s);
            amp_stmt_dump(prog->stmts + s);
        }
        fprintf(stderr, "          the   n_array   is: %d \n", prog->n_array);
        for (int a = 0; a < prog->n_array; a++)
        {
            fprintf(stderr, "          the  array[%d]  is: \n", a);
            amp_array_info_dump(prog->array + a);
        }
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_prog is NULL !!! \n");
    }
}

// dump amp_ppcg_kernel
void amp_ppcg_kernel_dump(struct amp_ppcg_kernel *kernel)
{
    if (kernel)
    {
        fprintf(stderr, "       the   amp_ppcg_kernel  is on the below :\n");
        fprintf(stderr, "          the       id        is: %d \n", kernel->id);
        fprintf(stderr, "          the      prog       is: \n");
        amp_prog_dump(kernel->prog);
        fprintf(stderr, "          the   size_expr     is: \n");
        isl_ast_expr_dump(kernel->size_expr);
        fprintf(stderr, "          the    context      is: \n");
        isl_set_dump(kernel->context);
        fprintf(stderr, "          the     core        is: \n");
        isl_union_set_dump(kernel->core);
        fprintf(stderr, "          the    arrays       is: \n");
        isl_union_set_dump(kernel->arrays);
        fprintf(stderr, "          the    contraction  is: \n");
        isl_union_pw_multi_aff_dump(kernel->contraction);
        fprintf(stderr, "          the expanded_domain is: \n");
        isl_union_set_dump(kernel->expanded_domain);
        fprintf(stderr, "          the     space       is: \n");
        isl_space_dump(kernel->space);
        fprintf(stderr, "          the    n_array      is: %d \n", kernel->n_array);
        for (int r = 0; r < kernel->n_array; r++)
        {
            fprintf(stderr, "          the     array[%d]       is: \n", r);
            amp_local_array_info_dump(kernel->array);
        }
        // fprintf(stderr, "          the     array       is: \n");
        // amp_local_array_info_dump(kernel->array);
        fprintf(stderr, "          the      n_var      is: %d \n", kernel->n_var);
        fprintf(stderr, "          the  thread_filter  is: \n");
        isl_union_set_dump(kernel->thread_filter);
        fprintf(stderr, "          the  copy_schedule  is: \n");
        isl_union_pw_multi_aff_dump(kernel->copy_schedule);
        fprintf(stderr, "          the copy_schedule_dim is: %d \n", kernel->copy_schedule_dim);
        fprintf(stderr, "          the       tree       is: \n");
        isl_ast_node_dump(kernel->tree);
    }
    else
    {
        fprintf(stderr, "    @WARN: \n       sorry, the amp_ppcg_kernel is NULL !!! \n");
    }
}

/* Return the name of the outer array (of structs) accessed by "access".
 */
static const char *get_outer_array_name(__isl_keep isl_map *access)
{
    isl_space *space;
    const char *name;

    space = isl_space_range(isl_map_get_space(access));
    while (space && isl_space_is_wrapping(space))
        space = isl_space_domain(isl_space_unwrap(space));
    name = isl_space_get_tuple_name(space, isl_dim_set);
    isl_space_free(space);

    return name;
}

/* Collect all references to the given array and store pointers to them
 * in array->refs.
 */
static isl_stat collect_references(amp_prog *prog, struct amp_array_info *array)
{
    // #define DEBUG_COLLECT_REFERRNCES
    int i;
    int n;

    n = 0;
    for (i = 0; i < prog->n_stmts; ++i)
    {
        struct amp_stmt *stmt = &prog->stmts[i];
        struct amp_stmt_access *access;

        for (access = stmt->accesses; access; access = access->next)
        {
            const char *name;
            name = get_outer_array_name(access->access);
            if (name && !strcmp(array->name, name))
                n++;
        }
    }

    array->refs = isl_alloc_array(prog->ctx, struct amp_stmt_access *, n);
    if (!array->refs)
        return isl_stat_error;
    array->n_ref = n;

    n = 0;
    for (i = 0; i < prog->n_stmts; ++i)
    {
        struct amp_stmt *stmt = &prog->stmts[i];
        struct amp_stmt_access *access;

        for (access = stmt->accesses; access; access = access->next)
        {
            const char *name;
            name = get_outer_array_name(access->access);
            if (!name || strcmp(array->name, name))
                continue;

            array->refs[n++] = access;
        }
    }

#ifdef DEBUG_COLLECT_REFERRNCES
    fprintf(stderr, "@DEBUG: \n       the array->refs of the collect_references function is:\n");
    for (int r = 0; r < array->n_ref; r++)
    {
        fprintf(stderr, "       array->refs[%d] :\n", r);
        amp_stmt_access_dump(array->refs[r]);
        fprintf(stderr, "\n");
    }
#endif // DEBUG_COLLECT_REFERRNCES
    return isl_stat_ok;
}

/* Compute and return the extent of "array", taking into account the set of
 * accessed elements.
 *
 * In particular, the extent in the outer dimension is taken
 * from "accessed", while the extents in the remaining dimensions
 * are taken from array->extent.
 *
 * The extent in the outer dimension cannot be taken from array->extent
 * because that may be unbounded.  Furthermore, even if it is bounded,
 * it may be larger than the piece of the array that is being accessed.
 */
static __isl_give isl_set *amp_compute_extent(struct pet_array *array, __isl_keep isl_set *accessed)
{
    int n_index;
    isl_id *id;
    isl_set *outer;
    isl_set *extent;

    extent = isl_set_copy(array->extent);

    n_index = isl_set_dim(accessed, isl_dim_set);
    if (n_index == 0)
        return extent;

    extent = isl_set_project_out(extent, isl_dim_set, 0, 1);
    outer = isl_set_copy(accessed);
    outer = isl_set_project_out(outer, isl_dim_set, 1, n_index - 1);
    extent = isl_set_flat_product(outer, extent);
    id = isl_set_get_tuple_id(accessed);
    extent = isl_set_set_tuple_id(extent, id);

    return extent;
}

/* Is the array "array" being extracted a read-only scalar?
 *
 * That is, is "array" a scalar that is never possibly written to.
 * An array containing structures is never considered to be a scalar.
 */
static int is_read_only_scalar(struct amp_array_info *array, amp_prog *prog)
{
    isl_set *space;
    isl_union_map *write;
    int empty;

    if (array->has_compound_element)
        return 0;
    if (array->n_index != 0)
        return 0;

    write = isl_union_map_copy(prog->may_write);
    space = isl_set_universe(isl_space_copy(array->space));
    write = isl_union_map_intersect_range(write, isl_union_set_from_set(space));
    empty = isl_union_map_is_empty(write);
    isl_union_map_free(write);

    return empty;
}

/* Is "array" only accessed as individual, fixed elements?
 * That is, does each access to "array" access a single, fixed element?
 */
static isl_bool only_fixed_element_accessed(struct amp_array_info *array)
{
    int i;

    for (i = 0; i < array->n_ref; ++i)
        if (!array->refs[i]->fixed_element)
            return isl_bool_false;

    return isl_bool_true;
}

/* Compute bounds on the host array "pa" based on the corresponding
 * accessed elements in "arrays"
 * and collect all references to the array.
 * Store the results in "info".
 *
 * If the array is zero-dimensional and does not contain structures,
 * i.e., if the array is a scalar, we check whether it is read-only.
 * We also check whether the array is accessed at all.
 */
static isl_stat extract_array_info(amp_prog *prog, struct amp_array_info *info, struct pet_array *pa, __isl_keep isl_union_set *arrays)
{
    int empty;
    const char *name;
    int n_index;
    isl_multi_pw_aff *bounds;
    isl_set *accessed, *extent;

    n_index = isl_set_dim(pa->extent, isl_dim_set);
    name = isl_set_get_tuple_name(pa->extent);

    info->space = isl_set_get_space(pa->extent);
    info->name = strdup(name);
    info->n_index = n_index;
    info->linearize = 0;

    info->type = strdup(pa->element_type);
    info->size = pa->element_size;
    info->local = pa->declared && !pa->exposed;
    info->has_compound_element = pa->element_is_record;
    info->read_only_scalar = is_read_only_scalar(info, prog);

    info->declared_extent = isl_set_copy(pa->extent);
    accessed = isl_union_set_extract_set(arrays, isl_space_copy(info->space));
    empty = isl_set_is_empty(accessed);
    extent = amp_compute_extent(pa, accessed);
    isl_set_free(accessed);
    info->extent = extent;
    if (empty < 0)
        return isl_stat_error;
    info->accessed = !empty;
    bounds = ppcg_size_from_extent(isl_set_copy(extent));
    bounds = isl_multi_pw_aff_gist(bounds, isl_set_copy(prog->context));
    if (!bounds)
        return isl_stat_error;
    if (!isl_multi_pw_aff_is_cst(bounds))
        info->linearize = 0;
    info->bound = bounds;

    if (collect_references(prog, info) < 0)
        return isl_stat_error;
    info->only_fixed_element = only_fixed_element_accessed(info);

    return isl_stat_ok;
}

/* Remove independence from the order constraints "order" on array "array".
 * Since the pairs of iterations in the filter relation of an independence
 * are guaranteed to be completely independent by the user, there is
 * no need to ensure that live ranges are ordered along those pairs.
 * We make an exception for local variables, though, as the independence
 * guarantee does not apply to those.
 *
 * The order constraints are used in two places.
 * Those on scalars are used in check_scalar_live_ranges to check if
 * we need to force the scalar to be private.  Any non-local scalar
 * should not be forced scalar if it only appears in independent loops.
 * Those on non-scalars are added to the coincidence constraints
 * in compute_schedule because we do not support any array expansion.
 * Accesses to non-local arrays should not prevent a loop from being
 * considered coincident so we should indeed remove those constraints
 * from the order constraints.
 */
static __isl_give isl_union_map *remove_independences(amp_prog *prog, struct amp_array_info *array, __isl_take isl_union_map *order)
{
    int i;

    for (i = 0; i < prog->scop->pet->n_independence; ++i)
    {
        struct pet_independence *pi = prog->scop->pet->independences[i];
        if (isl_union_set_contains(pi->local, array->space))
            continue;

        order = isl_union_map_subtract(order, isl_union_map_copy(pi->filter));
    }

    return order;
}

/* Can "array" be mapped to private memory?
 * That is, is it only accessed as individual elements with
 * constant index expressions?
 */
isl_bool amp_array_can_be_private(struct amp_array_info *array)
{
    if (!array)
        return isl_bool_error;
    return array->only_fixed_element;
}

/* For each array in "prog", store the (untagged) order dependences
 * derived from the array in array->dep_order.
 * In particular, consider all references that access the given array
 * and take the order dependences that have one of these references
 * as source.  (Since an order dependence relates two references to
 * the same array, the target of these order dependences will also
 * be one of these references.)
 * Additionally, store the union of these array->dep_order relations
 * for all arrays that cannot be mapped to private memory in prog->array_order.
 */
/* multiple definition with GPU */
void amp_collect_order_dependences(amp_prog *prog)
{
    int i;
    isl_space *space;
    isl_union_map *accesses;

    space = isl_union_map_get_space(prog->read);
    prog->array_order = isl_union_map_empty(space);

    accesses = isl_union_map_copy(prog->scop->tagged_reads);
    accesses = isl_union_map_union(accesses, isl_union_map_copy(prog->scop->tagged_may_writes));
    accesses = isl_union_map_universe(accesses);
    accesses = isl_union_map_apply_range(accesses, isl_union_map_copy(prog->to_outer));

    for (i = 0; i < prog->n_array; ++i)
    {
        struct amp_array_info *array = &prog->array[i];
        isl_set *set;
        isl_union_set *uset;
        isl_union_map *order;

        set = isl_set_universe(isl_space_copy(array->space));
        uset = isl_union_set_from_set(set);
        uset = isl_union_map_domain(isl_union_map_intersect_range(isl_union_map_copy(accesses), uset));
        order = isl_union_map_copy(prog->scop->tagged_dep_order);
        order = isl_union_map_intersect_domain(order, uset);
        order = isl_union_map_zip(order);
        order = isl_union_set_unwrap(isl_union_map_domain(order));
        order = remove_independences(prog, array, order);
        array->dep_order = order;

        if (amp_array_can_be_private(array))
            continue;

        prog->array_order = isl_union_map_union(prog->array_order, isl_union_map_copy(array->dep_order));
    }

    isl_union_map_free(accesses);
}

/* Construct a amp_array_info for each array referenced by amp->scop and
 * collect them in amp->array.
 *
 * The sizes are based on the extents and the set of possibly accessed
 * elements by "prog".
 * If there are any member accesses involved, then they are first mapped
 * to the outer arrays of structs.
 * Only extract amp_array_info entries for these outer arrays.
 *
 * If we are allowing live range reordering, then also set
 * the dep_order field.  Otherwise leave it NULL.
 */
static isl_stat collect_array_info(amp_prog *prog)
{
    int i;
    isl_stat r = isl_stat_ok;
    isl_union_set *arrays;

    prog->n_array = 0;
    prog->array = isl_calloc_array(prog->ctx, struct amp_array_info, prog->scop->pet->n_array);
    if (!prog->array)
        return isl_stat_error;

    arrays = isl_union_map_range(isl_union_map_copy(prog->read));
    arrays = isl_union_set_union(arrays, isl_union_map_range(isl_union_map_copy(prog->may_write)));

    arrays = isl_union_set_apply(arrays, isl_union_map_copy(prog->to_outer));

    arrays = isl_union_set_coalesce(arrays);

    for (i = 0; i < prog->scop->pet->n_array; ++i)
    {
        isl_bool field;

        field = isl_set_is_wrapping(prog->scop->pet->arrays[i]->extent);
        if (field < 0)
            break;
        if (field)
            continue;
        if (extract_array_info(prog, &prog->array[prog->n_array++], prog->scop->pet->arrays[i], arrays) < 0)
            r = isl_stat_error;
    }
    if (i < prog->scop->pet->n_array)
        r = isl_stat_error;

    isl_union_set_free(arrays);

    if (prog->scop->options->live_range_reordering)
        amp_collect_order_dependences(prog);

    return r;
}

// free array_info in amp_prog
static void free_array_info(amp_prog *prog)
{
    int i;

    for (i = 0; i < prog->n_array; ++i)
    {
        free(prog->array[i].type);
        free(prog->array[i].name);
        isl_multi_pw_aff_free(prog->array[i].bound);
        isl_ast_expr_free(prog->array[i].bound_expr);
        isl_space_free(prog->array[i].space);
        isl_set_free(prog->array[i].declared_extent);
        isl_set_free(prog->array[i].extent);
        isl_ast_expr_free(prog->array[i].declared_size);
        free(prog->array[i].refs);
        isl_union_map_free(prog->array[i].dep_order);
    }
    free(prog->array);
}

/* Check if a amp array is a scalar.  A scalar is a value that is not stored
 * as an array or through a pointer reference, but as a single data element.
 * At the moment, scalars are represented as zero-dimensional arrays.
 * Note that the single data element may be an entire structure.
 */
int amp_array_is_scalar(struct amp_array_info *array)
{
    return array->n_index == 0;
}

/* Is "array" a read-only scalar?
 */
int amp_array_is_read_only_scalar(struct amp_array_info *array)
{
    return array->read_only_scalar;
}

/* Has statement "stmt" been killed from "scop"?
 * That is, is the instance set of "scop" free from any
 * instances of "stmt"?
 */
static isl_bool is_stmt_killed(struct ppcg_scop *scop, struct pet_stmt *stmt)
{
    isl_space *space;
    isl_set *left;
    isl_bool empty;

    if (!scop || !stmt)
        return isl_bool_error;
    space = isl_set_get_space(stmt->domain);
    left = isl_union_set_extract_set(scop->domain, space);
    empty = isl_set_plain_is_empty(left);
    isl_set_free(left);

    return empty;
}

/* Does the index expression "index" of "expr" represent an access
 * to a single element?
 * That is, is "index" completely specified?
 *
 * If "expr" accesses elements from different spaces (i.e., fields
 * of a structure), then it does not access a single element.
 * Otherwise, if the single space of the access matches the space
 * of "index", then the index expression is completely specified
 * (no pointer to a lower-dimensional slice of the accessed array)
 * and a single element is being accessed.
 */
static isl_bool complete_index(__isl_keep pet_expr *expr, __isl_keep isl_multi_pw_aff *index)
{
    isl_union_map *read, *write, *all;
    isl_map *map;
    isl_space *space1, *space2;
    isl_bool complete;

    read = pet_expr_access_get_may_read(expr);
    write = pet_expr_access_get_may_write(expr);
    all = isl_union_map_union(read, write);
    if (!all)
        return isl_bool_error;
    if (isl_union_map_n_map(all) != 1)
    {
        isl_union_map_free(all);
        return isl_bool_false;
    }
    map = isl_map_from_union_map(all);
    space1 = isl_map_get_space(map);
    isl_map_free(map);
    space2 = isl_multi_pw_aff_get_space(index);
    complete = isl_space_tuple_is_equal(space1, isl_dim_out,
                                        space2, isl_dim_out);
    isl_space_free(space1);
    isl_space_free(space2);

    return complete;
}

static isl_bool accesses_fixed_element(__isl_keep pet_expr *expr)
{
    int i, n;
    isl_multi_pw_aff *index;
    isl_bool fixed = isl_bool_true;

    index = pet_expr_access_get_index(expr);
    if (index < 0)
        return isl_bool_error;
    n = isl_multi_pw_aff_dim(index, isl_dim_out);
    for (i = 0; i < n; ++i)
    {
        isl_pw_aff *pa;

        pa = isl_multi_pw_aff_get_pw_aff(index, 0);
        fixed = isl_pw_aff_n_piece(pa) == 1;
        if (fixed)
            fixed = isl_pw_aff_is_cst(pa);
        isl_pw_aff_free(pa);
        if (fixed < 0 || !fixed)
            break;
    }
    if (fixed >= 0 && fixed)
        fixed = complete_index(expr, index);
    isl_multi_pw_aff_free(index);

    return fixed;
}
/* Internal data structure for extract_access.
 * "next_access" points to the end of a linked list that is extended
 * by extract_access.
 * "single_expression" is set if the access expressions belong to
 * an expression statement (i.e., a statement without internal control).
 * "any_to_outer" maps all intermediate arrays to their outer arrays.
 */
struct ppcg_extract_access_data
{
    struct amp_stmt_access **next_access;
    int single_expression;
    isl_union_map *any_to_outer;
};

/* Given a tagged access relation to a single array "tagged", extract it
 * as a map, taking into account that the input may be empty.
 * If the access relation is empty, then it does not contain
 * any space information, so we try to recover it from the index
 * expression.
 * The space of the index expression is of the form I -> A,
 * with I the statement instances and A the array, or [I -> F] -> A,
 * with F the filters corresponding to arguments.
 * We first drop F, if present, obtaining I -> A.
 * Then we construct I -> R, with R the reference tag,
 * combine the two into I -> [R -> A] and uncurry to obtain
 * the final result [I -> R] -> A.
 * Note that the index expression may have a lower dimension
 * than that of the array, but this dimension is not used
 * if the access relation is empty.
 */
static __isl_give isl_map *extract_single_tagged_access(
    __isl_take isl_union_map *tagged, __isl_keep pet_expr *expr)
{
    int empty;
    isl_id *id;
    isl_space *space, *space2;
    isl_multi_pw_aff *index;

    empty = isl_union_map_is_empty(tagged);
    if (empty < 0)
        goto error;
    if (!empty)
        return isl_map_from_union_map(tagged);
    isl_union_map_free(tagged);

    index = pet_expr_access_get_index(expr);
    space = isl_multi_pw_aff_get_space(index);
    isl_multi_pw_aff_free(index);
    if (isl_space_domain_is_wrapping(space))
        space = isl_space_domain_factor_domain(space);
    space2 = isl_space_from_domain(isl_space_domain(isl_space_copy(space)));
    id = pet_expr_access_get_ref_id(expr);
    space2 = isl_space_set_tuple_id(space2, isl_dim_out, id);
    space = isl_space_range_product(space2, space);
    space = isl_space_uncurry(space);

    return isl_map_empty(space);
error:
    isl_union_map_free(tagged);
    return NULL;
}

/* Extract a amp_stmt_access from "expr", append it to the list
 * that ends in *data->next_access and update the end of the list.
 * If the access expression performs a write, then it is considered
 * exact only if it appears in a single expression statement and
 * if its may access relation is equal to its must access relation.
 *
 * The combined set of may accesses may be a union if member accesses
 * are involved, but the entire set is derived from a single reference and
 * therefore from a single index expression.  These accesses therefore
 * all map to the same outer array.
 */
static int extract_access(__isl_keep pet_expr *expr, void *user)
{
    struct ppcg_extract_access_data *data = user;
    isl_union_map *tagged;
    struct amp_stmt_access *access;
    isl_ctx *ctx = pet_expr_get_ctx(expr);
    isl_multi_pw_aff *index;

    access = isl_alloc_type(ctx, struct amp_stmt_access);
    if (!access)
        return -1;
    access->next = NULL;
    access->read = pet_expr_access_is_read(expr);
    access->write = pet_expr_access_is_write(expr);
    tagged = pet_expr_access_get_tagged_may_read(expr);
    tagged = isl_union_map_union(tagged, pet_expr_access_get_tagged_may_write(expr));
    tagged = isl_union_map_apply_range(tagged, isl_union_map_copy(data->any_to_outer));
    if (!access->write)
    {
        access->exact_write = 1;
    }
    else if (!data->single_expression)
    {
        access->exact_write = 0;
    }
    else
    {
        isl_union_map *must, *may;
        may = isl_union_map_copy(tagged);
        may = isl_union_map_domain_factor_domain(may);
        must = pet_expr_access_get_must_write(expr);
        access->exact_write = isl_union_map_is_equal(must, may);
        isl_union_map_free(must);
        isl_union_map_free(may);
    }
    index = pet_expr_access_get_index(expr);
    access->n_index = isl_multi_pw_aff_dim(index, isl_dim_out);
    isl_multi_pw_aff_free(index);
    access->ref_id = pet_expr_access_get_ref_id(expr);
    access->tagged_access = extract_single_tagged_access(tagged, expr);
    access->access = isl_map_copy(access->tagged_access);
    access->access = isl_map_domain_factor_domain(access->access);
    access->fixed_element = accesses_fixed_element(expr);

    *data->next_access = access;
    data->next_access = &(*data->next_access)->next;

    if (!access->access || access->fixed_element < 0)
        return -1;

    return 0;
}

/* Construct a linked list of amp_stmt_access objects,
 * one for each access expression in the statement body.
 * "any_to_outer" maps all intermediate arrays to their outer arrays.
 */
static int pet_stmt_extract_accesses(struct amp_stmt *stmt, __isl_keep isl_union_map *any_to_outer)
{
    struct ppcg_extract_access_data data;

    stmt->accesses = NULL;
    data.next_access = &stmt->accesses;
    data.single_expression = pet_tree_get_type(stmt->stmt->body) == pet_tree_expr;
    data.any_to_outer = any_to_outer;
    return pet_tree_foreach_access_expr(stmt->stmt->body, &extract_access, &data);
}

static void *free_stmts(struct amp_stmt *stmts, int n)
{
    int i;

    if (!stmts)
        return NULL;

    for (i = 0; i < n; ++i)
    {
        struct amp_stmt_access *access, *next;

        for (access = stmts[i].accesses; access; access = next)
        {
            next = access->next;
            isl_id_free(access->ref_id);
            isl_map_free(access->access);
            isl_map_free(access->tagged_access);
            free(access);
        }

        isl_id_free(stmts[i].id);
    }
    free(stmts);

    return NULL;
}

/* Return an array of amp_stmt representing the statements in "scop".
 * Do not collect array accesses for statements that have been killed.
 */
static struct amp_stmt *extract_stmts(isl_ctx *ctx, struct ppcg_scop *scop,
                                      __isl_keep isl_union_map *any_to_outer)
{
    int i;
    struct amp_stmt *stmts;

    stmts = isl_calloc_array(ctx, struct amp_stmt, scop->pet->n_stmt);
    if (!stmts)
        return NULL;

    for (i = 0; i < scop->pet->n_stmt; ++i)
    {
        struct amp_stmt *s = &stmts[i];
        isl_bool killed;

        s->id = isl_set_get_tuple_id(scop->pet->stmts[i]->domain);
        s->stmt = scop->pet->stmts[i];
        killed = is_stmt_killed(scop, scop->pet->stmts[i]);
        if (killed < 0)
            return free_stmts(stmts, i + 1);
        if (killed)
            continue;
        if (pet_stmt_extract_accesses(s, any_to_outer) < 0)
            return free_stmts(stmts, i + 1);
    }

    return stmts;
}

/* Compute the set of inner array elements that may have their values
 * preserved by "prog".  In particular, collect the array elements of
 * arrays that are not local to "prog" and remove those elements that
 * are definitely killed or definitely written by "prog".
 */
static __isl_give isl_union_set *compute_may_persist(amp_prog *prog)
{
    int i;
    isl_union_set *may_persist, *killed;
    isl_union_map *must_kill;

    may_persist = isl_union_set_empty(isl_set_get_space(prog->context));
    for (i = 0; i < prog->n_array; ++i)
    {
        isl_set *extent;

        if (prog->array[i].local)
            continue;

        extent = isl_set_copy(prog->array[i].extent);
        may_persist = isl_union_set_add_set(may_persist, extent);
    }

    may_persist = isl_union_set_intersect_params(may_persist, isl_set_copy(prog->context));
    may_persist = isl_union_set_apply(may_persist, isl_union_map_copy(prog->to_inner));
    must_kill = isl_union_map_copy(prog->tagged_must_kill);
    killed = isl_union_map_range(must_kill);
    must_kill = isl_union_map_copy(prog->must_write);
    killed = isl_union_set_union(killed, isl_union_map_range(must_kill));

    may_persist = isl_union_set_subtract(may_persist, killed);
    return may_persist;
}

void *amp_prog_free(amp_prog *prog)
{
    if (!prog)
        return NULL;
    free_array_info(prog);
    free_stmts(prog->stmts, prog->n_stmts);
    isl_union_map_free(prog->any_to_outer);
    isl_union_map_free(prog->to_outer);
    isl_union_map_free(prog->to_inner);
    isl_union_map_free(prog->read);
    isl_union_map_free(prog->may_write);
    isl_union_map_free(prog->must_write);
    isl_union_map_free(prog->tagged_must_kill);
    isl_union_map_free(prog->array_order);
    isl_union_set_free(prog->may_persist);
    isl_set_free(prog->context);
    free(prog);

    return NULL;
}

amp_prog *amp_prog_alloc(__isl_take isl_ctx *ctx, struct ppcg_scop *scop)
{
    amp_prog *prog;
    isl_space *space;
    isl_map *id;

    if (!scop)
        return NULL;

    prog = isl_calloc_type(ctx, amp_prog);
    if (!prog)
        return NULL;

    prog->ctx = ctx;
    prog->scop = scop;
    prog->kernel_id = 0;
    prog->context = isl_set_copy(scop->context);
    prog->n_stmts = scop->pet->n_stmt;
    prog->any_to_outer = pet_scop_compute_outer_to_any(scop->pet);
    prog->any_to_outer = isl_union_map_reverse(prog->any_to_outer);
    space = isl_union_map_get_space(prog->any_to_outer);
    space = isl_space_set_from_params(space);
    space = isl_space_add_dims(space, isl_dim_set, 1);
    space = isl_space_map_from_set(space);
    id = isl_map_identity(space);
    prog->any_to_outer = isl_union_map_add_map(prog->any_to_outer, id);
    prog->stmts = extract_stmts(ctx, scop, prog->any_to_outer);
    prog->read = isl_union_map_copy(scop->reads);
    prog->may_write = isl_union_map_copy(scop->may_writes);
    prog->must_write = isl_union_map_copy(scop->must_writes);
    prog->tagged_must_kill = isl_union_map_copy(scop->tagged_must_kills);
    prog->to_inner = pet_scop_compute_outer_to_inner(scop->pet);
    prog->to_outer = isl_union_map_copy(prog->to_inner);
    prog->to_outer = isl_union_map_reverse(prog->to_outer);

    if (!prog->stmts)
        return amp_prog_free(prog);

    if (collect_array_info(prog) < 0)
        return amp_prog_free(prog);

    prog->may_persist = compute_may_persist(prog);

    return prog;
}

/** 获得比当前精度更低的数据类型 **/
char *amp_get_lower_precision_type(char *type)
{
    // 生成更低精度的
    if (strcmp("double", type) == 0)
        return "float";
    else if (strcmp("float", type) == 0)
        return "int";
    else if (strcmp("int", type) == 0)
        return "short int";
    else
        return type;
}

/* Does "array" need to be allocated on the device?
 * If it is a read-only scalar, then it will be passed as an argument
 * to the kernel and therefore does not require any allocation.
 * If this device memory is not accessed at all, then it does not
 * need to be allocated either.
 */
static int amp_array_requires_allocation(struct amp_array_info *array)
{
    // if (amp_array_is_read_only_scalar(array))
    //     return 0;
    return 1;
}

/* Build AST expressions for the amp array sizes of all arrays in "prog"
 * that require allocation on the device using "build", as well as
 * for the original array sizes of all arrays that need to be declared
 * on the host.
 * "node" is freed in case of error.
 */
__isl_give isl_ast_node *amp_build_array_bounds(__isl_take isl_ast_node *node, amp_prog *prog, __isl_keep isl_ast_build *build)
{
    int i;

    for (i = 0; i < prog->n_array; ++i)
    {
        struct amp_array_info *array = &prog->array[i];
        isl_multi_pw_aff *size;
        isl_ast_expr *expr;

        // if (!amp_array_requires_allocation(array))
        //     continue;

        size = isl_multi_pw_aff_copy(array->bound);
        expr = ppcg_build_size_expr(size, build);
        array->bound_expr = expr;
        if (!expr)
            return isl_ast_node_free(node);
    }

    for (i = 0; i < prog->n_array; ++i)
    {
        struct amp_array_info *array = &prog->array[i];
        isl_set *extent;
        isl_multi_pw_aff *size;
        isl_ast_expr *expr;

        if (!array->declare_local)
            continue;
        extent = isl_set_copy(array->declared_extent);
        size = ppcg_size_from_extent(extent);
        expr = ppcg_build_size_expr(size, build);
        array->declared_size = expr;
        if (!expr)
            return isl_ast_node_free(node);
    }

    return node;
}

/* Is "node" a mark node with an identifier called "name"?
 */
static int is_marked(__isl_keep isl_schedule_node *node, const char *name)
{
    isl_id *mark;
    int has_name;

    if (!node)
        return -1;

    if (isl_schedule_node_get_type(node) != isl_schedule_node_mark)
        return 0;

    mark = isl_schedule_node_mark_get_id(node);
    if (!mark)
        return -1;

    has_name = !strcmp(isl_id_get_name(mark), name);
    isl_id_free(mark);

    return has_name;
}

/* Is "node" a mark node with an identifier called "kernel"?
 */
int amp_tree_node_is_kernel(__isl_keep isl_schedule_node *node)
{
    return is_marked(node, "amp_kernel");
}

/* Is "node" a mark node with an identifier called "amp_higher"?
 */
static int node_is_amp_higher(__isl_keep isl_schedule_node *node)
{
    return is_marked(node, "amp_higher");
}

/* Is "node" a mark node with an identifier called "amp_lower"?
 */
static int node_is_amp_lower(__isl_keep isl_schedule_node *node)
{
    return is_marked(node, "amp_lower");
}

/* Is "node" a mark node with an identifier called "shared"?
 */
static int node_is_shared(__isl_keep isl_schedule_node *node)
{
    return is_marked(node, "shared");
}

/* Is "node" a mark node with an identifier called "thread"?
 * in the amp_kernel, 'thread' mark means an atomic calculation, which without modification.
 */
static int node_is_thread(__isl_keep isl_schedule_node *node)
{
    return is_marked(node, "thread");
}

/* Should this array reference group be mapped to private, shared or global
 * memory?
 * If we have computed both a private and a shared tile, then
 * the tile with the smallest depth is used.  If both have the same depth,
 * then the private tile is used.
 */
enum ppcg_group_access_type amp_array_ref_group_type(struct amp_array_ref_group *group)
{
    if (group->shared_tile)
        return ppcg_access_shared;
    return ppcg_access_global;
}

/* Print the name of the local copy of a given group of array references.
 */
__isl_give isl_printer *amp_array_ref_group_print_name(
    struct amp_array_ref_group *group, __isl_take isl_printer *p)
{
    int global = 0;
    enum ppcg_group_access_type type;

    type = amp_array_ref_group_type(group);
    if (type == ppcg_access_shared)
        p = isl_printer_print_str(p, "amp_lower_");
    else
    {
        global = 1;
        p = isl_printer_print_str(p, "amp_lower_");
    }

    p = isl_printer_print_str(p, group->array->name);
    if (!global && group->local_array->n_group > 1)
    {
        p = isl_printer_print_str(p, "_");
        p = isl_printer_print_int(p, group->nr);
    }

    return p;
}

/* Insert a mark node with identifier "shared" in front of "node".
 */
static __isl_give isl_schedule_node *insert_shared(
    __isl_take isl_schedule_node *node)
{
    isl_ctx *ctx;
    isl_id *id;

    ctx = isl_schedule_node_get_ctx(node);
    id = isl_id_alloc(ctx, "shared", NULL);
    node = isl_schedule_node_insert_mark(node, id);

    return node;
}

/* Insert a mark node with identifier "amp_lower" in front of "node".
 */
static __isl_give isl_schedule_node *insert_amp_lower(
    __isl_take isl_schedule_node *node)
{
    isl_ctx *ctx;
    isl_id *id;

    ctx = isl_schedule_node_get_ctx(node);
    id = isl_id_alloc(ctx, "amp_lower", NULL);
    node = isl_schedule_node_insert_mark(node, id);

    return node;
}

/* Insert a mark node with identifier "amp_higher" in front of "node".
 */
static __isl_give isl_schedule_node *insert_amp_higher(
    __isl_take isl_schedule_node *node)
{
    isl_ctx *ctx;
    isl_id *id;

    ctx = isl_schedule_node_get_ctx(node);
    id = isl_id_alloc(ctx, "amp_higher", NULL);
    node = isl_schedule_node_insert_mark(node, id);

    return node;
}

/* Mark all dimensions in the current band node atomic.
 */
static __isl_give isl_schedule_node *atomic(__isl_take isl_schedule_node *node)
{
    return ppcg_set_schedule_node_type(node, isl_ast_loop_atomic);
}

/* Mark "node" atomic, if it is a band node.
 * Do the same for all ancestors.
 * Return a pointer to "node" (in the updated schedule tree).
 */
static __isl_give isl_schedule_node *atomic_ancestors(
    __isl_take isl_schedule_node *node)
{
    // #define DEBUG_ATOMIC_ANCESTORS

#ifdef DEBUG_ATOMIC_ANCESTORS
    fprintf(stderr, "@DEBUG: \n       at the start of the atomic_ancestors function,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ATOMIC_ANCESTORS
    int pos;

    if (!node)
        return NULL;
    if (!isl_schedule_node_has_parent(node))
        return node;

    pos = isl_schedule_node_get_child_position(node);
    node = isl_schedule_node_parent(node);
    if (isl_schedule_node_get_type(node) == isl_schedule_node_band)
        node = atomic(node);
    node = atomic_ancestors(node);
    node = isl_schedule_node_child(node, pos);

#ifdef DEBUG_ATOMIC_ANCESTORS
    fprintf(stderr, "@DEBUG: \n       at the end of the atomic_ancestors function,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ATOMIC_ANCESTORS

    return node;
}

/* Insert a "shared" mark in front of the "thread" mark
 * provided the linear branch between "node" and the "thread" mark
 * does not contain such a "shared" mark already.
 *
 * As a side effect, this function checks that the subtree at "node"
 * actually contains a "thread" mark and that there is no branching
 * in between "node" and this "thread" mark.
 */
__isl_give isl_schedule_node *amp_tree_insert_shared_before_thread(
    __isl_take isl_schedule_node *node)
{
    int depth0, depth;
    int any_shared = 0;

    if (!node)
        return NULL;

    depth0 = isl_schedule_node_get_tree_depth(node);

    for (;;)
    {
        int is_thread;
        int n;

        if (!any_shared)
        {
            any_shared = node_is_shared(node);
            if (any_shared < 0)
                return isl_schedule_node_free(node);
        }
        is_thread = node_is_thread(node);
        if (is_thread < 0)
            return isl_schedule_node_free(node);
        if (is_thread)
            break;
        n = isl_schedule_node_n_children(node);
        if (n == 0)
            isl_die(isl_schedule_node_get_ctx(node),
                    isl_error_invalid,
                    "no thread marker found",
                    return isl_schedule_node_free(node));
        if (n > 1)
            isl_die(isl_schedule_node_get_ctx(node),
                    isl_error_invalid,
                    "expecting single thread marker",
                    return isl_schedule_node_free(node));

        node = isl_schedule_node_child(node, 0);
    }

    if (!any_shared)
        node = insert_shared(node);
    depth = isl_schedule_node_get_tree_depth(node);
    node = isl_schedule_node_ancestor(node, depth - depth0);

    return node;
}

/* Assuming "node" is a filter node, does it correspond to the branch
 * that contains the "thread" mark, i.e., does it contain any elements
 * in "core"?
 */
static int node_is_core(__isl_keep isl_schedule_node *node,
                        __isl_keep isl_union_set *core)
{
    int disjoint;
    isl_union_set *filter;

    filter = isl_schedule_node_filter_get_filter(node);
    disjoint = isl_union_set_is_disjoint(filter, core);
    isl_union_set_free(filter);
    if (disjoint < 0)
        return -1;

    return !disjoint;
}

/* Move to the only child of "node" that has the "thread" mark as descendant,
 * where the branch containing this mark is identified by the domain elements
 * in "core".
 *
 * If "node" is not a sequence, then it only has one child and we move
 * to that single child.
 * Otherwise, we check each of the filters in the children, pick
 * the one that corresponds to "core" and return a pointer to the child
 * of the filter node.
 */
static __isl_give isl_schedule_node *core_child(
    __isl_take isl_schedule_node *node, __isl_keep isl_union_set *core)
{
    int i, n;

    if (isl_schedule_node_get_type(node) != isl_schedule_node_sequence)
        return isl_schedule_node_child(node, 0);

    n = isl_schedule_node_n_children(node);
    for (i = 0; i < n; ++i)
    {
        int is_core;

        node = isl_schedule_node_child(node, i);
        is_core = node_is_core(node, core);

        if (is_core < 0)
            return isl_schedule_node_free(node);
        if (is_core)
            return isl_schedule_node_child(node, 0);

        node = isl_schedule_node_parent(node);
    }

    isl_die(isl_schedule_node_get_ctx(node), isl_error_internal,
            "core child not found", return isl_schedule_node_free(node));
}

/* Move down the branch between "kernel" and "thread" until
 * the "shared" mark is reached, where the branch containing the "shared"
 * mark is identified by the domain elements in "core".
 */
__isl_give isl_schedule_node *amp_tree_move_down_to_shared(
    __isl_take isl_schedule_node *node, __isl_keep isl_union_set *core)
{
    int is_shared;

    while ((is_shared = node_is_shared(node)) == 0)
        node = core_child(node, core);
    if (is_shared < 0)
        node = isl_schedule_node_free(node);

    return node;
}

/* Move down the branch between "kernel" and "thread" until
 * the "thread" mark is reached, where the branch containing the "thread"
 * mark is identified by the domain elements in "core".
 */
__isl_give isl_schedule_node *amp_tree_move_down_to_thread(
    __isl_take isl_schedule_node *node, __isl_keep isl_union_set *core)
{
    int is_thread;

    while ((is_thread = node_is_thread(node)) == 0)
        node = core_child(node, core);
    if (is_thread < 0)
        node = isl_schedule_node_free(node);

    return node;
}

/* Move up the tree underneath the "thread" mark until
 * the "thread" mark is reached.
 */
__isl_give isl_schedule_node *amp_tree_move_up_to_thread(
    __isl_take isl_schedule_node *node)
{
    int is_thread;

    while ((is_thread = node_is_thread(node)) == 0)
        node = isl_schedule_node_parent(node);
    if (is_thread < 0)
        node = isl_schedule_node_free(node);

    return node;
}

/* Move up the tree underneath the "kernel" mark until
 * the "kernel" mark is reached.
 */
__isl_give isl_schedule_node *amp_tree_move_up_to_kernel(
    __isl_take isl_schedule_node *node)
{
    int is_kernel;

    while ((is_kernel = amp_tree_node_is_kernel(node)) == 0)
        node = isl_schedule_node_parent(node);
    if (is_kernel < 0)
        node = isl_schedule_node_free(node);

    return node;
}

/* Move down from the "kernel" mark (or at least a node with schedule
 * depth smaller than or equal to "depth") to a band node at schedule
 * depth "depth".  The "thread" mark is assumed to have a schedule
 * depth greater than or equal to "depth".  The branch containing the
 * "thread" mark is identified by the domain elements in "core".
 *
 * If the desired schedule depth is in the middle of band node,
 * then the band node is split into two pieces, the second piece
 * at the desired schedule depth.
 */
__isl_give isl_schedule_node *amp_tree_move_down_to_depth(
    __isl_take isl_schedule_node *node, int depth,
    __isl_keep isl_union_set *core)
{

    while (node && isl_schedule_node_get_schedule_depth(node) < depth)
    {
        if (isl_schedule_node_get_type(node) == isl_schedule_node_band)
        {
            int node_depth, node_dim;
            node_depth = isl_schedule_node_get_schedule_depth(node);
            node_dim = isl_schedule_node_band_n_member(node);
            if (node_depth + node_dim > depth)
                node = isl_schedule_node_band_split(node,
                                                    depth - node_depth);
        }
        node = core_child(node, core);
    }

    return node;
}

/* If group->n_ref == 1, then group->refs was set by
 * populate_array_references to point directly into
 * group->array->refs and should not be freed.
 * If group->n_ref > 1, then group->refs was set by join_groups
 * to point to a newly allocated array.
 */
struct amp_array_ref_group *amp_array_ref_group_free(
    struct amp_array_ref_group *group)
{
    if (!group)
        return NULL;

    isl_map_free(group->access);
    if (group->n_ref > 1)
        free(group->refs);
    free(group);
    return NULL;
}

struct amp_ppcg_kernel *amp_ppcg_kernel_free(struct amp_ppcg_kernel *kernel)
{
    int i, j;

    if (!kernel)
        return NULL;
    isl_set_free(kernel->context);
    isl_union_set_free(kernel->core);
    isl_union_set_free(kernel->arrays);
    isl_union_pw_multi_aff_free(kernel->contraction);
    isl_union_set_free(kernel->expanded_domain);
    isl_space_free(kernel->space);
    isl_ast_node_free(kernel->tree);

    isl_union_pw_multi_aff_free(kernel->copy_schedule);

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *array = &kernel->array[i];

        for (j = 0; j < array->n_group; ++j)
            amp_array_ref_group_free(array->groups[j]);
        free(array->groups);

        isl_multi_pw_aff_free(array->bound);
        isl_ast_expr_free(array->bound_expr);
    }
    free(kernel->array);

    for (i = 0; i < kernel->n_var; ++i)
    {
        free(kernel->var[i].name);
        isl_vec_free(kernel->var[i].size);
    }
    free(kernel->var);

    free(kernel);

    return NULL;
}

/* Create the array of gpu_local_array_info structures "array"
 * inside "kernel".  The number of elements in this array is
 * the same as the number of arrays in "prog".
 * Initialize the "array" field of each local array to point
 * to the corresponding array in "prog".
 */
static struct amp_ppcg_kernel *amp_ppcg_kernel_create_local_arrays(
    struct amp_ppcg_kernel *kernel, struct amp_prog *prog)
{
    int i;
    isl_ctx *ctx;

    if (!kernel)
        return NULL;

    ctx = isl_set_get_ctx(prog->context);
    kernel->array = isl_calloc_array(ctx, struct amp_local_array_info, prog->n_array);
    if (!kernel->array)
        return amp_ppcg_kernel_free(kernel);
    kernel->n_array = prog->n_array;

    for (i = 0; i < prog->n_array; ++i)
        kernel->array[i].array = &prog->array[i];

    return kernel;
}

/* Return the union of all read (read = 1) and/or write (write = 1)
 * access relations in the group.
 */
__isl_give isl_union_map *amp_array_ref_group_access_relation(
    struct amp_array_ref_group *group, int read, int write)
{
    int i;
    isl_union_map *access;

    access = isl_union_map_empty(isl_map_get_space(group->access));
    for (i = 0; i < group->n_ref; ++i)
    {
        isl_map *map_i;

        if (!((read && group->refs[i]->read) ||
              (write && group->refs[i]->write)))
            continue;
        map_i = isl_map_copy(group->refs[i]->access);
        access = isl_union_map_union(access, isl_union_map_from_map(map_i));
    }

    return access;
}

/* Return the union of all tagged access relations in the group.
 */
static __isl_give isl_union_map *amp_group_tagged_access_relation(
    struct amp_array_ref_group *group)
{
    int i;
    isl_union_map *access;

    access = isl_union_map_empty(isl_map_get_space(group->access));
    for (i = 0; i < group->n_ref; ++i)
    {
        isl_map *map_i;

        map_i = isl_map_copy(group->refs[i]->tagged_access);
        access = isl_union_map_union(access, isl_union_map_from_map(map_i));
    }

    return access;
}

/* Given a set of wrapped references "ref", return the corresponding
 * access relations based on the tagged access relations "tagged".
 *
 * The elements of "ref" are of the form
 *
 *	[D -> R]
 *
 * with D an iteration domains and R a reference.
 * The elements of "tagged" are of the form
 *
 *	[D -> R] -> A
 *
 * with A an array.
 *
 * Extend "tagged" to include the iteration domain in the range, i.e.,
 *
 *	[D -> R] -> [D -> A]
 *
 * apply the result to "ref" and then unwrap the resulting set
 * to obtain relations of the form
 *
 *	D -> A
 */
static __isl_give isl_union_map *wrapped_reference_to_access(
    __isl_take isl_union_set *ref, __isl_take isl_union_map *tagged)
{
    isl_union_map *tag2access;

    tag2access = isl_union_map_copy(tagged);
    tag2access = isl_union_map_universe(tag2access);
    tag2access = isl_union_set_unwrap(isl_union_map_domain(tag2access));
    tag2access = isl_union_map_domain_map(tag2access);
    tag2access = isl_union_map_range_product(tag2access, tagged);

    ref = isl_union_set_coalesce(ref);
    ref = isl_union_set_apply(ref, tag2access);

    return isl_union_set_unwrap(ref);
}

/* Given an access relation "access" from one or more array reference groups,
 * remove those reads if ("read" is 1) or writes (if "read" is 0)
 * that are only needed to communicate data within
 * the same iteration of "sched".
 * The domain of "sched" corresponds to the original statement instances,
 * i.e., those that appear in the domains of the access relations.
 * "tagged" contains all tagged access relations to all
 * the array reference groups accessed by "access" from statement
 * instances scheduled by "sched".
 *
 * If the access is a read then it is either an element of
 *
 *	live_in union (range flow)
 *
 * where live_in and flow may be overapproximations, or
 * it reads an uninitialized value (that is not live-in because
 * there is an intermediate kill) or it reads a value that was
 * written within the same (compound) statement instance.
 * If the access is a write then it is either an element of
 *
 *	live_out union (domain flow)
 *
 * or it writes a value that is never read (and is not live-out
 * because of an intermediate kill) or only
 * within the same (compound) statement instance.
 * In both cases, the access relation is also a subset of
 * the group access relation.
 *
 * The cases where an uninitialized value is read or a value is written
 * that is never read or where the dataflow occurs within a statement
 * instance are also considered local and may also be removed.
 *
 * Essentially, we compute the intersection of "access" with either
 *
 *	live_in union (range non-local-flow)
 *
 * or
 *
 *	live_out union (domain non-local-flow)
 *
 * We first construct a relation "local"
 *
 *	[[D -> R] -> [D' -> R']]
 *
 * of pairs of domain iterations accessing the reference group
 * and references in the group that are coscheduled by "sched".
 *
 * If this relation does not intersect the dataflow dependences,
 * then there is nothing we can possibly remove, unless the dataflow
 * dependences themselves only relate a subset of the accesses.
 * In particular, the accesses may not be involved in any dataflow
 * dependences, either because they are uninitialized reads/dead writes
 * or because the dataflow occurs inside a statement instance.
 *
 * Since the computation below may break up the access relation
 * into smaller pieces, we only perform the intersection with
 * the non-local dependent accesses if the local pairs
 * intersect the dataflow dependences.  Otherwise, we intersect
 * with the universe of the non-local dependent accesses.
 * This should at least remove accesses from statements that
 * do not participate in any dependences.
 *
 * In particular, we remove the "local" dataflow dependences from
 * the set of all dataflow dependences, or at least those
 * that may contribute to a domain/range that intersects
 * the domain of "access".
 * Note that if the potential dataflow dependences are an overapproximation
 * of the actual dataflow dependences, then the result remains an
 * overapproximation of the non-local dataflow dependences.
 * Copying to/from global memory is only needed for the references
 * in the domain/range of the result or for accesses that are live out/in
 * for the entire scop.
 *
 * We therefore map the domain/range of the "external" relation
 * to the corresponding access relation and take the union with
 * the live out/in relation.
 */
static __isl_give isl_union_map *amp_remove_local_accesses(
    struct amp_prog *prog, __isl_take isl_union_map *tagged,
    __isl_take isl_union_map *access, __isl_take isl_union_map *sched,
    int read)
{
    // #define DEBUG_AMP_REMOVE_LOCAL_ACCESSES

    int empty;
    isl_union_pw_multi_aff *tagger;
    isl_union_set *domain, *access_domain;
    isl_union_map *local, *external, *universe;
    isl_union_set *tag_set;
    struct ppcg_scop *ps = prog->scop;

    if (isl_union_map_is_empty(access))
    {
        isl_union_map_free(sched);
        isl_union_map_free(tagged);
        return access;
    }

#ifdef DEBUG_AMP_REMOVE_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       in start of amp_remove_local_accesses(no group), the access is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_REMOVE_LOCAL_ACCESSES

    // compute_live_out(ps);

    tagger = isl_union_pw_multi_aff_copy(ps->tagger);
    domain = isl_union_map_domain(isl_union_map_copy(tagged));
    tagger = isl_union_pw_multi_aff_intersect_domain(tagger, isl_union_set_copy(domain));
    sched = isl_union_map_preimage_domain_union_pw_multi_aff(sched, tagger);

    local = isl_union_map_apply_range(sched, isl_union_map_reverse(isl_union_map_copy(sched)));
    local = isl_union_map_intersect(local, isl_union_map_copy(ps->tagged_dep_flow));

    empty = isl_union_map_is_empty(local);

    external = isl_union_map_copy(ps->tagged_dep_flow);
    universe = isl_union_map_universe(isl_union_map_copy(access));
    access_domain = isl_union_map_domain(universe);
    domain = isl_union_set_universe(domain);
    universe = isl_union_set_unwrap(domain);
    universe = isl_union_map_intersect_domain(universe, access_domain);
    domain = isl_union_map_wrap(universe);
    if (read)
        external = isl_union_map_intersect_range(external, domain);
    else
        external = isl_union_map_intersect_domain(external, domain);
    external = isl_union_map_intersect_params(external, isl_set_copy(ps->context));
    external = isl_union_map_subtract(external, local);

#ifdef DEBUG_AMP_REMOVE_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       in the prog->scop, the live_in and live_out is on the below :\n");
    isl_union_map_dump(ps->live_in);
    isl_union_map_dump(ps->live_out);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_REMOVE_LOCAL_ACCESSES

    if (read)
    {
        tag_set = isl_union_map_range(external);
        external = wrapped_reference_to_access(tag_set, tagged);
        external = isl_union_map_union(external, isl_union_map_copy(ps->live_in));
    }
    else
    {
        tag_set = isl_union_map_domain(external);
        external = wrapped_reference_to_access(tag_set, tagged);
        external = isl_union_map_union(external, isl_union_map_copy(ps->live_out));
    }

    if (empty < 0)
        external = isl_union_map_free(external);
    else if (empty)
        external = isl_union_map_universe(external);

    access = isl_union_map_intersect(access, external);

    return access;
}

/* Given an access relation "access" from "group", remove those reads
 * if ("read" is 1) or writes (if "read" is 0) that are only needed to
 * communicate data within the same iteration of the schedule "prefix"
 * at the position where the copying of the group is inserted.
 * That is, the output dimension of "prefix"
 * is equal to tile->depth.
 * The domain of "prefix" corresponds to the original statement instances,
 * i.e., those that appear in the domains of the access relations.
 *
 * Extract the tagged access relation of "group" and
 * then call remove_local_accesses.
 */
static __isl_give isl_union_map *amp_remove_local_accesses_group(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_union_map *access, __isl_keep isl_union_map *prefix,
    int read)
{
    // #define DEBUG_AMP_REMOVE_LOCAL_ACCESSES_GROUP
    isl_union_map *sched, *tagged;

    if (isl_union_map_is_empty(access))
        return access;

#ifdef DEBUG_AMP_REMOVE_LOCAL_ACCESSES_GROUP
    fprintf(stderr, "@DEBUG: \n       in start of amp_remove_local_accesses_group, the access is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif
    tagged = amp_group_tagged_access_relation(group);
    sched = isl_union_map_copy(prefix);

#ifdef DEBUG_AMP_REMOVE_LOCAL_ACCESSES_GROUP
    fprintf(stderr, "@DEBUG: \n       the tagged is:\n");
    isl_union_map_dump(tagged);
    fprintf(stderr, "\n       the sched is:\n");
    isl_union_map_dump(sched);
    fprintf(stderr, "\n\n");
#endif

    return amp_remove_local_accesses(kernel->prog, tagged, access, sched, read);
}

/* Return the effective gpu_array_tile associated to "group" or
 * NULL if there is no such gpu_array_tile.
 */
struct amp_array_tile *amp_array_ref_group_tile(struct amp_array_ref_group *group)
{
    switch (amp_array_ref_group_type(group))
    {
    case ppcg_access_shared:
        return group->shared_tile;
    case ppcg_access_global:
        return NULL;
    }
}

/* Return a read ("read" is 1) or write access relation for "group"
 * with those accesses removed that are only needed to communicate data
 * within the subtree of the schedule rooted at "node".
 * Furthermore, include the prefix schedule at "node".
 * That is, return a relation of the form
 *
 *	S -> [D -> A]
 *
 * with D the outer schedule dimensions at "node".
 */
static __isl_give isl_union_map *anchored_non_local_accesses(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
    // #define DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    isl_union_map *access;
    isl_union_map *prefix;

#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       anchored_non_local_accesses, the node is:\n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n       the read's value is: %d \n\n", read);
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    prefix = isl_schedule_node_get_prefix_schedule_relation(node);
#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       the prefix from current node is:\n");
    isl_union_map_dump(prefix);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    prefix = isl_union_map_preimage_domain_union_pw_multi_aff(prefix, isl_union_pw_multi_aff_copy(kernel->contraction));
#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       after insert prefix(from current node) into the kernel->contraction,the result is:\n");
    isl_union_map_dump(prefix);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    access = amp_array_ref_group_access_relation(group, read, !read);
#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       all read (read = 1) and/or write (write = 1) access relations in the group, is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n       the read's value is: %d \n\n", read);
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    access = amp_remove_local_accesses_group(kernel, group, access, prefix, read);
#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       after amp_remove_local_accesses_group, the access is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    access = isl_union_map_range_product(prefix, access);
#ifdef DEBUG_ANCHORED_NON_LOCAL_ACCESSES
    fprintf(stderr, "@DEBUG: \n       after isl_union_map_range_product function, the returned access is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_ANCHORED_NON_LOCAL_ACCESSES

    return access;
}

/* Given an array reference group "group", create a mapping
 *
 *	read[D -> A] -> [D -> A]
 *
 * if "read" is set or
 *
 *	write[D -> A] -> [D -> A]
 *
 * if "read" is not set.
 * D corresponds to the outer tile->depth dimensions of
 * the kernel schedule.
 */
static __isl_give isl_multi_aff *create_from_access(isl_ctx *ctx, struct amp_array_ref_group *group, int read)
{
    isl_space *space;
    isl_id *id;

    space = isl_space_copy(group->array->space);
    space = isl_space_from_range(space);
    space = isl_space_wrap(space);
    space = isl_space_map_from_set(space);

    id = isl_id_alloc(ctx, read ? "read" : "write", group);
    space = isl_space_set_tuple_id(space, isl_dim_in, id);

    return isl_multi_aff_identity(space);
}

/* Add copy statements to the schedule tree of "node"
 * for reading from global memory to shared memory (if "read" is set) or
 * for writing back from shared memory to global memory
 * (if "read" is not set) for the array reference group "group" that
 * is mapped to shared memory.
 * On input, "node" points to the kernel node, and it is moved
 * back there on output.
 *
 * The copies are performed in the order of the corresponding shared
 * memory tile.
 * The copy statement instances include a reference to the outer
 * tile->depth dimensions of the kernel schedule for ease of
 * combining them with the group tiling.
 *
 * If we are performing a read from global memory to shared memory and
 * if the array involved is not a scalar, then we copy
 * the entire tile to shared memory.  This may result in some extra
 * elements getting copied, but it should lead to simpler code
 * (which means that fewer registers may be needed) and less divergence.
 *
 * Otherwise, we only copy the elements that will be read or have been written
 * in the kernel.
 *
 * That is, the extra schedule is of the form
 *
 *	type[D -> A] -> T
 *
 * where D corresponds to the outer tile->depth dimensions of
 * the kernel schedule, A to the global array and T is the corresponding
 * shared memory tile.
 *
 * The copying is inserted in the schedule tree through an extension
 * of the form
 *
 *	D -> type[D -> A]
 *
 * where the extra domain elements type[D -> A] are those accessed
 * by the group.  In the case of read from a non-scalar, this set
 * is replaced by the entire shared memory tile.
 *
 * If the "unroll_copy_shared" option is set, then the AST generator
 * is instructed to unroll the copying code.
 *
 * A filter is inserted on type[D -> A] to map the copy instances
 * to the threads.  In particular, the thread identifiers are
 * equated to the position inside the shared memory tile (T)
 * modulo the block size.
 * We try to align the innermost tile dimension with the innermost
 * thread identifier (x) as a heuristic to improve coalescing.
 * In particular, if the dimension of the tile is greater than
 * the dimension of the block, then the schedule mapping to the tile
 * is broken up into two pieces and the filter is applied to the inner part.
 * If, on the other hand, the dimension of the tile is smaller than
 * the dimension of the block, then the initial thread identifiers
 * are equated to zero and the remaining thread identifiers are
 * matched to the memory tile.
 *
 * The extension is inserted before the core computation in case of a read
 * and after the core computation in case of a write.
 * In the case of a read, we first need to make sure there is some
 * synchronization before the core computation such that we can put the read
 * from global memory to shared memory before that synchronization.
 * This ensures that all threads have finished copying into shared memory
 * before the shared memory is used.
 * We also need to make sure that there is a synchronization node after
 * the core computation to ensure that the next load into shared memory
 * only happens after all data has been used.  There is no need for
 * this synchronization if we are at the outer level since then there
 * won't be a next load.
 * In the case of a write, we need to make sure there is some synchronization
 * after the core computation such that we can put the write from shared
 * memory to global memory after that synchronization.
 * Unless we are at the outer level, we also need a synchronization node
 * after the write to ensure the data is saved to global memory
 * before the next iteration writes to the same shared memory.
 * It also makes sure the data has arrived in global memory before
 * it is read in a subsequent iteration.
 */
//  参考CUDA的add_copies_group_shared函数实现,同时又引入了和add_copies_group_private的对比

static __isl_give isl_schedule_node *amp_add_copies_group_shared(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
    // #define DEBUG_AMP_ADD_COPIES_GROUP

    struct amp_array_tile *tile;
    isl_union_map *access;
    isl_union_set *domain;
    isl_multi_aff *ma;
    isl_multi_aff *from_access;
    isl_multi_pw_aff *mpa;
    isl_multi_union_pw_aff *mupa;
    isl_schedule_node *graft;
    isl_union_set *filter;
    isl_space *space;
    int kernel_depth;
    int empty;

    // if (amp_array_is_scalar(group->array))
    // {
    //     return node;
    // }

    tile = amp_array_ref_group_tile(group);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_SHARED
    fprintf(stderr, "\n\n\n@DEBUG \n       the tile of group in amp_add_copies_group function is: \n");
    if (tile)
        amp_array_tile_dump(tile);
    fprintf(stderr, "\n\n\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       at start of the amp add copies group function, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n       the read is :%d \n\n", read);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    node = amp_tree_move_down_to_depth(node, tile->depth, kernel->core);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    access = anchored_non_local_accesses(kernel, group, node, read);
    empty = isl_union_map_is_empty(access);
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        if (empty < 0)
            return isl_schedule_node_free(node);
        return amp_tree_move_up_to_kernel(node);
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       after anchored_non_local_accesses, the access is: \n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       the from_access variable is: \n");
    isl_multi_aff_dump(from_access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    if (tile->tiling)
    {
        ma = isl_multi_aff_copy(tile->tiling);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       the ma(tile->tiling) is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

        ma = isl_multi_aff_pullback_multi_aff(ma, isl_multi_aff_copy(from_access));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       plug in from_access in ma, the ma is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    }

    else
    {
        ma = isl_multi_aff_copy(from_access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       the ma(from_access) is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    }

    mpa = isl_multi_pw_aff_from_multi_aff(ma);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       the mpa is: \n");
    isl_multi_pw_aff_dump(mpa);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);
    domain = isl_union_map_range(access);

    // 只对读生效,注释掉就正确了,后面再深究
    // if (read && !amp_array_is_scalar(group->array)) {
    // 	isl_map *map;
    // 	isl_union_set_free(domain);
    // 	map = group_tile(group);
    // 	domain = isl_union_set_from_set(isl_map_wrap(map));
    // }

    domain = isl_union_set_preimage_multi_aff(domain, isl_multi_aff_copy(from_access));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       the domain is: \n");
    isl_union_set_dump(domain);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    // 基本上从下面开始就没办法修改了
    access = isl_union_set_wrapped_domain_map(domain);
    access = isl_union_map_reverse(access);
    access = isl_union_map_coalesce(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       before the isl_schedule_node_from_extension(access), the access is: \n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    graft = isl_schedule_node_from_extension(access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       after the isl_schedule_node_from_extension(access), the graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, isl_multi_union_pw_aff_copy(mupa));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       after the isl_schedule_node_insert_partial_schedule(graft, mupa), the graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n       the mupa(isl_multi_union_pw_aff) is : \n");
    isl_multi_union_pw_aff_dump(mupa);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    if (kernel->options->unroll_copy_shared)
        graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    while (graft && isl_schedule_node_has_parent(graft))
        graft = isl_schedule_node_parent(graft);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       the final graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       before insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    if (read)
    {
        node = amp_tree_move_down_to_shared(node, kernel->core);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_before(node, graft);
    }
    else
    {
        node = amp_tree_move_down_to_shared(node, kernel->core);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_after(node, graft);
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    fprintf(stderr, "@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    return node;
}

static __isl_give isl_schedule_node *amp_add_copies_group_global(struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
                                                                 __isl_take isl_schedule_node *node, int read)
{
    // #define DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    isl_union_map *access;
    isl_union_set *domain;
    isl_multi_aff *ma;
    isl_multi_aff *from_access;
    isl_multi_pw_aff *mpa;
    isl_multi_union_pw_aff *mupa;
    isl_union_pw_multi_aff *contraction;
    isl_schedule_node *graft;
    isl_union_set *filter;
    isl_space *space;
    int kernel_depth;
    int empty;
    isl_printer *p;
    char *local_name;
    isl_multi_aff *tiling;

    if (!amp_array_is_scalar(group->array))
    {
        return node;
    }

    // tile = amp_array_ref_group_tile(group);
    // #ifdef DEBUG_AMP_ADD_COPIES_GROUP
    //     fprintf(stderr, "\n\n\n@DEBUG \n       the tile of group in amp_add_copies_group function is: \n");
    //     amp_array_tile_dump(tile);
    //     fprintf(stderr, "\n\n\n\n");
    // #endif // DEBUG_AMP_ADD_COPIES_GROUP

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       at start of the amp add copies group function, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n       the read is :%d \n\n", read);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    // node = amp_tree_move_down_to_depth(node, tile->depth, kernel->core);

    // #ifdef DEBUG_AMP_ADD_COPIES_GROUP
    //     fprintf(stderr, "@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    //     isl_schedule_node_dump(node);
    //     fprintf(stderr, "\n\n");
    // #endif // DEBUG_AMP_ADD_COPIES_GROUP

    access = anchored_non_local_accesses(kernel, group, node, read);
    empty = isl_union_map_is_empty(access);
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        if (empty < 0)
            return isl_schedule_node_free(node);
        return amp_tree_move_up_to_kernel(node);
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       after anchored_non_local_accesses, the access is: \n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       the from_access variable is: \n");
    isl_multi_aff_dump(from_access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    space = isl_map_get_space(group->access);
    space = isl_space_from_range(isl_space_range(space));
    space = isl_space_add_dims(space, isl_dim_in, 0);

    tiling = isl_multi_aff_range_map(isl_space_copy(space));

    p = isl_printer_to_str(isl_multi_aff_get_ctx(tiling));
    p = amp_array_ref_group_print_name(group, p);
    local_name = isl_printer_get_str(p);
    isl_printer_free(p);
    tiling = isl_multi_aff_set_tuple_name(tiling, isl_dim_out, local_name);
    free(local_name);

    if (tiling)
    {
        ma = isl_multi_aff_copy(tiling);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       the ma(tiling) is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

        ma = isl_multi_aff_pullback_multi_aff(ma, isl_multi_aff_copy(from_access));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       plug in from_access in ma, the ma is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    }
    else
    {
        ma = isl_multi_aff_copy(from_access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       the ma(from_access) is: \n");
        isl_multi_aff_dump(ma);
        fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    }

    mpa = isl_multi_pw_aff_from_multi_aff(ma);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       the mpa is: \n");
    isl_multi_pw_aff_dump(mpa);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);
    domain = isl_union_map_range(access);

    // 只对读生效,注释掉就正确了,后面再深究
    // if (read && !amp_array_is_scalar(group->array)) {
    // 	isl_map *map;
    // 	isl_union_set_free(domain);
    // 	map = group_tile(group);
    // 	domain = isl_union_set_from_set(isl_map_wrap(map));
    // }

    domain = isl_union_set_preimage_multi_aff(domain, isl_multi_aff_copy(from_access));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       the domain is: \n");
    isl_union_set_dump(domain);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    // 基本上从下面开始就没办法修改了
    access = isl_union_set_wrapped_domain_map(domain);
    access = isl_union_map_reverse(access);
    access = isl_union_map_coalesce(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       before the isl_schedule_node_from_extension(access), the access is: \n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    graft = isl_schedule_node_from_extension(access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       after the isl_schedule_node_from_extension(access), the graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, isl_multi_union_pw_aff_copy(mupa));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       after the isl_schedule_node_insert_partial_schedule(graft, mupa), the graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n       the mupa(isl_multi_union_pw_aff) is : \n");
    isl_multi_union_pw_aff_dump(mupa);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    // if (kernel->options->unroll_copy_shared)
    //     graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    while (graft && isl_schedule_node_has_parent(graft))
        graft = isl_schedule_node_parent(graft);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       the final graft is: \n");
    isl_schedule_node_dump(graft);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       before insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    if (read)
    {
        node = amp_tree_move_down_to_shared(node, kernel->core);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_before(node, graft);
    }
    else
    {
        node = amp_tree_move_down_to_shared(node, kernel->core);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_after(node, graft);
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL
    fprintf(stderr, "@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_GLOBAL

    return node;
}

/* Extract the set of parameter values and outer schedule dimensions
 * for which any statement instance
 * in the kernel inserted at "node" needs to be executed.
 * Intersect the set of parameter values derived from the host schedule
 * relation with the context of "prog".
 */
static __isl_give isl_set *extract_context(__isl_keep isl_schedule_node *node,
                                           struct amp_prog *prog)
{
    isl_union_map *schedule;
    isl_union_set *schedule_domain;
    isl_set *context;
    int empty;

    schedule = isl_schedule_node_get_prefix_schedule_relation(node);
    schedule_domain = isl_union_map_range(schedule);
    empty = isl_union_set_is_empty(schedule_domain);
    if (empty < 0)
    {
        isl_union_set_free(schedule_domain);
        return NULL;
    }
    if (empty)
    {
        int depth;
        isl_space *space;

        space = isl_union_set_get_space(schedule_domain);
        isl_union_set_free(schedule_domain);
        space = isl_space_set_from_params(space);
        depth = isl_schedule_node_get_schedule_depth(node);
        space = isl_space_add_dims(space, isl_dim_set, depth);
        context = isl_set_empty(space);
    }
    else
    {
        context = isl_set_from_union_set(schedule_domain);
    }
    context = isl_set_intersect_params(context,
                                       isl_set_copy(prog->context));

    return context;
}

/* Fill up the groups array with singleton groups, i.e., one group
 * per reference, initializing the array, access, write, n_ref and refs fields.
 * In particular the access field is initialized to the scheduled
 * access relation of the array reference.
 *
 * Return the number of elements initialized, i.e., the number of
 * active references in the current kernel.
 */
static int populate_array_references(struct amp_local_array_info *local,
                                     struct amp_array_ref_group **groups, struct amp_group_data *data)
{
    int i;
    int n;
    isl_ctx *ctx = isl_union_map_get_ctx(data->copy_sched);

    n = 0;
    for (i = 0; i < local->array->n_ref; ++i)
    {
        isl_union_map *umap;
        isl_map *map;
        struct amp_array_ref_group *group;
        struct amp_stmt_access *access = local->array->refs[i];

        map = isl_map_copy(access->access);
        umap = isl_union_map_from_map(map);
        umap = isl_union_map_apply_domain(umap, isl_union_map_copy(data->copy_sched));

        if (isl_union_map_is_empty(umap))
        {
            isl_union_map_free(umap);
            continue;
        }

        map = isl_map_from_union_map(umap);
        map = isl_map_detect_equalities(map);

        group = isl_calloc_type(ctx, struct amp_array_ref_group);
        if (!group)
        {
            isl_map_free(map);
            return -1;
        }
        group->local_array = local;
        group->array = local->array;
        group->access = map;
        group->write = access->write;
        group->exact_write = access->exact_write;
        group->slice = access->n_index < local->array->n_index;
        group->refs = &local->array->refs[i];
        group->n_ref = 1;

        groups[n++] = group;
    }

    return n;
}

/* Combine the given two groups into a single group, containing
 * the references of both groups.
 */
static struct amp_array_ref_group *join_groups(
    struct amp_array_ref_group *group1,
    struct amp_array_ref_group *group2)
{
    int i;
    isl_ctx *ctx;
    struct amp_array_ref_group *group;

    if (!group1 || !group2)
        return NULL;

    ctx = isl_map_get_ctx(group1->access);
    group = isl_calloc_type(ctx, struct amp_array_ref_group);
    if (!group)
        return NULL;
    group->local_array = group1->local_array;
    group->array = group1->array;
    group->access = isl_map_union(isl_map_copy(group1->access),
                                  isl_map_copy(group2->access));
    group->write = group1->write || group2->write;
    group->exact_write = group1->exact_write && group2->exact_write;
    group->slice = group1->slice || group2->slice;
    group->n_ref = group1->n_ref + group2->n_ref;
    group->refs = isl_alloc_array(ctx, struct amp_stmt_access *,
                                  group->n_ref);
    if (!group->refs)
        return amp_array_ref_group_free(group);
    for (i = 0; i < group1->n_ref; ++i)
        group->refs[i] = group1->refs[i];
    for (i = 0; i < group2->n_ref; ++i)
        group->refs[group1->n_ref + i] = group2->refs[i];

    return group;
}

/* Combine the given two groups into a single group and free
 * the original two groups.
 */
static struct amp_array_ref_group *join_groups_and_free(
    struct amp_array_ref_group *group1,
    struct amp_array_ref_group *group2)
{
    struct amp_array_ref_group *group;

    group = join_groups(group1, group2);
    amp_array_ref_group_free(group1);
    amp_array_ref_group_free(group2);
    return group;
}

/* Combine all groups in "groups" into a single group and return
 * the new number of groups (1 or 0 if there were no groups to start with).
 */
static int join_all_groups(int n, struct amp_array_ref_group **groups)
{
    int i;

    for (i = n - 1; i > 0; --i)
    {
        groups[0] = join_groups_and_free(groups[0], groups[i]);
        groups[i] = NULL;
        n--;
    }

    return n;
}

/* Set array->n_group and array->groups to n and groups.
 *
 * Additionally, set the "nr" field of each group.
 */
static void set_array_groups(struct amp_local_array_info *array,
                             int n, struct amp_array_ref_group **groups)
{
    int i;

    array->n_group = n;
    array->groups = groups;

    for (i = 0; i < n; ++i)
        groups[i]->nr = i;
}

/* Check if the access relations of group1 and group2 overlap within
 * copy_sched.
 */
static int accesses_overlap(struct amp_array_ref_group *group1,
                            struct amp_array_ref_group *group2)
{
    int disjoint;

    disjoint = isl_map_is_disjoint(group1->access, group2->access);
    if (disjoint < 0)
        return -1;

    return !disjoint;
}

/* Construct a map from domain_space to domain_space that increments
 * the dimension at position "pos" and leaves all other dimensions
 * constant.
 */
static __isl_give isl_map *next(__isl_take isl_space *domain_space, int pos)
{
    isl_space *space;
    isl_aff *aff;
    isl_multi_aff *next;

    space = isl_space_map_from_set(domain_space);
    next = isl_multi_aff_identity(space);
    aff = isl_multi_aff_get_aff(next, pos);
    aff = isl_aff_add_constant_si(aff, 1);
    next = isl_multi_aff_set_aff(next, pos, aff);

    return isl_map_from_multi_aff(next);
}

/* Check if the given access is coalesced (or if there is no point
 * in trying to coalesce the access by mapping the array to shared memory).
 * That is, check whether incrementing the dimension that will get
 * wrapped over the last thread index results in incrementing
 * the last array index.
 *
 * If no two consecutive array elements are ever accessed by "access",
 * then mapping the corresponding array to shared memory will not
 * improve coalescing.  In fact, the copying will likely be performed
 * by a single thread.  Consider the access as coalesced such that
 * the caller will not try and map the array to shared memory just
 * to improve coalescing.
 *
 * This function is only called for access relations without reuse and
 * kernels with at least one thread identifier.
 */
static int access_is_coalesced(struct amp_group_data *data,
                               __isl_keep isl_union_map *access)
{
    // #define DEBUG_ACCESS_IS_COALESCED

    int dim;
    isl_space *space;
    isl_set *accessed;
    isl_map *access_map;
    isl_map *next_thread_x;
    isl_map *next_element;
    isl_map *map;
    int coalesced, empty;

#ifdef DEBUG_ACCESS_IS_COALESCED
    fprintf(stderr, "@DEBUG: \n       in access_is_coalesced function, the data is: \n");
    amp_group_data_dump(data);
    fprintf(stderr, "        the access is: \n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n");
#endif // DEBUG_ACCESS_IS_COALESCED

    access = isl_union_map_copy(access);
    access = isl_union_map_apply_domain(access, isl_union_map_copy(data->full_sched));
    access_map = isl_map_from_union_map(access);

    space = isl_map_get_space(access_map);
    space = isl_space_range(space);
    dim = isl_space_dim(space, isl_dim_set);
    if (dim == 0)
        next_element = isl_map_empty(isl_space_map_from_set(space));
    else
        next_element = next(space, dim - 1);

    accessed = isl_map_range(isl_map_copy(access_map));
    map = isl_map_copy(next_element);
    map = isl_map_intersect_domain(map, isl_set_copy(accessed));
    map = isl_map_intersect_range(map, accessed);
    empty = isl_map_is_empty(map);
    isl_map_free(map);

    if (empty < 0 || empty)
    {
        isl_map_free(next_element);
        isl_map_free(access_map);
        return empty;
    }

    space = isl_map_get_space(access_map);
    space = isl_space_domain(space);
    next_thread_x = next(space, data->thread_depth + data->n_thread - 1);

    map = isl_map_apply_domain(next_thread_x, isl_map_copy(access_map));
    map = isl_map_apply_range(map, access_map);

    coalesced = isl_map_is_subset(map, next_element);

    isl_map_free(next_element);
    isl_map_free(map);

    return coalesced;
}

/* Replace the host schedule dimensions in the access relation "access"
 * by parameters, so that they are treated as fixed when checking for reuse
 * (within a kernel) or whether two consecutive elements are accessed
 * (within a kernel).
 */
static __isl_give isl_union_map *localize_access(struct amp_group_data *data,
                                                 __isl_take isl_union_map *access)
{
    // #define DEBUG_LOCALIZE_ACCESS

    int n;
    isl_space *space;
    isl_set *param;
    isl_union_map *umap;
    isl_id_list *ids;

    umap = isl_union_map_copy(data->host_sched);
    space = isl_union_map_get_space(umap);

#ifdef DEBUG_LOCALIZE_ACCESS
    fprintf(stderr, "@DEBUG: \n       in start 1 of localize_access function, the space is:\n");
    isl_space_dump(space);
    fprintf(stderr, "\n       the umap is:\n");
    isl_union_map_dump(umap);
    fprintf(stderr, "\n\n");
#endif // DEBUG_LOCALIZE_ACCESS

    n = data->kernel_depth;
    ids = ppcg_scop_generate_names(data->scop, n, "__ppcg_host_");
    param = parametrization(space, n, 0, ids);

#ifdef DEBUG_LOCALIZE_ACCESS
    fprintf(stderr, "@DEBUG: \n       in start 2 of localize_access function, the param is:\n");
    isl_set_dump(param);
    fprintf(stderr, "\n       the ids is:\n");
    isl_id_list_dump(ids);
    fprintf(stderr, "\n       the umap is:\n");
    isl_union_map_dump(umap);
    fprintf(stderr, "\n\n");
#endif // DEBUG_LOCALIZE_ACCESS

    isl_id_list_free(ids);
    umap = isl_union_map_intersect_range(umap, isl_union_set_from_set(param));

#ifdef DEBUG_LOCALIZE_ACCESS
    fprintf(stderr, "@DEBUG: \n       in middle of localize_access function, the umap is:\n");
    isl_union_map_dump(umap);
    fprintf(stderr, "\n");
#endif // DEBUG_LOCALIZE_ACCESS

    access = isl_union_map_intersect_domain(access, isl_union_map_domain(umap));

#ifdef DEBUG_LOCALIZE_ACCESS
    fprintf(stderr, "@DEBUG: \n       in end of  localize_access function, the final access is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "\n\n");
#endif // DEBUG_LOCALIZE_ACCESS

    return access;
}

/* Map the domain of "access" to the outer data->shared_depth
 * schedule dimensions.  When data->shared_depth is equal to
 * data->thread_depth, this result is already available in group->access.
 */
static __isl_give isl_map *shared_access(struct amp_array_ref_group *group,
                                         __isl_keep isl_union_map *access, struct amp_group_data *data)
{
    isl_union_map *shared;
    if (data->shared_depth == data->thread_depth)
        return isl_map_copy(group->access);

    shared = isl_union_map_copy(access);
    shared = isl_union_map_apply_domain(shared, isl_union_map_copy(data->shared_sched));

    return isl_map_from_union_map(shared);
}

/* Given an array access "access", check if for any index i there is
 * a shift a(p) and a stride g such that
 *
 *	a(p) + i = 0 mod g
 *
 * If so, record the information in tile->bound[i]->stride and
 * tile->bound[i]->shift.
 * Otherwise, set tile->bound[i]->stride to 1 (and tile->bound[i]->shift to 0).
 * Return isl_bool_true if any non-trivial stride was found.
 *
 * Note that the stride info returned by isl_map_get_range_stride_info
 * is of the form
 *
 *	i = o(p) + g n
 *
 * a(p) can therefore be taken to be equal to -o(p).
 */
static isl_bool detect_strides(struct amp_array_tile *tile,
                               __isl_keep isl_map *access)
{
    int i;
    isl_bool has_strides = isl_bool_false;

    for (i = 0; i < tile->n; ++i)
    {
        struct amp_array_bound *bound = &tile->bound[i];
        isl_stride_info *si;

        si = isl_map_get_range_stride_info(access, i);
        bound->stride = isl_stride_info_get_stride(si);
        bound->shift = isl_aff_neg(isl_stride_info_get_offset(si));
        isl_stride_info_free(si);

        if (!has_strides)
            has_strides = isl_val_gt_si(bound->stride, 1);
        if (has_strides < 0)
            return isl_bool_error;
    }

    return has_strides;
}

/* Given an array access "access", remove the strides based
 * on the information in tile->bound[i]->stride and tile->bound[i]->shift.
 *
 * In particular let the access be A[a] and
 * let the shifts s_i(p) and the strides g_i be such that
 *
 *  S(p) + a = 0 mod G
 *
 * Replace the access by
 *
 *  A[(a + S(p))/G]
 *
 * First collect the shifts s_i into an isl_multi_aff and
 * the strides into the scaling function A[i] -> A[G i].
 * Then add the shifts to the original access and
 * take the preimage over the scaling.
 */
static __isl_give isl_map *remove_strides(__isl_take isl_map *access,
                                          struct amp_array_tile *tile)
{
    int i;
    isl_space *space;
    isl_multi_aff *shift, *scale;
    isl_multi_val *stride;

    space = isl_map_get_space(access);
    shift = isl_multi_aff_zero(isl_space_copy(space));
    space = isl_space_range(space);
    stride = isl_multi_val_zero(isl_space_copy(space));
    scale = isl_multi_aff_identity(isl_space_map_from_set(space));
    for (i = 0; i < tile->n; ++i)
    {
        struct amp_array_bound *bound = &tile->bound[i];
        isl_aff *shift_i;
        isl_val *stride_i;

        shift_i = isl_aff_copy(bound->shift);
        stride_i = isl_val_copy(bound->stride);
        shift = isl_multi_aff_set_aff(shift, i, shift_i);
        stride = isl_multi_val_set_val(stride, i, stride_i);
    }
    scale = isl_multi_aff_scale_multi_val(scale, stride);

    access = isl_map_sum(access, isl_map_from_multi_aff(shift));
    access = isl_map_preimage_range_multi_aff(access, scale);

    return access;
}

/* Check if we can find a memory tile for the given array
 * based on the given accesses, and if so, put the results in "tile".
 *
 * We project the accesses on each index in turn and look for a parametric
 * offset such that the size is constant, after removing
 * any stride that may appear in the accesses.
 *
 * tile->depth is initialized to the input dimension of the computed bounds.
 */
static isl_bool can_tile(__isl_keep isl_map *access,
                         struct amp_array_tile *tile)
{
    // #define DEBUG_CAN_TILE

    int i;
    isl_bool has_strides, valid;
    isl_fixed_box *box;
    isl_multi_aff *offset;
    isl_multi_val *size;

#ifdef DEBUG_CAN_TILE
    fprintf(stderr, "@DEBUG: \n       the access is: \n");
    isl_map_dump(access);
    fprintf(stderr, "       the tile is: \n");
    amp_array_tile_dump(tile);
#endif // DEBUG_CAN_TILE

    if (!tile)
        return isl_bool_error;

    isl_map_free(isl_map_detect_equalities(isl_map_copy(access)));

    has_strides = detect_strides(tile, access);
#ifdef DEBUG_CAN_TILE
    fprintf(stderr, "@DEBUG: \n       the access is: \n");
    isl_map_dump(access);
    fprintf(stderr, "       the has_strides is: %d \n\n", has_strides);
#endif // DEBUG_CAN_TILE
    if (has_strides < 0)
        return isl_bool_error;

    tile->depth = isl_map_dim(access, isl_dim_in);
#ifdef DEBUG_CAN_TILE
    fprintf(stderr, "@DEBUG: \n       the tile->depth is: %d \n\n", tile->depth);
#endif // DEBUG_CAN_TILE

    access = isl_map_copy(access);
    if (has_strides)
        access = remove_strides(access, tile);

    box = isl_map_get_range_simple_fixed_box_hull(access);
#ifdef DEBUG_CAN_TILE
    fprintf(stderr, "@DEBUG: \n       the box is: \n");
    isl_fixed_box_dump(box);
    fprintf(stderr, "\n");
#endif // DEBUG_CAN_TILE

    // isl_map_free(access);

    valid = isl_fixed_box_is_valid(box);
#ifdef DEBUG_CAN_TILE
    fprintf(stderr, "@DEBUG: \n       the valid is: %d \n\n", valid);
#endif // DEBUG_CAN_TILE

    if (valid >= 0 && valid)
    {
        offset = isl_fixed_box_get_offset(box);
        size = isl_fixed_box_get_size(box);
#ifdef DEBUG_CAN_TILE
        fprintf(stderr, "@DEBUG: \n       in the can tile function, the offset is: \n");
        isl_multi_aff_dump(offset);
        fprintf(stderr, "@DEBUG: \n       in the can tile function, the  size  is: \n");
        isl_multi_val_dump(size);
#endif // DEBUG_CAN_TILE
        for (i = 0; i < tile->n; ++i)
        {
            tile->bound[i].size = isl_multi_val_get_val(size, i);
            tile->bound[i].lb = isl_multi_aff_get_aff(offset, i);
        }
        isl_multi_aff_free(offset);
        isl_multi_val_free(size);
    }
    else if (!valid)
    {
        box = isl_map_get_range_lattice_tile(access);
#ifdef DEBUG_CAN_TILE
        fprintf(stderr, "@DEBUG: \n       get_range_lattice_tile is: \n");
        isl_fixed_box_dump(box);
#endif // DEBUG_CAN_TILE

        offset = isl_fixed_box_get_offset(box);
        size = isl_fixed_box_get_size(box);
        for (i = 0; i < tile->n; ++i)
        {
            tile->bound[i].size = isl_multi_val_get_val(size, i);
            tile->bound[i].lb = isl_multi_aff_get_aff(offset, i);
        }
        isl_multi_aff_free(offset);
        isl_multi_val_free(size);
    }
    isl_fixed_box_free(box);

    return valid;
}

struct amp_array_tile *amp_array_tile_free(struct amp_array_tile *tile)
{
    int j;

    if (!tile)
        return NULL;

    for (j = 0; j < tile->n; ++j)
    {
        isl_val_free(tile->bound[j].size);
        isl_val_free(tile->bound[j].stride);
        isl_aff_free(tile->bound[j].lb);
        isl_aff_free(tile->bound[j].shift);
    }
    free(tile->bound);
    isl_multi_aff_free(tile->tiling);
    free(tile);

    return NULL;
}

/* Create a gpu_array_tile for an array of dimension "n_index".
 */
static struct amp_array_tile *amp_array_tile_create(isl_ctx *ctx, int n_index)
{
    // #define DEBUG_AMP_ARRAY_TILE

    int i;
    struct amp_array_tile *tile;

    tile = isl_calloc_type(ctx, struct amp_array_tile);
    if (!tile)
        return NULL;

    tile->ctx = ctx;
    tile->bound = isl_alloc_array(ctx, struct amp_array_bound, n_index);
    if (!tile->bound)
        return amp_array_tile_free(tile);

    tile->n = n_index;

    for (i = 0; i < n_index; ++i)
    {
        tile->bound[i].size = NULL;
        tile->bound[i].lb = NULL;
        tile->bound[i].stride = NULL;
        tile->bound[i].shift = NULL;
    }
#ifdef DEBUG_AMP_ARRAY_TILE
    fprintf(stderr, "@DEBUG: \n       in amp_array_tile_create over, the tile information are: \n");
    fprintf(stderr, "          the tile->n is %d, the tile->depth is %d . \n\n", tile->n, tile->depth);
    fprintf(stderr, "          the tile->n is %d, the tile->depth is %d, the tile->tiling is:\n", tile->n, tile->depth);
    isl_multi_aff_dump(tile->tiling);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ARRAY_TILE
    return tile;
}

/* Report that the array reference group with the given access relation
 * is not mapped to shared memory in the given kernel because
 * it does not exhibit any reuse and is considered to be coalesced.
 */
static void report_no_reuse_and_coalesced(struct amp_ppcg_kernel *kernel,
                                          __isl_keep isl_union_map *access)
{
    isl_ctx *ctx;
    isl_printer *p;

    ctx = isl_union_map_get_ctx(access);
    p = isl_printer_to_file(ctx, stdout);
    p = isl_printer_print_str(p, "Array reference group ");
    p = isl_printer_print_union_map(p, access);
    p = isl_printer_print_str(p,
                              " not considered for mapping to shared memory in kernel");
    p = isl_printer_print_int(p, kernel->id);
    p = isl_printer_print_str(p,
                              " because it exhibits no reuse and is considered to be coalesced");
    p = isl_printer_end_line(p);
    isl_printer_free(p);
}

/* Return the lowest depth between data->kernel_depth and data->thread_depth
 * at which every array element accessed through "acc" is accessed
 * by a single thread.  The input dimension of "acc" is
 * data->thread_depth + data->n_thread, where the final data->n_thread
 * dimensions are those that will be mapped to threads.
 * If the values for these dimensions are uniquely determined
 * by the array index and a given number of outer dimensions, then
 * there is only one thread accessing that array element within those
 * outer dimensions.
 *
 * The input space of "acc" is first split up, such that it has the form
 *
 *	[O -> T] -> A
 *
 * with O the outer dimensions, T the dimensions that will be mapped to threads
 * and A the array index.
 *
 * Then the positions of T and A are interchanged to simplify the test
 * whether T uniquely depends on O and A.
 * In particular, the above access relation is first combined with
 *
 *	[O -> T] -> T
 *
 * to form
 *
 *	[O -> T] -> [A -> T]
 *
 * from which
 *
 *	O -> [A -> T]
 *
 * is extracted, which is then uncurried to
 *
 *	[O -> A] -> T
 *
 * Finally, the final dimensions of O are projected out one by one
 * until T is no longer uniquely determined by A and the remaining
 * dimensions in O.  The value returned is that of the last dimension
 * that was successfully projected out.
 * Note that there is no need to test whether [O -> A] -> T itself
 * is single-valued as that was already tested in access_is_bijective.
 */
static int compute_accessed_by_single_thread_depth(struct amp_group_data *data,
                                                   __isl_keep isl_map *acc)
{
    int i;
    isl_space *space;
    isl_map *map;
    isl_bool sv;

    if (data->thread_depth == data->kernel_depth)
        return data->thread_depth;

    acc = isl_map_copy(acc);

    space = isl_map_get_space(acc);
    space = isl_space_params(space);
    space = isl_space_set_from_params(space);
    space = isl_space_add_dims(space, isl_dim_set, data->thread_depth);
    space = isl_space_from_domain(space);
    space = isl_space_add_dims(space, isl_dim_out, data->n_thread);
    space = isl_space_wrap(space);
    map = isl_set_flatten_map(isl_set_universe(space));
    acc = isl_map_apply_range(map, acc);

    space = isl_space_domain(isl_map_get_space(acc));
    map = isl_map_range_map(isl_map_universe(isl_space_unwrap(space)));
    acc = isl_map_range_product(acc, map);
    acc = isl_map_domain_factor_domain(acc);
    acc = isl_map_uncurry(acc);

    for (i = data->thread_depth - 1; i >= data->kernel_depth; --i)
    {
        acc = isl_map_project_out(acc, isl_dim_in, i, 1);
        sv = isl_map_is_single_valued(acc);
        if (sv < 0)
            goto error;
        if (!sv)
            break;
    }

    isl_map_free(acc);

    return ++i;
error:
    isl_map_free(acc);
    return -1;
}

/* Given an access relation in terms of at least data->thread_depth initial
 * dimensions of the computed schedule, check if it is bijective for
 * fixed values of the first data->thread_depth dimensions.
 * We perform this check by equating these dimensions to parameters.
 */
static int access_is_bijective(struct amp_group_data *data,
                               __isl_keep isl_map *access)
{
    int res;
    int dim;
    isl_set *par;
    isl_space *space;
    isl_id_list *ids;

    access = isl_map_copy(access);
    space = isl_space_params(isl_map_get_space(access));
    ids = ppcg_scop_generate_names(data->scop, data->thread_depth, "s");
    dim = isl_map_dim(access, isl_dim_in);
    par = parametrization(space, dim, 0, ids);
    isl_id_list_free(ids);
    access = isl_map_intersect_domain(access, par);
    res = isl_map_is_bijective(access);
    isl_map_free(access);

    return res;
}

/* Given an access relation in terms of the data->thread_depth initial
 * dimensions of the computed schedule and the thread identifiers
 * (as parameters), check if the use of the corresponding private tile
 * requires unrolling.
 *
 * If we are creating a private tile because we are forced to,
 * then no unrolling is required.
 * Otherwise we check if "access" is bijective and unrolling
 * is required if it is not.  Note that the access relation
 * has already been determined to be bijective before the introduction
 * of the thread identifiers and the removal of the schedule dimensions
 * that are mapped to these threads.  If the access relation is no longer
 * bijective, then this means that more than one value of one of those
 * schedule dimensions is mapped to the same thread and therefore
 * unrolling is required.
 */
static int check_requires_unroll(struct amp_group_data *data,
                                 __isl_keep isl_map *access, int force_private)
{
    int bijective;

    if (force_private)
        return 0;
    bijective = access_is_bijective(data, access);
    if (bijective < 0)
        return -1;
    return !bijective;
}

/* Compute the private and/or shared memory tiles for the array
 * reference group "group" of array "array".
 * Return isl_stat_ok on success and isl_stat_error on error.
 *
 * If the array is a read-only scalar or if the user requested
 * not to use shared or private memory, then we do not need to do anything.
 *
 * If any reference in the reference group accesses more than one element,
 * then we would have to make sure that the layout in shared memory
 * is the same as that in global memory.  Since we do not handle this yet
 * (and it may not even be possible), we refuse to map to private or
 * shared memory in such cases.
 *
 * If the array group involves any may writes (that are not must writes),
 * then we would have to make sure that we load the data into shared/private
 * memory first in case the data is not written by the kernel
 * (but still written back out to global memory).
 * Since we don't have any such mechanism at the moment, we don't
 * compute shared/private tiles for groups involving may writes.
 *
 * We only try to compute a shared memory tile if there is any reuse
 * or if the access is not coalesced.
 * Reuse and coalescing are checked within the given kernel.
 *
 * For computing a private memory tile, we also require that there is
 * some reuse.  Moreover, we require that the access is private
 * to the thread.  That is, we check that any given array element
 * is only accessed by a single thread.
 * We compute an access relation that maps the outer
 * data->thread_depth + data->n_thread schedule dimensions.
 * The latter data->n_thread will be mapped to thread identifiers.
 * We actually check that those iterators that will be wrapped
 * partition the array space.  This check is stricter than necessary
 * since several iterations may be mapped onto the same thread
 * and then they could be allowed to access the same memory elements,
 * but our check does not allow this situation.
 *
 * For private memory tiles, the number of schedule dimensions that
 * affect the offset is computed and stored in tile->depth, with
 * a lower bound of data->kernel_depth.  If this depth is smaller
 * than the minimal depth that still ensures that every element
 * is accessed by a single thread, then the depth is raised
 * to this minimal depth.
 * The fields of the tile are then adjusted to only refer to the tile->depth
 * outer schedule dimensions.
 *
 * We also check that the index expression only depends on parallel
 * loops.  That way, we can move those loops innermost and unroll them.
 * Again, we use a test that is stricter than necessary.
 * We actually check whether the index expression only depends
 * on the iterators that are wrapped over the threads.
 * These are necessarily parallel, but there may be more parallel loops.
 *
 * Combining the injectivity of the first test with the single-valuedness
 * of the second test, we simply test for bijectivity.
 *
 * If the use of the private tile requires unrolling, but some
 * of the other arrays are forcibly mapped to private memory,
 * then we do not allow the use of this private tile since
 * we cannot move the schedule dimensions that need to be unrolled down
 * without performing some kind of expansion on those arrays
 * that are forcibly mapped to private memory.
 *
 * If the array is marked force_private, then we bypass all checks
 * and assume we can (and should) use registers only.
 *
 * If it turns out we can (or have to) use registers, we compute
 * the private memory tile size using can_tile, after introducing a dependence
 * on the thread indices.
 */
static isl_stat compute_group_bounds_core(struct amp_ppcg_kernel *kernel,
                                          struct amp_array_ref_group *group, struct amp_group_data *data)
{
    // #define DEBUG_COMPUTE_GROUP_BOUNDS_CORE

    isl_ctx *ctx = isl_space_get_ctx(group->array->space);
    isl_union_map *access, *local;
    int n_index = group->array->n_index;
    int no_reuse, coalesced;
    isl_map *acc;
    int force_private = 0;
    int use_shared = 1;
    int use_private = 0;
    isl_stat r = isl_stat_ok;
    isl_bool ok;
    int requires_unroll;
    int unique_depth;

    if (!use_shared && !use_private)
        return isl_stat_ok;
    // if (amp_array_is_read_only_scalar(group->array))
    //     return isl_stat_ok;
    if (!group->exact_write)
        return isl_stat_ok;
    if (group->slice)
        return isl_stat_ok;

    access = amp_array_ref_group_access_relation(group, 1, 1);
    local = localize_access(data, isl_union_map_copy(access));
    no_reuse = isl_union_map_is_injective(local);

#ifdef DEBUG_COMPUTE_GROUP_BOUNDS_CORE
    fprintf(stderr, "@DEBUG: \n       in start of compute_group_bounds_core, the data is: \n");
    // amp_group_data_dump(data);
    fprintf(stderr, "        the group is:        ");
    // amp_array_ref_group_dump(group);
    fprintf(stderr, "        the access(read and write) is:\n");
    isl_union_map_dump(access);
    fprintf(stderr, "        the local is:\n");
    isl_union_map_dump(local);
    fprintf(stderr, "        the no_reuser value is: %d \n\n", no_reuse);
#endif // DEBUG_COMPUTE_GROUP_BOUNDS_CORE

    if (no_reuse < 0)
        r = isl_stat_error;
    if (no_reuse)
        coalesced = access_is_coalesced(data, local);
    isl_union_map_free(local);

    if (r >= 0 && kernel->options->debug->verbose && use_shared && no_reuse && coalesced)
        report_no_reuse_and_coalesced(kernel, access);

#ifdef DEBUG_COMPUTE_GROUP_BOUNDS_CORE
    fprintf(stderr, "@DEBUG: \n       the coalesced is: %d \n", coalesced);
#endif // DEBUG_COMPUTE_GROUP_BOUNDS_CORE
    // if (!no_reuse || !coalesced)
    if (1)
    {
        group->shared_tile = amp_array_tile_create(ctx, group->array->n_index);
        acc = shared_access(group, access, data);
#ifdef DEBUG_COMPUTE_GROUP_BOUNDS_CORE
        fprintf(stderr, "@TEMP_DEBUG: \n       the acc is : \n");
        isl_map_dump(acc);
#endif

        ok = can_tile(acc, group->shared_tile);

        if (ok < 0)
            r = isl_stat_error;
        // else if (!ok)
        // {
        //     struct amp_array_tile *tile = group->shared_tile;
        //     struct amp_array_info *array = group->array;
        //     for (int i = 0; i < tile->n; i++)
        //     {
        //         fprintf(stderr, "\n\n  hbobiuhfeiopfh  \n\n");
        //         isl_multi_pw_aff_get_at(array->bound, i);

        //         isl_val *v = isl_multi_val_get_val(isl_set_get_plain_multi_val_if_fixed(array->declared_extent), i);
        //         if (isl_val_is_one(tile->bound[i].size))
        //             tile->bound[i].size = v;

        //         isl_val_dump(tile->bound[i].size);
        //         fprintf(stderr, "\n\n 123456789987654321 \n\n");
        //     }
        // }

        // #ifdef DEBUG_COMPUTE_GROUP_BOUNDS_CORE
        //         if (group->shared_tile && group) {
        //             fprintf(stderr, "@DEBUG: \n       the shared_tile create over!    the ok of can_tile function is: %d \n       the shared_tile is : \n", ok);
        //             amp_array_tile_dump(group->shared_tile);
        //             fprintf(stderr, "\n       the group is : \n");
        //             amp_array_ref_group_dump(group);
        //         } else {
        //             fprintf(stderr, "@ERROR: \n       the created shared_tile is NULL!!! \n");
        //             fprintf(stderr, "          the group->array->n_index is %d 、the acc is:\n", group->array->n_index);
        //             isl_map_dump(acc);
        //             fprintf(stderr, "          the ok of can_tile function is: %d\n", ok);
        //             fprintf(stderr, "\n\n");
        //         }
        // #endif // DEBUG_COMPUTE_GROUP_BOUNDS_CORE
        isl_map_free(acc);
    }

    if (r < 0 || (!force_private && (!use_private || no_reuse)))
    {
        isl_union_map_free(access);
        return r;
    }

    // access = isl_union_map_apply_domain(access, isl_union_map_copy(data->thread_sched));

    // acc = isl_map_from_union_map(access);

    // if (!force_private && !access_is_bijective(data, acc))
    // {
    //     isl_map_free(acc);
    //     return isl_stat_ok;
    // }

    // unique_depth = compute_accessed_by_single_thread_depth(data, acc);

    // acc = isl_map_intersect_domain(acc, isl_set_copy(data->privatization));
    // acc = isl_map_project_out(acc, isl_dim_in, data->thread_depth, data->n_thread);
    // requires_unroll = check_requires_unroll(data, acc, force_private);
    // if (unique_depth < 0 || requires_unroll < 0 || (requires_unroll))
    // {
    //     isl_map_free(acc);
    //     return requires_unroll < 0 ? isl_stat_error : isl_stat_ok;
    // }

    // group->private_tile = amp_array_tile_create(ctx, n_index);
    // group->private_tile->requires_unroll = requires_unroll;
    // ok = can_tile(acc, group->private_tile);
    // if (ok >= 0 && !ok)
    //     group->private_tile = amp_array_tile_free(group->private_tile);
    // isl_map_free(acc);
    // if (ok < 0)
    //     return isl_stat_error;

    // if (group->private_tile)
    // {
    //     struct amp_array_tile *tile = group->private_tile;
    //     int tile_depth = compute_tile_depth(data, tile);
    //     if (tile_depth < unique_depth)
    //         tile_depth = unique_depth;
    //     if (tile_adjust_depth(tile, tile_depth) < 0)
    //         return isl_stat_error;
    // }

    // if (force_private && !group->private_tile)
    //     isl_die(ctx, isl_error_internal,
    //             "unable to map array reference group to registers",
    //             return isl_stat_error);

    return isl_stat_ok;
}

/* Compute the number of outer schedule tile dimensions that affect
 * the offset of "tile".
 * If there is no such dimension, then return the index
 * of the first kernel dimension, i.e., data->kernel_depth.
 */
static int compute_tile_depth(struct amp_group_data *data,
                              struct amp_array_tile *tile)
{
    int i, j;

    for (j = tile->depth - 1; j >= data->kernel_depth; --j)
    {
        for (i = 0; i < tile->n; ++i)
        {
            isl_aff *lb;
            isl_aff *shift;

            lb = tile->bound[i].lb;
            if (isl_aff_involves_dims(lb, isl_dim_in, j, 1))
                break;

            shift = tile->bound[i].shift;
            if (!shift)
                continue;
            if (isl_aff_involves_dims(shift, isl_dim_in, j, 1))
                break;
        }
        if (i < tile->n)
            break;
    }

    return ++j;
}

/* Adjust the fields of "tile" to reflect the new input dimension "depth".
 * The dimension beyond "depth" are assumed not to affect the tile,
 * so they can simply be dropped.
 */
static int tile_adjust_depth(struct amp_array_tile *tile, int depth)
{
    int i;

    if (tile->depth == depth)
        return 0;

    for (i = 0; i < tile->n; ++i)
    {
        tile->bound[i].lb = isl_aff_drop_dims(tile->bound[i].lb,
                                              isl_dim_in, depth, tile->depth - depth);
        if (!tile->bound[i].lb)
            return -1;
        if (!tile->bound[i].shift)
            continue;
        tile->bound[i].shift = isl_aff_drop_dims(tile->bound[i].shift,
                                                 isl_dim_in, depth, tile->depth - depth);
        if (!tile->bound[i].shift)
            return -1;
    }

    tile->depth = depth;

    return 0;
}

/* Determine the number of schedule dimensions that affect the offset of the
 * shared or private tile "tile" and store the result in tile->depth, with
 * a lower bound of data->kernel_depth.
 * Also adjust the fields of the tile to only refer to the tile->depth
 * outer schedule dimensions.
 */
static isl_stat tile_set_depth(struct amp_group_data *data,
                               struct amp_array_tile *tile)
{
    if (tile_adjust_depth(tile, compute_tile_depth(data, tile)) < 0)
        return isl_stat_error;

    return isl_stat_ok;
}

/* Determine the number of schedule dimensions that affect the offset of the
 * shared tile and store the minimum of the private and shared tile depth
 * in group->min_depth, with a lower bound of data->kernel_depth.
 * If there is no tile defined on the array reference group,
 * then set group->min_depth to data->thread_depth.
 */
static int set_depth(struct amp_group_data *data,
                     struct amp_array_ref_group *group)
{
    group->min_depth = data->thread_depth;

    if (group->shared_tile)
    {
        if (tile_set_depth(data, group->shared_tile) < 0)
            return -1;
        if (group->shared_tile->depth < group->min_depth)
            group->min_depth = group->shared_tile->depth;
    }

    return 0;
}

/* Compute the private and/or shared memory tiles for the array
 * reference group "group" of array "array" and set the tile depth.
 * Return 0 on success and -1 on error.
 */
static int compute_group_bounds(struct amp_ppcg_kernel *kernel,
                                struct amp_array_ref_group *group, struct amp_group_data *data)
{

    if (!group)
        return -1;
    if (compute_group_bounds_core(kernel, group, data) < 0)
        return -1;
    if (set_depth(data, group) < 0)
        return -1;

    return 0;
}

/* If two groups have overlapping access relations (as determined by
 * the "overlap" function) and if one of them involves a write,
 * then merge the two groups into one.
 * If "compute_bounds" is set, then call compute_group_bounds
 * on the merged groups.
 * If any group is merged into the current group, then its access
 * relation may have changed or it may have been turned into a write.
 * The combined group might therefore overlap with groups that
 * the original group did not overlap with.  The groups therefore
 * need to be checked again.
 *
 * Return the updated number of groups.
 * Return -1 on error.
 */
static int group_writes(struct amp_ppcg_kernel *kernel,
                        int n, struct amp_array_ref_group **groups,
                        int (*overlap)(struct amp_array_ref_group *group1,
                                       struct amp_array_ref_group *group2),
                        int compute_bounds,
                        struct amp_group_data *data)
{
    int i, j;
    int any_merge;

    for (i = 0; i < n; i += !any_merge)
    {
        any_merge = 0;
        for (j = n - 1; j > i; --j)
        {
            if (!groups[i]->write && !groups[j]->write)
                continue;

            if (!overlap(groups[i], groups[j]))
                continue;

            any_merge = 1;
            groups[i] = join_groups_and_free(groups[i], groups[j]);
            if (j != n - 1)
                groups[j] = groups[n - 1];
            groups[n - 1] = NULL;
            n--;

            if (!groups[i])
                return -1;
            if (compute_bounds &&
                compute_group_bounds(kernel, groups[i], data) < 0)
                return -1;
        }
    }

    return n;
}

/* If two groups have overlapping access relations (within the innermost
 * loop) and if one of them involves a write, then merge the two groups
 * into one.
 *
 * Return the updated number of groups.
 */
static int group_overlapping_writes(struct amp_ppcg_kernel *kernel,
                                    int n, struct amp_array_ref_group **groups,
                                    struct amp_group_data *data)
{
    return group_writes(kernel, n, groups, &accesses_overlap, 0, data);
}

/* Check if the access relations of group1 and group2 overlap within
 * the outermost min(group1->min_depth, group2->min_depth) loops.
 */
static int depth_accesses_overlap(struct amp_array_ref_group *group1,
                                  struct amp_array_ref_group *group2)
{
    int depth;
    int dim;
    int empty;
    isl_map *map_i, *map_j, *map;

    depth = group1->min_depth;
    if (group2->min_depth < depth)
        depth = group2->min_depth;
    map_i = isl_map_copy(group1->access);
    dim = isl_map_dim(map_i, isl_dim_in);
    map_i = isl_map_eliminate(map_i, isl_dim_in, depth, dim - depth);
    map_j = isl_map_copy(group2->access);
    map_j = isl_map_eliminate(map_j, isl_dim_in, depth, dim - depth);
    map = isl_map_intersect(map_i, map_j);
    empty = isl_map_is_empty(map);
    isl_map_free(map);

    return !empty;
}

/* Compute the size of the tile specified by "tile"
 * in number of elements and return the result.
 */
static __isl_give isl_val *amp_array_tile_size(struct amp_array_tile *tile)
{
    int i;
    isl_val *size;

    if (!tile)
        return NULL;

    size = isl_val_one(tile->ctx);

    for (i = 0; i < tile->n; ++i)
        size = isl_val_mul(size, isl_val_copy(tile->bound[i].size));

    return size;
}

/* Is the size of the tile specified by "tile" smaller than the sum of
 * the sizes of the tiles specified by "tile1" and "tile2"?
 */
static int smaller_tile(struct amp_array_tile *tile,
                        struct amp_array_tile *tile1, struct amp_array_tile *tile2)
{
    int smaller;
    isl_val *size, *size1, *size2;

    size = amp_array_tile_size(tile);
    size1 = amp_array_tile_size(tile1);
    size2 = amp_array_tile_size(tile2);

    size = isl_val_sub(size, size1);
    size = isl_val_sub(size, size2);
    smaller = isl_val_is_neg(size);

    isl_val_free(size);

    return smaller;
}

/* If two groups have overlapping access relations (within the outer
 * depth loops) and if one of them involves a write,
 * then merge the two groups into one.
 *
 * Return the updated number of groups.
 */
static int group_depth_overlapping_writes(struct amp_ppcg_kernel *kernel,
                                          int n, struct amp_array_ref_group **groups, struct amp_group_data *data)
{
    return group_writes(kernel, n, groups, &depth_accesses_overlap, 1, data);
}

/* Given an initial grouping of array references and shared memory tiles
 * for each group that allows for a shared memory tile, merge two groups
 * if both have a shared memory tile, the merged group also has
 * a shared memory tile and the size of the tile for the merge group
 * is smaller than the sum of the tile sizes of the individual groups.
 * If any group is merged into the current group, then it may become
 * profitable to combine it with groups that were considered before
 * the merge.  The groups are therefore checked again after a merge.
 *
 * If merging two groups decreases the depth of the tile of
 * one or both of the two groups, then we need to check for overlapping
 * writes again.
 *
 * Return the number of groups after merging.
 * Return -1 on error.
 */
static int group_common_shared_memory_tile(struct amp_ppcg_kernel *kernel,
                                           struct amp_array_info *array, int n,
                                           struct amp_array_ref_group **groups, struct amp_group_data *data)
{
    int i, j;
    int recompute_overlap = 0;
    int any_merge;

    for (i = 0; i < n; i += !any_merge)
    {
        any_merge = 0;
        if (!groups[i]->shared_tile)
            continue;
        for (j = n - 1; j > i; --j)
        {
            struct amp_array_ref_group *group;

            if (!groups[j]->shared_tile)
                continue;

            if (!depth_accesses_overlap(groups[i], groups[j]))
                continue;

            group = join_groups(groups[i], groups[j]);
            if (compute_group_bounds(kernel, group, data) < 0)
            {
                amp_array_ref_group_free(group);
                return -1;
            }
            if (!group->shared_tile || !smaller_tile(group->shared_tile, groups[i]->shared_tile, groups[j]->shared_tile))
            {
                amp_array_ref_group_free(group);
                continue;
            }

            any_merge = 1;
            if (group->min_depth < groups[i]->min_depth ||
                group->min_depth < groups[j]->min_depth)
                recompute_overlap = 1;
            amp_array_ref_group_free(groups[i]);
            amp_array_ref_group_free(groups[j]);
            groups[i] = group;
            if (j != n - 1)
                groups[j] = groups[n - 1];
            n--;
        }
    }

    if (recompute_overlap)
        n = group_depth_overlapping_writes(kernel, n, groups, data);
    return n;
}

/* Group array references that should be considered together when
 * deciding whether to access them from private, shared or global memory.
 * Return -1 on error.
 *
 * In particular, if two array references overlap and if one of them
 * is a write, then the two references are grouped together.
 * We first perform an initial grouping based only on the access relation.
 * After computing shared and private memory tiles, we check for
 * overlapping writes again, but this time taking into account
 * the depth of the effective tile.
 *
 * Furthermore, if two groups admit a shared memory tile and if the
 * combination of the two also admits a shared memory tile, we merge
 * the two groups.
 *
 * If the array contains structures, then we compute a single
 * reference group without trying to find any tiles
 * since we do not map such arrays to private or shared
 * memory.  The only exception is when those arrays of structures
 * are required to be mapped to private memory.
 */
static int amp_group_array_references(struct amp_ppcg_kernel *kernel,
                                      struct amp_local_array_info *local, struct amp_group_data *data)
{
    // #define DEBUG_AMP_GROUP_ARRAY_REFERENCES

    int i;
    int n;
    isl_ctx *ctx = isl_union_map_get_ctx(data->shared_sched);
    struct amp_array_ref_group **groups;

    groups = isl_calloc_array(ctx, struct amp_array_ref_group *, local->array->n_ref);
    if (!groups)
        return -1;

    n = populate_array_references(local, groups, data);

    if (local->array->has_compound_element)
    {
        n = join_all_groups(n, groups);
        set_array_groups(local, n, groups);
        return 0;
    }

    n = group_overlapping_writes(kernel, n, groups, data);

    for (i = 0; i < n; ++i)
        if (compute_group_bounds(kernel, groups[i], data) < 0)
            n = -1;

    n = group_depth_overlapping_writes(kernel, n, groups, data);

    n = group_common_shared_memory_tile(kernel, local->array, n, groups, data);

    set_array_groups(local, n, groups);

#ifdef DEBUG_AMP_GROUP_ARRAY_REFERENCES
    fprintf(stderr, "@DEBUG: \n       the group information in in end of amp_group_array_references is on the below:\n");
    for (i = 0; i < n; ++i)
    {
        fprintf(stderr, "           in amp_group_array_references structs, the index is %d.           the group is:\n", i);
        amp_array_ref_group_dump(groups[i]);
        fprintf(stderr, "\n\n");
    }
#endif // DEBUG_AMP_GROUP_ARRAY_REFERENCES

    if (n >= 0)
        return 0;

    // for (i = 0; i < local->array->n_ref; ++i)
    //     amp_array_ref_group_free(groups[i]);
    return -1;
}

/* Expand the domain of the schedule "s" by plugging in
 * the contraction "contraction" and return the result.
 */
static __isl_give isl_union_map *expand(__isl_take isl_union_map *s,
                                        __isl_keep isl_union_pw_multi_aff *contraction)
{
    contraction = isl_union_pw_multi_aff_copy(contraction);
    s = isl_union_map_preimage_domain_union_pw_multi_aff(s, contraction);
    return s;
}

/* Return the prefix schedule at "node" as a relation
 * between domain elements and schedule dimensions after detecting
 * equalities in this relation.
 */
static __isl_give isl_union_map *prefix_with_equalities(
    __isl_keep isl_schedule_node *node)
{
    isl_union_map *schedule;

    schedule = isl_schedule_node_get_prefix_schedule_relation(node);
    schedule = isl_union_map_detect_equalities(schedule);

    return schedule;
}

/* Group references of all arrays in "kernel".
 * "node" points to the kernel mark.
 * The mapping to shared memory in computed at the "shared" mark.
 *
 * We first extract all required schedule information into
 * a gpu_group_data structure and then consider each array
 * in turn.
 */
int amp_group_references(struct amp_ppcg_kernel *kernel,
                         __isl_keep isl_schedule_node *node)
{
    // #define DEBUG_AMP_GROUP_REFERENCES

    int i;
    int r = 0;
    isl_union_pw_multi_aff *contraction;
    struct amp_group_data data;

#ifdef DEBUG_AMP_GROUP_REFERENCES
    fprintf(stderr, "@DEBUG: \n       in amp_group_references, node is :\n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_GROUP_REFERENCES

    // check_can_be_private_live_ranges(kernel, node);

    data.scop = kernel->prog->scop;

    data.kernel_depth = isl_schedule_node_get_schedule_depth(node);
    data.host_sched = isl_schedule_node_get_prefix_schedule_relation(node);

    node = isl_schedule_node_copy(node);
    node = amp_tree_move_down_to_shared(node, kernel->core);
    data.shared_depth = isl_schedule_node_get_schedule_depth(node);
    data.shared_sched = prefix_with_equalities(node);

#ifdef DEBUG_AMP_GROUP_REFERENCES
    fprintf(stderr, "@DEBUG: \n       data.shared_sched is :\n");
    isl_union_map_dump(data.shared_sched);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_GROUP_REFERENCES
    node = amp_tree_move_down_to_thread(node, kernel->core);
    node = isl_schedule_node_child(node, 0);
    data.thread_depth = isl_schedule_node_get_schedule_depth(node);
    data.n_thread = isl_schedule_node_band_n_member(node);
    if (data.thread_depth == data.shared_depth)
        data.copy_sched = isl_union_map_copy(data.shared_sched);
    else
        data.copy_sched = prefix_with_equalities(node);
    data.thread_sched = isl_union_map_copy(data.copy_sched);
    data.thread_sched = isl_union_map_flat_range_product(data.thread_sched, isl_schedule_node_band_get_partial_schedule_union_map(node));
    data.thread_sched = isl_union_map_detect_equalities(data.thread_sched);

    contraction = isl_union_pw_multi_aff_copy(kernel->contraction);
    data.host_sched = expand(data.host_sched, contraction);
    data.shared_sched = expand(data.shared_sched, contraction);
    if (data.thread_depth == data.shared_depth)
    {
        isl_union_map_free(data.copy_sched);
        data.copy_sched = isl_union_map_copy(data.shared_sched);
    }
    else
    {
        data.copy_sched = expand(data.copy_sched, contraction);
    }
    data.thread_sched = expand(data.thread_sched, contraction);
    isl_union_pw_multi_aff_free(contraction);

    node = isl_schedule_node_child(node, 0);
    // data.full_sched = isl_union_map_copy(data.shared_sched);
    data.full_sched = isl_union_map_copy(data.thread_sched);
    data.full_sched = isl_union_map_flat_range_product(data.full_sched, isl_schedule_node_get_subtree_schedule_union_map(node));
    isl_schedule_node_free(node);

    // compute_privatization(&data, kernel);

    for (i = 0; i < kernel->n_array; ++i)
    {
        r = amp_group_array_references(kernel, &kernel->array[i], &data);
        if (r < 0)
            break;
    }

    isl_union_map_free(data.host_sched);
    isl_union_map_free(data.shared_sched);
    isl_union_map_free(data.copy_sched);
    isl_union_map_free(data.thread_sched);
    isl_union_map_free(data.full_sched);
    // isl_set_free(data.privatization);

    return r;
}

/* Check whether the array reference group "group" is mapped to
 * private or shared memory and, if so,
 * add copy statements to the schedule tree of "node"
 * for reading from global memory to private or shared memory
 * (if "read" is set) or for writing back from private or shared memory
 * to global memory (if "read" is not set) for this group.
 * On input, "node" points to the kernel node, and it is moved
 * back there on output.
 */
static __isl_give isl_schedule_node *amp_add_copies_group(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
    // #define DEBUG_AMP_ADD_COPIES_GROUP

    enum ppcg_group_access_type type;

    type = amp_array_ref_group_type(group);
    if (type == ppcg_access_shared)
    {
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        fprintf(stderr, "@DEBUG: \n       the shared (amp_add_copies_group) : group->array is:\n");
        amp_array_info_dump(group->array);
#endif // DEBUG_AMP_ADD_COPIES_GROUP

        return amp_add_copies_group_shared(kernel, group, node, read);
        // } else if (type == ppcg_access_global) {
        //     return amp_add_copies_group_shared(kernel, group, node, read);
    }

    return node;
}

/* For each array reference group that is mapped to private or shared memory,
 * add copy statements to the schedule tree of "node"
 * for reading from global memory to private or shared memory
 * and for writing back.
 * On input, "node" points to the kernel node, and it is moved
 * back there on output.
 */
static __isl_give isl_schedule_node *amp_add_copies(struct amp_ppcg_kernel *kernel,
                                                    __isl_take isl_schedule_node *node)
{
    int i, j;

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *array = &kernel->array[i];

        for (j = 0; j < array->n_group; ++j)
        {
            struct amp_array_ref_group *group = array->groups[j];

            node = amp_add_copies_group(kernel, group, node, 1);
            if (!node)
                return NULL;
            node = amp_add_copies_group(kernel, group, node, 0);
            if (!node)
                return NULL;
        }
    }

    return node;
}

static void create_amp_kernel_var(isl_ctx *ctx, struct amp_array_ref_group *group,
                                  struct amp_ppcg_kernel_var *var)
{
    int j;
    struct amp_array_tile *tile;
    isl_printer *p;

    var->array = group->array;

    var->type = amp_array_ref_group_type(group);
    tile = amp_array_ref_group_tile(group);

    p = isl_printer_to_str(ctx);
    p = amp_array_ref_group_print_name(group, p);
    var->name = isl_printer_get_str(p);
    isl_printer_free(p);

    var->size = isl_vec_alloc(ctx, group->array->n_index);

    for (j = 0; j < group->array->n_index; ++j)
    {
        var->size = isl_vec_set_element_val(var->size, j, isl_val_copy(tile->bound[j].size));
    }
}

static isl_stat create_amp_kernel_vars(struct amp_ppcg_kernel *kernel)
{
    int i, j, n;

    n = 0;
    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *array = &kernel->array[i];

        for (j = 0; j < array->n_group; ++j)
        {
            // struct amp_array_ref_group *group = array->groups[j];
            // enum ppcg_group_access_type type;

            // type = amp_array_ref_group_type(group);
            // if (type != ppcg_access_global)
            // ++n;
            ++n;
        }
    }

    kernel->var = isl_calloc_array(kernel->ctx, struct amp_ppcg_kernel_var, n);
    if (!kernel->var)
        return isl_stat_error;
    kernel->n_var = n;

    n = 0;
    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *array = &kernel->array[i];

        for (j = 0; j < array->n_group; ++j)
        {
            struct amp_array_ref_group *group = array->groups[j];
            enum ppcg_group_access_type type;

            type = amp_array_ref_group_type(group);
            // if (type == ppcg_access_global)
            //     continue;
            create_amp_kernel_var(kernel->ctx, group, &kernel->var[n]);
            ++n;
        }
    }

    return isl_stat_ok;
}

/* Return the set of outer array elements accessed by
 * by the statement instances in "domain" in "prog".
 * The instances in "domain" are those that appear
 * in the domains of the access relations in "prog".
 */
static __isl_give isl_union_set *accessed_by_domain(
    __isl_take isl_union_set *domain, struct amp_prog *prog)
{
    isl_union_map *access;
    isl_union_set *arrays;

    access = isl_union_map_union(isl_union_map_copy(prog->read), isl_union_map_copy(prog->may_write));
    access = isl_union_map_intersect_domain(access, domain);
    arrays = isl_union_map_range(access);
    arrays = isl_union_set_apply(arrays, isl_union_map_copy(prog->to_outer));

    return arrays;
}

/* Wrapper around ppcg_kernel_free for use as a isl_id_set_free_user callback.
 */
static void amp_ppcg_kernel_free_wrap(void *user)
{
    struct amp_ppcg_kernel *kernel = user;

    amp_ppcg_kernel_free(kernel);
}

/* Group the domain elements into a single space, named kernelX,
 * with X the kernel sequence number "kernel_id".
 */
static __isl_give isl_schedule_node *group_statements(
    __isl_take isl_schedule_node *node, int kernel_id)
{
    char buffer[20];
    isl_id *id;

    if (!node)
        return NULL;

    snprintf(buffer, sizeof(buffer), "kernel%d", kernel_id);
    id = isl_id_alloc(isl_schedule_node_get_ctx(node), buffer, NULL);
    return isl_schedule_node_group(node, id);
}

/* Replace "pa" by the zero function defined over the universe domain
 * in the space of "pa".
 */
static __isl_give isl_pw_aff *set_universally_zero(__isl_take isl_pw_aff *pa)
{
    isl_space *space;
    isl_aff *zero;

    space = isl_space_domain(isl_pw_aff_get_space(pa));
    isl_pw_aff_free(pa);
    zero = isl_aff_zero_on_domain(isl_local_space_from_space(space));

    return isl_pw_aff_from_aff(zero);
}

/* The sizes of the arrays on the host that have been computed by
 * extract_array_info may depend on the parameters.  Use the extra
 * constraints on the parameters that are valid at "host_domain"
 * to simplify these expressions and store the results in kernel->array.
 *
 * We only need these localized bounds for arrays that are accessed
 * by the current kernel.  If we have found at least one reference group
 * then the array is accessed by the kernel.
 *
 * The resulting sizes may be functions that are nowhere defined
 * in case the access function cannot possibly access anything inside
 * the kernel for some reason.  If so, they are replaced by the zero
 * function.  Since the access function cannot actually access anything,
 * there is no harm in printing the array sizes as zero.
 */
static void localize_bounds(struct amp_ppcg_kernel *kernel,
                            __isl_keep isl_set *host_domain)
{
    int i, j;
    isl_set *context;

    context = isl_set_copy(host_domain);
    context = isl_set_params(context);

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *local = &kernel->array[i];
        isl_multi_pw_aff *bound;
        int n_index;

        if (local->n_group == 0)
            continue;

        n_index = local->array->n_index;
        bound = isl_multi_pw_aff_copy(local->array->bound);

        for (j = 0; j < n_index; ++j)
        {
            isl_pw_aff *pwaff;
            int empty;

            pwaff = isl_multi_pw_aff_get_pw_aff(bound, j);
            pwaff = isl_pw_aff_gist(pwaff, isl_set_copy(context));
            empty = isl_pw_aff_is_empty(pwaff);
            if (empty < 0)
                pwaff = isl_pw_aff_free(pwaff);
            else if (empty)
                pwaff = set_universally_zero(pwaff);
            bound = isl_multi_pw_aff_set_pw_aff(bound, j, pwaff);
        }

        local->n_index = n_index;
        local->bound = bound;
    }
    isl_set_free(context);
}

/* If max_shared_memory is not set to infinity (-1), then make
 * sure that the total amount of shared memory required by the
 * array reference groups mapped to shared memory by "kernel"
 * is no larger than this maximum.
 *
 * We apply a greedy approach and discard (keep in global memory)
 * those groups that would result in a total memory size that
 * is larger than the maximum.
 *
 * This function should be called after any function that may
 * affect the decision on whether to place a reference group
 * in private, shared or global memory.
 */
static void check_shared_memory_bound(struct amp_ppcg_kernel *kernel)
{
    // #define DEBUG_CHECK_SHRED_MEMORY_BOUND

    int i, j;
    isl_val *left, *size;

    // left = isl_val_int_from_si(kernel->ctx,kernel->options->max_shared_memory);
    left = isl_val_int_from_si(kernel->ctx, 40960000);

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *local = &kernel->array[i];

        for (j = 0; j < local->n_group; ++j)
        {
            struct amp_array_ref_group *group;
            enum ppcg_group_access_type type;

            group = local->groups[j];
            type = amp_array_ref_group_type(group);
            if (type != ppcg_access_shared)
                continue;

            size = amp_array_tile_size(group->shared_tile);
            size = isl_val_mul_ui(size, local->array->size);
#ifdef DEBUG_CHECK_SHRED_MEMORY_BOUND
            fprintf(stderr, "\n the size and left is: \n");
            isl_val_dump(size);
            isl_val_dump(left);
#endif // DEBUG_CHECK_SHRED_MEMORY_BOUND

            if (isl_val_le(size, left))
            {
                left = isl_val_sub(left, size);
                continue;
            }
            isl_val_free(size);

            // 这里也尝试新增了注释
            // fprintf(stderr, "\n\n@WARN\n      check_shared_memory_bound meets errors!\n\n");
            group->shared_tile = amp_array_tile_free(group->shared_tile);
        }
    }

    isl_val_free(left);
}

/* Mark all arrays of "kernel" that have an array reference group
 * that is not mapped to private or shared memory as
 * accessing the corresponding global device memory.
 */
static void mark_data_copy_arrays(struct amp_ppcg_kernel *kernel)
{
    int i, j;

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *local = &kernel->array[i];

        if (local->global)
            continue;
        for (j = 0; j < local->n_group; ++j)
        {
            if (amp_array_ref_group_tile(local->groups[j]))
                continue;

            local->global = 1;
            local->array->global = 1;
            break;
        }
    }
}

/* Given a description of an array tile "tile" and the "space"
 *
 *	{ D -> A }
 *
 * where D represents the first tile->depth schedule dimensions
 * and A represents the array, construct an isl_multi_aff
 *
 *	{ [D[i] -> A[a]] -> A'[a'] }
 *
 * with A' a scaled down copy of A according to the shifts and strides
 * in "tile".  In particular,
 *
 *	a' = (a + shift(i))/stride
 *
 * "insert_array" represents
 *
 *	{ [D -> A] -> D }
 *
 * and is used to insert A into the domain of functions that only
 * reference D.
 */
static __isl_give isl_multi_aff *strided_tile(
    struct amp_array_tile *tile, __isl_keep isl_space *space,
    __isl_keep isl_multi_aff *insert_array)
{
    int i;
    isl_ctx *ctx;
    isl_multi_aff *shift;
    isl_multi_val *stride;
    isl_space *space2;
    isl_local_space *ls;
    isl_multi_aff *tiling;

    ctx = isl_space_get_ctx(space);
    space2 = isl_space_domain(isl_space_copy(space));
    ls = isl_local_space_from_space(space2);
    space2 = isl_space_range(isl_space_copy(space));
    stride = isl_multi_val_zero(space2);
    shift = isl_multi_aff_zero(isl_space_copy(space));

    for (i = 0; i < tile->n; ++i)
    {
        struct amp_array_bound *bound = &tile->bound[i];
        isl_val *stride_i;
        isl_aff *shift_i;

        stride_i = isl_val_copy(bound->stride);
        shift_i = isl_aff_copy(bound->shift);

        stride = isl_multi_val_set_val(stride, i, stride_i);
        shift = isl_multi_aff_set_aff(shift, i, shift_i);
    }
    isl_local_space_free(ls);

    shift = isl_multi_aff_pullback_multi_aff(shift,
                                             isl_multi_aff_copy(insert_array));

    tiling = isl_multi_aff_range_map(isl_space_copy(space));
    tiling = isl_multi_aff_add(tiling, shift);
    tiling = isl_multi_aff_scale_down_multi_val(tiling, stride);

    return tiling;
}

/* Compute a tiling for the array reference group "group".
 *
 * The tiling is of the form
 *
 *	{ [D[i] -> A[a]] -> T[t] }
 *
 * where D represents the first tile->depth schedule dimensions,
 * A represents the global array and T represents the shared or
 * private memory tile.  The name of T is the name of the local
 * array.
 *
 * If there is any stride in the accesses, then the mapping is
 *
 *	t = (a + shift(i))/stride - lb(i)
 *
 * otherwise, it is simply
 *
 *	t = a - lb(i)
 */
static void amp_array_ref_group_compute_tiling(struct amp_array_ref_group *group)
{
    // #define DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING

    int i;
    struct amp_array_tile *tile;
    isl_space *space;
    isl_multi_aff *tiling, *lb, *insert_array;
    isl_printer *p;
    char *local_name;

    tile = amp_array_ref_group_tile(group);
    if (!tile)
    {
        // fprintf(stderr, "@WARN_INFO: \n       in the amp_array_ref_group_compute_tiling        function, the tile of group is null !!! please notice! the group is:         \n");
        fprintf(stderr, "@WARN_INFO: \n       in the         amp_array_ref_group_compute_tiling function, the tile of group is null         !!! please notice! the group is: \n");
        amp_array_ref_group_dump(group);
        // amp_array_tile_dump(tile);
        fprintf(stderr, "\n\n");
        return;
    }
#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       the start of amp_array_ref_group_compute_tiling, the tile is:\n");
    amp_array_tile_dump(tile);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING

    space = isl_map_get_space(group->access);
    space = isl_space_from_range(isl_space_range(space));
    space = isl_space_add_dims(space, isl_dim_in, tile->depth);
    insert_array = isl_multi_aff_domain_map(isl_space_copy(space));
#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       the insert_array is:\n");
    isl_multi_aff_dump(insert_array);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING

    for (i = 0; i < tile->n; ++i)
        if (tile->bound[i].shift)
            break;

    if (i < tile->n)
        tiling = strided_tile(tile, space, insert_array);
    else
        tiling = isl_multi_aff_range_map(isl_space_copy(space));

#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       when i(%d) and tile->n(%d), tiling is:\n", i, tile->n);
    isl_multi_aff_dump(tiling);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING

    lb = isl_multi_aff_zero(space);
    for (i = 0; i < tile->n; ++i)
    {
        isl_aff *lb_i = isl_aff_copy(tile->bound[i].lb);
        lb = isl_multi_aff_set_aff(lb, i, lb_i);
    }
    lb = isl_multi_aff_pullback_multi_aff(lb, insert_array);

#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       the final lb is: \n");
    isl_multi_aff_dump(lb);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    if (lb)
        tiling = isl_multi_aff_sub(tiling, lb);

#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       the tiling - lb is: \n");
    isl_multi_aff_dump(tiling);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING

    p = isl_printer_to_str(isl_multi_aff_get_ctx(tiling));
    p = amp_array_ref_group_print_name(group, p);
    local_name = isl_printer_get_str(p);
    isl_printer_free(p);
    tiling = isl_multi_aff_set_tuple_name(tiling, isl_dim_out, local_name);
    free(local_name);

#ifdef DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    fprintf(stderr, "@DEBUG: \n       the end of amp_array_ref_group_compute_tiling, the tile->tiling is:\n");
    isl_multi_aff_dump(tiling);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_ARRAY_REFGROUP_COMPUTE_TILING
    tile->tiling = tiling;
}

/* Compute a tiling for all the array reference groups in "kernel".
 */
static void compute_group_tilings(struct amp_ppcg_kernel *kernel)
{
    // #define DEBUG_COMPUTE_GROUP_TILINGS

    int i, j;

    for (i = 0; i < kernel->n_array; ++i)
    {
        struct amp_local_array_info *array = &kernel->array[i];

        for (j = 0; j < array->n_group; ++j)
        {
            amp_array_ref_group_compute_tiling(array->groups[j]);

#ifdef DEBUG_COMPUTE_GROUP_TILINGS
            fprintf(stderr, "@DEBUG: \n       in end of compute_group_tilings,the n_array is %d,         the n_group is %d! \n", kernel->n_array, array->n_group);
            fprintf(stderr, "       the kernel->array[%d]->groups[%d]->array->name is %s .\n", i, j, array->groups[j]->array->name);
            // fprintf(stderr, "       the kernel->array[%d]->groups[%d]->shared_tile is: \n", i, j);
            // amp_array_tile_dump(array->groups[j]->shared_tile);
            // fprintf(stderr, "\n\n");
#endif // DEBUG_COMPUTE_GROUP_TILINGS
        }
    }
}

/* Create a ppcg_kernel representing the domain instances that reach "node"
 * and insert a mark node pointing to the ppcg_kernel before "node".
 * The band that "node" points to is the band that needs to be mapped
 * to block identifiers.  The band that needs to be mapped to thread
 * identifiers should be marked by a "thread" mark by the caller.
 * The linear branch between the current node and the "thread" mark
 * may also have a "shared" mark.  If present, the mapping to shared
 * memory is computed at that point.
 * Both marks are removed by this function.
 * If "scale" is set, then the band that "node" points to is scaled
 * by "sizes".
 *
 * Mark all outer band nodes as atomic to ensure each kernel is only
 * scheduled once.
 * If the domain elements that reach "node" live in more than one space,
 * then group the domain elements into a single space, named kernelX,
 * with X the kernel sequence number.
 *
 * Insert a guard node governing the kernel node to ensure that
 * no kernels with zero blocks are launched.
 *
 * Insert a context node describing the block and thread
 * identifiers inside the kernel mark.
 * The context node needs to be inserted after the effective block size
 * has been determined such that the bounds on the thread identifiers
 * would reflect the effective block size.
 * Insert a filter node inside the context node mapping the statement
 * instances to block identifiers.  In particular, the block identifiers
 * are equated to the partial schedule of band that was marked for mapping
 * to blocks modulo the grid size.
 * Insert a filter node inside the "thread" mark mapping the statement
 * instances to thread identifiers.  In particular, the thread identifiers
 * are equated to the partial schedule of band that was marked for mapping
 * to threads modulo the block size.
 *
 * Compute array reference groups for all arrays, set the local
 * array bounds based on the set of domain instances that reach
 * the kernel node, check the total amount of shared memory used
 * and compute all group tilings.
 * The array reference groups are computed after the block filter
 * has been inserted because it affects the mapping to shared or
 * private memory.  This computation also requires the thread filter
 * (in the ppcg_kernel object), but this thread filter should not
 * have been added to the schedule tree yet since the computation
 * requires the schedule of the band that needs to be mapped to
 * threads before the privatization is applied.
 *
 * If any array reference group requires the band mapped to threads
 * to be unrolled, then we perform the required unrolling.
 *
 * We save a copy of the schedule that may influence the mappings
 * to shared or private memory in kernel->copy_schedule.
 *
 * Finally, we add synchronization and copy statements to the schedule tree,
 * remove the "thread" mark and create representations for the local
 * variables in the kernel.
 *
 * We keep a copy of the isl_id that points to the kernel to ensure
 * that the kernel does not get destroyed if the schedule node
 * is freed due to some error condition.
 */
__isl_give isl_schedule_node *amp_create_kernel(struct amp_prog *prog, __isl_take isl_schedule_node *node) // __isl_keep isl_multi_val *sizes)
{
    // #define DEBUG_AMP_CREATE_KERNEL

    struct amp_ppcg_kernel *kernel;
    isl_id *id;
    isl_schedule_node *node_thread;
    isl_union_map *host_schedule;
    isl_union_pw_multi_aff *contraction;
    isl_set *host_domain;
    isl_union_set *domain, *expanded;
    int single_statement;

    // isl_schedule_node_dump(node);
    node = amp_tree_insert_shared_before_thread(node);
#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       after insert shared mark, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL
    if (!node)
        return NULL;

    kernel = isl_calloc_type(prog->ctx, struct amp_ppcg_kernel);
    kernel = amp_ppcg_kernel_create_local_arrays(kernel, prog);
    if (!kernel)
        return isl_schedule_node_free(node);

    domain = isl_schedule_node_get_domain(node);
    single_statement = isl_union_set_n_set(domain) == 1;
#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       amp_create_kernel get the domain is: \n");
    isl_union_set_dump(domain);
    fprintf(stderr, "\n       the sigle_statementt is :%d \n\n", single_statement);
#endif // DEBUG_AMP_CREATE_KERNEL

    kernel->ctx = prog->ctx;
    kernel->prog = prog;
    kernel->options = prog->scop->options;
    kernel->context = extract_context(node, prog);
    kernel->core = isl_union_set_universe(isl_union_set_copy(domain));
    contraction = isl_schedule_node_get_subtree_contraction(node);
    kernel->contraction = isl_union_pw_multi_aff_copy(contraction);
    expanded = isl_union_set_copy(domain);
    expanded = isl_union_set_preimage_union_pw_multi_aff(expanded, contraction);
    kernel->expanded_domain = isl_union_set_copy(expanded);
    kernel->arrays = accessed_by_domain(expanded, prog);

    kernel->id = prog->kernel_id++;

    host_schedule = isl_schedule_node_get_prefix_schedule_union_map(node);
    host_domain = isl_set_from_union_set(isl_union_map_range(host_schedule));

    // 插入mark结点
    node = atomic_ancestors(node);
    node = isl_schedule_node_child(node, 0);
    id = isl_id_alloc(prog->ctx, "amp_kernel", kernel);
    id = isl_id_set_free_user(id, &amp_ppcg_kernel_free_wrap);
    node = isl_schedule_node_insert_mark(node, isl_id_copy(id));

#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       after insert amp_kernel, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    // if (!single_statement)
    //     node = group_statements(node, kernel->id);

    // node = isl_schedule_node_child(node, 0);
    // node = insert_context(kernel, node);
    // node = isl_schedule_node_child(node, 0);
    // node = isl_schedule_node_insert_filter(node,
    // 			    isl_union_set_copy(kernel->block_filter))

    node = amp_tree_move_up_to_kernel(node);

    if (amp_group_references(kernel, node) < 0)
        node = isl_schedule_node_free(node);
    localize_bounds(kernel, host_domain);
    isl_set_free(host_domain);

    // check_shared_memory_bound(kernel);
    mark_data_copy_arrays(kernel);
    compute_group_tilings(kernel);

    node = amp_tree_move_down_to_thread(node, kernel->core);
    kernel->copy_schedule_dim = isl_schedule_node_get_schedule_depth(node);
    kernel->copy_schedule = isl_schedule_node_get_prefix_schedule_union_pw_multi_aff(node);
    contraction = isl_union_pw_multi_aff_copy(kernel->contraction);
    kernel->copy_schedule = isl_union_pw_multi_aff_pullback_union_pw_multi_aff(kernel->copy_schedule, contraction);
    // node = amp_tree_move_up_to_kernel(node);

    node = amp_tree_move_up_to_kernel(node);
#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       before the amp add copies function, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    node = amp_add_copies(kernel, node);

#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       after the amp add copies function and before delete thread and shared mark, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    // node = isl_schedule_node_parent(node);
    node = amp_tree_move_down_to_shared(node, kernel->core);
    node = isl_schedule_node_delete(node);

    node = amp_tree_move_down_to_thread(node, kernel->core);
    node = isl_schedule_node_delete(node);

    node = amp_tree_move_up_to_kernel(node);

#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       after delete thread and shared mark, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    if (create_amp_kernel_vars(kernel) < 0)
        node = isl_schedule_node_free(node);

    if (!single_statement)
        node = isl_schedule_node_parent(node);
        // node = isl_schedule_node_parent(node);

#ifdef DEBUG_AMP_CREATE_KERNEL
    fprintf(stderr, "@DEBUG: \n       after the amp add create kernrl function, the node is: \n");
    isl_schedule_node_dump(node);
    fprintf(stderr, "\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    return node;
}

/**
 * @brief   为了更好的划分amp划分后的filter设计的结构体.
 *          left代表划分后上面(对应树中的左分支)的filter(isl_basic_set),lower代表下面(对应树中的右分支)的filter.
 *          flag代表对应的filter是否存在,true即为存在.
 * @note
 *      注意：left非空并不是就一定存在该迭代空间，其存在与否是通过其flag控制的。
 *      当rate = 0,或者取整后的partition_val = 0时，left_flag应该为false
 *      当rate = 100,right_flag应该为false.(目前这种情况直接返回原始调度了,不会进行迭代空间划分.所以不会出现这种情况)
 */
struct amp_basic_set
{
    // the left basic set and its flag
    isl_basic_set *left;
    isl_bool left_flag;

    // the right basic set and its flag
    isl_basic_set *right;
    isl_bool right_flag;

    // split need data
    isl_ctx *ctx;
    int rate;
    isl_size split_dim;
};

/**
 * @brief 初始化amp_basic_set
 */
__isl_give struct amp_basic_set *amp_basic_set_init(__isl_keep isl_ctx *ctx)
{
    if (!ctx)
        goto error;
    struct amp_basic_set *amp_bset = isl_calloc_type(ctx, struct amp_basic_set);
    amp_bset->left = NULL;
    amp_bset->right = NULL;
    amp_bset->left_flag = isl_bool_error;
    amp_bset->right_flag = isl_bool_error;

    amp_bset->ctx = ctx;
    amp_bset->rate = 50;
    amp_bset->split_dim = -1;

    return amp_bset;

error:
    return NULL;
}

void amp_basic_set_free(__isl_take struct amp_basic_set *amp_bset)
{
    if (!amp_bset)
        fprintf(stderr, "\n@WARN:       amp_basic_set is NULL, maybe it is freed early! \n");
    else
    {
        isl_basic_set_free(amp_bset->left);
        isl_basic_set_free(amp_bset->right);
        free(amp_bset);
    }
}

/**
 * @brief 计算自动混合精度依据rate进行划分,rate是高精度所占百分比的数值，rate = 50,即高精度占比是50/100(50%)
 */
__isl_give isl_aff *amp_get_partition_aff(__isl_keep isl_ctx *ctx, __isl_take isl_aff *x, __isl_take isl_aff *y, int rate)
{
    // #define DEBUG_AMP_GET_PARTITION_AFF

    // x,y,rate 是否有异常值
    if (!x || !y || rate < 0 || rate > 100)
        goto error;

    isl_aff *length = isl_aff_add(isl_aff_copy(y), isl_aff_copy(x));
#ifdef DEBUG_AMP_GET_PARTITION_AFF
    fprintf(stderr, "@DEBUG: \n       区间长度是： \n");
    isl_aff_dump(length);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_GET_PARTITION_AFF

    isl_val *rate_val = isl_val_div_ui(isl_val_int_from_ui(ctx, rate), 100);
#ifdef DEBUG_AMP_GET_PARTITION_AFF
    fprintf(stderr, "@DEBUG: \n       高精度所占比例是： \n");
    isl_val_dump(rate_val);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_GET_PARTITION_AFF

    isl_aff *scale_down = isl_aff_scale_val(length, rate_val);
#ifdef DEBUG_AMP_GET_PARTITION_AFF
    fprintf(stderr, "@DEBUG: \n       length * rate_val 后是： \n");
    isl_aff_dump(scale_down);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_GET_PARTITION_AFF

    isl_aff *floor = isl_aff_floor(scale_down);
#ifdef DEBUG_AMP_GET_PARTITION_AFF
    fprintf(stderr, "@DEBUG: \n       向下取整后的结果是： \n");
    isl_aff_dump(floor);
    fprintf(stderr, "\n");
#endif // DEBUG_AMP_GET_PARTITION_AFF

    isl_aff_free(x);
    isl_aff_free(y);

    return floor;
error:
    isl_aff_free(x);
    isl_aff_free(y);

    return NULL;
}

/*
 * check_loop_index_constraint_on_bound_pair function's  data
 */
struct check_date
{
    isl_size **constraint_array;
    isl_size dims;

    isl_size current_dim;
};

/**
 * 检查循环索引之间的依赖关系
 */
static isl_stat check_loop_index_constraint_on_bound_pair(__isl_take isl_constraint *lower,
                                                          __isl_take isl_constraint *upper,
                                                          __isl_take isl_basic_set *bset,
                                                          void *user)
{
    // #define DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
    struct check_date *data = (struct check_date *)user;
    isl_size nvar;

    nvar = isl_basic_set_dim(bset, isl_dim_set);
#ifdef DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       current_dim is : %d\n", data->current_dim);
    fprintf(stderr, "       lower constraint and upper constraint is :\n");
    isl_constraint_dump(lower);
    isl_constraint_dump(upper);
    fprintf(stderr, "\n");
    fprintf(stderr, "@DEBUG: \n       isl_basic_set(Remove the current dimension) is :\n");
    isl_basic_set_dump(bset);
    fprintf(stderr, "       n var of isl_basic_set is : %d \n\n", nvar);
#endif // DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
    if (nvar < 0)
        goto error;

    for (isl_size i = 0; i < data->current_dim; i++)
    {
        if (isl_constraint_involves_dims(lower, isl_dim_set, i, 1) == 1)
        {
#ifdef DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
            fprintf(stderr, "@DEBUG: \n       The lower bound of %d dimension depends on the %d dimension .\n\n", data->current_dim, i);
#endif // DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
            data->constraint_array[data->current_dim][i] = 1;
        }

        if (isl_constraint_involves_dims(upper, isl_dim_set, i, 1) == 1)
        {
#ifdef DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
            fprintf(stderr, "@DEBUG: \n       The upper bound of %d dimension depends on the %d dimension .\n\n", data->current_dim, i);
#endif // DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
            data->constraint_array[i][data->current_dim] = 1;
        }
    }

#ifdef DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR
    // 打印出来数组。
    fprintf(stderr, "@DEBUG: \n       the loop index constraint array is:\n");
    for (int i = 0; i < data->dims; i++)
    {
        // 打印表头
        if (i == 0)
        {
            fprintf(stderr, "dim\\dim\t");
            for (int j = 0; j < data->dims; j++)
            {
                fprintf(stderr, "%d \t", j);
            }
            fprintf(stderr, "\n");
        }
        // 打印列号
        fprintf(stderr, "%d \t", i);
        // 打印内容
        for (int j = 0; j < data->dims; j++)
        {
            fprintf(stderr, "%d \t", data->constraint_array[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif // DEBUG_CHECK_LOOP_INDEX_CONSTRAINT_ON_BOUND_PAIR

    isl_constraint_free(upper);
    isl_constraint_free(lower);
    isl_basic_set_free(bset);

    return isl_stat_ok;
error:
    isl_constraint_free(upper);
    isl_constraint_free(lower);
    isl_basic_set_free(bset);

    return isl_stat_error;
}

/**
 * @brief   通过struct check_date中的数据，依据划分策略，获得划分的维度.
 * @note    目前的迭代空间划分的维度选择的策略是：循环索引不被内层循环索引所依赖的最外层循环对应的维度.
 */
isl_size get_split_dim_by_loop_index_constraint_relation_array(struct check_date *data)
{
    for (isl_size i = 0; i < data->dims; i++)
    {
        isl_bool upper_flag = isl_bool_true;
        isl_bool lower_flag = isl_bool_true;
        for (isl_size j = i + 1; j < data->dims; j++)
        {
            // 先判断其后面循环的上界是否依赖于当前的第i维度
            if (data->constraint_array[i][j] == 1)
            {
                upper_flag = isl_bool_false;
            }
            // 再判断其后面循环的下界是否依赖于当前的第i维度
            if (data->constraint_array[j][i] == 1)
            {
                lower_flag = isl_bool_false;
            }
        }
        if (upper_flag && lower_flag)
        {
            return i;
        }
    }
    return -1;
}

/**
 * @brief   通过该维度的上下界的约束（lower和upper）进行迭代空间的划分
 * @note    划分的方法是：
 */
static isl_stat split_on_bound_pair(__isl_take isl_constraint *lower,
                                    __isl_take isl_constraint *upper, __isl_take isl_basic_set *bset,
                                    void *user)
{
    // #define DEBUG_SPLIT_ON_BOUND_PAIR
    struct amp_basic_set *amp_bset = (struct amp_basic_set *)user;
    isl_aff *aff_x, *aff_y;

#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       split_on_bound_pair 接收到的一些参数有： \n");
    isl_basic_set_dump(bset);
    fprintf(stderr, "\n");
    isl_constraint_dump(upper);
    isl_constraint_dump(lower);
    fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

    // 获取当前维度的左边约束(循环下界)的仿射表达式：[] -> { S[] -> [(-x + t)]}.示例中当前维度对应的循环索引为t，x为下界值
    aff_x = isl_constraint_get_aff(lower);
    // 获取当前维度的右边约束(循环上界)的仿射表达式: [] -> { S[] -> [(y - t)]}.示例中当前维度对应的循环索引为t，y为上界值
    aff_y = isl_constraint_get_aff(upper);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       切分维度的左约束的仿射表达式是： \n");
    isl_aff_dump(aff_x);
    fprintf(stderr, "\n       切分维度的右约束的仿射表达式是： \n");
    isl_aff_dump(aff_y);
    fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

    // 获取左右边约束对应维度的系数，如果不一致，则缩放使其一致。
    isl_val *upper_val = isl_constraint_get_coefficient_val(upper, isl_dim_set, amp_bset->split_dim);
    // 对upper_val取反，使得系数的符号一致(均为+)
    upper_val = isl_val_neg(upper_val);
    isl_val *lower_val = isl_constraint_get_coefficient_val(lower, isl_dim_set, amp_bset->split_dim);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       计算amp划分值前,upper和lower在划分维度(%d)的系数是： \n", amp_bset->split_dim);
    isl_val_dump(upper_val);
    isl_val_dump(lower_val);
    fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR
    if (!isl_val_eq(upper_val, lower_val))
    {
        isl_val *scale_val = isl_val_div(lower_val, upper_val);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       统一系数时的缩放因子是： \n");
        isl_val_dump(scale_val);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR
        aff_y = isl_aff_scale_val(aff_y, scale_val);
    }
    else
    {
        isl_val_free(upper_val);
        isl_val_free(lower_val);
    }
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       统一系数后, 切分维度的左约束的仿射表达式是： \n");
    isl_aff_dump(aff_x);
    fprintf(stderr, "\n       统一系数后, 切分维度的右约束的仿射表达式是： \n");
    isl_aff_dump(aff_y);
    fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

    /**
     * @brief 获取amp_partition_val(根据rate划分出来高精度计算所需的长度值，也叫amp划分值)
     *           即求floor( (y-x) * rate/100 ),数学上则是[ (y-x) * rate/100 ]
     * @note 注意这里是向下取整
     */
    isl_aff *amp_partition_val = amp_get_partition_aff(amp_bset->ctx, isl_aff_copy(aff_x), isl_aff_copy(aff_y), amp_bset->rate);
    if (!amp_partition_val)
        goto error;
    assert(amp_partition_val);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
    fprintf(stderr, "@DEBUG: \n       amp_get_partition_aff 返回的值是： \n");
    isl_aff_dump(amp_partition_val);
    fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

    // 如果amp_partition_val为0，那么左分支（第一个迭代空间）应该为空
    if (isl_aff_plain_is_zero(amp_partition_val))
    {
        amp_bset->left_flag = isl_bool_false;
        amp_bset->right_flag = isl_bool_true;
    }
    // 如果amp_partition_val不为0，那么左分支和左分支两个迭代空间均为非空
    else
    {
        amp_bset->left_flag = isl_bool_true;
        amp_bset->right_flag = isl_bool_true;
        /**
         * @brief 先获取新的左右分支需要用到的仿射表达式.
         *
         * @note 显然，左分支需要求新y,右分支需要求新x.但是两者的边界其实只差一个1.
         *       所以先求新y,然后就顺理成章得到了新x.
         */
        isl_aff *aff_left = isl_aff_copy(aff_x);
        isl_aff *aff_right = isl_aff_copy(aff_y);
        isl_aff_free(aff_x);
        isl_aff_free(aff_y);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       得到amp划分值后,拼接前的左分支的仿射表达式是： \n");
        isl_aff_dump(aff_left);
        fprintf(stderr, "       得到amp划分值后,拼接前的右分支的仿射表达式是： \n");
        isl_aff_dump(aff_right);
        fprintf(stderr, "       amp_partition_val 是： \n");
        isl_aff_dump(amp_partition_val);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

        // 获取(拼接)新的左分支的新仿射表达式（右约束，新y）：[] -> { S[] -> [(  amp_partition_val + x - t )]}
        aff_left = isl_aff_sub(isl_aff_copy(amp_partition_val), aff_left);
        /**
         * @brief 获取(拼接)新的右分支的新仿射表达式（左约束,新x）：[] -> { S[] -> [( -x - amp_partition_val - 1 + t )]}
         *            这里是先用 (-x + t) - amp_partition_val, 然后再判断两个空间是否有重合的元素(取反策略肯定必有)。如果有重合元素，再减 1 求得最终结果。
         * @note  注意分成了两步，第一步是获得(拼接)了aff_right = [] -> { S[] -> [( -amp_partition_val - x + t )]}，第二步只需减1即可(去除切平面或者切线上的点)
         * @note  也即：切线或者切平面上的点,只会出现在左分支上(第一个迭代空间上),不会出现在右分支上.
         */
        // // 先获取新的右分支的仿射表达式（左约束,新x）：[] -> { S[] -> [( -amp_partition_val - x + t )]}
        // // aff_right = isl_aff_sub(isl_aff_copy(amp_partition_val), aff_right);
        /**
         * @brief   为了方便，也可以直接对aff_left取反(但该方法会出现重复).
         * @note    也即：切线或者切平面上的点会同时出现在左右两个迭代空间造成重复,所以还需要第二步来去除(第二个迭代空间中所包含的)切平面或者切线上的点.
         */
        aff_right = isl_aff_neg(isl_aff_copy(aff_left));
        isl_aff_free(amp_partition_val);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       拼接后,去除切平面或者切线上的点前.左分支的仿射表达式是： \n");
        isl_aff_dump(aff_left);
        fprintf(stderr, "       拼接后,去除切平面或者切线上的点前.右分支的仿射表达式是： \n");
        isl_aff_dump(aff_right);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR
       // 如果左右分支包含了相同的元素，则右分支去掉这个元素(取反策略必有)
        isl_set *common_set = isl_aff_eq_set(isl_aff_copy(aff_left), isl_aff_copy(aff_right));
        if (!isl_set_is_empty(common_set))
        {
            // 初始化来获取与aff_right同空间的常数为1的仿射表达式
            isl_space *space = isl_aff_get_domain_space(aff_right);
            isl_val *coeff_val = isl_constraint_get_coefficient_val(lower, isl_dim_set, amp_bset->split_dim);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
            fprintf(stderr, "@DEBUG: \n       去除切平面上的点时, 公共的集合和减去的值依次是： \n");
            isl_set_dump(common_set);
            isl_val_dump(coeff_val);
            fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR
            isl_aff *one = isl_aff_val_on_domain_space(space, coeff_val);
            // 减1 得到最终的仿射表达式
            aff_right = isl_aff_sub(aff_right, one);
        }
        isl_set_free(common_set);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       切分的维度是： %d .\n", amp_bset->split_dim);
        fprintf(stderr, "       新的左分支的仿射表达式是： \n");
        isl_aff_dump(aff_left);
        fprintf(stderr, "       新的右分支的仿射表达式是： \n");
        isl_aff_dump(aff_right);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

        // 将左右的新的仿射表达式转换成约束
        isl_constraint *c_left = isl_inequality_from_aff(aff_left);
        isl_constraint *c_right = isl_inequality_from_aff(aff_right);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       新的左分支的右约束是：  \n");
        isl_constraint_dump(c_left);
        fprintf(stderr, "\n      新的右分支的左约束是：  \n");
        isl_constraint_dump(c_right);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR

        // 将新约束添加到左右的filter中
        amp_bset->left = isl_basic_set_add_constraint(amp_bset->left, c_left);
        amp_bset->right = isl_basic_set_add_constraint(amp_bset->right, c_right);
#ifdef DEBUG_SPLIT_ON_BOUND_PAIR
        fprintf(stderr, "@DEBUG: \n       左filter是:  \n");
        isl_basic_set_dump(amp_bset->left);
        fprintf(stderr, "\n      右filter是:  \n");
        isl_basic_set_dump(amp_bset->right);
        fprintf(stderr, "\n");
#endif // DEBUG_SPLIT_ON_BOUND_PAIR
    }

    isl_constraint_free(upper);
    isl_constraint_free(lower);
    isl_basic_set_free(bset);

    return isl_stat_ok;
error:
    isl_constraint_free(upper);
    isl_constraint_free(lower);
    isl_basic_set_free(bset);
    isl_aff_free(amp_partition_val);
    isl_aff_free(aff_x);
    isl_aff_free(aff_y);

    return isl_stat_error;
}

/**
 * @brief 根据划分比例，将基本语句basic_set的迭代空间划分成两个迭代空间，并返回
 */
static __isl_give struct amp_basic_set *partition_by_rate(__isl_keep isl_ctx *ctx, __isl_take isl_basic_set *basic_set, int rate)
{
    // #define DEBUG_PARTITION_BY_RATE
    struct amp_basic_set *amp_bset = amp_basic_set_init(ctx);

    // 检查basic_set不能为空
    if (!basic_set)
    {
        goto error;
    }

    // 复制basic_set,初始化amp_bset的左右分支的两个basic_set
    amp_bset->left = isl_basic_set_copy(basic_set);
    amp_bset->right = isl_basic_set_copy(basic_set);

    // 获取basic_set的维度数,也即当前语句是在S_0[c0,c1,c2]是在几层循环里
    isl_size bset_dims = isl_basic_set_dim(basic_set, isl_dim_set);
    // 声明一个 basic_set_dims * basic_set_dims的二维数组，用于存储循环索引间的依赖关系
    isl_size **loop_index_constraint_relation_array = isl_alloc_array(ctx, isl_size *, bset_dims);
    /**
     * @brief   初始化数组为0
     * @note    对角线必须是0，表示自己对自己不存在依赖关系.
     */
    for (isl_size i = 0; i < bset_dims; i++)
    {
        loop_index_constraint_relation_array[i] = isl_alloc_array(ctx, isl_size, bset_dims);
        for (isl_size j = 0; j < bset_dims; j++)
        {
            loop_index_constraint_relation_array[i][j] = 0;
        }
    }

    // 申请check_date
    struct check_date *check_date = isl_alloc_array(ctx, struct check_date, 1);
    check_date->constraint_array = loop_index_constraint_relation_array;
    check_date->dims = bset_dims;

    /**
     * @brief   遍历每一个维度的上界和下界的约束，检查当前维度的索引是否对前面的循环索引存在依赖关系.
     * @note    如果检查异常结束（也即除去该维度后的isl_basic_set为空）则提示一个ERROR.
     */
    for (isl_size i = 1; i < bset_dims; i++)
    {
        check_date->current_dim = i;
        if (isl_basic_set_foreach_bound_pair(basic_set, isl_dim_set, i, &check_loop_index_constraint_on_bound_pair, check_date) < 0)
        {
            fprintf(stderr, "\n\033[31m@ERROR:\n       check_loop_index_constraint_on_bound_pair meets some ERRORS!!!  \033[0m\n\n");
        }
    }

#ifdef DEBUG_PARTITION_BY_RATE
    // 打印出来数组。
    fprintf(stderr, "@DEBUG: \n       the loop index constraint array is:\n");
    for (int i = 0; i < bset_dims; i++)
    {
        // 打印表头
        if (i == 0)
        {
            fprintf(stderr, "dim\\dim\t");
            for (int j = 0; j < bset_dims; j++)
            {
                fprintf(stderr, "%d \t", j);
            }
            fprintf(stderr, "\n");
        }
        // 打印列号
        fprintf(stderr, "%d \t", i);
        // 打印内容
        for (int j = 0; j < bset_dims; j++)
        {
            fprintf(stderr, "%d \t", loop_index_constraint_relation_array[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif // DEBUG_PARTITION_BY_RATE

    // 获取应该划分的维度
    isl_size split_dim = get_split_dim_by_loop_index_constraint_relation_array(check_date);
    assert(split_dim != -1);
#ifdef DEBUG_PARTITION_BY_RATE
    fprintf(stderr, "@DEBUG: \n       the split dim should be %d .\n\n", split_dim);
#endif // DEBUG_PARTITION_BY_RATE

    // 保存split时需要用到的信息
    amp_bset->rate = rate;
    amp_bset->split_dim = split_dim;
    // 划分split_dim对应的维度
    if (isl_basic_set_foreach_bound_pair(basic_set, isl_dim_set, split_dim, &split_on_bound_pair, amp_bset) < 0)
        fprintf(stderr, "\n\033[31m@ERROR:\n       split_on_bound_pair function meets some ERRORS!!!  \033[0m\n\n");
#ifdef DEBUG_PARTITION_BY_RATE
    fprintf(stderr, "@DEBUG: \n       对当前isl_basic_set进行迭代空间划分后,左迭代空间是:  \n");
    isl_basic_set_dump(amp_bset->left);
    fprintf(stderr, "       对当前isl_basic_set进行迭代空间划分后,右迭代空间是:  \n");
    isl_basic_set_dump(amp_bset->right);
    fprintf(stderr, "\n");
#endif // DEBUG_PARTITION_BY_RATE

    isl_basic_set_free(basic_set);

    return amp_bset;
error:
    isl_basic_set_free(basic_set);
    amp_basic_set_free(amp_bset);

    return NULL;
}

/**
 * @brief 存储迭代空间划分后的左右两个迭代空间(isl_union_set)以及划分的比例
 */
struct amp_domain
{
    isl_ctx *ctx;
    int rate;

    /**
     * left和right分别代表划分后的‘左’和‘右’两个迭代空间.
     */
    isl_union_set *left;
    isl_union_set *right;
};

/**
 * This function is called for each set in a union_set.
 * If the dimension of the set is 0, we store the set in the both left and right.
 * else the dimension of the set > 0, we split the set, and then store it(which is splited) in left and right.
 * @note if the dimension of the set < 0, that means errors.
 */
static isl_stat repartition_set(__isl_take isl_set *set, void *user)
{
    struct amp_domain *amp_domain = (struct amp_domain *)user;
    isl_size dims;

    // 获取该语句(isl_set)的维度数
    dims = isl_set_n_dim(set);
    if (dims < 0)
        goto error;
    // 如果dims为0，说明语句就是简单的赋值，且不在循环体中，则将其同时放到左右两个子迭代空间中
    else if (dims == 0)
    {
        // 如果右分支的迭代空间为NULL,说明还未放入,直接转换即可
        if (!amp_domain->right)
        {
            // 如果左分支（第一个迭代空间）存在.
            if (amp_domain->rate > 0)
                amp_domain->left = isl_union_set_from_set(set);
            amp_domain->right = isl_union_set_from_set(set);
        }
        else
        {
            // 如果左分支（第一个迭代空间）存在.
            if (amp_domain->rate > 0)
                amp_domain->left = isl_union_set_add_set(amp_domain->left, set);
            amp_domain->right = isl_union_set_add_set(amp_domain->right, set);
        }
    }
    // 如果dims大于0，说明语句循环体中，则将其切分后放到左右两个子迭代空间中
    else
    {
        /**
         * 获取basic_set_list及其大小,遍历所有的basic_set并对其进行划分
         * @note 由于之前已经压缩过，对所有语句size均应为1).
         */
        isl_basic_set_list *bset_list = isl_set_get_basic_set_list(set);
        isl_size sizes = isl_basic_set_list_n_basic_set(bset_list);
        for (isl_size i = 0; i < sizes; i++)
        {
            isl_basic_set *bset = isl_basic_set_list_get_basic_set(bset_list, i);
            // 根据rate将当前的一个basic_set划分成left和right两个，也即将该语句的迭代空间分成2个
            struct amp_basic_set *amp_bset = partition_by_rate(amp_domain->ctx, bset, amp_domain->rate);
            /** 根据情况进行判断，添加该语句的迭代空间到amp_domain中. */
            if (!amp_domain->right)
            {
                // 如果左分支（第一个迭代空间）存在.
                if (amp_bset->left_flag)
                    amp_domain->left = isl_union_set_from_set(isl_set_from_basic_set(isl_basic_set_copy(amp_bset->left)));
                // 如果右分支（第二个迭代空间）存在，目前一定是存在的.
                if (amp_bset->right_flag)
                    amp_domain->right = isl_union_set_from_set(isl_set_from_basic_set(isl_basic_set_copy(amp_bset->right)));
            }
            else
            {
                // 如果左分支（第一个迭代空间）存在.
                if (amp_bset->left_flag)
                    amp_domain->left = isl_union_set_add_set(amp_domain->left, isl_set_from_basic_set(isl_basic_set_copy(amp_bset->left)));
                // 如果右分支（第二个迭代空间）存在，目前一定是存在的.
                if (amp_bset->right_flag)
                    amp_domain->right = isl_union_set_add_set(amp_domain->right, isl_set_from_basic_set(isl_basic_set_copy(amp_bset->right)));
            }
            amp_basic_set_free(amp_bset);
        }
    }
    isl_set_free(set);

    return isl_stat_ok;
error:
    isl_set_free(set);

    return isl_stat_error;
}

/**
 * @brief AMP重新调度，首先是拆分isl_union_set变成上下两个部分合在一起的isl_union_set_list;然后插入不同的mark标签结点;再然后根据插入的mark结点，生成对应的kernel,并删去多余的mark结点.
 *           最后在将修改好的调度树，转换成调度(isl_schedule)返回,如果在这里遇到问题，则返回最原始的调度。
 */
__isl_give isl_schedule *amp_reschedule(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule *sched, int rate)
{
    // #define DEBUG_AMP_RESCHEDULE
    isl_schedule *schedule;
    isl_schedule_node *node;
    isl_union_set *domain;

    // 切分的次数
    int division_number = 1;

    // 获取调度树的根结点,以便后续在这里插入新的调度
    node = isl_schedule_get_root(sched);
    // 获取domain(isl_union_set)
    domain = isl_schedule_get_domain(sched);

    for (int d = 0; d < division_number; d++)
    {
        // 初始化amp_domian()
        struct amp_domain *amp_domain;
        amp_domain = isl_alloc_array(ctx, struct amp_domain, 1);
        amp_domain->ctx = ctx;
        amp_domain->rate = rate;
        amp_domain->left = NULL;
        amp_domain->right = NULL;

        // 遍历所有的语句，对其迭代空间进行划分.
        if (isl_union_set_foreach_set(domain, repartition_set, amp_domain) < 0)
        {
            fprintf(stderr, "\n\033[31m@ERROR:\n       There are some errors because isl_union_set_foreach_set() < 0 !!! \033[0m\n\n");
            goto error;
        }
        assert(amp_domain);

        // 将两个迭代空间拼接成filters(isl_union_set_list)
        isl_union_set_list *filters = isl_union_set_list_alloc(ctx, 0);
        // 如果存在第一个迭代空间，就加入
        if (amp_domain->left)
            filters = isl_union_set_list_add(filters, amp_domain->left);
        filters = isl_union_set_list_add(filters, amp_domain->right);

        // 将filters插入到根节点
        node = isl_schedule_node_child(node, 0);
        node = isl_schedule_node_insert_sequence(node, filters);
        if (!node)
        {
            fprintf(stderr, "\n\033[31m@ERROR:\n       There are some errors because insert sequence filed !!!  Now will return the original schedule !!! \033[0m\n\n");
            goto error;
        }
        // 如果存在第一个迭代空间
        if (amp_domain->left)
            node = isl_schedule_node_child(node, 1);
        // 如果不存在第一个迭代空间
        else
            node = isl_schedule_node_child(node, 0);
        node = isl_schedule_node_child(node, 0);

        // 插入mark(thread)结点
        isl_id *id = isl_id_alloc(ctx, "thread", NULL);
        node = isl_schedule_node_insert_mark(node, id);
        // 返回到插入前的位置（filter的位置）
        // node = isl_schedule_node_parent(node);
#ifdef DEBUG_AMP_RESCHEDULE
        fprintf(stderr, "@DEBUG: \n       after insert mark('thread'), the node is : \n");
        isl_schedule_node_dump(node);
        fprintf(stderr, "\n\n");
#endif

        // 插入mark(amp_lower)结点
        node = insert_amp_lower(node);
        // 返回到插入前的位置（filter的位置）
        node = isl_schedule_node_parent(node);
#ifdef DEBUG_AMP_RESCHEDULE
        fprintf(stderr, "@DEBUG: \n       after insert mark('amp_lower'), the node is : \n");
        isl_schedule_node_dump(node);
        fprintf(stderr, "\n\n");
#endif
        // 如果要增加多次切分，这里需要更新node和domain,但是更建议把这里做成一个函数，或许会更好一点，因为按理说division_number=2,最终的结果应该是4个filter.这个目前先不考虑，后续再设计实现吧
    }

    // 根据amp_lower标记,生成amp_kernel（自动混合精度计算核心）
    node = amp_create_kernel(prog, node);
    if (!node)
    {
        fprintf(stderr, "\n\033[31m@ERROR:\n       There are some errors because the node (the amp_create_kernel function returned) is NULL, Now will return the original schedule !!! \033[0m\n\n");
        goto error;
    }
    assert(node);

    // 获取新的带有amp_kernel的调度
    schedule = isl_schedule_node_get_schedule(node);
#ifdef DEBUG_AMP_RESCHEDULE
    fprintf(stderr, "@DEBUG: \n       automatic_mixed_precision_reschedule generate's schedule  is : \n");
    isl_schedule_dump(schedule);
    fprintf(stderr, "\n\n");
#endif

    // 释放内存
    isl_schedule_free(sched);
    isl_schedule_node_free(node);
    isl_union_set_free(domain);

    return schedule;
error:
    isl_schedule_node_free(node);
    isl_union_set_free(domain);

    return sched;
}

// Compute a new schedule based on the sched
__isl_give isl_schedule *amp_schedule_again(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule *sched)
{
    int rate;
    isl_schedule *schedule;

    // 检查参数
    if ((!sched) || (!prog))
        goto error;

    // 如果不进行自动混合精度,则直接返回原始的调度即可.--这个检查是多余的,主要是检查一下,确保参数的正确性
    if (!prog->scop->options->automatic_mixed_precision)
    {
        return sched;
    }

    // 获取混合精度的比例
    rate = (int)prog->scop->options->automatic_mixed_precision_rate;
    // 如果比例不合法(小于0或者大于99),则不进行混合精度,将原始调度返回.
    if (rate < 0 || rate > 99)
    {
        fprintf(stderr, "\n\033[31m@WARNING:\n       automatic mixed precision rate < 0 / rate > 99 , that is incorrect. \033[0m\n\n");
        return sched;
    }

    schedule = amp_reschedule(ctx, prog, sched, rate);
    if (!schedule)
    {
        fprintf(stderr, "\n\033[31m@ERROR:\n       There are some errors because the schedule (the amp_reschedule function returned) is NULL, Now will return the original schedule !!! \033[0m\n\n");
        goto error;
    }
    assert(schedule);

    return schedule;
error:
    return NULL;
}
