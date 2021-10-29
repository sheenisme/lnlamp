#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include <isl/polynomial.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
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
    info->linearize = prog->scop->options->linearize_device_arrays;

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
        info->linearize = 1;
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

static __isl_give isl_printer *print_amp_macros(__isl_take isl_printer *p)
{
    const char *macros =
        "#define ampCheckReturn(ret) \\\n"
        "  if (ret != AMP_SUCCESS) {\\\n"
        "    fprintf(stderr, \"AMP error: %s\\n\", "
        "amp_error_string(ret)); \\\n"
        "    fflush(stderr); \\\n"
        "    assert(ret == AMP_SUCCESS);\\\n  }\n";

    p = isl_printer_print_str(p, macros);

    p = isl_printer_start_line(p);
    p = isl_printer_end_line(p);

    return p;
}

/* Set the names of the macros that may appear in a printed isl AST.
 */
__isl_give isl_printer *amp_print_macros(__isl_take isl_printer *p)
{
    p = print_amp_macros(p);

    return p;
}

/** 获得比当前精度更低的数据类型 **/
char *amp_get_lower_precision_type(char *type)
{
    // 生成更低精度的
    if (strcmp("double", type) == 0)
        return "float";
    else
        return "bfloat";
}

/* Does "array" need to be allocated on the device?
 * If it is a read-only scalar, then it will be passed as an argument
 * to the kernel and therefore does not require any allocation.
 * If this device memory is not accessed at all, then it does not
 * need to be allocated either.
 */
int amp_array_requires_allocation(struct amp_array_info *array)
{
    if (amp_array_is_read_only_scalar(array))
        return 0;
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

        if (!amp_array_requires_allocation(array))
            continue;

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

/* Print a declaration for the amp array corresponding to "array" on "p".
 */
__isl_give isl_printer *declare_amp_lower_precision_array(__isl_take isl_printer *p, struct amp_array_info *array)
{
    int i;

    p = isl_printer_start_line(p);
    p = isl_printer_print_str(p, amp_get_lower_precision_type(array->type)); // 换成更低精度的类型
    p = isl_printer_print_str(p, " ");
    if (array->n_index > 1)
        p = isl_printer_print_str(p, "(");
    p = isl_printer_print_str(p, "*amp_lower_");
    p = isl_printer_print_str(p, array->name);
    if (array->n_index > 1)
    {
        p = isl_printer_print_str(p, ")");
        for (i = 1; i < array->n_index; i++)
        {
            isl_ast_expr *bound;
            bound = isl_ast_expr_get_op_arg(array->bound_expr, 1 + i);
            p = isl_printer_print_str(p, "[");
            p = isl_printer_print_ast_expr(p, bound);
            p = isl_printer_print_str(p, "]");
            isl_ast_expr_free(bound);
        }
    }
    p = isl_printer_print_str(p, ";");
    p = isl_printer_end_line(p);

    return p;
}

__isl_give isl_printer *declare_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog)
{
    int i;

    for (i = 0; i < prog->n_array; ++i)
    {
        if (!amp_array_requires_allocation(&prog->array[i]))
            continue;

        p = declare_amp_lower_precision_array(p, &prog->array[i]);
    }
    p = isl_printer_start_line(p);
    p = isl_printer_end_line(p);

    return p;
}

/* Print an expression for the size of "array" in bytes.
 */
__isl_give isl_printer *amp_array_info_print_size(__isl_take isl_printer *prn, struct amp_array_info *array)
{
    int i;

    for (i = 0; i < array->n_index; ++i)
    {
        isl_ast_expr *bound;

        prn = isl_printer_print_str(prn, "(");
        bound = isl_ast_expr_get_op_arg(array->bound_expr, 1 + i);
        prn = isl_printer_print_ast_expr(prn, bound);
        isl_ast_expr_free(bound);
        prn = isl_printer_print_str(prn, ") * ");
    }
    prn = isl_printer_print_str(prn, "sizeof(");
    prn = isl_printer_print_str(prn, amp_get_lower_precision_type(array->type)); // 换成更低精度的
    prn = isl_printer_print_str(prn, ")");

    return prn;
}

__isl_give isl_printer *allocate_amp_lower_precision_arrays(__isl_take isl_printer *p, amp_prog *prog)
{
    int i;

    for (i = 0; i < prog->n_array; ++i)
    {
        struct amp_array_info *array = &prog->array[i];

        if (!amp_array_requires_allocation(&prog->array[i]))
            continue;
        p = ppcg_ast_expr_print_macros(array->bound_expr, p);
        p = isl_printer_start_line(p);
        p = isl_printer_print_str(p, "ampCheckReturn(ampMalloc((void **) &amp_lower_");
        p = isl_printer_print_str(p, prog->array[i].name);
        p = isl_printer_print_str(p, ", ");
        p = amp_array_info_print_size(p, &prog->array[i]);
        p = isl_printer_print_str(p, "));");
        p = isl_printer_end_line(p);
    }
    p = isl_printer_start_line(p);
    p = isl_printer_end_line(p);

    return p;
}

/**
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * dew
 * 
 */

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

#ifdef DEBUG_ATOMIC_ANCESTORS
    printf("@DEBUG: \n       at the start of the atomic_ancestors function,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
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
    printf("@DEBUG: \n       at the end of the atomic_ancestors function,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_ATOMIC_ANCESTORS

    return node;
}

// find the suitable position for lower precision calculate insert, in the node.
static __isl_give isl_schedule_node *find_amp_lower_precision_calculate_insert_position_of_node(
    __isl_take isl_schedule_node *node)
{
#define DEBUG_FIND_AMP_LOWER_PRECISION_POSITION

    // isl_schedule_node *node = atomic_ancestors(orig_node);
    int depth = isl_schedule_node_get_schedule_depth(node);

#ifdef DEBUG_FIND_AMP_LOWER_PRECISION_POSITION
    printf("@DEBUG: \n       at the start of the find_amp_lower_precision_calculate_insert_position_of_node function,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n the depth of the node is: %d \n\n", depth);
#endif // DEBUG_FIND_AMP_LOWER_PRECISION_POSITION

    if (!node)
        return NULL;
    if (!isl_schedule_node_has_parent(node))
        return node;

    // for (int i = 0; i < depth; i++)
    // {
    //     node = isl_schedule_node_get_child(node, i);
    //     if (isl_schedule_node_get_type(node) == isl_schedule_node_band)
    //     {
    //         break;
    //     }
    // }
    node = isl_schedule_node_get_child(isl_schedule_node_get_child(node, 1), 0);

    // isl_schedule_node_free(orig_node);

    // node = isl_schedule_node_child(node, pos);

#ifdef DEBUG_FIND_AMP_LOWER_PRECISION_POSITION
    printf("@DEBUG: \n       at the end of the find amp lower location function,the retuened node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_FIND_AMP_LOWER_PRECISION_POSITION

    return node;
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
    kernel->array = isl_calloc_array(ctx,
                                     struct amp_local_array_info, prog->n_array);
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
        access = isl_union_map_union(access,
                                     isl_union_map_from_map(map_i));
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
        access = isl_union_map_union(access,
                                     isl_union_map_from_map(map_i));
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
    int empty;
    isl_union_pw_multi_aff *tagger;
    isl_union_set *domain, *access_domain;
    isl_union_map *local, *external, *universe;
    isl_union_set *tag_set;

    if (isl_union_map_is_empty(access))
    {
        isl_union_map_free(sched);
        isl_union_map_free(tagged);
        return access;
    }

    tagger = isl_union_pw_multi_aff_copy(prog->scop->tagger);
    domain = isl_union_map_domain(isl_union_map_copy(tagged));
    tagger = isl_union_pw_multi_aff_intersect_domain(tagger,
                                                     isl_union_set_copy(domain));
    sched = isl_union_map_preimage_domain_union_pw_multi_aff(sched, tagger);

    local = isl_union_map_apply_range(sched,
                                      isl_union_map_reverse(isl_union_map_copy(sched)));
    local = isl_union_map_intersect(local,
                                    isl_union_map_copy(prog->scop->tagged_dep_flow));

    empty = isl_union_map_is_empty(local);

    external = isl_union_map_copy(prog->scop->tagged_dep_flow);
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
    external = isl_union_map_intersect_params(external,
                                              isl_set_copy(prog->scop->context));
    external = isl_union_map_subtract(external, local);

    if (read)
    {
        tag_set = isl_union_map_range(external);
        external = wrapped_reference_to_access(tag_set, tagged);
        external = isl_union_map_union(external,
                                       isl_union_map_copy(prog->scop->live_in));
    }
    else
    {
        tag_set = isl_union_map_domain(external);
        external = wrapped_reference_to_access(tag_set, tagged);
        external = isl_union_map_union(external,
                                       isl_union_map_copy(prog->scop->live_out));
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
    isl_union_map *sched, *tagged;

    if (isl_union_map_is_empty(access))
        return access;

    tagged = amp_group_tagged_access_relation(group);
    sched = isl_union_map_copy(prefix);

    return amp_remove_local_accesses(kernel->prog, tagged, access, sched, read);
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
        if (isl_schedule_node_get_type(node) ==
            isl_schedule_node_band)
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
    isl_union_map *access;
    isl_union_map *prefix;

    prefix = isl_schedule_node_get_prefix_schedule_relation(node);
    prefix = isl_union_map_preimage_domain_union_pw_multi_aff(prefix,
                                                              isl_union_pw_multi_aff_copy(kernel->contraction));
    access = amp_array_ref_group_access_relation(group, read, !read);
    access = amp_remove_local_accesses_group(kernel, group, access, prefix,
                                             read);
    access = isl_union_map_range_product(prefix, access);

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
static __isl_give isl_multi_aff *create_from_access(isl_ctx *ctx,
                                                    struct amp_array_ref_group *group, int read)
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

static __isl_give isl_schedule_node *amp_add_copies_group(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
#define DEBUG_AMP_ADD_COPIES_GROUP

    // struct amp_array_tile *tile;
    isl_union_map *access;
    isl_union_set *domain;
    isl_multi_aff *ma;
    isl_multi_aff *from_access;
    isl_multi_pw_aff *mpa;
    isl_multi_union_pw_aff *mupa;
    isl_schedule_node *graft;
    isl_union_set *filter;
    int skip;
    int kernel_depth;
    int empty;

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       at start of the amp add copies group function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n       the read is :%d \n\n", read);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    // tile = gpu_array_ref_group_tile(group);
    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    access = anchored_non_local_accesses(kernel, group, node, read);
    empty = isl_union_map_is_empty(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the access[1] is: \n");
    isl_union_map_dump(access);
    printf("\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        // if (empty < 0)
        //     return isl_schedule_node_free(node);
        // return isl_schedule_node_parent(node);
        return node;
    }

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the from access variable is: \n");
    isl_multi_aff_dump(from_access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    ma = isl_multi_aff_copy(from_access);
    mpa = isl_multi_pw_aff_from_multi_aff(ma);
    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);
    domain = isl_union_map_range(access);

    // if (read && !amp_array_is_scalar(group->array))
    // {
    //     isl_map *map;
    //     isl_union_set_free(domain);
    //     map = group_tile(group);
    //     domain = isl_union_set_from_set(isl_map_wrap(map));
    // }

    domain = isl_union_set_preimage_multi_aff(domain, from_access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the domain is: \n");
    isl_union_set_dump(domain);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    access = isl_union_set_wrapped_domain_map(domain);
    access = isl_union_map_reverse(access);
    access = isl_union_map_coalesce(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       before the graft = isl_schedule_node_from_extension(access);the access is: \n");
    isl_union_map_dump(access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    graft = isl_schedule_node_from_extension(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after the graft = isl_schedule_node_from_extension(access);the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, mupa);
    /*
        if (kernel->options->unroll_copy_shared)
            graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    if (tile->n > kernel->n_block && kernel->n_block > 0)
    {
        graft = isl_schedule_node_band_split(graft,
                                             tile->n - kernel->n_block);
        graft = isl_schedule_node_child(graft, 0);
    }
    if (tile->n < kernel->n_block)
        skip = kernel->n_block - tile->n;
    else
        skip = 0;
    filter = set_schedule_modulo(graft, kernel->thread_ids,
                                 kernel->block_dim);
    if (!kernel->options->wrap)
        graft = snap_band_to_sizes(graft, kernel->block_dim + skip,
                                   kernel->options);
    if (tile->n > kernel->n_block && kernel->n_block > 0)
        graft = isl_schedule_node_parent(graft);
    graft = isl_schedule_node_insert_filter(graft, filter);
    */

    while (graft && isl_schedule_node_has_parent(graft))
        graft = isl_schedule_node_parent(graft);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the final graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    // node = find_amp_lower_precision_calculate_insert_position_of_node(node);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       before insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    if (read)
    {
        // node = amp_tree_move_left_to_sync(node, kernel);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_before(node, graft);
    }
    else
    {
        // node = amp_tree_move_right_to_sync(node, kernel);
        // node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);
        node = isl_schedule_node_graft_after(node, graft);
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

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
    isl_ctx *ctx = isl_union_map_get_ctx(data->full_sched);

    n = 0;
    for (i = 0; i < local->array->n_ref; ++i)
    {
        isl_union_map *umap;
        isl_map *map;
        struct amp_array_ref_group *group;
        struct amp_stmt_access *access = local->array->refs[i];

        map = isl_map_copy(access->access);
        umap = isl_union_map_from_map(map);

        // umap = isl_union_map_apply_domain(umap, isl_union_map_copy(data->copy_sched));

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
    // isl_ctx *ctx = isl_space_get_ctx(group->array->space);
    // isl_union_map *access, *local;
    // int n_index = group->array->n_index;
    // int no_reuse, coalesced;
    // isl_map *acc;
    int force_private = group->local_array->force_private;
    // isl_stat r = isl_stat_ok;
    // isl_bool ok;
    // int requires_unroll;
    // int unique_depth;

    if (amp_array_is_read_only_scalar(group->array))
        return isl_stat_ok;
    if (!force_private && !group->exact_write)
        return isl_stat_ok;
    if (group->slice)
        return isl_stat_ok;

    // access = amp_array_ref_group_access_relation(group, 1, 1);
    // local = localize_access(data, isl_union_map_copy(access));
    // no_reuse = isl_union_map_is_injective(local);
    // if (no_reuse < 0)
    //     r = isl_stat_error;
    // if (use_shared && no_reuse)
    //     coalesced = access_is_coalesced(data, local);
    // isl_union_map_free(local);

    // if (r >= 0 && kernel->options->debug->verbose &&
    //     use_shared && no_reuse && coalesced)
    //     report_no_reuse_and_coalesced(kernel, access);

    // if (use_shared && (!no_reuse || !coalesced))
    // {
    //     group->shared_tile = gpu_array_tile_create(ctx,
    //                                                group->array->n_index);
    //     acc = shared_access(group, access, data);
    //     ok = can_tile(acc, group->shared_tile);
    //     if (ok < 0)
    //         r = isl_stat_error;
    //     else if (!ok)
    //         group->shared_tile =
    //             gpu_array_tile_free(group->shared_tile);
    //     isl_map_free(acc);
    // }

    // if (r < 0 || (!force_private && (!use_private || no_reuse)))
    // {
    //     isl_union_map_free(access);
    //     return r;
    // }

    // access = isl_union_map_apply_domain(access,
    //                                     isl_union_map_copy(data->thread_sched));

    // acc = isl_map_from_union_map(access);

    // if (!force_private && !access_is_bijective(data, acc))
    // {
    //     isl_map_free(acc);
    //     return isl_stat_ok;
    // }

    // unique_depth = compute_accessed_by_single_thread_depth(data, acc);

    // acc = isl_map_intersect_domain(acc, isl_set_copy(data->privatization));
    // acc = isl_map_project_out(acc, isl_dim_in, data->thread_depth,
    //                           data->n_thread);
    // requires_unroll = check_requires_unroll(data, acc, force_private);
    // if (unique_depth < 0 || requires_unroll < 0 ||
    //     (requires_unroll && kernel->any_force_private))
    // {
    //     isl_map_free(acc);
    //     return requires_unroll < 0 ? isl_stat_error : isl_stat_ok;
    // }

    // group->private_tile = gpu_array_tile_create(ctx, n_index);
    // group->private_tile->requires_unroll = requires_unroll;
    // ok = can_tile(acc, group->private_tile);
    // if (ok >= 0 && !ok)
    //     group->private_tile = gpu_array_tile_free(group->private_tile);
    // isl_map_free(acc);
    // if (ok < 0)
    //     return isl_stat_error;

    // if (group->private_tile)
    // {
    //     struct gpu_array_tile *tile = group->private_tile;
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
    int i;
    int n;
    isl_ctx *ctx = isl_union_map_get_ctx(data->full_sched);
    struct amp_array_ref_group **groups;

    groups = isl_calloc_array(ctx, struct amp_array_ref_group *,
                              local->array->n_ref);
    if (!groups)
        return -1;

    n = populate_array_references(local, groups, data);

    if (local->array->has_compound_element && !local->force_private)
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

    // n = group_common_shared_memory_tile(kernel, local->array,
    //                                     n, groups, data);

    set_array_groups(local, n, groups);

    if (n >= 0)
        return 0;

    for (i = 0; i < local->array->n_ref; ++i)
        amp_array_ref_group_free(groups[i]);
    return -1;
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
    int i;
    int r = 0;
    isl_union_pw_multi_aff *contraction;
    struct amp_group_data data;

    // check_can_be_private_live_ranges(kernel, node);

    data.scop = kernel->prog->scop;

    // node = isl_schedule_node_copy(node);

    // node = gpu_tree_move_down_to_thread(node, kernel->core);
    // node = isl_schedule_node_child(node, 0);
    // data.thread_depth = isl_schedule_node_get_schedule_depth(node);
    // data.n_thread = isl_schedule_node_band_n_member(node);
    // if (data.thread_depth == data.shared_depth)
    //     data.copy_sched = isl_union_map_copy(data.shared_sched);
    // else
    //     data.copy_sched = prefix_with_equalities(node);
    // data.thread_sched = isl_union_map_copy(data.copy_sched);
    // data.thread_sched = isl_union_map_flat_range_product(data.thread_sched,
    //                                                      isl_schedule_node_band_get_partial_schedule_union_map(node));
    // data.thread_sched = isl_union_map_detect_equalities(data.thread_sched);

    // contraction = isl_union_pw_multi_aff_copy(kernel->contraction);
    // data.host_sched = expand(data.host_sched, contraction);
    // data.shared_sched = expand(data.shared_sched, contraction);
    // if (data.thread_depth == data.shared_depth)
    // {
    //     isl_union_map_free(data.copy_sched);
    //     data.copy_sched = isl_union_map_copy(data.shared_sched);
    // }
    // else
    // {
    //     data.copy_sched = expand(data.copy_sched, contraction);
    // }
    // data.thread_sched = expand(data.thread_sched, contraction);
    // isl_union_pw_multi_aff_free(contraction);

    // node = isl_schedule_node_child(node, 0);
    data.full_sched = isl_union_map_copy(isl_schedule_node_get_subtree_schedule_union_map(node));
    // data.full_sched = isl_union_map_flat_range_product(data.full_sched, isl_schedule_node_get_subtree_schedule_union_map(node));
    // isl_schedule_node_free(node);

    // compute_privatization(&data, kernel);

    for (i = 0; i < kernel->n_array; ++i)
    {
        r = amp_group_array_references(kernel, &kernel->array[i], &data);
        if (r < 0)
            break;
    }

    isl_union_map_free(data.full_sched);

    return r;
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

/* Print the name of the local copy of a given group of array references.
 */
__isl_give isl_printer *amp_array_ref_group_print_name(
    struct amp_array_ref_group *group, __isl_take isl_printer *p)
{
    int global = 0;

    p = isl_printer_print_str(p, "amp_lower_");

    global = 1;
    p = isl_printer_print_str(p, group->array->name);
    if (!global && group->local_array->n_group > 1)
    {
        p = isl_printer_print_str(p, "_");
        p = isl_printer_print_int(p, group->nr);
    }

    return p;
}

static void create_amp_kernel_var(isl_ctx *ctx, struct amp_array_ref_group *group,
                                  struct amp_ppcg_kernel_var *var)
{
    int j;
    isl_printer *p;

    var->array = group->array;

    p = isl_printer_to_str(ctx);
    p = amp_array_ref_group_print_name(group, p);
    var->name = isl_printer_get_str(p);
    isl_printer_free(p);

    var->size = isl_vec_alloc(ctx, group->array->n_index);
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
            struct amp_array_ref_group *group = array->groups[j];
            if (group)
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
__isl_give isl_schedule_node *amp_create_kernel(struct amp_prog *prog,
                                                __isl_take isl_schedule_node *node) // __isl_keep isl_multi_val *sizes)
{
#define DEBUG_AMP_CREATE_KERNEL
    struct amp_ppcg_kernel *kernel;
    isl_id *id;
    // isl_schedule_node *node_thread;
    // isl_union_map *host_schedule;
    isl_union_pw_multi_aff *contraction;
    //  isl_set *host_domain;
    isl_union_set *domain, *expanded;
    int single_statement;

    kernel = isl_calloc_type(prog->ctx, struct amp_ppcg_kernel);
    kernel = amp_ppcg_kernel_create_local_arrays(kernel, prog);
    if (!kernel)
        return isl_schedule_node_free(node);

    domain = isl_schedule_node_get_domain(node);
    single_statement = isl_union_set_n_set(domain) == 1;
#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       amp_create_kernel get the domain is: \n");
    isl_union_set_dump(domain);
    printf("\n       the sigle_statementt is :%d \n\n", single_statement);
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
    // node = isl_schedule_node_parent(node);

    // if (!single_statement)
    //     node = group_statements(node, id);

    kernel->copy_schedule_dim = isl_schedule_node_get_schedule_depth(node);
    kernel->copy_schedule = isl_schedule_node_get_prefix_schedule_union_pw_multi_aff(node);
    contraction = isl_union_pw_multi_aff_copy(kernel->contraction);
    kernel->copy_schedule = isl_union_pw_multi_aff_pullback_union_pw_multi_aff(kernel->copy_schedule, contraction);

    if (amp_group_references(kernel, node) < 0)
        node = isl_schedule_node_free(node);
        // node = isl_schedule_node_parent(node);
#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       before the find_amp_lower_calculate_insert_position_of_node function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    node = atomic_ancestors(node);
    node = find_amp_lower_precision_calculate_insert_position_of_node(node);
    id = isl_id_alloc(prog->ctx, "amp_kernel", kernel);
    node = isl_schedule_node_insert_mark(node, isl_id_copy(id));

#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       after the find_amp_lower_calculate_insert_position_of_node and isl_schedule_node_insert_mark function, the node & id is: \n");
    isl_schedule_node_dump(node);
    isl_id_dump(id);
    printf("\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       before the amp add copies function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    node = amp_add_copies(kernel, node);

#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       after the amp add copies function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    // node = isl_schedule_node_parent(node);

    if (create_amp_kernel_vars(kernel) < 0)
        node = isl_schedule_node_free(node);

    if (!single_statement)
        node = isl_schedule_node_parent(node);
        // node = isl_schedule_node_parent(node);

#ifdef DEBUG_AMP_CREATE_KERNEL
    printf("@DEBUG: \n       after the amp add create kernrl function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_CREATE_KERNEL

    return node;
}

/**
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * -----------------------------------
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

// 检查rate的合法性
isl_bool amp_check_rate(__isl_keep isl_constraint_list *cons_list, unsigned long rate)
{
    // #define DEBUG

    // d is the abs value of the sub of v0 and v1
    unsigned long d;
    isl_size constraint_list_dim = isl_constraint_list_n_constraint(cons_list);
    if (constraint_list_dim)
    {
        // isl_constraint_list -> isl_constraint
        isl_constraint *c0 = isl_constraint_list_get_constraint(cons_list, 0);
        isl_constraint *c1 = isl_constraint_list_get_constraint(cons_list, 1);
        isl_val *v0 = isl_constraint_get_constant_val(c0);
        isl_val *v1 = isl_constraint_get_constant_val(c1);
        isl_val *v = isl_val_sub(v1, v0);
        v = isl_val_abs(v);
        v = isl_val_add_ui(v, 1);
        d = (unsigned long)isl_val_get_num_si(v);

#ifdef DEBUG
        printf("@DEBUG: \n        d rate is: %ld 、%ld \n\n", d, rate);
#endif // DEBUG
        isl_constraint_free(c0);
        isl_constraint_free(c1);
        isl_val_free(v);
        if (!(d % rate))
        {
            return isl_bool_true;
        }
    }
    printf("\n\033[31m@ERROR:\n       the number of the rate is incorrect! this time the rate is %ld, \n       but the abs of the sub of v0 and v1 is %ld !!!\n\n\033[0m", rate, d);
    return isl_bool_error;
}

// amp修改调度
__isl_give isl_schedule *get_amp_schedule(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule *sched, unsigned long rate)
{
    // #define DEBUG_GET_AMP_SCHEDULE
    isl_schedule *schedule;
    isl_schedule_node *node;
    isl_union_set_list *uset_list;
    isl_union_set *domain, *uset1, *uset2;
    isl_set *set;
    isl_basic_set_list *bset_list;
    isl_basic_set *bset, *bset_new;
    isl_constraint_list *cons_list;
    isl_constraint *c, *c1, *c2, *c3;
    isl_size basic_set_list_dim, constraint_list_dim;
    // unsigned long rate = 2;

    // 获取调度树的 node 结点
    // node = isl_schedule_node_from_domain(isl_union_set_copy(domain)); 这里是不对的，具体区别还不确定。
    node = isl_schedule_get_root(sched);
    // 获取domain(isl_schedule -> isl_union_set)
    domain = isl_schedule_get_domain(sched);
    // 释放之前的调度
    isl_schedule_free(sched);
    // isl_union_set -> isl_set
    set = isl_set_from_union_set(domain);
    // isl_set -> isl_basic_set_list
    bset_list = isl_set_get_basic_set_list(set);
    basic_set_list_dim = isl_basic_set_list_n_basic_set(bset_list);
    // 初始化isl_union_set_list
    uset_list = isl_union_set_list_alloc(ctx, 2 * (int)basic_set_list_dim);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       domain、node、set 、basci_set_list、basic_set_list_dim is:\n");
    isl_union_set_dump(domain);
    isl_schedule_node_dump(node);
    isl_set_dump(set);
    isl_basic_set_list_dump(bset_list);
    printf("basic_set_list_dim is :%d \n\n", basic_set_list_dim);
#endif // DEBUG_GET_AMP_SCHEDULE
    isl_set_free(set);
    for (isl_size i = 0; i < basic_set_list_dim; ++i)
    {
        // isl_basic_set_list -> isl_basic_set
        bset = isl_basic_set_list_get_basic_set(bset_list, i);
        bset_new = isl_basic_set_list_get_basic_set(bset_list, i);
        isl_basic_set_list_free(bset_list);
        // isl_basic_set -> isl_constraint_list
        cons_list = isl_basic_set_get_constraint_list(bset);
        // isl_constraint_list *cons_list_new = isl_basic_set_get_constraint_list(bset);
        // isl_basic_set_free(bset);
        constraint_list_dim = isl_constraint_list_n_constraint(cons_list);

#ifdef DEBUG_GET_AMP_SCHEDULE
        printf("@DEBUG: \n       basci_set、isl_constraint_list、constraint_list_dim is:\n");
        isl_basic_set_dump(bset);
        isl_constraint_list_dump(cons_list);
        printf("constraint_list_dim is:%d \n\n", constraint_list_dim);
#endif
        // 依次获取，修改约束列表
        for (isl_size j = 0; j < constraint_list_dim; ++j)
        {
            // isl_constraint_list -> isl_constraint
            c = isl_constraint_list_get_constraint(cons_list, j);
            // isl_constraint_dump(c);

            if (j == 0) // 如果是amp的前面的一个原始计算
            {
                if (amp_check_rate(cons_list, rate) == isl_bool_error)
                {
                    printf("检查rate的合法性,遇到错误!\n");
                    goto error;
                }

                // create new val in new constraint
                isl_val *v = isl_constraint_get_constant_val(c);
                // calcute v1
                c1 = isl_constraint_list_get_constraint(cons_list, j + 1);
                isl_val *v1 = isl_constraint_get_constant_val(c1);
                v1 = isl_val_add_ui(v1, 1);
                v1 = isl_val_div_ui(v1, rate);
                // sub v1 in v
                v = isl_val_sub(v, v1);
                c2 = isl_constraint_set_constant_val(isl_constraint_copy(c), v);
                // cons_list_new = isl_constraint_list_drop(cons_list_new, j, 1);
                // cons_list_new = isl_constraint_list_add(cons_list_new, isl_constraint_copy(c2));
                bset_new = isl_basic_set_add_constraint(bset_new, c2);
                isl_constraint_free(c1);

#ifdef DEBUG_GET_AMP_SCHEDULE
                printf("@DEBUG: \n       created the new val 、isl_constraint、isl_constraint_list、isl_basic_set is:\n");
                // isl_val_dump(v);
                // isl_constraint_dump(c2);
                // isl_constraint_list_dump(cons_list_new);
                isl_basic_set_dump(bset_new);
                printf(" when j == %d \n\n", j);
#endif
            }
            else if (j == 1) // 如果是amp的前面的一个原始计算
            {
                // create new val of origion
                isl_val *v = isl_constraint_get_constant_val(c);
                v = isl_val_add_ui(v, 1);
                v = isl_val_div_ui(v, rate);
                v = isl_val_sub_ui(v, 1);
                c3 = isl_constraint_set_constant_val(isl_constraint_copy(c), v);
                // cons_list = isl_constraint_list_drop(cons_list, j, 1);
                // cons_list = isl_constraint_list_add(cons_list, isl_constraint_copy(c3));
                bset = isl_basic_set_add_constraint(bset, c3);

#ifdef DEBUG_GET_AMP_SCHEDULE
                printf("@DEBUG: \n       new val and revised isl_constraint、isl_constraint_list、isl_basic_set is:\n");
                // isl_val_dump(v);
                // isl_constraint_dump(c3);
                // isl_constraint_list_dump(cons_list);
                isl_basic_set_dump(bset);
                printf(" when j == %d \n\n", j);
#endif
            }
            isl_constraint_free(c);
        }
        // 释放掉constraint list
        isl_constraint_list_free(cons_list);

        // isl_basic_set -> isl_set
        uset1 = isl_union_set_from_basic_set(bset);
        uset2 = isl_union_set_from_basic_set(bset_new);
        uset_list = isl_union_set_list_add(uset_list, uset1);
        uset_list = isl_union_set_list_add(uset_list, uset2);
#ifdef DEBUG_GET_AMP_SCHEDULE
        printf("@DEBUG: \n       the two isl_union_set and the isl_union_set_list is:\n");
        // isl_union_set_dump(uset1);
        // isl_union_set_dump(uset2);
        isl_union_set_list_dump(uset_list);
        printf("when basic_set_list_dim = %d\n\n", basic_set_list_dim);
#endif
    }

    node = isl_schedule_node_child(node, 0);
    node = isl_schedule_node_insert_sequence(node, uset_list);
    // amp 添加数据转换
    // node = amp_add_data_conversion(ctx, prog, node);

    // amp_create_kernel
    node = amp_create_kernel(prog, node);

    schedule = isl_schedule_node_get_schedule(node);
    isl_schedule_node_free(node);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       amp generate's schedule  is:  \n");
    isl_schedule_dump(sched);
    printf("\n\n");
#endif

    return schedule;
error:
    isl_constraint_free(c);
    isl_constraint_list_free(cons_list);
    isl_basic_set_free(bset);
    isl_basic_set_free(bset_new);
    isl_union_set_list_free(uset_list);
    isl_union_set_free(uset1);
    isl_union_set_free(uset2);
    isl_schedule_node_free(node);
    return NULL;
}

// Compute a new schedule based on the schd
__isl_give isl_schedule *amp_schedule_again(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule *schedule)
{
    unsigned long rate;

    if ((!schedule) || (!prog))
        return NULL;

    // 如果不进行自动混合精度，这个检查是多余的，主要是检查确保参数的正确性
    if (!prog->scop->options->automatic_mixed_precision)
    {
        return NULL;
    }

    // 获取混合精度的比例
    rate = (unsigned long)prog->scop->options->automatic_mixed_precision_rate;

    // 如果比例 <= 0, 则不进行混合精度，将原始调度返回。
    if (rate <= 0)
    {
        printf("\n\033[31m@ERROR:\n       automatic mixed precision rate is 0 , that is incorrect.\n\n\033[0m");
        return schedule;
    }

    schedule = get_amp_schedule(ctx, prog, schedule, rate);

    // if (prog->scop->options->tile)
    //     printf("tile 时候，amp再进行一次调度！\n");

    return schedule;
}