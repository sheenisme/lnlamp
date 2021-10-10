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
static int is_read_only_scalar(struct amp_array_info *array, struct amp_prog *prog)
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
static isl_stat extract_array_info(struct amp_prog *prog, struct amp_array_info *info, struct pet_array *pa, __isl_keep isl_union_set *arrays)
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
static __isl_give isl_union_map *remove_independences(struct amp_prog *prog, struct amp_array_info *array, __isl_take isl_union_map *order)
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

static void free_array_info(struct amp_prog *prog)
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
static __isl_give isl_union_set *compute_may_persist(struct amp_prog *prog)
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

void *amp_prog_free(struct amp_prog *prog)
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
__isl_give isl_ast_node *amp_build_array_bounds(__isl_take isl_ast_node *node, struct amp_prog *prog, __isl_keep isl_ast_build *build)
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