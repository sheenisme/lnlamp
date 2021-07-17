/*
 * Copyright 2010-2011 INRIA Saclay
 * Copyright 2021      Sven Verdoolaege
 *
 * Use of this software is governed by the MIT license
 *
 * Written by Sven Verdoolaege, INRIA Saclay - Ile-de-France,
 * Parc Club Orsay Universite, ZAC des vignes, 4 rue Jacques Monod,
 * 91893 Orsay, France
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include <isl/aff.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/constraint.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include "grouping.h"
#include "schedule.h"

/* Add parameters with identifiers "ids" to "set".
 */
static __isl_give isl_set *add_params(__isl_take isl_set *set,
	__isl_keep isl_id_list *ids)
{
	int i, n;
	unsigned nparam;

	n = isl_id_list_n_id(ids);

	nparam = isl_set_dim(set, isl_dim_param);
	set = isl_set_add_dims(set, isl_dim_param, n);

	for (i = 0; i < n; ++i) {
		isl_id *id;

		id = isl_id_list_get_id(ids, i);
		set = isl_set_set_dim_id(set, isl_dim_param, nparam + i, id);
	}

	return set;
}

/* Equate the dimensions of "set" starting at "first" to
 * freshly created parameters with identifiers "ids".
 * The number of equated dimensions is equal to the number of elements in "ids".
 */
static __isl_give isl_set *parametrize(__isl_take isl_set *set,
	int first, __isl_keep isl_id_list *ids)
{
	int i, n;
	unsigned nparam;

	nparam = isl_set_dim(set, isl_dim_param);

	set = add_params(set, ids);

	n = isl_id_list_n_id(ids);
	for (i = 0; i < n; ++i)
		set = isl_set_equate(set, isl_dim_param, nparam + i,
					isl_dim_set, first + i);

	return set;
}

/* Given a parameter space "space", create a set of dimension "len"
 * of which the dimensions starting at "first" are equated to
 * freshly created parameters with identifiers "ids".
 */
__isl_give isl_set *parametrization(__isl_take isl_space *space,
	int len, int first, __isl_keep isl_id_list *ids)
{
	isl_set *set;

	space = isl_space_set_from_params(space);
	space = isl_space_add_dims(space, isl_dim_set, len);
	set = isl_set_universe(space);

	return parametrize(set, first, ids);
}

/* Load and return a schedule from a file called "filename".
 */
static __isl_give isl_schedule *load_schedule(isl_ctx *ctx,
	const char *filename)
{
	FILE *file;
	isl_schedule *schedule;

	file = fopen(filename, "r");
	if (!file) {
		fprintf(stderr, "Unable to open '%s' for reading\n", filename);
		return NULL;
	}
	schedule = isl_schedule_read_from_file(ctx, file);
	fclose(file);

	return schedule;
}

/* Save the schedule "schedule" to a file called "filename".
 * The schedule is printed in block style.
 */
static void save_schedule(__isl_keep isl_schedule *schedule,
	const char *filename)
{
	FILE *file;
	isl_ctx *ctx;
	isl_printer *p;

	if (!schedule)
		return;

	file = fopen(filename, "w");
	if (!file) {
		fprintf(stderr, "Unable to open '%s' for writing\n", filename);
		return;
	}
	ctx = isl_schedule_get_ctx(schedule);
	p = isl_printer_to_file(ctx, file);
	p = isl_printer_set_yaml_style(p, ISL_YAML_STYLE_BLOCK);
	p = isl_printer_print_schedule(p, schedule);
	isl_printer_free(p);
	fclose(file);
}

/* Compute a schedule on the domain of "sc" that respects the schedule
 * constraints in "sc", without trying to combine groups of statements.
 */
__isl_give isl_schedule *ppcg_compute_non_grouping_schedule(
	__isl_take isl_schedule_constraints *sc, struct ppcg_options *options)
{
	if (options->debug->dump_schedule_constraints)
		isl_schedule_constraints_dump(sc);
	return isl_schedule_constraints_compute_schedule(sc);
}

/* Compute a schedule on the domain of "sc" that respects the schedule
 * constraints in "sc".
 *
 * "schedule" is a known correct schedule that is used to combine
 * groups of statements if options->group_chains is set.
 */
__isl_give isl_schedule *ppcg_compute_schedule(
	__isl_take isl_schedule_constraints *sc,
	__isl_keep isl_schedule *schedule, struct ppcg_options *options)
{
	if (options->group_chains)
		return ppcg_compute_grouping_schedule(sc, schedule, options);
	return ppcg_compute_non_grouping_schedule(sc, options);
}

/* Obtain a schedule, either by reading it form a file
 * or by computing it using "compute".
 * Also take care of saving the computed schedule and/or
 * dumping the obtained schedule if requested by the user.
 */
__isl_give isl_schedule *ppcg_get_schedule(isl_ctx *ctx,
	struct ppcg_options *options,
	__isl_give isl_schedule *(*compute)(void *user), void *user)
{
	isl_schedule *schedule;

	if (options->load_schedule_file) {
		schedule = load_schedule(ctx, options->load_schedule_file);
	} else {
		schedule = compute(user);
		if (options->save_schedule_file)
			save_schedule(schedule, options->save_schedule_file);
	}
	if (options->debug->dump_schedule)
		isl_schedule_dump(schedule);

	return schedule;
}

/* Mark all dimensions in the band node "node" to be of "type".
 */
__isl_give isl_schedule_node *ppcg_set_schedule_node_type(
	__isl_take isl_schedule_node *node, enum isl_ast_loop_type type)
{
	int i, n;

	n = isl_schedule_node_band_n_member(node);
	for (i = 0; i < n; ++i)
		node = isl_schedule_node_band_member_set_ast_loop_type(node, i,
							type);

	return node;
}

/* Try and shift the given band node to the origin.
 *
 * In particular, obtain the set of schedule values and
 * compute the element-wise minimal value, which may
 * depend on the parameters.
 * If this results in any piecewise expressions,
 * then do not perform any shifting as that may
 * very well make the resulting code more complicated.
 * Otherwise shift the band by the opposite of this minimal value.
 */
static __isl_give isl_schedule_node *shift_to_origin(
	__isl_take isl_schedule_node *node)
{
	isl_bool is_multi_aff;
	isl_space *space;
	isl_union_set *domain, *range;
	isl_union_map *partial;
	isl_set *min_domain;
	isl_set *set;
	isl_multi_pw_aff *min;
	isl_multi_aff *min_ma, *shift;
	isl_multi_union_pw_aff *mupa;

	partial = isl_schedule_node_band_get_partial_schedule_union_map(node);
	domain = isl_schedule_node_get_domain(node);
	range = isl_union_set_apply(isl_union_set_copy(domain), partial);
	space = isl_schedule_node_band_get_space(node);
	set = isl_union_set_extract_set(range, space);
	isl_union_set_free(range);
	domain = isl_union_set_universe(domain);

	min = isl_set_min_multi_pw_aff(set);
	min_domain = isl_multi_pw_aff_domain(isl_multi_pw_aff_copy(min));
	min = isl_multi_pw_aff_gist(min, min_domain);
	is_multi_aff = isl_multi_pw_aff_isa_multi_aff(min);
	if (is_multi_aff < 0 || !is_multi_aff) {
		isl_union_set_free(domain);
		isl_multi_pw_aff_free(min);
		if (is_multi_aff < 0)
			return isl_schedule_node_free(node);
		return node;
	}

	min_ma = isl_multi_pw_aff_as_multi_aff(min);
	shift = isl_multi_aff_neg(min_ma);
	mupa = isl_multi_union_pw_aff_multi_aff_on_domain(domain, shift);
	node = isl_schedule_node_band_shift(node, mupa);

	return node;
}

/* Tile "node" with tile sizes "sizes", but first try and shift the band
 * to the origin.
 * Tiling a band that does not start at the origin is likely
 * to result in initial partial tiles.
 */
__isl_give isl_schedule_node *ppcg_tile(__isl_take isl_schedule_node *node,
	__isl_take isl_multi_val *sizes)
{
	node = shift_to_origin(node);
	return isl_schedule_node_band_tile(node, sizes);
}
