/*
 * Copyright 2012 INRIA Paris-Rocquencourt
 * Copyright 2012 Ecole Normale Superieure
 *
 * Use of this software is governed by the MIT license
 *
 * Written by Tobias Grosser, INRIA Paris-Rocquencourt,
 * Domaine de Voluceau, Rocquenqourt, B.P. 105,
 * 78153 Le Chesnay Cedex France
 * and Sven Verdoolaege,
 * Ecole Normale Superieure, 45 rue d'Ulm, 75230 Paris, France
 */

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include <isl/aff.h>
#include <isl/ctx.h>
#include <isl/flow.h>
#include <isl/map.h>
#include <isl/ast_build.h>
#include <isl/schedule.h>
#include <isl/schedule_node.h>
#include <pet.h>

#include "ppcg.h"
#include "ppcg_options.h"
#include "cpu.h"
#include "print.h"
#include "schedule.h"
#include "util.h"

// 引入混合精度的头文件
#include "amp.h"

/* Representation of a statement inside a generated AST.
 *
 * "stmt" refers to the original statement.
 * "ref2expr" maps the reference identifier of each access in
 * the statement to an AST expression that should be printed
 * at the place of the access.
 */
struct ppcg_stmt {
	struct pet_stmt *stmt;

	isl_id_to_ast_expr *ref2expr;
};

/* Internal data structure for at_domain.
 *
 * "prog" represents the entire scop.
 * "kernel" points to the kernel to which the current schedule node
 * belongs.  It is set by before_mark and reset by after_mark.
 * It may be NULL if we are outside any kernel.
 */
struct ppcg_at_domain_data
{
	struct amp_prog *prog;
	struct amp_ppcg_kernel *kernel;
};

/* Internal data structure for the index and AST expression transformation
 * callbacks for pet_stmt_build_ast_exprs.
 *
 * "kernel" is the kernel for which are computing AST expressions and
 * may be NULL if we are not inside a kernel.
 * "accesses" is the list of gpu_stmt_access in the statement.
 * "iterator_map" expresses the statement iterators in terms of
 * the AST loop iterators.
 * "sched2copy" expresses the outer copy_schedule_dim dimensions of
 * the kernel schedule in terms of the AST loop iterators and
 * may be NULL if we are not inside a kernel.
 *
 * The following fields are set in transform_index and used in transform_expr.
 * "array" is the array that is being accessed.
 * "global" is set if the global array is accessed (rather than
 * shared/private memory).
 * "local_array" refers to information on the array specialized
 * to the current kernel.
 */
struct ppcg_transform_data
{
	struct amp_ppcg_kernel *kernel;
	struct amp_stmt_access *accesses;
	isl_pw_multi_aff *iterator_map;
	isl_pw_multi_aff *sched2copy;

	struct amp_array_info *array;
	int global;
	struct amp_local_array_info *local_array;
};

// enum ppcg_group_access_type
// {
// 	ppcg_access_global,
// 	ppcg_access_shared,
// 	ppcg_access_private
// };

enum ppcg_kernel_stmt_type
{
	ppcg_kernel_copy,
	ppcg_kernel_domain
};

/* Representation of special statements, in particular copy statements
 * and __syncthreads statements, inside a kernel.
 *
 * type represents the kind of statement
 *
 *
 * for ppcg_kernel_copy statements we have
 *
 * read is set if the statement should copy data from global memory
 * to shared memory or registers.
 *
 * index expresses an access to the array element that needs to be copied
 * local_index expresses the corresponding element in the tile
 *
 * array refers to the original array being copied
 * local_array is a pointer to the appropriate element in the "array"
 *	array of the ppcg_kernel to which this copy access belongs
 *
 *
 * for ppcg_kernel_domain statements we have
 *
 * stmt is the corresponding input statement
 *
 * n_access is the number of accesses in stmt
 * access is an array of local information about the accesses
 */
struct ppcg_kernel_stmt
{
	enum ppcg_kernel_stmt_type type;

	union
	{
		struct
		{
			int read;
			isl_ast_expr *index;
			isl_ast_expr *local_index;
			struct amp_array_info *array;
			struct amp_local_array_info *local_array;
		} c;
		struct
		{
			struct amp_stmt *stmt;
			isl_id_to_ast_expr *ref2expr;
		} d;
	} u;
};

static void ppcg_stmt_free(void *user)
{
	struct ppcg_stmt *stmt = user;

	if (!stmt)
		return;

	isl_id_to_ast_expr_free(stmt->ref2expr);

	free(stmt);
}

/* Derive the output file name from the input file name.
 * 'input' is the entire path of the input file. The output
 * is the file name plus the additional extension.
 *
 * We will basically replace everything after the last point
 * with '.ppcg.c'. This means file.c becomes file.ppcg.c
 */
static FILE *get_output_file(const char *input, const char *output)
{
	char name[PATH_MAX];
	const char *ext;
	const char ppcg_marker[] = ".ppcg";
	int len;
	FILE *file;

	len = ppcg_extract_base_name(name, input);

	strcpy(name + len, ppcg_marker);
	ext = strrchr(input, '.');
	strcpy(name + len + sizeof(ppcg_marker) - 1, ext ? ext : ".c");

	if (!output)
		output = name;

	file = fopen(output, "w");
	if (!file) {
		fprintf(stderr, "Unable to open '%s' for writing\n", output);
		return NULL;
	}

	return file;
}

/* Derive the output file name from the input file name.
 * 'input' is the entire path of the input file. The output
 * is the file name plus the additional extension.
 *
 * We will basically replace everything after the last point
 * with '.ppcg.c'. This means file.c becomes file.ppcg.c
 */
static FILE *get_output_file_with_amp(const char *input, const char *output, struct ppcg_options *options)
{
	char name[PATH_MAX];
	const char *ext;
	int len;
	FILE *file;

	len = ppcg_extract_base_name(name, input);

	if (options->automatic_mixed_precision)
	{
		char ppcg_marker[] = ".amp_ppcg";
		strcpy(name + len, ppcg_marker);
		ext = strrchr(input, '.');
		strcpy(name + len + sizeof(ppcg_marker) - 1, ext ? ext : ".c");
	}
	else
	{
		char ppcg_marker[] = ".ppcg";
		strcpy(name + len, ppcg_marker);
		ext = strrchr(input, '.');
		strcpy(name + len + sizeof(ppcg_marker) - 1, ext ? ext : ".c");
	}

	if (!output)
		output = name;

	file = fopen(output, "w");
	if (!file)
	{
		fprintf(stderr, "Unable to open '%s' for writing\n", output);
		return NULL;
	}

	if (options->automatic_mixed_precision)
	{
		// 在文件开始加入AMP相关的头文件
		fprintf(file, "#include <assert.h>\n");
		fprintf(file, "#include <stdio.h>\n");
		fprintf(file, "#include \"amp_utilities.h\"\n\n");
	}

	return file;
}

/* Data used to annotate for nodes in the ast.
 */
struct ast_node_userinfo {
	/* The for node is an openmp parallel for node. */
	int is_openmp;
};

/* Information used while building the ast.
 */
struct ast_build_userinfo {
	/* The current ppcg scop. */
	struct ppcg_scop *scop;

	/* Are we currently in a parallel for loop? */
	int in_parallel_for;

	/* The contraction of the entire schedule tree. */
	isl_union_pw_multi_aff *contraction;
};

/* Check if the current scheduling dimension is parallel.
 *
 * We check for parallelism by verifying that the loop does not carry any
 * dependences.
 *
 * If any expansion nodes are present in the schedule tree,
 * then they are assumed to be situated near the leaves of the schedule tree,
 * underneath any node that may result in a for loop.
 * In particular, these expansions may have been introduced
 * by the call to isl_schedule_expand inside ppcg_compute_grouping_schedule.
 * The dependence relations are formulated in terms of the expanded
 * domains, while, by assumption, the partial schedule returned
 * by isl_ast_build_get_schedule refers to the contracted domains.
 * Plug in the contraction such that the schedule would also
 * refer to the expanded domains.
 * Note that if the schedule tree does not contain any expansions,
 * then the contraction is an identity function.
 *
 * If the live_range_reordering option is set, then this currently
 * includes the order dependences.  In principle, non-zero order dependences
 * could be allowed, but this would require privatization and/or expansion.
 *
 * Parallelism test: if the distance is zero in all outer dimensions, then it
 * has to be zero in the current dimension as well.
 * Implementation: first, translate dependences into time space, then force
 * outer dimensions to be equal.  If the distance is zero in the current
 * dimension, then the loop is parallel.
 * The distance is zero in the current dimension if it is a subset of a map
 * with equal values for the current dimension.
 */
static int ast_schedule_dim_is_parallel(__isl_keep isl_ast_build *build,
	struct ast_build_userinfo *build_info)
{
	struct ppcg_scop *scop = build_info->scop;
	isl_union_map *schedule, *deps;
	isl_map *schedule_deps, *test;
	isl_space *schedule_space;
	unsigned i, dimension, is_parallel;

	schedule = isl_ast_build_get_schedule(build);
	schedule = isl_union_map_preimage_domain_union_pw_multi_aff(schedule,
		isl_union_pw_multi_aff_copy(build_info->contraction));
	schedule_space = isl_ast_build_get_schedule_space(build);

	dimension = isl_space_dim(schedule_space, isl_dim_out) - 1;

	deps = isl_union_map_copy(scop->dep_flow);
	deps = isl_union_map_union(deps, isl_union_map_copy(scop->dep_false));
	if (scop->options->live_range_reordering) {
		isl_union_map *order = isl_union_map_copy(scop->dep_order);
		deps = isl_union_map_union(deps, order);
	}
	deps = isl_union_map_apply_range(deps, isl_union_map_copy(schedule));
	deps = isl_union_map_apply_domain(deps, schedule);

	if (isl_union_map_is_empty(deps)) {
		isl_union_map_free(deps);
		isl_space_free(schedule_space);
		return 1;
	}

	schedule_deps = isl_map_from_union_map(deps);

	for (i = 0; i < dimension; i++)
		schedule_deps = isl_map_equate(schedule_deps, isl_dim_out, i,
					       isl_dim_in, i);

	test = isl_map_universe(isl_map_get_space(schedule_deps));
	test = isl_map_equate(test, isl_dim_out, dimension, isl_dim_in,
			      dimension);
	is_parallel = isl_map_is_subset(schedule_deps, test);

	isl_space_free(schedule_space);
	isl_map_free(test);
	isl_map_free(schedule_deps);

	return is_parallel;
}

/* Mark a for node openmp parallel, if it is the outermost parallel for node.
 */
static void mark_openmp_parallel(__isl_keep isl_ast_build *build,
	struct ast_build_userinfo *build_info,
	struct ast_node_userinfo *node_info)
{
	if (build_info->in_parallel_for)
		return;

	if (ast_schedule_dim_is_parallel(build, build_info)) {
		build_info->in_parallel_for = 1;
		node_info->is_openmp = 1;
	}
}

/* Allocate an ast_node_info structure and initialize it with default values.
 */
static struct ast_node_userinfo *allocate_ast_node_userinfo()
{
	struct ast_node_userinfo *node_info;
	node_info = (struct ast_node_userinfo *)
		malloc(sizeof(struct ast_node_userinfo));
	node_info->is_openmp = 0;
	return node_info;
}

/* Free an ast_node_info structure.
 */
static void free_ast_node_userinfo(void *ptr)
{
	struct ast_node_userinfo *info;
	info = (struct ast_node_userinfo *) ptr;
	free(info);
}

/* This method is executed before the construction of a for node. It creates
 * an isl_id that is used to annotate the subsequently generated ast for nodes.
 *
 * In this function we also run the following analyses:
 *
 * 	- Detection of openmp parallel loops
 */
static __isl_give isl_id *ast_build_before_for(
	__isl_keep isl_ast_build *build, void *user)
{
	isl_id *id;
	struct ast_build_userinfo *build_info;
	struct ast_node_userinfo *node_info;

	build_info = (struct ast_build_userinfo *) user;
	node_info = allocate_ast_node_userinfo();
	id = isl_id_alloc(isl_ast_build_get_ctx(build), "", node_info);
	id = isl_id_set_free_user(id, free_ast_node_userinfo);

	mark_openmp_parallel(build, build_info, node_info);

	return id;
}

/* This method is executed after the construction of a for node.
 *
 * It performs the following actions:
 *
 * 	- Reset the 'in_parallel_for' flag, as soon as we leave a for node,
 * 	  that is marked as openmp parallel.
 *
 */
static __isl_give isl_ast_node *ast_build_after_for(
	__isl_take isl_ast_node *node, __isl_keep isl_ast_build *build,
	void *user)
{
	isl_id *id;
	struct ast_build_userinfo *build_info;
	struct ast_node_userinfo *info;

	id = isl_ast_node_get_annotation(node);
	info = isl_id_get_user(id);

	if (info && info->is_openmp) {
		build_info = (struct ast_build_userinfo *) user;
		build_info->in_parallel_for = 0;
	}

	isl_id_free(id);

	return node;
}

/* Find the element in scop->stmts that has the given "id".
 */
static struct pet_stmt *find_stmt(struct ppcg_scop *scop, __isl_keep isl_id *id)
{
	int i;

	for (i = 0; i < scop->pet->n_stmt; ++i) {
		struct pet_stmt *stmt = scop->pet->stmts[i];
		isl_id *id_i;

		id_i = isl_set_get_tuple_id(stmt->domain);
		isl_id_free(id_i);

		if (id_i == id)
			return stmt;
	}

	isl_die(isl_id_get_ctx(id), isl_error_internal,
		"statement not found", return NULL);
}

/* Print a user statement in the generated AST.
 * The ppcg_stmt has been attached to the node in at_each_domain.
 */
static __isl_give isl_printer *print_user(__isl_take isl_printer *p,
	__isl_take isl_ast_print_options *print_options,
	__isl_keep isl_ast_node *node, void *user)
{
	struct ppcg_stmt *stmt;
	isl_id *id;

	id = isl_ast_node_get_annotation(node);
	stmt = isl_id_get_user(id);
	isl_id_free(id);

	p = pet_stmt_print_body(stmt->stmt, p, stmt->ref2expr);

	isl_ast_print_options_free(print_options);

	return p;
}

/* Print an access to the element in the private/shared memory copy
 * described by "stmt".  The index of the copy is recorded in
 * stmt->local_index as an access to the array.
 */
static __isl_give isl_printer *stmt_print_local_index(__isl_take isl_printer *p,
													  struct ppcg_kernel_stmt *stmt)
{
	return isl_printer_print_ast_expr(p, stmt->u.c.local_index);
}

/* Print an access to the element in the global memory copy
 * described by "stmt".  The index of the copy is recorded in
 * stmt->index as an access to the array.
 */
static __isl_give isl_printer *stmt_print_global_index(
	__isl_take isl_printer *p, struct ppcg_kernel_stmt *stmt)
{
	struct amp_array_info *array = stmt->u.c.array;
	isl_ast_expr *index;

	if (amp_array_is_scalar(array))
	{
		if (!amp_array_is_read_only_scalar(array))
			p = isl_printer_print_str(p, "*");
		p = isl_printer_print_str(p, array->name);
		return p;
	}

	index = isl_ast_expr_copy(stmt->u.c.index);

	p = isl_printer_print_ast_expr(p, index);
	isl_ast_expr_free(index);

	return p;
}

/* Print a copy statement.
 *
 * A read copy statement is printed as
 *
 *	local = global;
 *
 * while a write copy statement is printed as
 *
 *	global = local;
 */
static __isl_give isl_printer *ppcg_kernel_print_copy(__isl_take isl_printer *p,
													  struct ppcg_kernel_stmt *stmt)
{
	p = isl_printer_start_line(p);
	if (stmt->u.c.read)
	{
		p = stmt_print_local_index(p, stmt);
		p = isl_printer_print_str(p, " = (float)");
		p = stmt_print_global_index(p, stmt);
	}
	else
	{
		p = stmt_print_global_index(p, stmt);
		p = isl_printer_print_str(p, " = (double)");
		p = stmt_print_local_index(p, stmt);
	}
	p = isl_printer_print_str(p, ";");
	p = isl_printer_end_line(p);

	return p;
}

static __isl_give isl_printer *ppcg_kernel_print_domain(__isl_take isl_printer *p,
														struct ppcg_kernel_stmt *stmt)
{
	return pet_stmt_print_body(stmt->u.d.stmt->stmt, p, stmt->u.d.ref2expr);
}

/* This function is called for each user statement in the AST,
 * i.e., for each kernel body statement, copy statement or sync statement.
 */
static __isl_give isl_printer *print_kernel_stmt(__isl_take isl_printer *p,
												 __isl_take isl_ast_print_options *print_options,
												 __isl_keep isl_ast_node *node, void *user)
{
	isl_id *id;
	struct ppcg_kernel_stmt *stmt;

	id = isl_ast_node_get_annotation(node);
	stmt = isl_id_get_user(id);
	isl_id_free(id);

	isl_ast_print_options_free(print_options);

	switch (stmt->type)
	{
	case ppcg_kernel_copy:
		return ppcg_kernel_print_copy(p, stmt);
	case ppcg_kernel_domain:
		return ppcg_kernel_print_domain(p, stmt);
	}

	return p;
}

/* Print a user statement in the generated AST.
 * The ppcg_stmt has been attached to the node in at_each_domain.
 */
static __isl_give isl_printer *print_user_with_amp(__isl_take isl_printer *p,
												   __isl_take isl_ast_print_options *print_options,
												   __isl_keep isl_ast_node *node, void *user)
{
	// #define DEBUG_PRINT_USER_WITH_AMP

	isl_id *id;
	int is_user, is_amp_kernel;
	struct amp_ppcg_kernel *kernel;
	struct ppcg_kernel_stmt *stmt;
	struct amp_prog *prog;

	isl_ast_print_options_free(print_options);

	prog = (struct amp_prog *)user;

	id = isl_ast_node_get_annotation(node);
	is_user = !strcmp(isl_id_get_name(id), "user");
	is_amp_kernel = !strcmp(isl_id_get_name(id), "amp_kernel");

#ifdef DEBUG_PRINT_USER_WITH_AMP
	printf("@DEBUG: \n       in print_user_with_amp, the id is : \n       ");
	isl_id_dump(id);
	printf("\n       the name is %s \n\n", isl_id_get_name(id));
	// printf("\n\n");
#endif // DEBUG_PRINT_USER_WITH_AMP

	// if (is_amp_kernel)
	// {
	// 	kernel = isl_id_get_user(id);
	// 	isl_id_free(id);
	// 	isl_ctx *ctx = isl_ast_node_get_ctx(kernel->tree);
	// 	isl_ast_print_options *print_options;

	// 	// p = print_kernel_vars(p, kernel);
	// 	// p = isl_printer_end_line(p);
	// 	print_options = isl_ast_print_options_alloc(ctx);
	// 	print_options = isl_ast_print_options_set_print_user(print_options, &print_kernel_stmt, NULL);
	// 	p = isl_ast_node_print(kernel->tree, p, print_options);
	// 	return p;
	// }
	// else if (is_user)
	if (is_user)
	{
		stmt = isl_id_get_user(id);
		isl_id_free(id);
		return ppcg_kernel_print_domain(p, stmt);
	}
	else
	{
		stmt = isl_id_get_user(id);
		isl_id_free(id);
	}

#ifdef DEBUG_PRINT_USER_WITH_AMP
	printf("@DEBUG: \n       in print_user_with_amp, the id is : \n       ");
	isl_id_dump(id);
	printf("\n\n");
#endif // DEBUG_PRINT_USER_WITH_AMP
	if (stmt)
	{
		switch (stmt->type)
		{
		case ppcg_kernel_copy:
			return ppcg_kernel_print_copy(p, stmt);
		case ppcg_kernel_domain:
			return ppcg_kernel_print_domain(p, stmt);
		}
	}

	return p;
}

/* Print a for loop node as an openmp parallel loop.
 *
 * To print an openmp parallel loop we print a normal for loop, but add
 * "#pragma openmp parallel for" in front.
 *
 * Variables that are declared within the body of this for loop are
 * automatically openmp 'private'. Iterators declared outside of the
 * for loop are automatically openmp 'shared'. As ppcg declares all iterators
 * at the position where they are assigned, there is no need to explicitly mark
 * variables. Their automatically assigned type is already correct.
 *
 * This function only generates valid OpenMP code, if the ast was generated
 * with the 'atomic-bounds' option enabled.
 *
 */
static __isl_give isl_printer *print_for_with_openmp(
	__isl_keep isl_ast_node *node, __isl_take isl_printer *p,
	__isl_take isl_ast_print_options *print_options)
{
	p = isl_printer_start_line(p);
	p = isl_printer_print_str(p, "#pragma omp parallel for");
	p = isl_printer_end_line(p);

	p = isl_ast_node_for_print(node, p, print_options);

	return p;
}

/* Print a for node.
 *
 * Depending on how the node is annotated, we either print a normal
 * for node or an openmp parallel for node.
 */
static __isl_give isl_printer *print_for(__isl_take isl_printer *p,
	__isl_take isl_ast_print_options *print_options,
	__isl_keep isl_ast_node *node, void *user)
{
	isl_id *id;
	int openmp;

	openmp = 0;
	id = isl_ast_node_get_annotation(node);

	if (id) {
		struct ast_node_userinfo *info;

		info = (struct ast_node_userinfo *) isl_id_get_user(id);
		if (info && info->is_openmp)
			openmp = 1;
	}

	if (openmp)
		p = print_for_with_openmp(node, p, print_options);
	else
		p = isl_ast_node_for_print(node, p, print_options);

	isl_id_free(id);

	return p;
}

/* Index transformation callback for pet_stmt_build_ast_exprs.
 *
 * "index" expresses the array indices in terms of statement iterators
 * "iterator_map" expresses the statement iterators in terms of
 * AST loop iterators.
 *
 * The result expresses the array indices in terms of
 * AST loop iterators.
 */
static __isl_give isl_multi_pw_aff *pullback_index(
	__isl_take isl_multi_pw_aff *index, __isl_keep isl_id *id, void *user)
{
	isl_pw_multi_aff *iterator_map = user;

	iterator_map = isl_pw_multi_aff_copy(iterator_map);
	return isl_multi_pw_aff_pullback_pw_multi_aff(index, iterator_map);
}

/* Transform the accesses in the statement associated to the domain
 * called by "node" to refer to the AST loop iterators, construct
 * corresponding AST expressions using "build",
 * collect them in a ppcg_stmt and annotate the node with the ppcg_stmt.
 */
static __isl_give isl_ast_node *at_each_domain(__isl_take isl_ast_node *node,
	__isl_keep isl_ast_build *build, void *user)
{
	struct ppcg_scop *scop = user;
	isl_ast_expr *expr, *arg;
	isl_ctx *ctx;
	isl_id *id;
	isl_map *map;
	isl_pw_multi_aff *iterator_map;
	struct ppcg_stmt *stmt;

	ctx = isl_ast_node_get_ctx(node);
	stmt = isl_calloc_type(ctx, struct ppcg_stmt);
	if (!stmt)
		goto error;

	expr = isl_ast_node_user_get_expr(node);
	arg = isl_ast_expr_get_op_arg(expr, 0);
	isl_ast_expr_free(expr);
	id = isl_ast_expr_get_id(arg);
	isl_ast_expr_free(arg);
	stmt->stmt = find_stmt(scop, id);
	isl_id_free(id);
	if (!stmt->stmt)
		goto error;

	map = isl_map_from_union_map(isl_ast_build_get_schedule(build));
	map = isl_map_reverse(map);
	iterator_map = isl_pw_multi_aff_from_map(map);
	stmt->ref2expr = pet_stmt_build_ast_exprs(stmt->stmt, build, &pullback_index, iterator_map, NULL, NULL);
	isl_pw_multi_aff_free(iterator_map);

	id = isl_id_alloc(isl_ast_node_get_ctx(node), NULL, stmt);
	id = isl_id_set_free_user(id, &ppcg_stmt_free);
	return isl_ast_node_set_annotation(node, id);
error:
	ppcg_stmt_free(stmt);
	return isl_ast_node_free(node);
}

/* Find the element in gen->stmt that has the given "id".
 * Return NULL if no such gpu_stmt can be found.
 */
static struct amp_stmt *find_amp_stmt(struct amp_prog *prog, __isl_keep isl_id *id)
{
	int i;

	for (i = 0; i < prog->n_stmts; ++i)
	{
		if (id == prog->stmts[i].id)
			break;
	}

	return i < prog->n_stmts ? &prog->stmts[i] : NULL;
}

/* Given a mapping "iterator_map" from the AST schedule to a domain,
 * return the corresponding mapping from the AST schedule
 * to the outer kernel->copy_schedule_dim dimensions of
 * the schedule computed by PPCG for this kernel.
 *
 * Note that kernel->copy_schedule_dim is at least as large as
 * the largest depth of any array reference group associated to the kernel.
 * This is needed as the returned schedule is used to extract a mapping
 * to the outer tile->depth dimensions in transform_index.
 */
static __isl_give isl_pw_multi_aff *compute_sched_to_copy(
	struct amp_ppcg_kernel *kernel, __isl_take isl_pw_multi_aff *iterator_map)
{
	isl_union_pw_multi_aff *upma;
	isl_pw_multi_aff *pma;
	isl_space *space;

	space = isl_space_range(isl_pw_multi_aff_get_space(iterator_map));
	space = isl_space_from_domain(space);
	space = isl_space_add_dims(space, isl_dim_out,
							   kernel->copy_schedule_dim);

	upma = isl_union_pw_multi_aff_copy(kernel->copy_schedule);
	pma = isl_union_pw_multi_aff_extract_pw_multi_aff(upma, space);
	isl_union_pw_multi_aff_free(upma);

	return isl_pw_multi_aff_pullback_pw_multi_aff(pma, iterator_map);
}

/* Return the gpu_stmt_access in the list "accesses" that corresponds
 * to "ref_id".
 */
static struct amp_stmt_access *find_access(struct amp_stmt_access *accesses,
										   __isl_keep isl_id *ref_id)
{
	struct amp_stmt_access *access;

	for (access = accesses; access; access = access->next)
		if (access->ref_id == ref_id)
			return access;

	return NULL;
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

/* Return the index of the array called "name" in the list of arrays.
 */
static int find_array_index(struct amp_ppcg_kernel *kernel, const char *name)
{
	int i;

	for (i = 0; i < kernel->n_array; ++i)
		if (!strcmp(name, kernel->array[i].array->name))
			return i;

	return -1;
}

/* Return a pointer to the gpu_array_ref_group in "local"
 * that contains the reference "access".
 * Return NULL if no such group can be found.
 */
static struct amp_array_ref_group *find_ref_group(
	struct amp_local_array_info *local, struct amp_stmt_access *access)
{
	int i, j;

	for (i = 0; i < local->n_group; ++i)
	{
		struct amp_array_ref_group *group = local->groups[i];

		for (j = 0; j < group->n_ref; ++j)
			if (group->refs[j] == access)
				return group;
	}

	return NULL;
}

/* Given an index expression "index" of the form
 *
 *	L -> F(A),
 *
 * with F(A) either A or some subfield of A and L the AST loop iterators,
 * and a tiling "tiling" of the form
 *
 *	[L -> A] -> T
 *
 * apply the tiling to the outer array in the index expression to obtain
 *
 *	L -> T(A)
 *
 * If F(A) is some subfield of A, then separate the member access
 * into the base index expression and the field index expression,
 * apply the tiling to the base index expression and combine the result
 * with the field index expression.
 *
 * If F(A) is A, then modify index to keep track of the iterators
 *
 *	L -> [L -> A]
 *
 * and combine the result with the tiling to obtain a tiled index expression
 * in terms of the AST loop iterators
 *
 *	L -> T
 */
static __isl_give isl_multi_pw_aff *tile_outer(
	__isl_take isl_multi_pw_aff *index, __isl_take isl_multi_pw_aff *tiling)
{
	isl_bool is_wrapping;
	isl_space *space;
	isl_multi_pw_aff *mpa;

	is_wrapping = isl_multi_pw_aff_range_is_wrapping(index);
	if (is_wrapping < 0)
		goto error;
	if (is_wrapping)
	{
		isl_multi_pw_aff *field;

		field = isl_multi_pw_aff_copy(index);
		field = isl_multi_pw_aff_range_factor_range(field);
		index = isl_multi_pw_aff_range_factor_domain(index);
		index = tile_outer(index, tiling);
		return isl_multi_pw_aff_range_product(index, field);
	}

	space = isl_space_domain(isl_multi_pw_aff_get_space(index));
	space = isl_space_map_from_set(space);
	mpa = isl_multi_pw_aff_identity(space);
	index = isl_multi_pw_aff_range_product(mpa, index);
	index = isl_multi_pw_aff_pullback_multi_pw_aff(tiling, index);

	return index;
error:
	isl_multi_pw_aff_free(index);
	isl_multi_pw_aff_free(tiling);
	return NULL;
}

/* Index transformation callback for pet_stmt_build_ast_exprs.
 *
 * "index" expresses the array indices in terms of statement iterators
 *
 * We first reformulate "index" in terms of the AST loop iterators.
 * Then we check if we are accessing the global array or
 * a shared/private copy.  In particular, if we are not inside a kernel
 * then we must be accessing a global array.
 * In the former case, we simply return
 * the updated index.  If "index" is an affine expression rather
 * than an array access, then we also return the updated index here.
 *
 * If no reference groups have been computed for the array,
 * then we can only be accessing the global array.
 *
 * Otherwise, we apply the tiling to the index.
 * This tiling is of the form
 *
 *	[D -> A] -> T
 *
 * where D corresponds to the outer tile->depth dimensions of
 * the kernel schedule.
 * The index is of the form
 *
 *	L -> A
 *
 * We update the tiling to refer to the AST loop iterators
 *
 *	[L -> A] -> T
 *
 * and combine it with the index to obtain a tiled index expression in terms
 * of the AST loop iterators
 *
 *	L -> T
 *
 * Note that while the tiling applies directly to an outer array.
 * the index may refer to some subfield of this outer array.
 * In such cases, the result will refer to the same subfield of the tile.
 * That is, an index expression of the form  L -> F(A) will be transformed
 * into an index expression of the form L -> F(T).
 */
static __isl_give isl_multi_pw_aff *transform_index(
	__isl_take isl_multi_pw_aff *index, __isl_keep isl_id *ref_id,
	void *user)
{
	struct ppcg_transform_data *data = user;
	struct amp_stmt_access *access;
	struct amp_array_ref_group *group;
	struct amp_array_tile *tile;
	isl_pw_multi_aff *iterator_map;
	int i;
	int dim;
	const char *name;
	isl_space *space;
	isl_multi_pw_aff *tiling;
	isl_pw_multi_aff *pma;
	isl_pw_multi_aff *sched2depth;

	data->array = NULL;

	iterator_map = isl_pw_multi_aff_copy(data->iterator_map);
	index = isl_multi_pw_aff_pullback_pw_multi_aff(index, iterator_map);

	if (!data->kernel)
		return index;

	access = find_access(data->accesses, ref_id);
	if (!access)
		return index;
	if (!isl_map_has_tuple_name(access->access, isl_dim_out))
		return index;

	name = get_outer_array_name(access->access);
	if (!name)
		return isl_multi_pw_aff_free(index);
	i = find_array_index(data->kernel, name);
	if (i < 0)
		isl_die(isl_multi_pw_aff_get_ctx(index), isl_error_internal,
				"cannot find array",
				return isl_multi_pw_aff_free(index));
	data->local_array = &data->kernel->array[i];
	data->array = data->local_array->array;

	group = find_ref_group(data->local_array, access);
	if (!group)
	{
		data->global = 1;
		return index;
	}

	tile = amp_array_ref_group_tile(group);
	data->global = !tile;
	if (!tile)
		return index;

	space = isl_space_domain(isl_multi_aff_get_space(tile->tiling));
	space = isl_space_range(isl_space_unwrap(space));
	space = isl_space_map_from_set(space);
	pma = isl_pw_multi_aff_identity(space);
	sched2depth = isl_pw_multi_aff_copy(data->sched2copy);
	dim = isl_pw_multi_aff_dim(sched2depth, isl_dim_out);
	sched2depth = isl_pw_multi_aff_drop_dims(sched2depth, isl_dim_out, tile->depth, dim - tile->depth);
	pma = isl_pw_multi_aff_product(sched2depth, pma);
	tiling = isl_multi_pw_aff_from_multi_aff(isl_multi_aff_copy(tile->tiling));
	tiling = isl_multi_pw_aff_pullback_pw_multi_aff(tiling, pma);

	index = tile_outer(index, tiling);

	return index;
}
/* Dereference "expr" by adding an index [0].
 * The original "expr" is assumed not to have any indices.
 *
 * If "expr" is a member access, then the dereferencing needs
 * to be applied to the structure argument of this member access.
 */
static __isl_give isl_ast_expr *dereference(__isl_take isl_ast_expr *expr)
{
	isl_ctx *ctx;
	isl_ast_expr *arg0, *res;
	isl_ast_expr_list *list;

	arg0 = isl_ast_expr_get_op_arg(expr, 0);
	if (!arg0)
		return isl_ast_expr_free(expr);
	if (isl_ast_expr_get_type(arg0) == isl_ast_expr_op &&
		isl_ast_expr_get_op_type(arg0) == isl_ast_op_member)
	{
		isl_ast_expr *arg;

		arg = isl_ast_expr_get_op_arg(arg0, 0);
		arg = dereference(arg);
		arg0 = isl_ast_expr_set_op_arg(arg0, 0, arg);
		expr = isl_ast_expr_set_op_arg(expr, 0, arg0);

		return expr;
	}
	isl_ast_expr_free(arg0);

	ctx = isl_ast_expr_get_ctx(expr);
	res = isl_ast_expr_from_val(isl_val_zero(ctx));
	list = isl_ast_expr_list_from_ast_expr(res);
	res = isl_ast_expr_get_op_arg(expr, 0);
	res = isl_ast_expr_access(res, list);
	isl_ast_expr_free(expr);

	return res;
}

/* Linearize the index expression "expr" based on the array bounds
 * of "array".
 *
 * That is, transform expression
 *
 *	A[i_0][i_1]...[i_n]
 *
 * to
 *
 *	A[(..((i_0 * b_1 + i_1) ... ) * b_n + i_n]
 *
 * where b_0, b_1, ..., b_n are the bounds on the array.
 *
 * If the base of "expr" is a member access, then the linearization needs
 * to be applied to the structure argument of this member access.
 *
 * In the base case, if "expr" has no arguments (other than the name of
 * the array), then we are passing an entire array to a function.
 * In this case, there is nothing to linearize.
 * Note that at this point an expression with no arguments can
 * only be an entire array because the scalar case and
 * the case of single struct are handled by the caller.
 *
 * If the number of specified index expressions in "expr"
 * is smaller than the dimension of the accessed array,
 * then the missing i_j also do not appear in the linearized expression.
 * Furthermore, since such an expression does not refer to a single
 * element while the default linearized expression would refer to
 * a single element, we return the expression
 *
 *	A + (..((i_0 * b_1 + i_1) ... ) * b_l + i_l)
 *
 * instead.  Note that because of the special case handling above,
 * we can assume here that there is at least one index expression.
 */
__isl_give isl_ast_expr *amp_local_array_info_linearize_index(
	struct amp_local_array_info *array, __isl_take isl_ast_expr *expr)
{
	int i, n;
	isl_ast_expr *arg0;
	isl_ast_expr *res;
	isl_ast_expr_list *list;

	arg0 = isl_ast_expr_get_op_arg(expr, 0);
	if (isl_ast_expr_get_type(arg0) == isl_ast_expr_op &&
		isl_ast_expr_get_op_type(arg0) == isl_ast_op_member)
	{
		isl_ast_expr *arg;

		arg = isl_ast_expr_get_op_arg(arg0, 0);
		arg = amp_local_array_info_linearize_index(array, arg);
		arg0 = isl_ast_expr_set_op_arg(arg0, 0, arg);
		expr = isl_ast_expr_set_op_arg(expr, 0, arg0);

		return expr;
	}
	isl_ast_expr_free(arg0);

	if (isl_ast_expr_get_op_n_arg(expr) == 1)
		return expr;

	n = isl_ast_expr_get_op_n_arg(expr);
	res = isl_ast_expr_get_op_arg(expr, 1);
	for (i = 1; i < array->n_index; ++i)
	{
		isl_ast_expr *expr_i;

		expr_i = isl_ast_expr_get_op_arg(array->bound_expr, 1 + i);
		res = isl_ast_expr_mul(res, expr_i);

		if (i + 1 >= n)
			continue;
		expr_i = isl_ast_expr_get_op_arg(expr, i + 1);
		res = isl_ast_expr_add(res, expr_i);
	}

	if (1 + array->n_index > n)
	{
		res = isl_ast_expr_add(isl_ast_expr_get_op_arg(expr, 0), res);
	}
	else
	{
		list = isl_ast_expr_list_from_ast_expr(res);
		res = isl_ast_expr_get_op_arg(expr, 0);
		res = isl_ast_expr_access(res, list);
	}

	isl_ast_expr_free(expr);

	return res;
}

/* AST expression transformation callback for pet_stmt_build_ast_exprs.
 *
 * If the AST expression refers to an array that is not accessed
 * at all, then this means the value of the expression is not used,
 * so we might as well print zero (NULL pointer) instead.
 *
 * If the AST expression refers to a global scalar that is not
 * a read-only scalar, then its address was passed to the kernel and
 * we need to dereference it.
 *
 * If the AST expression refers to an access to a global array,
 * then we linearize the access exploiting the bounds in data->local_array.
 */
static __isl_give isl_ast_expr *transform_expr(__isl_take isl_ast_expr *expr,
											   __isl_keep isl_id *id, void *user)
{
	struct ppcg_transform_data *data = user;

	if (!data->array)
		return expr;
	if (!data->array->accessed)
	{
		isl_ctx *ctx;

		ctx = isl_ast_expr_get_ctx(expr);
		isl_ast_expr_free(expr);
		return isl_ast_expr_from_val(isl_val_zero(ctx));
	}
	if (amp_array_is_read_only_scalar(data->array))
		return expr;
	if (!data->global)
		return expr;
	if (data->array->n_index == 0)
		return dereference(expr);
	if (!data->array->linearize)
		return expr;

	return amp_local_array_info_linearize_index(data->local_array, expr);
}

static void ppcg_kernel_stmt_free(void *user)
{
	struct ppcg_kernel_stmt *stmt = user;

	if (!stmt)
		return;

	switch (stmt->type)
	{
	case ppcg_kernel_copy:
		isl_ast_expr_free(stmt->u.c.index);
		isl_ast_expr_free(stmt->u.c.local_index);
		break;
	case ppcg_kernel_domain:
		isl_id_to_ast_expr_free(stmt->u.d.ref2expr);
		break;
	}

	free(stmt);
}

/* This function is called for each instance of a user statement
 * in the kernel "kernel", identified by "gpu_stmt".
 * "kernel" may be NULL if we are not inside a kernel.
 *
 * We attach a struct ppcg_kernel_stmt to the "node", containing
 * a computed AST expression for each access, through an annotation
 * with name "user".
 * These AST expressions are computed from iterator_map,
 * which expresses the domain
 * elements in terms of the generated loops, and sched2copy,
 * which expresses the outer copy_schedule_dim dimensions of
 * the kernel schedule computed by PPCG in terms of the generated loops.
 */
static __isl_give isl_ast_node *create_domain_leaf(
	struct amp_ppcg_kernel *kernel, __isl_take isl_ast_node *node,
	__isl_keep isl_ast_build *build, struct amp_stmt *amp_stmt)
{
	struct ppcg_transform_data data;
	struct ppcg_kernel_stmt *stmt;
	isl_ctx *ctx;
	isl_id *id;
	isl_pw_multi_aff *sched2copy;
	isl_map *map;
	isl_pw_multi_aff *iterator_map;
	isl_union_map *schedule;

	if (!node)
		return NULL;
	ctx = isl_ast_node_get_ctx(node);

	stmt = isl_calloc_type(ctx, struct ppcg_kernel_stmt);
	if (!stmt)
		return isl_ast_node_free(node);

	schedule = isl_ast_build_get_schedule(build);
	map = isl_map_reverse(isl_map_from_union_map(schedule));
	iterator_map = isl_pw_multi_aff_from_map(map);
	if (kernel)
		sched2copy = compute_sched_to_copy(kernel, isl_pw_multi_aff_copy(iterator_map));
	else
		sched2copy = NULL;

	stmt->type = ppcg_kernel_domain;
	stmt->u.d.stmt = amp_stmt;

	data.kernel = kernel;
	data.accesses = stmt->u.d.stmt->accesses;
	data.iterator_map = iterator_map;
	data.sched2copy = sched2copy;
	stmt->u.d.ref2expr = pet_stmt_build_ast_exprs(stmt->u.d.stmt->stmt,
												  build, &transform_index, &data,
												  &transform_expr, &data);

	isl_pw_multi_aff_free(iterator_map);
	isl_pw_multi_aff_free(sched2copy);

	id = isl_id_alloc(ctx, "user", stmt);
	id = isl_id_set_free_user(id, &ppcg_kernel_stmt_free);
	if (!id)
		ppcg_kernel_stmt_free(stmt);
	return isl_ast_node_set_annotation(node, id);
}

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
/* This function is called for each statement node in the AST
 * for copying to or from shared/private memory.
 * Attach a pointer to a ppcg_kernel_stmt representing the copy
 * statement to the node.
 * The statement name is "read" or "write", depending on whether we are
 * reading from global memory or writing to global memory.
 *
 * The schedule is of the form
 *
 *	type[D -> A] -> L
 *
 * where D corresponds to the outer tile->depth dimensions of
 * the kernel schedule, A to the global array and L to the outer
 * generated AST schedule.
 * We compute the inverse and strip off the type, resulting in
 *
 *	L -> [D -> A]
 *
 * We combine this mapping with on the one hand the projection
 *
 *	[D -> A] -> A
 *
 * and on the other hand the group tiling
 *
 *	[D -> A] -> T
 *
 * resulting in
 *
 *	L -> A		and 	L -> T
 *
 * and store the corresponding expressions in stmt->index and stmt->local_index,
 * where stmt points to the ppcg_kernel_stmt that is attached to the node.
 * stmt->index is linearized if the global memory array is linearized.
 */
static __isl_give isl_ast_node *create_access_leaf(struct amp_ppcg_kernel *kernel,
												   struct amp_array_ref_group *group, __isl_take isl_ast_node *node,
												   __isl_keep isl_ast_build *build)
{
	// #define DEBUG_CREATE_ACCESS_LEAF

	struct ppcg_kernel_stmt *stmt;
	struct amp_array_tile *tile;
	isl_id *id;
	isl_ast_expr *expr;
	isl_space *space;
	isl_map *access;
	isl_pw_multi_aff *pma, *pma2;
	const char *type;
	isl_ctx *ctx;

	if (kernel == NULL)
		printf("\n\033[31m@ERROR:\n       the amp_ppcg_kernel in create_access_leaf,is NULL!!!  \n\n\033[0m");
	if (!node)
		return NULL;
	ctx = isl_ast_node_get_ctx(node);

#ifdef DEBUG_CREATE_ACCESS_LEAF
	printf("@DEBUG: \n       in create_access_leaf, the node is: \n");
	isl_ast_node_dump(node);
	printf("\n       the build is: \n");
	isl_ast_build_dump(build);
	printf("\n\n");
#endif // DEBUG_CREATE_ACCESS_LEAF

	stmt = isl_calloc_type(ctx, struct ppcg_kernel_stmt);
	if (!stmt)
		return isl_ast_node_free(node);

	access = isl_map_from_union_map(isl_ast_build_get_schedule(build));
	type = isl_map_get_tuple_name(access, isl_dim_in);
	stmt->u.c.read = type && !strcmp(type, "read");
	access = isl_map_reverse(access);
	pma = isl_pw_multi_aff_from_map(access);
	pma = isl_pw_multi_aff_reset_tuple_id(pma, isl_dim_out);

	space = isl_space_range(isl_pw_multi_aff_get_space(pma));
	space = isl_space_unwrap(space);
	pma2 = isl_pw_multi_aff_range_map(space);
	pma2 = isl_pw_multi_aff_pullback_pw_multi_aff(pma2, isl_pw_multi_aff_copy(pma));
	expr = isl_ast_build_access_from_pw_multi_aff(build, pma2);
	if (group->array->linearize)
		expr = amp_local_array_info_linearize_index(group->local_array, expr);
	stmt->u.c.index = expr;

	tile = amp_array_ref_group_tile(group);
	if (tile->tiling == NULL)
	{

		printf("\n\033[31m@ERROR:\n       the tile->tiling in create_access_leaf,is NULL!!!  \n\n\033[0m");
		isl_multi_aff *ma = isl_multi_aff_copy(create_from_access(kernel->ctx, group, 1));
		pma2 = isl_pw_multi_aff_from_multi_aff(isl_multi_aff_copy(ma));
	}
	else
	{
		pma2 = isl_pw_multi_aff_from_multi_aff(isl_multi_aff_copy(tile->tiling));
	}
	// from_access = create_from_access(kernel->ctx, group, read);
	// pma2 = isl_pw_multi_aff_from_multi_aff(isl_multi_aff_copy(tile->tiling));
	pma2 = isl_pw_multi_aff_pullback_pw_multi_aff(pma2, pma);
	expr = isl_ast_build_access_from_pw_multi_aff(build, pma2);
	stmt->u.c.local_index = expr;

	stmt->u.c.array = group->array;
	stmt->u.c.local_array = group->local_array;
	stmt->type = ppcg_kernel_copy;

	id = isl_id_alloc(kernel->ctx, "copy", stmt);
	id = isl_id_set_free_user(id, &ppcg_kernel_stmt_free);
	if (!id)
		ppcg_kernel_stmt_free(stmt);
	return isl_ast_node_set_annotation(node, id);
}

/* Transform the accesses in the statement associated to the domain
 * called by "node" to refer to the AST loop iterators, construct
 * corresponding AST expressions using "build",
 * collect them in a ppcg_stmt and annotate the node with the ppcg_stmt.
 */
static __isl_give isl_ast_node *at_each_domain_with_amp(__isl_take isl_ast_node *node, __isl_keep isl_ast_build *build, void *user)
{
	// #define DEBUG_AT_EACH_DOMAIN_WITH_AMP

	struct ppcg_at_domain_data *data = user;
	amp_prog *prog = data->prog;
	struct ppcg_scop *scop = data->prog->scop;
	isl_ast_expr *expr, *arg;
	isl_ctx *ctx;
	isl_id *id;
	isl_map *map;
	isl_pw_multi_aff *iterator_map;
	// struct ppcg_stmt *stmt;
	struct amp_stmt *amp_stmt;
	const char *name;
	void *p;

#ifdef DEBUG_AT_EACH_DOMAIN_WITH_AMP
	printf("@DEBUG: \n       at start of the at_each_domain_with_amp, the node is: \n");
	isl_ast_node_dump(node);
	printf("\n       the build is :\n");
	isl_ast_build_dump(build);
	printf("\n\n");
#endif // DEBUG_AT_EACH_DOMAIN_WITH_AMP

	ctx = isl_ast_node_get_ctx(node);
	amp_stmt = isl_calloc_type(ctx, struct amp_stmt);
	if (!amp_stmt)
		goto error;

	expr = isl_ast_node_user_get_expr(node);
	arg = isl_ast_expr_get_op_arg(expr, 0);
	isl_ast_expr_free(expr);
	id = isl_ast_expr_get_id(arg);
	name = isl_id_get_name(id);
	p = isl_id_get_user(id);
	isl_ast_expr_free(arg);

#ifdef DEBUG_AT_EACH_DOMAIN_WITH_AMP
	printf("@DEBUG: \n       the id、name、user of arg is: \n");
	isl_id_dump(id);
	printf("\n       the name is %s :\n", name);
	printf("\n       the user(p) is %p :\n", p);
	printf("\n\n");
#endif // DEBUG_AT_EACH_DOMAIN_WITH_AMP

	amp_stmt = find_amp_stmt(data->prog, id);

	if (amp_stmt)
		return create_domain_leaf(data->kernel, node, build, amp_stmt);

	if (!strcmp(name, "read") || !strcmp(name, "write"))
	{
		struct amp_array_ref_group *group = p;
		/** Build AST expressions for the amp array sizes of all arrays in "prog" **/
		node = amp_build_array_bounds(node, prog, build);
		return create_access_leaf(data->kernel, group, node, build);
	}

	printf("\n\033[31m@ERROR:\n       the at_each_domain_with_amp function meets an unexpected errors.  \n\n\033[0m");
	return isl_ast_node_set_annotation(node, id);
	// // 下面是原始的
	// stmt->stmt = find_stmt(scop, id);
	// isl_id_free(id);
	// if (!stmt->stmt)
	// 	goto error;

	// map = isl_map_from_union_map(isl_ast_build_get_schedule(build));
	// map = isl_map_reverse(map);
	// iterator_map = isl_pw_multi_aff_from_map(map);
	// stmt->ref2expr = pet_stmt_build_ast_exprs(stmt->stmt, build, &pullback_index, iterator_map, NULL, NULL);
	// isl_pw_multi_aff_free(iterator_map);

	// id = isl_id_alloc(isl_ast_node_get_ctx(node), NULL, stmt);
	// id = isl_id_set_free_user(id, &ppcg_stmt_free);

	// return isl_ast_node_set_annotation(node, id);
error:
	// ppcg_stmt_free(stmt);
	return isl_ast_node_free(node);
}

/* Build access AST expressions for the localized array sizes using "build".
 * Store the result in local->bound_expr.
 * Only do this for arrays for which localized bounds have been computed.
 */
static isl_stat build_local_array_sizes(struct amp_ppcg_kernel *kernel,
										__isl_keep isl_ast_build *build)
{
	int i;

	for (i = 0; i < kernel->n_array; ++i)
	{
		struct amp_local_array_info *local = &kernel->array[i];
		isl_multi_pw_aff *size;

		if (local->n_group == 0)
			continue;
		size = isl_multi_pw_aff_copy(local->bound);
		local->bound_expr = ppcg_build_size_expr(size, build);
		kernel->size_expr = ppcg_build_size_expr(size, build);
		if (!local->bound_expr)
			return isl_stat_error;
	}

	return isl_stat_ok;
}

/* Build access AST expressions for the effective grid size and
 * the localized array sizes using "build".
 */
static isl_stat build_amp_array_sizes(struct amp_ppcg_kernel *kernel,
									  __isl_keep isl_ast_build *build)
{
	if (build_local_array_sizes(kernel, build) < 0)
		return isl_stat_error;
	return isl_stat_ok;
}

/* This function is called before the AST generator starts traversing
 * the schedule subtree of a node with mark "mark".
 *
 * If the mark is called "kernel", store the kernel pointer in data->kernel
 * for use in at_domain and build AST expressions for the grid size and
 * the localized array sizes.
 */
static isl_stat before_mark_with_amp(__isl_keep isl_id *mark,
									 __isl_keep isl_ast_build *build, void *user)
{
	struct ppcg_at_domain_data *data = user;

	if (!mark)
		return isl_stat_error;
	if (!strcmp(isl_id_get_name(mark), "amp_kernel"))
	{
		data->kernel = isl_id_get_user(mark);
		if (build_amp_array_sizes(data->kernel, build) < 0)
			return isl_stat_error;
	}
	return isl_stat_ok;
}

/* This function is called after the AST generator has finished traversing
 * the schedule subtree of a mark node.  "node" points to the corresponding
 * mark AST node.
 *
 * If the mark is called "kernel", then replace "node" by a user node
 * that "calls" the kernel, representing the launch of the kernel.
 * The original "node" is stored inside the kernel object so that
 * it can be used to print the device code.
 * Note that this assumes that a kernel is only launched once.
 * Also clear data->kernel.
 */
static __isl_give isl_ast_node *after_mark_with_amp(__isl_take isl_ast_node *node,
													__isl_keep isl_ast_build *build, void *user)
{
	isl_ctx *ctx;
	isl_id *id;
	isl_ast_expr *expr;
	isl_ast_expr_list *list;
	struct amp_ppcg_kernel *kernel;
	struct ppcg_at_domain_data *data = user;

	// ctx = isl_ast_node_get_ctx(node);
	id = isl_ast_node_mark_get_id(node);
	if (!id)
		return isl_ast_node_free(node);
	// // if (strcmp(isl_id_get_name(id), "amp_kernel") || !data->kernel)
	// if (strcmp(isl_id_get_name(id), "amp_kernel"))
	// {
	// 	isl_id_free(id);
	// 	return node;
	// }
	// kernel = data->kernel;
	// data->kernel = NULL;
	// kernel->space = isl_ast_build_get_schedule_space(build);
	// kernel->tree = isl_ast_node_mark_get_node(node);
	// isl_ast_node_free(node);

	// expr = isl_ast_expr_from_id(isl_id_copy(id));
	// list = isl_ast_expr_list_alloc(ctx, 0);
	// expr = isl_ast_expr_call(expr, list);
	// node = isl_ast_node_alloc_user(expr);
	// node = isl_ast_node_set_annotation(node, id);

	isl_id_free(id);
	return node;
}

/* Set *depth (initialized to 0 by the caller) to the maximum
 * of the schedule depths of the leaf nodes for which this function is called.
 */
static isl_bool update_depth(__isl_keep isl_schedule_node *node, void *user)
{
	int *depth = user;
	int node_depth;

	if (isl_schedule_node_get_type(node) != isl_schedule_node_leaf)
		return isl_bool_true;
	node_depth = isl_schedule_node_get_schedule_depth(node);
	if (node_depth > *depth)
		*depth = node_depth;

	return isl_bool_false;
}

/* This function is called for each node in a CPU AST.
 * In case of a user node, print the macro definitions required
 * for printing the AST expressions in the annotation, if any.
 * For other nodes, return true such that descendants are also
 * visited.
 *
 * In particular, print the macro definitions needed for the substitutions
 * of the original user statements.
 */
static isl_bool at_node(__isl_keep isl_ast_node *node, void *user)
{
	struct ppcg_stmt *stmt;
	isl_id *id;
	isl_printer **p = user;

	if (isl_ast_node_get_type(node) != isl_ast_node_user)
		return isl_bool_true;

	id = isl_ast_node_get_annotation(node);
	stmt = isl_id_get_user(id);
	isl_id_free(id);

	if (!stmt)
		return isl_bool_error;

	*p = ppcg_print_body_macros(*p, stmt->ref2expr);
	if (!*p)
		return isl_bool_error;

	return isl_bool_false;
}

/* This function is called for each node in a CPU AST.
 * In case of a user node, print the macro definitions required
 * for printing the AST expressions in the annotation, if any.
 * For other nodes, return true such that descendants are also
 * visited.
 *
 * In particular, print the macro definitions needed for the substitutions
 * of the original user statements.
 */
static isl_bool at_node_with_amp(__isl_keep isl_ast_node *node, void *user)
{
	// #define DEBUG_AT_NODE_WITH_AMP

	const char *name;
	int is_kernel;
	struct amp_ppcg_kernel *kernel;
	struct ppcg_kernel_stmt *stmt;
	// struct ppcg_stmt *stmt;
	isl_id *id;
	isl_printer **p = user;

#ifdef DEBUG_AT_NODE_WITH_AMP
	printf("@DEBUG: \n       at_node_with_amp function,the ast node is: \n");
	isl_ast_node_dump(node);
	printf("\n\n");
#endif // DEBUG_AT_NODE_WITH_AMP

	if (isl_ast_node_get_type(node) != isl_ast_node_user)
		return isl_bool_true;

	id = isl_ast_node_get_annotation(node);
	if (!id)
		return isl_bool_false;

	name = isl_id_get_name(id);
	if (!name)
		return isl_bool_error;

	is_kernel = !strcmp(name, "amp_kernel");
	kernel = is_kernel ? isl_id_get_user(id) : NULL;
	stmt = is_kernel ? NULL : isl_id_get_user(id);
	isl_id_free(id);

	if ((is_kernel && !kernel) || (!is_kernel && !stmt))
		return isl_bool_error;

	if (is_kernel)
	{

		// *p = ppcg_ast_expr_print_macros(kernel->tree, *p);
		return isl_bool_true;
	}
	if (stmt->type == ppcg_kernel_copy)
	{
		*p = ppcg_ast_expr_print_macros(stmt->u.c.index, *p);
		*p = ppcg_ast_expr_print_macros(stmt->u.c.local_index, *p);
	}
	else if (stmt->type == ppcg_kernel_domain)
	{
		*p = ppcg_print_body_macros(*p, stmt->u.d.ref2expr);
	}
	if (!*p)
		return isl_bool_error;

	return isl_bool_false;

	// stmt = isl_id_get_user(id);
	// isl_id_free(id);

	// if (!stmt)
	// 	return isl_bool_error;

	// *p = ppcg_print_body_macros(*p, stmt->ref2expr);
	// if (!*p)
	// 	return isl_bool_error;

	// return isl_bool_false;
}

/* Print the required macros for the CPU AST "node" to "p",
 * including those needed for the user statements inside the AST.
 */
static __isl_give isl_printer *cpu_print_macros_with_amp(__isl_take isl_printer *p,
														 __isl_keep isl_ast_node *node)
{
	if (isl_ast_node_foreach_descendant_top_down(node, &at_node_with_amp, &p) < 0)
		return isl_printer_free(p);
	p = ppcg_print_macros(p, node);
	return p;
}

/* Print the required macros for the CPU AST "node" to "p",
 * including those needed for the user statements inside the AST.
 */
static __isl_give isl_printer *cpu_print_macros(__isl_take isl_printer *p,
												__isl_keep isl_ast_node *node)
{
	if (isl_ast_node_foreach_descendant_top_down(node, &at_node, &p) < 0)
		return isl_printer_free(p);
	p = ppcg_print_macros(p, node);
	return p;
}
/* Initialize the fields of "build_info".
 *
 * Initially, the AST generation is not inside any parallel for loop.
 *
 * The contraction of the entire schedule tree is extracted
 * right underneath the root node.
 */
static isl_stat init_build_info(struct ast_build_userinfo *build_info,
								struct ppcg_scop *scop, __isl_keep isl_schedule *schedule)
{
	isl_schedule_node *node = isl_schedule_get_root(schedule);
	node = isl_schedule_node_child(node, 0);

	build_info->scop = scop;
	build_info->in_parallel_for = 0;
	build_info->contraction =
		isl_schedule_node_get_subtree_contraction(node);

	isl_schedule_node_free(node);

	return isl_stat_non_null(build_info->contraction);
}

/* Clear all memory allocated by "build_info".
 */
static void clear_build_info(struct ast_build_userinfo *build_info)
{
	isl_union_pw_multi_aff_free(build_info->contraction);
}

/* Code generate the scop 'scop' using "schedule"
 * and print the corresponding C code to 'p'.
 */
static __isl_give isl_printer *print_scop(struct ppcg_scop *scop,
										  __isl_take isl_schedule *schedule, __isl_take isl_printer *p,
										  struct ppcg_options *options)
{
	isl_ctx *ctx = isl_printer_get_ctx(p);
	isl_ast_build *build;
	isl_ast_print_options *print_options;
	isl_ast_node *tree;
	isl_id_list *iterators;
	struct ast_build_userinfo build_info;
	int depth;

	depth = 0;
	if (isl_schedule_foreach_schedule_node_top_down(schedule, &update_depth,
													&depth) < 0)
		goto error;

	build = isl_ast_build_alloc(ctx);
	iterators = ppcg_scop_generate_names(scop, depth, "c");
	build = isl_ast_build_set_iterators(build, iterators);
	build = isl_ast_build_set_at_each_domain(build, &at_each_domain, scop);

	if (options->openmp)
	{
		if (init_build_info(&build_info, scop, schedule) < 0)
			build = isl_ast_build_free(build);

		build = isl_ast_build_set_before_each_for(build,
												  &ast_build_before_for,
												  &build_info);
		build = isl_ast_build_set_after_each_for(build,
												 &ast_build_after_for,
												 &build_info);
	}

	tree = isl_ast_build_node_from_schedule(build, schedule);
	isl_ast_build_free(build);

	if (options->openmp)
		clear_build_info(&build_info);

	print_options = isl_ast_print_options_alloc(ctx);
	print_options = isl_ast_print_options_set_print_user(print_options,
														 &print_user, NULL);

	print_options = isl_ast_print_options_set_print_for(print_options,
														&print_for, NULL);

	p = cpu_print_macros(p, tree);

	p = isl_ast_node_print(tree, p, print_options);

	isl_ast_node_free(tree);

	return p;
error:
	isl_schedule_free(schedule);
	isl_printer_free(p);
	return NULL;
}

/* Code generate the scop 'scop' using "schedule"
 * and print the corresponding C code to 'p'.
 */
/** 待完善 **/
static __isl_give isl_printer *print_scop_with_amp(__isl_take isl_schedule *schedule, __isl_take isl_printer *p, struct ppcg_options *options, amp_prog *prog)
{
	// #define DEBUG_PRINT_SCOP_WITH_AMP

	struct ppcg_at_domain_data data;
	struct ppcg_scop *scop = prog->scop;
	isl_ctx *ctx = isl_printer_get_ctx(p);
	isl_ast_build *build;
	isl_ast_print_options *print_options;
	isl_ast_node *tree;
	isl_id_list *iterators;
	struct ast_build_userinfo build_info;
	int depth;

#ifdef DEBUG_PRINT_SCOP_WITH_AMP
	printf("\n@DEBUG: \n       at start of cpu.c-print scop,the schedule is: \n");
	isl_schedule_dump(schedule);
	printf("\n\n");
#endif // DEBUG_PRINT_SCOP_WITH_AMP

	data.prog = prog;
	data.kernel = NULL;

	depth = 0;
	if (isl_schedule_foreach_schedule_node_top_down(schedule, &update_depth, &depth) < 0)
		goto error;

	build = isl_ast_build_alloc(ctx);
	iterators = ppcg_scop_generate_names(scop, depth, "c");
	build = isl_ast_build_set_iterators(build, iterators);

	if (options->automatic_mixed_precision)
	{
		/**-- 修改了哈 --**/
		build = isl_ast_build_set_at_each_domain(build, &at_each_domain_with_amp, &data);
		/**-- 仿照gpu，新增的两行代码 */
		build = isl_ast_build_set_before_each_mark(build, &before_mark_with_amp, &data);
		build = isl_ast_build_set_after_each_mark(build, &after_mark_with_amp, &data);
	}
	else
	{
		build = isl_ast_build_set_at_each_domain(build, &at_each_domain, scop);
	}

	if (options->openmp)
	{
		if (init_build_info(&build_info, scop, schedule) < 0)
			build = isl_ast_build_free(build);

		build = isl_ast_build_set_before_each_for(build, &ast_build_before_for, &build_info);
		build = isl_ast_build_set_after_each_for(build, &ast_build_after_for, &build_info);
	}

	tree = isl_ast_build_node_from_schedule(build, schedule);
	isl_ast_build_free(build);

	if (options->openmp)
		clear_build_info(&build_info);

	print_options = isl_ast_print_options_alloc(ctx);
	print_options = isl_ast_print_options_set_print_user(print_options, &print_user_with_amp, prog);
	print_options = isl_ast_print_options_set_print_for(print_options, &print_for, NULL);

	p = cpu_print_macros_with_amp(p, tree);
	if (options->automatic_mixed_precision)
	{
		// 加入混合精度的宏定义
		p = amp_print_macros(p);
		// 打印低精度数组
		p = declare_amp_lower_precision_arrays(p, prog);
		p = allocate_amp_lower_precision_arrays(p, prog);
	}

	p = isl_ast_node_print(tree, p, print_options);

	isl_ast_node_free(tree);

	return p;
error:
	isl_schedule_free(schedule);
	isl_printer_free(p);
	return NULL;
}

/* Tile the band node "node" with tile sizes "sizes" and
 * mark all members of the resulting tile node as "atomic".
 */
static __isl_give isl_schedule_node *tile(__isl_take isl_schedule_node *node,
	__isl_take isl_multi_val *sizes)
{
	node = ppcg_tile(node, sizes);
	node = ppcg_set_schedule_node_type(node, isl_ast_loop_atomic);

	return node;
}

/* Tile "node", if it is a band node with at least 2 members.
 * The tile sizes are set from the "tile_size" option.
 */
static __isl_give isl_schedule_node *tile_band(
	__isl_take isl_schedule_node *node, void *user)
{
	struct ppcg_scop *scop = user;
	int n;
	isl_space *space;
	isl_multi_val *sizes;

	if (isl_schedule_node_get_type(node) != isl_schedule_node_band)
		return node;

	n = isl_schedule_node_band_n_member(node);
	if (n <= 1)
		return node;

	space = isl_schedule_node_band_get_space(node);
	sizes = ppcg_multi_val_from_int(space, scop->options->tile_size);

	return tile(node, sizes);
}

/* Construct schedule constraints from the dependences in ps
 * for the purpose of computing a schedule for a CPU.
 *
 * The proximity constraints are set to the flow dependences.
 *
 * If live-range reordering is allowed then the conditional validity
 * constraints are set to the order dependences with the flow dependences
 * as condition.  That is, a live-range (flow dependence) will be either
 * local to an iteration of a band or all adjacent order dependences
 * will be respected by the band.
 * The validity constraints are set to the union of the flow dependences
 * and the forced dependences, while the coincidence constraints
 * are set to the union of the flow dependences, the forced dependences and
 * the order dependences.
 *
 * If live-range reordering is not allowed, then both the validity
 * and the coincidence constraints are set to the union of the flow
 * dependences and the false dependences.
 *
 * Note that the coincidence constraints are only set when the "openmp"
 * options is set.  Even though the way openmp pragmas are introduced
 * does not rely on the coincident property of the schedule band members,
 * the coincidence constraints do affect the way the schedule is constructed,
 * such that more schedule dimensions should be detected as parallel
 * by ast_schedule_dim_is_parallel.
 * Since the order dependences are also taken into account by
 * ast_schedule_dim_is_parallel, they are also added to
 * the coincidence constraints.  If the openmp handling learns
 * how to privatize some memory, then the corresponding order
 * dependences can be removed from the coincidence constraints.
 */
static __isl_give isl_schedule_constraints *construct_cpu_schedule_constraints(
	struct ppcg_scop *ps)
{
	isl_schedule_constraints *sc;
	isl_union_map *validity, *coincidence;

	sc = isl_schedule_constraints_on_domain(isl_union_set_copy(ps->domain));
	if (ps->options->live_range_reordering) {
		sc = isl_schedule_constraints_set_conditional_validity(sc,
				isl_union_map_copy(ps->tagged_dep_flow),
				isl_union_map_copy(ps->tagged_dep_order));
		validity = isl_union_map_copy(ps->dep_flow);
		validity = isl_union_map_union(validity,
				isl_union_map_copy(ps->dep_forced));
		if (ps->options->openmp) {
			coincidence = isl_union_map_copy(validity);
			coincidence = isl_union_map_union(coincidence,
					isl_union_map_copy(ps->dep_order));
		}
	} else {
		validity = isl_union_map_copy(ps->dep_flow);
		validity = isl_union_map_union(validity,
				isl_union_map_copy(ps->dep_false));
		if (ps->options->openmp)
			coincidence = isl_union_map_copy(validity);
	}
	if (ps->options->openmp)
		sc = isl_schedule_constraints_set_coincidence(sc, coincidence);
	sc = isl_schedule_constraints_set_validity(sc, validity);
	sc = isl_schedule_constraints_set_proximity(sc,
					isl_union_map_copy(ps->dep_flow));

	return sc;
}

/* Compute a schedule for the scop "ps".
 *
 * First derive the appropriate schedule constraints from the dependences
 * in "ps" and then compute a schedule from those schedule constraints,
 * possibly grouping statement instances based on the input schedule.
 */
static __isl_give isl_schedule *compute_cpu_schedule(struct ppcg_scop *ps)
{
	isl_schedule_constraints *sc;
	isl_schedule *schedule;

	if (!ps)
		return NULL;

	sc = construct_cpu_schedule_constraints(ps);

	schedule = ppcg_compute_schedule(sc, ps->schedule, ps->options);

	return schedule;
}

/* Compute a new schedule to the scop "ps" if the reschedule option is set.
 * Otherwise, return a copy of the original schedule.
 */
static __isl_give isl_schedule *optionally_compute_schedule(void *user)
{
	struct ppcg_scop *ps = user;

	if (!ps)
		return NULL;
	if (!ps->options->reschedule)
		return isl_schedule_copy(ps->schedule);
	return compute_cpu_schedule(ps);
}

/* Compute a schedule based on the dependences in "ps" and
 * tile it if requested by the user.
 */
static __isl_give isl_schedule *get_schedule(struct ppcg_scop *ps,
	struct ppcg_options *options)
{
	isl_ctx *ctx;
	isl_schedule *schedule;

	if (!ps)
		return NULL;

	ctx = isl_union_set_get_ctx(ps->domain);
	schedule = ppcg_get_schedule(ctx, options,
				    &optionally_compute_schedule, ps);
	if (ps->options->tile)
		schedule = isl_schedule_map_schedule_node_bottom_up(schedule,
							&tile_band, ps);

	return schedule;
}

/* Generate CPU code for the scop "ps" using "schedule" and
 * print the corresponding C code to "p", including variable declarations.
 */
static __isl_give isl_printer *print_cpu_with_schedule(
	__isl_take isl_printer *p, struct ppcg_scop *ps,
	__isl_take isl_schedule *schedule, struct ppcg_options *options)
{
	int hidden;
	isl_set *context;

	p = isl_printer_start_line(p);
	p = isl_printer_print_str(p, "/* ppcg generated CPU code */");
	p = isl_printer_end_line(p);

	p = isl_printer_start_line(p);
	p = isl_printer_end_line(p);

	p = ppcg_set_macro_names(p);
	p = ppcg_print_exposed_declarations(p, ps);
	hidden = ppcg_scop_any_hidden_declarations(ps);
	if (hidden) {
		p = ppcg_start_block(p);
		p = ppcg_print_hidden_declarations(p, ps);
	}

	context = isl_set_copy(ps->context);
	context = isl_set_from_params(context);
	schedule = isl_schedule_insert_context(schedule, context);
	if (options->debug->dump_final_schedule)
		isl_schedule_dump(schedule);
	p = print_scop(ps, schedule, p, options);
	if (hidden)
		p = ppcg_end_block(p);

	return p;
}

/* Generate CPU code for the scop "ps" and print the corresponding C code
 * to "p", including variable declarations.
 */
__isl_give isl_printer *print_cpu(__isl_take isl_printer *p,
	struct ppcg_scop *ps, struct ppcg_options *options)
{
	isl_schedule *schedule;

	schedule = isl_schedule_copy(ps->schedule);
	return print_cpu_with_schedule(p, ps, schedule, options);
}

/* Generate CPU code for the scop "ps" using "schedule" and
 * print the corresponding C code to "p", including variable declarations.
 */
/** 自动混合精度 打印代码**/
/** 待修改 ***/
static __isl_give isl_printer *print_cpu_with_amp(__isl_take isl_printer *p, __isl_take isl_schedule *schedule, struct ppcg_options *options, amp_prog *prog)
{
	struct ppcg_scop *ps = prog->scop;
	int hidden;
	isl_set *context;

	p = isl_printer_start_line(p);
	p = isl_printer_print_str(p, "/* ppcg generated CPU code with AMP */");
	p = isl_printer_end_line(p);

	p = isl_printer_start_line(p);
	p = isl_printer_end_line(p);

	p = ppcg_set_macro_names(p);
	p = ppcg_print_exposed_declarations(p, ps);
	hidden = ppcg_scop_any_hidden_declarations(ps);
	if (hidden)
	{
		p = ppcg_start_block(p);
		p = ppcg_print_hidden_declarations(p, ps);
	}

	context = isl_set_copy(ps->context);
	context = isl_set_from_params(context);
	schedule = isl_schedule_insert_context(schedule, context);
	if (options->debug->dump_final_schedule)
		isl_schedule_dump(schedule);

	// 打印AST ?
	p = print_scop_with_amp(schedule, p, options, prog);

	if (hidden)
		p = ppcg_end_block(p);

	// 用完释放掉amp_prog
	amp_prog_free(prog);
	return p;
}

/* Generate CPU code for "scop" and print it to "p".
 *
 * First obtain a schedule for "scop" and then print code for "scop"
 * using that schedule.
 * 
 * 待进一步修改和完善。
 */
static __isl_give isl_printer *generate(__isl_take isl_printer *p,
	struct ppcg_scop *scop, struct ppcg_options *options)
{
#define DEBUG_GENERATE
// 调试显示参数
#ifdef DEBUG_GENERATE
	printf("\n@DEBUG: \n       automatic mixed precision paramaters are on the below:\n");
	printf("              the amp is   : %d ( 1==on, 0==off ) \n", options->automatic_mixed_precision);
	printf("              the amp rate : %d ( e.g: 2 means - the higher precision accounts for 1/2 )\n", options->automatic_mixed_precision_rate);
	printf("\n\n");
#endif // DEBUG_GENERATE

	isl_schedule *schedule;
	// 如果进行自动混合精度
	if (options->automatic_mixed_precision)
	{
		if (!scop)
			return isl_printer_free(p);

		// 这里先进行PPCG的调度
		schedule = get_schedule(scop, options);
#ifdef DEBUG_GENERATE
		printf("@DEBUG: \n       ppcg calcu schedule is: \n");
		isl_schedule_dump(schedule);
		printf("\n\n");
#endif // DEBUG_GENERATE

		isl_ctx *ctx = isl_printer_get_ctx(p);
		amp_prog *prog = amp_prog_alloc(ctx, scop);
		if (!prog)
		{
			printf("\n\033[31m@ERROR:\n       There are some errors because the alloced amp_prog is NULL, Now will print cpu code with the ppcg calcuted schedule !!! \n\n\033[0m");
			return print_cpu_with_schedule(p, scop, schedule, options);
		}

		// amp 再调度
		isl_schedule *reschedule = amp_schedule_again(ctx, prog, isl_schedule_copy(schedule));
#ifdef DEBUG_GENERATE
		printf("@DEBUG: \n       amp again calcu schedule is: \n");
		isl_schedule_dump(reschedule);
		printf("\n\n");
#endif // DEBUG_GENERATE
		if (!reschedule || (reschedule == schedule))
		{
			printf("\n\033[31m@ERROR:\n       There are some errors because the schedule calcuted again by amp is NULL or original schedule, Now will print cpu code with the ppcg calcuted schedule !!! \n\n\033[0m");
			amp_prog_free(prog);
			return print_cpu_with_schedule(p, scop, schedule, options);
		}

		isl_schedule_free(schedule);

		return print_cpu_with_amp(p, reschedule, options, prog);
	}
	// 不进行自动混合精度
	schedule = get_schedule(scop, options);

	return print_cpu_with_schedule(p, scop, schedule, options);
}

/* Wrapper around generate for use as a ppcg_transform callback.
 */
static __isl_give isl_printer *print_cpu_wrap(__isl_take isl_printer *p,
	struct ppcg_scop *scop, void *user)
{
	struct ppcg_options *options = user;

	return generate(p, scop, options);
}

/* Transform the code in the file called "input" by replacing
 * all scops by corresponding CPU code and write the results to a file
 * called "output".
 */
int generate_cpu(isl_ctx *ctx, struct ppcg_options *options,
	const char *input, const char *output)
{
	FILE *output_file;
	int r;

	/**
	 * 之前是get_output_file(),现在修改成了get_output_file_with_amp()
	 * 多传进去了options选项，方便判断是否进行自动混合精度（automatic mixed precision，amp）
	 */
	output_file = get_output_file_with_amp(input, output, options);
	if (!output_file)
		return -1;

	r = ppcg_transform(ctx, input, output_file, options,
					&print_cpu_wrap, options);

	fclose(output_file);

	return r;
}
