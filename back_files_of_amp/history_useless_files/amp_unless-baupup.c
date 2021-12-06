/** 这一个注释，记录了全部的专业词汇的翻译。
 * schedule: 调度，表示了程序的执行次序。
 * coincident：重合，是一个长度为 n 的数组。表示一个调度维度是否满足对应的依赖距离为零的意义上的重合约束。
 * permutable：可置换的，即循环是否是可置换的。
 * anchored：锚定的、固定的，如果节点取决于其在调度树中的位置，则设置锚定。 特别是，如果 AST 构建选项包括一个隔离选项
 * 
 * 
 * 
 */

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

/* Construct the string "<a>_<b>".
 */
static char *concat(isl_ctx *ctx, const char *a, const char *b)
{
    isl_printer *p;
    char *s;

    p = isl_printer_to_str(ctx);
    p = isl_printer_print_str(p, a);
    p = isl_printer_print_str(p, "_");
    p = isl_printer_print_str(p, b);
    s = isl_printer_get_str(p);
    isl_printer_free(p);

    return s;
}

/* For each array in "prog" of which an element appears in "accessed" and
 * that is not a read only scalar, create a zero-dimensional universe set
 * of which the tuple id has name "<prefix>_<name of array>" and a user
 * pointer pointing to the array (gpu_array_info).
 *
 * If the array is local to "prog", then make sure it will be declared
 * in the host code.
 *
 * Return the list of these universe sets.
 */
static __isl_give isl_union_set_list *create_data_conversion_filters(amp_prog *prog, const char *prefix, __isl_take isl_union_set *accessed)
{
    int i;
    isl_ctx *ctx;
    isl_union_set_list *filters;

    ctx = prog->ctx;
    filters = isl_union_set_list_alloc(ctx, 0);
    for (i = 0; i < prog->n_array; ++i)
    {
        struct amp_array_info *array = &prog->array[i];
        isl_space *space;
        isl_set *accessed_i;
        int empty;
        char *name;
        isl_id *id;
        isl_union_set *uset;

        if (amp_array_is_read_only_scalar(array))
            continue;

        space = isl_space_copy(array->space);
        accessed_i = isl_union_set_extract_set(accessed, space);
        empty = isl_set_plain_is_empty(accessed_i);
        isl_set_free(accessed_i);
        if (empty < 0)
        {
            filters = isl_union_set_list_free(filters);
            break;
        }
        if (empty)
            continue;

        array->global = 1;
        if (array->local)
            array->declare_local = 1;

        name = concat(ctx, prefix, array->name);
        id = name ? isl_id_alloc(ctx, name, array) : NULL;
        free(name);
        space = isl_space_set_alloc(ctx, 0, 0);
        space = isl_space_set_tuple_id(space, isl_dim_set, id);
        uset = isl_union_set_from_set(isl_set_universe(space));

        filters = isl_union_set_list_add(filters, uset);
    }
    isl_union_set_free(accessed);

    return filters;
}

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
    printf("\n\033[31m@ERROR:\n       the number of the rate is incorrect! this time the rate is %ld, \n       but the abs of the sub of v0 and v1 is %ld !!!\033[0m\n\n", rate, d);
    return isl_bool_error;
}

// amp修改调度
__isl_give isl_schedule *get_amp_schedule(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule *sched, unsigned long rate)
{
#define DEBUG_GET_AMP_SCHEDULE
    isl_schedule *schedule;
    isl_schedule_node *node;
    isl_union_set_list *uset_list, *uset_list_1, *uset_list_2;
    isl_union_set *domain, *uset1, *uset2;
    isl_basic_set_list *bset_list, *bset_list_new;
    isl_basic_set *bset, *bset_new;
    isl_constraint_list *cons_list;
    isl_constraint *c, *c1, *c2, *c3;
    isl_size basic_set_list_dim, constraint_list_dim;

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       the node of the get_amp_schedule is:\n");
    isl_schedule_dump(sched);
    printf("\n\n");
#endif // DEBUG_GET_AMP_SCHEDULE

    // 获取调度树的根结点
    // node = isl_schedule_node_from_domain(isl_union_set_copy(domain)); 这里是不对的，具体区别还不确定。
    node = isl_schedule_get_root(sched);
    // 获取domain(isl_schedule -> isl_union_set)
    domain = isl_schedule_get_domain(sched);
    // 释放之前的调度
    isl_schedule_free(sched);

    //  isl_union_set -> isl_basic_set_list
    bset_list = isl_union_set_get_basic_set_list(domain);
    // 再得到一个副本
    // bset_list_new = isl_union_set_get_basic_set_list(domain);
    // 获取isl_basic_set_list 的维度
    basic_set_list_dim = isl_basic_set_list_n_basic_set(bset_list);
    // 初始化最终的 isl_union_set_list
    uset_list = isl_union_set_list_alloc(ctx, 2);
    uset_list_1 = isl_union_set_list_alloc(ctx, basic_set_list_dim);
    uset_list_2 = isl_union_set_list_alloc(ctx, basic_set_list_dim);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       domain、node、set 、basci_set_list、basic_set_list_dim is:\n");
    isl_union_set_dump(domain);
    isl_schedule_node_dump(node);
    isl_basic_set_list_dump(bset_list);
    printf("basic_set_list_dim is :%d \n\n", basic_set_list_dim);
#endif // DEBUG_GET_AMP_SCHEDULE
    isl_union_set_free(domain);
    for (isl_size i = 0; i < basic_set_list_dim; ++i)
    {
        // isl_basic_set_list -> isl_basic_set
        bset = isl_basic_set_list_get_basic_set(bset_list, i);
        bset_new = isl_basic_set_list_get_basic_set(bset_list, i);
        // isl_basic_set -> isl_constraint_list
        cons_list = isl_basic_set_get_constraint_list(bset);
        // 获取isl_constraint_list的数量（维度）
        constraint_list_dim = isl_constraint_list_n_constraint(cons_list);

#ifdef DEBUG_GET_AMP_SCHEDULE
        printf("@DEBUG: \n       bset_list、basci_set、isl_constraint_list、constraint_list_dim is:\n");
        isl_basic_set_list_dump(bset_list);
        isl_basic_set_dump(bset);
        isl_constraint_list_dump(cons_list);
        printf("when i is = %d, constraint_list_dim is:%d \n\n", i, constraint_list_dim);
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
#ifdef DEBUG_GET_AMP_SCHEDULE
        printf("@DEBUG: \n       the two bset is:\n");
        isl_basic_set_dump(bset);
        isl_basic_set_dump(bset_new);
        printf("when i is = %d\n\n", i);
#endif

        // 释放掉constraint list
        isl_constraint_list_free(cons_list);

        // // 更新原始的isl_bacic_set_list
        // bset_list = isl_basic_set_list_drop(bset_list, i, 1);
        // bset_list = isl_basic_set_list_insert(bset_list, i, bset);
        // // 更新新生成的isl_bacic_set_list
        // bset_list_new = isl_basic_set_list_drop(bset_list_new, i, 1);
        // bset_list_new = isl_basic_set_list_insert(bset_list_new, i, bset_new);

        uset_list_1 = isl_union_set_list_insert(uset_list_1, i, isl_union_set_from_basic_set(bset));
        uset_list_2 = isl_union_set_list_insert(uset_list_2, i, isl_union_set_from_basic_set(bset_new));

        // #ifdef DEBUG_GET_AMP_SCHEDULE
        //         printf("@DEBUG: \n       the two isl_basic_set_list is:\n");
        //         isl_basic_set_list_dump(bset_list);
        //         isl_basic_set_list_dump(bset_list_new);
        //         printf("when i is = %d\n\n", i);
        // #endif
    }

    printf("\n     the two union set list is: \n");
    isl_union_set_list_dump(uset_list_1);
    isl_union_set_list_dump(uset_list_2);

    // intersect应该是去交集
    // uset1 = isl_union_set_from_basic_set(isl_basic_set_list_intersect(bset_list));
    // uset2 = isl_union_set_from_basic_set(isl_basic_set_list_intersect(bset_list_new));
    // isl_basic_set_list_intersect(bset_list);
    uset1 = isl_union_set_list_union(uset_list_1);
    uset2 = isl_union_set_list_union(uset_list_2);

    printf("\n\n   the two union set is: \n");
    isl_union_set_dump(uset1);
    isl_union_set_dump(uset2);

    // isl_set *set1 = isl_set_from_basic_set(bset_list);
    // isl_set *set2 = isl_set_from_basic_set(bset_list_new);
    // isl_set *set = isl_set_union_disjoint(set1, set2);

    // isl_set_dump(set);
    // #ifdef DEBUG_GET_AMP_SCHEDULE
    //     printf("@DEBUG: \n       the two uset is:\n");
    //     isl_union_set_dump(uset1);
    //     isl_union_set_dump(uset2);
    //     printf("\n\n");
    // #endif
    uset_list = isl_union_set_list_add(uset_list, uset1);
    uset_list = isl_union_set_list_add(uset_list, uset2);
    node = isl_schedule_node_child(node, 0);
    node = isl_schedule_node_insert_sequence(node, uset_list);
    // amp 添加数据转换,相关代码见amp_unless-baupup.c
    // node = amp_add_data_conversion(ctx, prog, node);
    // node = isl_schedule_node_child(node, 1);
    // node = isl_schedule_node_child(node, 0);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       after insert lower calcu and before amp_create_kernel, the node info is\n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif

    // amp_create_kernel
    node = isl_schedule_node_child(node, 1);
    node = isl_schedule_node_child(node, 0);
    isl_id *id = isl_id_alloc(ctx, "thread", NULL);
    node = isl_schedule_node_insert_mark(node, id);
    node = isl_schedule_node_parent(node);
    node = amp_create_kernel(prog, node);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       after amp_create_kernel, the node is\n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif

    schedule = isl_schedule_node_get_schedule(node);
    isl_schedule_node_free(node);

#ifdef DEBUG_GET_AMP_SCHEDULE
    printf("@DEBUG: \n       amp generate's schedule  is:  \n");
    isl_schedule_dump(schedule);
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

/* Return the set of parameter values for which the array has a positive
 * size in all dimensions.
 * If the sizes are only valid for some parameter values, then those
 * constraints are also taken into account.
 */
__isl_give isl_set *amp_array_positive_size_guard(struct amp_array_info *array)
{
    int i;
    isl_space *space;
    isl_set *guard;

    if (!array)
        return NULL;

    space = isl_space_params(isl_space_copy(array->space));
    guard = isl_set_universe(space);

    for (i = 0; i < array->n_index; ++i)
    {
        isl_pw_aff *bound;
        isl_set *guard_i, *zero;

        bound = isl_multi_pw_aff_get_pw_aff(array->bound, i);
        guard_i = isl_pw_aff_nonneg_set(isl_pw_aff_copy(bound));
        zero = isl_pw_aff_zero_set(bound);
        guard_i = isl_set_subtract(guard_i, zero);
        guard = isl_set_intersect(guard, guard_i);
    }

    return guard;
}

/* Make sure that code for the statements in "filters" that
 * copy arrays to or from the device is only generated when
 * the size of the corresponding array is positive.
 * That is, add a set node underneath "graft" with "filters" as children
 * and for each child add a guard that the selects the parameter
 * values for which the corresponding array has a positive size.
 * The array is available in the user pointer of the statement identifier.
 * "depth" is the schedule depth of the position where "graft"
 * will be added.
 */
static __isl_give isl_schedule_node *insert_positive_size_guards(
    __isl_take isl_schedule_node *graft,
    __isl_take isl_union_set_list *filters, int depth)
{
    int i, n;

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_set(graft, filters);
    n = isl_schedule_node_n_children(graft);
    for (i = 0; i < n; ++i)
    {
        isl_union_set *filter;
        isl_set *domain, *guard;
        isl_id *id;
        struct amp_array_info *array;

        graft = isl_schedule_node_child(graft, i);
        filter = isl_schedule_node_filter_get_filter(graft);
        domain = isl_set_from_union_set(filter);
        id = isl_set_get_tuple_id(domain);
        array = isl_id_get_user(id);
        isl_id_free(id);
        isl_set_free(domain);
        guard = amp_array_positive_size_guard(array);
        guard = isl_set_from_params(guard);
        guard = isl_set_add_dims(guard, isl_dim_set, depth);
        graft = isl_schedule_node_child(graft, 0);
        graft = isl_schedule_node_insert_guard(graft, guard);
        graft = isl_schedule_node_parent(graft);
        graft = isl_schedule_node_parent(graft);
    }
    graft = isl_schedule_node_parent(graft);

    return graft;
}

/* Create a graft for copying arrays to or from the device,
 * whenever the size of the array is strictly positive.
 * Each statement is called "<prefix>_<name of array>" and
 * the identifier has a user pointer pointing to the array.
 * The graft will be added at the position specified by "node".
 * "copy" contains the array elements that need to be copied.
 * Only arrays of which some elements need to be copied
 * will have a corresponding statement in the graph.
 * Note though that each such statement will copy the entire array.
 */
static __isl_give isl_schedule_node *create_data_conversion(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_keep isl_schedule_node *node, const char *prefix,
                                                            __isl_take isl_union_set *copy)
{
    int depth;
    isl_space *space;
    isl_union_set *all, *domain;
    isl_union_set_list *filters;
    isl_union_map *extension;
    isl_schedule_node *graft;

    depth = isl_schedule_node_get_schedule_depth(node);
    filters = create_data_conversion_filters(prog, prefix, copy);
    all = isl_union_set_list_union(isl_union_set_list_copy(filters));

    space = depth < 0 ? NULL : isl_space_set_alloc(ctx, 0, depth);
    domain = isl_union_set_from_set(isl_set_universe(space));
    extension = isl_union_map_from_domain_and_range(domain, all);
    graft = isl_schedule_node_from_extension(extension);

    if (!filters)
        return isl_schedule_node_free(graft);
    if (isl_union_set_list_n_union_set(filters) == 0)
    {
        isl_union_set_list_free(filters);
        return graft;
    }

    return insert_positive_size_guards(graft, filters, depth);
    // return graft;
}

// amp添加数据转换
__isl_give isl_schedule_node *amp_add_data_conversion(__isl_keep isl_ctx *ctx, amp_prog *prog, __isl_take isl_schedule_node *node)
{
    // #define DEBUG_AMP_ADD_DATA_CONVERSION

    isl_union_map *may_write, *must_write, *copy_out, *not_written;
    isl_union_map *read, *copy_in;
    isl_schedule_node *graft;

    read = isl_union_map_copy(prog->read);
    may_write = isl_union_map_copy(prog->may_write);
    must_write = isl_union_map_copy(prog->must_write);
    not_written = isl_union_map_subtract(may_write, must_write);
    copy_in = isl_union_map_union(read, not_written);
    copy_in = isl_union_map_apply_range(copy_in, isl_union_map_copy(prog->to_outer));
    copy_out = isl_union_map_copy(may_write);

    // 获取低精度所在的filter对应的node结点，注意这里的1和0的含义
    // isl_schedule_node *lower_filter_node = isl_schedule_node_get_child(isl_schedule_node_get_child(node, 1), 0);
    isl_schedule_node *lower_filter_node = find_amp_lower_precision_calculate_insert_position_of_node(node);

#ifdef DEBUG_AMP_ADD_DATA_CONVERSION
    printf("@DEBUG: \n       the lower precision calculate filter node is :\n");
    isl_schedule_node_dump(lower_filter_node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_DATA_CONVERSION

    // 添加原始数据到低精度数据转换的filter结点
    graft = create_data_conversion(ctx, prog, node, "to_lower_precision_data", isl_union_map_range(copy_in));
    lower_filter_node = isl_schedule_node_graft_before(lower_filter_node, graft);

#ifdef DEBUG_AMP_ADD_DATA_CONVERSION
    printf("@DEBUG: \n       before the lower precision calculate filter node add date coversion's result node is :\n");
    isl_schedule_node_dump(lower_filter_node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_DATA_CONVERSION

    // 添加低精度数据到原始数据转换的filter结点
    graft = create_data_conversion(ctx, prog, node, "from_lower_precision_data", isl_union_map_range(copy_out));
    lower_filter_node = isl_schedule_node_graft_after(lower_filter_node, graft);

#ifdef DEBUG_AMP_ADD_DATA_CONVERSION
    printf("@DEBUG: \n       after the lower precision calculate filter node add date coversion's result node is :\n");
    isl_schedule_node_dump(lower_filter_node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_DATA_CONVERSION

    isl_schedule_node_free(node);

    return lower_filter_node;
}

/*
// 将strRes中的t替换为s，替换成功返回1，否则返回0。
int strReplace(char strRes[], char from[], char to[])
{
    int i, flag = 0;
    char *p, *q, *ts;
    for (i = 0; strRes[i]; ++i)
    {
        if (strRes[i] == from[0])
        {
            p = strRes + i;
            q = from;
            while (*q && (*p++ == *q++))
                ;
            if (*q == '\0')
            {
                ts = (char *)malloc(strlen(strRes) + 1);
                strcpy(ts, p);
                strRes[i] = '\0';
                strcat(strRes, to);
                strcat(strRes, ts);
                free(ts);
                flag = 1;
            }
        }
    }
    return flag;
}

// 增加混合精度的混合因子到context中

__isl_give isl_schedule *amp_add_mixing_factor_context(isl_ctx *ctx, __isl_take isl_set *context, __isl_keep isl_schedule *sched)
{
    isl_schedule *schedule, *sched_1, *sched_2;
    isl_union_set *uset;

    if (!sched && !context)
        return NULL;

    char *str_1 = isl_schedule_to_str(isl_schedule_copy(sched));
    char *str_2 = isl_schedule_to_str(isl_schedule_copy(sched));
    isl_schedule_free(sched);

#ifdef DEBUG
    printf("@DEBUG: \n       str_1 is: %s\n", str_1);
    printf("@DEBUG: \n       str_2 is: %s\n", str_2);
#endif // DEBUG

    if (!strReplace(str_1, "31", "15"))
        printf("字符串替换出现异常！！！\n");
    sched_1 = isl_schedule_read_from_str(ctx, str_1);

    if (!strReplace(str_2, "0 <= t", "16 <= t"))
        printf("字符串替换出现异常！！！\n");
    sched_2 = isl_schedule_read_from_str(ctx, str_2);

#ifdef DEBUG
    printf("@DEBUG: \n       schedule is: %s\n", isl_schedule_to_str(sched_1));
    printf("@DEBUG: \n       schedule is: %s\n", isl_schedule_to_str(sched_2));
#endif // DEBUG

    schedule = isl_schedule_sequence(sched_1, sched_2);
    //isl_schedule_free(sched_1);
    //isl_schedule_free(sched_2);

    return schedule;
}
*/

static __isl_give isl_schedule_node *amp_add_copies_group(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
#define DEBUG_AMP_ADD_COPIES_GROUP

    struct amp_array_tile *tile;
    isl_union_map *access, *access2;
    isl_union_set *domain, *domain2;
    isl_multi_aff *ma;
    isl_multi_aff *from_access;
    isl_multi_pw_aff *mpa;
    isl_multi_union_pw_aff *mupa;
    isl_schedule_node *graft;
    isl_union_set *filter;
    isl_space *space;
    int skip;
    int kernel_depth;
    int empty;

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       at start of the amp add copies group function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n       the read is :%d \n\n", read);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    access = anchored_non_local_accesses(kernel, group, node, read);
    access2 = anchored_non_local_accesses(kernel, group, node, !read);
    empty = isl_union_map_is_empty(access);
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        if (empty < 0)
            return isl_schedule_node_free(node);
        return amp_tree_move_up_to_kernel(node);
        // return node;
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after anchored_non_local_accesses ,the two access is: \n");
    isl_union_map_dump(access);
    printf("\n");
    isl_union_map_dump(access2);
    printf("\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the from_access variable is: \n");
    isl_multi_aff_dump(from_access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
       /*
    tile = amp_array_ref_group_tile(group);
    ma = isl_multi_aff_copy(tile->tiling);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the ma is: \n");
    isl_multi_aff_dump(ma);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    ma = isl_multi_aff_pullback_multi_aff(ma, isl_multi_aff_copy(from_access));
    mpa = isl_multi_pw_aff_from_multi_aff(ma);
    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);
    domain = isl_union_map_range(access);
    domain2 = isl_union_map_range(access2);
    if (read && !amp_array_is_scalar(group->array))
    {
        isl_map *map;
        isl_union_set_free(domain);
        isl_union_set_free(domain2);
        map = group_tile(group);
        domain = isl_union_set_from_set(isl_map_wrap(isl_map_copy(map)));
        domain2 = isl_union_set_from_set(isl_map_wrap(map));

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
        printf("@DEBUG: \n       after group_tile, the domain and domian2 is: \n");
        isl_union_set_dump(domain);
        isl_union_set_dump(domain2);
        printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    }

    domain = isl_union_set_preimage_multi_aff(domain, isl_multi_aff_copy(from_access));
    domain2 = isl_union_set_preimage_multi_aff(domain2, from_access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       the two domain is: \n");
    isl_union_set_dump(domain);
    isl_union_set_dump(domain2);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    // if (read)
    // {
    //     access = isl_union_set_wrapped_domain_map(domain);
    // }
    // else
    // {
    //     access = isl_union_set_wrapped_domain_map(domain2);
    // }
    */
    space = isl_space_domain(isl_multi_aff_get_space(from_access));
    access = isl_union_map_preimage_range_multi_aff(access, from_access);

    // filter = isl_union_set_copy(kernel->thread_filter);
    // contraction = isl_union_pw_multi_aff_copy(kernel->contraction);
    // filter = isl_union_set_preimage_union_pw_multi_aff(filter, contraction);
    // filter = isl_union_set_apply(filter, isl_union_map_copy(access));
    // filter = isl_union_set_detect_equalities(filter);
    // filter = isl_union_set_coalesce(filter);
    space = isl_space_map_from_set(space);
    mpa = isl_multi_pw_aff_identity(space);
    mpa = isl_multi_pw_aff_range_factor_range(mpa);
    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);

    domain = isl_union_map_range(access);
    // 基本上从下面开始就没办法修改了
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
    printf("@DEBUG: \n       after the isl_schedule_node_from_extension(access), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, mupa);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       after the isl_schedule_node_insert_partial_schedule(graft, mupa), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n       the mupa(isl_multi_union_pw_aff) is : \n");
    isl_multi_union_pw_aff_dump(mupa);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP
    if (kernel->options->unroll_copy_shared)
        graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    /*
    if (tile->n > kernel->n_block && kernel->n_block > 0)
    {
        graft = isl_schedule_node_band_split(graft, tile->n - kernel->n_block);
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
    printf("@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP

    return node;
}

static __isl_give isl_schedule_node *amp_add_copies_group_read(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
#define DEBUG_AMP_ADD_COPIES_GROUP_READ

    struct amp_array_tile *tile;
    isl_union_map *access;
    isl_union_set *domain;
    isl_multi_aff *from_access;
    isl_multi_pw_aff *mpa;
    isl_multi_union_pw_aff *mupa;
    isl_schedule_node *graft;
    isl_space *space;
    int kernel_depth;
    int empty;

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       at start of the amp_add_copies_group_read function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n       the read is :%d \n\n", read);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);

    // #ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    //     printf("@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    //     isl_schedule_node_dump(node);
    //     printf("\n\n");
    // #endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    access = anchored_non_local_accesses(kernel, group, node, read);
    empty = isl_union_map_is_empty(access);
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        if (empty < 0)
            return isl_schedule_node_free(node);
        return amp_tree_move_up_to_kernel(node);
        // return node;
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       after anchored_non_local_accesses ,the access is: \n");
    isl_union_map_dump(access);
    printf("\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       the from_access variable is: \n");
    isl_multi_aff_dump(from_access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    space = isl_space_domain(isl_multi_aff_get_space(from_access));
    access = isl_union_map_preimage_range_multi_aff(access, from_access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       the space, the  isl_union_map_preimage_range_multi_aff(access, from_access) is: \n");
    isl_space_dump(space);
    printf("\n");
    isl_union_map_dump(access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    space = isl_space_map_from_set(space);
    mpa = isl_multi_pw_aff_identity(space);
    mpa = isl_multi_pw_aff_range_factor_range(mpa);
    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       the mupa is: \n");
    isl_multi_union_pw_aff_dump(mupa);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    domain = isl_union_map_range(access);
    // 基本上从下面开始就没办法修改了
    access = isl_union_set_wrapped_domain_map(domain);
    access = isl_union_map_reverse(access);
    access = isl_union_map_coalesce(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       before the isl_schedule_node_from_extension(access), the access is: \n");
    isl_union_map_dump(access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    graft = isl_schedule_node_from_extension(access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       after the isl_schedule_node_from_extension(access), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, mupa);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       after the isl_schedule_node_insert_partial_schedule(graft, mupa), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n       the mupa(isl_multi_union_pw_aff) is : \n");
    isl_multi_union_pw_aff_dump(mupa);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ
    if (kernel->options->unroll_copy_shared)
        graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    while (graft && isl_schedule_node_has_parent(graft))
        graft = isl_schedule_node_parent(graft);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       the final graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    // node = find_amp_lower_precision_calculate_insert_position_of_node(node);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       before insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

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

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_READ
    printf("@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_READ

    return node;
}

static __isl_give isl_schedule_node *amp_add_copies_group_write(
    struct amp_ppcg_kernel *kernel, struct amp_array_ref_group *group,
    __isl_take isl_schedule_node *node, int read)
{
#define DEBUG_AMP_ADD_COPIES_GROUP_WRITE

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
    int skip;
    int kernel_depth;
    int empty;

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       at the amp_add_copies_group_write function, the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n       the read is :%d \n\n", read);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    kernel_depth = isl_schedule_node_get_schedule_depth(node);
    node = amp_tree_move_down_to_depth(node, kernel_depth, kernel->core);

    // #ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    //     printf("@DEBUG: \n       after amp tree move down to depth, the node is: \n");
    //     isl_schedule_node_dump(node);
    //     printf("\n\n");
    // #endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    access = anchored_non_local_accesses(kernel, group, node, read);
    empty = isl_union_map_is_empty(access);
    if (empty < 0 || empty)
    {
        isl_union_map_free(access);
        if (empty < 0)
            return isl_schedule_node_free(node);
        return amp_tree_move_up_to_kernel(node);
        // return node;
    }

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       after anchored_non_local_accesses ,the two access is: \n");
    isl_union_map_dump(access);
    printf("\n       the empty number is :%d \n\n", empty);
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    group->array->global = 1;
    group->local_array->global = 1;

    from_access = create_from_access(kernel->ctx, group, read);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       the from_access variable is: \n");
    isl_multi_aff_dump(from_access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    tile = amp_array_ref_group_tile(group);
    ma = isl_multi_aff_copy(tile->tiling);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       the ma is: \n");
    isl_multi_aff_dump(ma);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    ma = isl_multi_aff_pullback_multi_aff(ma, isl_multi_aff_copy(from_access));
    mpa = isl_multi_pw_aff_from_multi_aff(ma);
    mupa = isl_multi_union_pw_aff_from_multi_pw_aff(mpa);

    domain = isl_union_map_range(access);

    domain = isl_union_set_preimage_multi_aff(domain, isl_multi_aff_copy(from_access));
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       the two domain is: \n");
    isl_union_set_dump(domain);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    // 基本上从下面开始就没办法修改了
    access = isl_union_set_wrapped_domain_map(domain);
    access = isl_union_map_reverse(access);
    access = isl_union_map_coalesce(access);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       before the graft = isl_schedule_node_from_extension(access);the access is: \n");
    isl_union_map_dump(access);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    graft = isl_schedule_node_from_extension(access);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       after the isl_schedule_node_from_extension(access), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    graft = isl_schedule_node_child(graft, 0);
    graft = isl_schedule_node_insert_partial_schedule(graft, mupa);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       after the isl_schedule_node_insert_partial_schedule(graft, mupa), the graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n       the mupa(isl_multi_union_pw_aff) is : \n");
    isl_multi_union_pw_aff_dump(mupa);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    if (kernel->options->unroll_copy_shared)
        graft = ppcg_set_schedule_node_type(graft, isl_ast_loop_unroll);

    while (graft && isl_schedule_node_has_parent(graft))
        graft = isl_schedule_node_parent(graft);

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       the final graft is: \n");
    isl_schedule_node_dump(graft);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    // node = find_amp_lower_precision_calculate_insert_position_of_node(node);
#ifdef DEBUG_AMP_ADD_COPIES_GROUP
    printf("@DEBUG: \n       before insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

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

#ifdef DEBUG_AMP_ADD_COPIES_GROUP_WRITE
    printf("@DEBUG: \n       after insert the final graft,the node is: \n");
    isl_schedule_node_dump(node);
    printf("\n\n");
#endif // DEBUG_AMP_ADD_COPIES_GROUP_WRITE

    return node;
}

// find the suitable position for lower precision calculate insert, in the node.
__isl_give isl_schedule_node *find_amp_lower_precision_calculate_insert_position_of_node(__isl_take isl_schedule_node *node)
{
    // #define DEBUG_FIND_AMP_LOWER_PRECISION_POSITION

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

/**
 * @brief 计算自动混合精度依据rate进行划分,划分之后的划分值.即,将[x,y]划分成[x,amp_partition_val - 1]和[amp_partition_val,y]两个区间,这里只返回amp_partition_val的值.
 *            同时,需要注意的是,如果 rate > (y - x),返回的其实是x的值,也就是说rate > 区间长度时,返回的是区间起点.
 * @note  这里要求x<=0、y>0、rate>1
 */
__isl_give isl_val *amp_get_partition_val(__isl_take isl_val *x, __isl_take isl_val *y, int rate)
{
#define DEBUG_AMP_GET_PARTITION_VAL

    isl_val *diff;

    // 如果' x>0 或 y<=0 或 rate<=1 '
    // if ((isl_val_cmp_si(x, 0) > 0) || (isl_val_cmp_si(y, 0) <= 0) || rate <= 1)
    // {
    //     printf("\n\033[31m@ERROR:\n       There are some errors because x>0 或 y<0 或 rate <=1, Now will return NULL !!! \033[0m\n\n");
    //     goto error;
    // }
    if (rate <= 1)
    {
        printf("\n\033[31m@ERROR:\n       There are some errors because rate <=1, Now will return NULL !!! \033[0m\n\n");
        goto error;
    }

    // 先求区间长度值
    diff = isl_val_sub(isl_val_abs(isl_val_copy(y)), isl_val_abs(isl_val_copy(x)));
    diff = isl_val_add_ui(diff, 1);

#ifdef DEBUG_AMP_GET_PARTITION_VAL
    printf("@DEBUG: \n       the length of |y| - |x| is : \n");
    isl_val_dump(diff);
    printf("\n\n");
#endif // DEBUG_AMP_GET_PARTITION_VAL

    if (!isl_val_cmp_si(diff, 0))
    {
        printf("\n\033[31m@ERROR:\n       There are some errors because |y| - |x| <= 0, Now will return NULL !!! \033[0m\n\n");
        goto error;
    }

    // 计算amp_partition_val(diff的最后即是)
    // diff = isl_val_set_si(diff, (int)(isl_val_get_d(diff) / rate));
    isl_val *rt = isl_val_copy(diff);
    diff = isl_val_div(isl_val_copy(diff), isl_val_set_si(rt, rate));
    diff = isl_val_floor(diff);
    diff = isl_val_add(diff, isl_val_abs(isl_val_copy(x)));

#ifdef DEBUG_AMP_GET_PARTITION_VAL
    printf("@DEBUG: \n       the partition_val is : \n");
    isl_val_dump(diff);
    printf("\n\n");
#endif // DEBUG_AMP_GET_PARTITION_VAL

    if ((isl_val_cmp_si(diff, 0) < 0) || !isl_val_is_int(diff))
    {
        printf("\n\033[31m@ERROR:\n       There are some errors because amp_partition_val < 0, Now will return NULL !!! \033[0m\n\n");
        goto error;
    }

    isl_val_free(x);
    isl_val_free(y);

    /*******
     * 
     * 
     * 
     * 
     * 
     * 
     *
    // 计算结果的值
    // isl_aff *rs = isl_aff_div(isl_aff_sub(aff_y, isl_aff_copy(aff_x)), rt);
    // rs = isl_aff_add(rs, aff_x);
    // rs = isl_aff_floor(rs);

    // printf(" isl_aff *rs is : \n");
    // isl_aff_dump(rs);

    // isl_val *val = isl_aff_get_constant_val(a3);
    // val = isl_val_floor(val);
    // isl_val_dump(val);

    // isl_val_dump(isl_aff_get_coefficient_val(a3, isl_dim_set, 0));
    // isl_aff_dump(a3);

    // #ifdef DEBUG_AMP_GET_SINGLE_STATEMENT_CONSTRAINTS
    //             printf("@DEBUG: \n       the v_x and v_y is : \n");
    //             isl_val_dump(v_x);
    //             isl_val_dump(v_y);
    //             printf("\n\n");
    // #endif // DEBUG_AMP_GET_SINGLE_STATEMENT_CONSTRAINTS

    // 对应rate时,amp对前后的划分的中间的划分值amp_partition_val
    // isl_val *amp_partition_val = amp_get_partition_val(isl_val_copy(v_x), v_y, rate);

    // if (!amp_partition_val)
    //     goto error;

    // // 更新upper的约束.如果|amp_amp_partition_val| == |v_x|,则说明upper应该为空,否则,upper的约束范围应该为[v_x,amp_partition_val - 1]
    // if (isl_val_abs_eq(amp_partition_val, v_x))
    // {
    //     amp_bset->flag = isl_bool_false;
    //     isl_constraint_free(con_back);
    //     isl_basic_set_free(amp_bset->upper);
    // }
    // else
    // {
    //     amp_bset->flag = isl_bool_true;
    //     con_back = isl_constraint_set_constant_val(con_back, isl_val_sub_ui(isl_val_copy(amp_partition_val), 1));
    //     amp_bset->upper = isl_basic_set_add_constraint(amp_bset->upper, con_back);
    // }
    // isl_val_free(v_x);

    // // 更新lower的约束,lower的约束范围应该为[amp_partition_val - 1, v_y]
    // con_front = isl_constraint_set_constant_val(con_front, isl_val_neg(amp_partition_val));

    // isl_aff *a = isl_constraint_get_aff(con_back);
    // printf(" isl_aff isl_constraint_get_aff(con_back) is : \n");
    // isl_aff_dump(a);

    // // con_back = isl_constraint_set_coefficient_val(con_back, isl_dim_set, 1, isl_aff_get_constant_val(rs));
    // // con_back = isl_constraint_set_constant_val(con_back, isl_aff_get_constant_val(rs));
    // printf(" isl_constraint con_back is : \n");
    // isl_constraint_dump(con_back);
    // printf(" isl_constraint isl_inequality_from_aff(rs) is : \n");
    // isl_constraint_dump(isl_inequality_from_aff(rs));
    // amp_bset->upper = isl_basic_set_add_constraint(amp_bset->upper, con_back);
    // // amp_bset->upper = isl_basic_set_add_constraint(amp_bset->upper, isl_inequality_from_aff(isl_aff_copy(rs)));
    // isl_basic_set_dump(amp_bset->upper);
    // amp_bset->flag = isl_bool_true;

    // // con_front = isl_constraint_set_coefficient_val(con_front, isl_dim_set, 0, isl_aff_get_constant_val(rs));
    // // isl_val_dump(isl_constraint_get_constant_val());
    // // isl_constraint_dump(isl_inequality_from_aff(rs));
    // // amp_bset->lower = isl_basic_set_drop_constraints_involving_dims(amp_bset->lower, isl_dim_set, 1, 1);
    // con_front = isl_constraint_set_constant_val(con_front, isl_aff_get_constant_val(rs));
    // // amp_bset->lower = isl_basic_set_add_constraint(amp_bset->lower, con_front);
    // isl_constraint_dump(con_front);
    // isl_constraint_dump(isl_inequality_from_aff(rs));
    // // amp_bset->lower = isl_basic_set_add_constraint(amp_bset->lower, isl_inequality_from_aff(isl_aff_copy(rs)));
    // // isl_basic_set_dump(amp_bset->lower);
    */

    return diff;
error:
    isl_val_free(x);
    isl_val_free(y);
    isl_val_free(diff);

    return NULL;
}

/**
 * @brief 目前修改完约束之后,往isl_union_set转换的时候,发现isl_basic_set_list不能直接转换回isl_union_set.于是在这部分,换了策略.
 * 将一个isl_basic_set转换成amp_basic_set,其中分别包含了前后两部分的isl_basic_set,将前后两部分的各合并成一个isl_union_set_list,然后对其求并集,获得上下两部分的filter(isl_union_set).
 */
__isl_give struct amp_union_set_list *amp_get_filters(__isl_keep isl_ctx *ctx, __isl_take isl_union_set *domain, int rate)
{
    // #define DEBUG_AMP_GET_FILTERS

    isl_basic_set_list *basic_set_list;
    isl_size basic_set_list_dim;
    isl_basic_set *bset;
    struct amp_basic_set *amp_bset;

    // 目前的往回转换策略用到的
    struct amp_union_set_list *amp_filters = isl_calloc_type(ctx, struct amp_union_set_list);
    isl_union_set *uset_left_filter, *uset_right_filter;
    isl_union_set_list *uset_list_left, *uset_list_right;

    if (!domain)
        goto error;

    //  isl_union_set(domain) -> isl_basic_set_list
    basic_set_list = isl_union_set_get_basic_set_list(domain);
#ifdef DEBUG_AMP_GET_FILTERS
    printf("@DEBUG: \n       (isl_union_set(domain) -> isl_basic_set_list) result is: \n");
    isl_basic_set_list_dump(basic_set_list);
    printf("\n\n");
#endif //  DEBUG_AMP_GET_FILTERS

    // 获取isl_basic_set_list的dim
    basic_set_list_dim = isl_basic_set_list_n_basic_set(basic_set_list);
    // 初始化需要返回的filters(isl_union_set_list)
    amp_filters->filters = isl_union_set_list_alloc(ctx, 2);
    // 初始化高低精度filter的union_set_list
    uset_list_left = isl_union_set_list_alloc(ctx, basic_set_list_dim);
    uset_list_right = isl_union_set_list_alloc(ctx, basic_set_list_dim);

    for (isl_size i = 0; i < basic_set_list_dim; i++)
    {
        // isl_basic_set_list获取对应dim=i时的isl_basic_set
        bset = isl_basic_set_list_get_basic_set(basic_set_list, i);
        if (!isl_basic_set_is_bounded(bset))
            continue;
#ifdef DEBUG_AMP_GET_FILTERS
        printf("@DEBUG: \n       对应dim=%d时的isl_basic_set is: \n", i);
        isl_basic_set_dump(bset);
        printf("\n\n");
#endif //  DEBUG_AMP_GET_FILTERS

        // 获取新的约束,目前默认的修改的deepth(循环的dim)是0
        amp_bset = amp_get_single_statement_constraints(ctx, bset, 0, rate);
        if (!amp_bset)
        {
            printf("\n\033[31m@ERROR:\n       There are some errors because the get_single_statement_constraints function returned amp_bset is NULL , Now will return the original domain in filters way!!! \033[0m\n\n");
            goto error;
        }
        // 逐语句将amp_basic_set的左右分支的isl_basic_set合并成isl_union_set_list
        if (amp_bset->left_flag)
            uset_list_left = isl_union_set_list_insert(uset_list_left, i, isl_union_set_from_basic_set(amp_bset->left));
        if (amp_bset->right_flag)
            uset_list_right = isl_union_set_list_insert(uset_list_right, i, isl_union_set_from_basic_set(amp_bset->right));
    }
    // 对isl_union_set_list求并集,获得左右分支的filter(isl_union_set)
    if (amp_bset->left_flag)
        uset_left_filter = isl_union_set_list_union(uset_list_left);
    if (amp_bset->right_flag)
        uset_right_filter = isl_union_set_list_union(uset_list_right);

    // 左右分支的filter合并成amp_filters
    if (amp_bset->left_flag)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_left_filter);
    if (amp_bset->right_flag)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_right_filter);

    // 判断类型
    if (amp_bset->left_flag && amp_bset->right_flag)
        amp_filters->type = amp_filters_both;
    else if (amp_bset->left_flag && !amp_bset->right_flag)
        amp_filters->type = amp_filters_only_left;
    else if (!amp_bset->left_flag && amp_bset->right_flag)
        amp_filters->type = amp_filters_only_right;
    else
        goto error;

    // 释放内存
    isl_union_set_free(domain);
    isl_basic_set_list_free(basic_set_list);

    return amp_filters;
error:
    isl_union_set_free(domain);
    isl_basic_set_list_free(basic_set_list);
    isl_union_set_list_free(uset_list_left);
    isl_union_set_list_free(uset_list_right);
    isl_union_set_free(uset_left_filter);
    isl_union_set_free(uset_right_filter);
    isl_union_set_list_free(amp_filters->filters);

    return NULL;
}

/*********************
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
__isl_give isl_set *isl_basic_set_list_union(__isl_take isl_basic_set_list *list);

/**
 * @brief 目前修改完约束之后,往isl_set转换的时候,发现isl_basic_set_list不能直接调用（）转换回isl_set.于是在这部分,换了策略.
 * 将一个isl_basic_set转换成amp_basic_set,其中分别包含了前后两部分的isl_basic_set,将前后两部分的各合并成一个isl_union_set_list,然后对其求并集,获得上下两部分的filter(isl_union_set).
 */
__isl_give struct amp_union_set_list *amp_get_filters_22(__isl_keep isl_ctx *ctx, __isl_take isl_union_set *domain, int rate)
{
    // #define DEBUG_AMP_GET_FILTERS

    // isl_basic_set_list *basic_set_list;
    // isl_size basic_set_list_dim;
    // isl_basic_set *bset;
    // struct amp_basic_set *amp_bset;

    // // 目前的往回转换策略用到的
    isl_union_set *union_set_left, *union_set_right;
    struct amp_basic_set *amp_bset = amp_basic_set_init(ctx);
    struct amp_union_set_list *amp_filters = isl_calloc_type(ctx, struct amp_union_set_list);
    int zero_deepth_flag = 0;

    isl_set *set, *set_copy;

    int div_dim = 1;
    if (!domain)
        goto error;

    isl_set_list *set_list = isl_union_set_get_set_list(domain);
    isl_size set_deepth = isl_union_set_n_set(domain);

    printf("domain is : \n");
    isl_union_set_dump(domain);
    // 初始化需要返回的filters(isl_union_set_list)
    amp_filters->filters = isl_union_set_list_alloc(ctx, 2);

    isl_set_list *set_list_left, *set_list_right;
    union_set_left = isl_union_set_copy(domain);
    union_set_left = isl_union_set_empty(isl_union_set_get_space(domain));
    union_set_right = isl_union_set_copy(union_set_left);

    for (isl_size set_dim = 0; set_dim < set_deepth; set_dim++)
    {

        set = isl_set_list_get_set(set_list, set_dim);
        set_copy = isl_set_copy(set);

        printf("when set_dim = %d, the set is:\n", set_dim);
        isl_set_dump(set);

        isl_basic_set_list *basic_set_list = isl_set_get_basic_set_list(set);
        isl_basic_set_list *basic_set_list_copy = isl_set_get_basic_set_list(set);
        isl_size basic_set_deepth = isl_basic_set_list_n_basic_set(basic_set_list);

        // // 判断isl_basic_set_list的deepth是否为0
        // if (basic_set_deepth == 0)
        //     zero_deepth_flag = 1;
        // else
        //     zero_deepth_flag = 0;

        for (isl_size basic_set_dim = 0; basic_set_dim < basic_set_deepth; basic_set_dim++)
        {
            set_list_left = isl_set_list_alloc(ctx, set_deepth);
            set_list_right = isl_set_list_alloc(ctx, set_deepth);

            // if (zero_deepth_flag == 1)
            //     break;

            isl_basic_set *basic_set = isl_basic_set_list_get_basic_set(basic_set_list, basic_set_dim);

            printf("when basic_set_dim = %d, the basic_set is:\n", basic_set_dim);
            isl_basic_set_dump(basic_set);

            isl_size con_deepth = isl_basic_set_n_constraint(basic_set);
            if (con_deepth < 1)
                continue;
            if (!amp_basic_set_reset(amp_bset))
                printf("amp_basic_set reset meets errors");
            amp_bset = amp_get_single_statement_constraints(ctx, basic_set, 1, rate);
            if (!amp_bset)
            {
                printf("\n\033[31m@ERROR:\n       There are some errors because the get_single_statement_constraints function returned amp_bset is NULL , Now will return the original domain in filters way!!! \033[0m\n\n");
                goto error;
            }

            printf("1. basic_set_list and copy is: \n");
            isl_basic_set_list_dump(basic_set_list);
            isl_basic_set_list_dump(basic_set_list_copy);

            // // 逐语句将amp_basic_set的左右分支的isl_basic_set合并成isl_union_set_list
            // if (amp_bset->left_flag)
            // {
            //     set_list_left = isl_set_list_add(set_list_left, isl_set_from_basic_set(amp_bset->left));
            // }
            // if (amp_bset->right_flag)
            // {
            //     set_list_right = isl_set_list_add(set_list_right, isl_set_from_basic_set(amp_bset->right));
            // }

            // // 逐语句将amp_basic_set的左右分支的isl_basic_set合并成isl_union_set_list
            // if (amp_bset->left_flag)
            // {
            //     basic_set_list = isl_basic_set_list_drop(basic_set_list, basic_set_dim, 1);
            //     basic_set_list = isl_basic_set_list_insert(basic_set_list, basic_set_dim, amp_bset->left);
            // }
            // if (amp_bset->right_flag)
            // {
            //     basic_set_list_copy = isl_basic_set_list_drop(basic_set_list_copy, basic_set_dim, 1);
            //     basic_set_list_copy = isl_basic_set_list_insert(basic_set_list_copy, basic_set_dim, amp_bset->right);
            // }

            // printf("2. basic_set_list and copy is: \n");
            // isl_basic_set_list_dump(basic_set_list);
            // isl_basic_set_list_dump(basic_set_list_copy);

            // set = isl_basic_set_list_union(basic_set_list);
            // set_copy = isl_basic_set_list_union(basic_set_list_copy);
            // printf("1. set and copy is : \n");
            // isl_set_dump(set);
            // isl_set_dump(set_copy);

            // set_list_left = isl_set_list_add(set_list_left, set);
            // printf("set_list_left is: \n");
            // isl_set_list_dump(set_list_left);

            // set_list_right = isl_set_list_add(set_list_right, set_copy);
            // printf("set_list_right is: \n");
            // isl_set_list_dump(set_list_right);

            // set = isl_set_list_union(set_list_left);
            // printf("unioned set is: \n");
            // isl_set_dump(set);
        }
        // // 对isl_union_set_list求并集,获得左右分支的filter(isl_union_set)
        // if (amp_bset->left_flag)
        //     set = isl_set_list_union(set_list_left);
        // if (amp_bset->right_flag)
        //     set_copy = isl_set_list_union(set_list_right);

        // printf(" 2. set and copy is : \n");
        // isl_set_dump(set);
        // isl_set_dump(set_copy);

        // printf("1. left and right union set is: \n");
        // isl_union_set_dump(union_set_left);
        // isl_union_set_dump(union_set_right);

        // // union_set_left = isl_union_set_add_set(union_set_left, isl_set_list_union(set_list_left));
        // // union_set_right = isl_union_set_add_set(union_set_right, isl_set_list_union(set_list_right));
        // union_set_left = isl_union_set_add_set(union_set_left, set);
        // union_set_right = isl_union_set_add_set(union_set_right, set_copy);

        // printf("2. left and right union set is: \n");
        // isl_union_set_dump(union_set_left);
        // isl_union_set_dump(union_set_right);
    }

    printf("hrt! \n");
    // 左右分支的filter合并成amp_filters
    if (amp_bset->left_flag)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, union_set_left);
    // else if (amp_bset->left_flag == isl_bool_false)
    //     amp_filters->filters = isl_union_set_list_add(amp_filters->filters, union_set_left);

    if (amp_bset->right_flag)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, union_set_right);

    // 判断类型
    if (amp_bset->left_flag && amp_bset->right_flag)
        amp_filters->type = amp_filters_both;
    else if (amp_bset->left_flag && !amp_bset->right_flag)
        amp_filters->type = amp_filters_only_left;
    else if (!amp_bset->left_flag && amp_bset->right_flag)
        amp_filters->type = amp_filters_only_right;

    return amp_filters;

error:

    return NULL;

    //     //  isl_union_set(domain) -> isl_basic_set_list
    //     basic_set_list = isl_union_set_get_basic_set_list(domain);
    // #ifdef DEBUG_AMP_GET_FILTERS
    //     printf("@DEBUG: \n       (isl_union_set(domain) -> isl_basic_set_list) result is: \n");
    //     isl_basic_set_list_dump(basic_set_list);
    //     printf("\n\n");
    // #endif //  DEBUG_AMP_GET_FILTERS

    //     // 获取isl_basic_set_list的dim
    //     basic_set_list_dim = isl_basic_set_list_n_basic_set(basic_set_list);
    //     // 初始化需要返回的filters(isl_union_set_list)
    //     amp_filters->filters = isl_union_set_list_alloc(ctx, 2);
    //     // 初始化高低精度filter的union_set_list
    //     uset_list_left = isl_union_set_list_alloc(ctx, basic_set_list_dim);
    //     uset_list_right = isl_union_set_list_alloc(ctx, basic_set_list_dim);

    //     for (isl_size i = 0; i < basic_set_list_dim; i++)
    //     {
    //         // isl_basic_set_list获取对应dim=i时的isl_basic_set
    //         bset = isl_basic_set_list_get_basic_set(basic_set_list, i);
    //         if (!isl_basic_set_is_bounded(bset))
    //             continue;
    // #ifdef DEBUG_AMP_GET_FILTERS
    //         printf("@DEBUG: \n       对应dim=%d时的isl_basic_set is: \n", i);
    //         isl_basic_set_dump(bset);
    //         printf("\n\n");
    // #endif //  DEBUG_AMP_GET_FILTERS

    //         // 获取新的约束,目前默认的修改的deepth(循环的dim)是0
    //         amp_bset = amp_get_single_statement_constraints(ctx, bset, 0, rate);
    //         if (!amp_bset)
    //         {
    //             printf("\n\033[31m@ERROR:\n       There are some errors because the get_single_statement_constraints function returned amp_bset is NULL , Now will return the original domain in filters way!!! \033[0m\n\n");
    //             goto error;
    //         }
    //         // 逐语句将amp_basic_set的左右分支的isl_basic_set合并成isl_union_set_list
    //         if (amp_bset->left_flag)
    //             uset_list_left = isl_union_set_list_insert(uset_list_left, i, isl_union_set_from_basic_set(amp_bset->left));
    //         if (amp_bset->right_flag)
    //             uset_list_right = isl_union_set_list_insert(uset_list_right, i, isl_union_set_from_basic_set(amp_bset->right));
    //     }
    //     // 对isl_union_set_list求并集,获得左右分支的filter(isl_union_set)
    //     if (amp_bset->left_flag)
    //         uset_left_filter = isl_union_set_list_union(uset_list_left);
    //     if (amp_bset->right_flag)
    //         uset_right_filter = isl_union_set_list_union(uset_list_right);

    //     // 左右分支的filter合并成amp_filters
    //     if (amp_bset->left_flag)
    //         amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_left_filter);
    //     if (amp_bset->right_flag)
    //         amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_right_filter);

    //     // 判断类型
    //     if (amp_bset->left_flag && amp_bset->right_flag)
    //         amp_filters->type = amp_filters_both;
    //     else if (amp_bset->left_flag && !amp_bset->right_flag)
    //         amp_filters->type = amp_filters_only_left;
    //     else if (!amp_bset->left_flag && amp_bset->right_flag)
    //         amp_filters->type = amp_filters_only_right;
    //     else
    //         goto error;

    //     // 释放内存
    //     isl_union_set_free(domain);
    //     isl_basic_set_list_free(basic_set_list);

    //     return amp_filters;
    // error:
    //     isl_union_set_free(domain);
    //     isl_basic_set_list_free(basic_set_list);
    //     isl_union_set_list_free(uset_list_left);
    //     isl_union_set_list_free(uset_list_right);
    //     isl_union_set_free(uset_left_filter);
    //     isl_union_set_free(uset_right_filter);
    //     isl_union_set_list_free(amp_filters->filters);

    //     return NULL;
}

/*****
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
 * 
 * 
 * 
 * 
 * 很好的正确版
 */
/**
 * @brief 目前修改完约束之后,往isl_union_set转换的时候,发现isl_basic_set_list不能直接转换回isl_union_set.于是在这部分,换了策略.
 * 将一个isl_basic_set转换成amp_basic_set,其中分别包含了前后两部分的isl_basic_set,将前后两部分的各合并成一个isl_union_set_list,然后对其求并集,获得上下两部分的filter(isl_union_set).
 */
__isl_give struct amp_union_set_list *amp_get_filters(__isl_keep isl_ctx *ctx, __isl_take isl_union_set *domain, int rate)
{
#define DEBUG_AMP_GET_FILTERS

    // 定义要拆分的循环的dim为1（目前拆第一层循环）
    int separate_dim = 1;
    isl_basic_set_list *basic_set_list;
    isl_size basic_set_list_deepth;
    isl_basic_set *bset;
    // 目前的往回转换策略用到的
    isl_union_set *uset_left_filter, *uset_right_filter;
    isl_union_set_list *uset_list_left, *uset_list_right;
    struct amp_basic_set *amp_bset;
    // 为返回的结构体申请内存
    struct amp_union_set_list *amp_filters = isl_calloc_type(ctx, struct amp_union_set_list);

    if (!domain)
        goto error;

    //  isl_union_set(domain) -> isl_basic_set_list
    basic_set_list = isl_union_set_get_basic_set_list(domain);
    // 获取isl_basic_set_list的deepth
    basic_set_list_deepth = isl_basic_set_list_n_basic_set(basic_set_list);
#ifdef DEBUG_AMP_GET_FILTERS
    printf("@DEBUG: \n       (isl_union_set(domain) -> isl_basic_set_list) result is: \n");
    isl_basic_set_list_dump(basic_set_list);
    printf("\n       the deepth of it is:  %d \n\n", basic_set_list_deepth);
#endif //  DEBUG_AMP_GET_FILTERS

    // 初始化(为其申请内存，设置大小)需要返回的filters(isl_union_set_list)
    amp_filters->filters = isl_union_set_list_alloc(ctx, 0);
    // 初始化左右分支filter的union_set_list（往回转换策略用）
    uset_list_left = isl_union_set_list_alloc(ctx, 0);
    uset_list_right = isl_union_set_list_alloc(ctx, 0);

    /**
     * @brief 如果第一个循环迭代空间外有不含循环约束的语句实例，则默认该语句只给左分支（右分支跳过 + 1）.且，这种情况下不能返回只有右分支的结果（因为左分支不为空）.
     *        如果是后面的迭代空间中间包含有不含循环约束的语句实例，则该语句同时包含进左右分支当中.
     * @note  skip代表左右分支拆分后插入到各自的isl_union_set_list（uset_list_left, uset_list_right）时需要跳过的大小.
     */
    int skip_left = 0;
    int skip_right = 0;
    isl_bool outside_of_first_cycle_flag = isl_bool_true;
    for (isl_size basic_set_list_dim = 0; basic_set_list_dim < basic_set_list_deepth; basic_set_list_dim++)
    {
        // isl_basic_set_list获取对应dim时的isl_basic_set
        bset = isl_basic_set_list_get_basic_set(basic_set_list, basic_set_list_dim);
        isl_size con_deepth = isl_basic_set_n_constraint(bset);
#ifdef DEBUG_AMP_GET_FILTERS
        printf("@DEBUG: \n       对应 basic_set_list_dim = %d 时的 isl_basic_set is: \n", basic_set_list_dim);
        isl_basic_set_dump(bset);
        printf("\n       对应约束的深度是： %d . \n\n", con_deepth);
#endif //  DEBUG_AMP_GET_FILTERS

        // 是否是第一个循环迭代空间外有不含循环约束的语句实例
        if (con_deepth < 1 && outside_of_first_cycle_flag)
        {
            printf("\n\n (con_deepth < 1 && outside_of_first_cycle_flag) is true!!\n\n");
            uset_list_left = isl_union_set_list_insert(uset_list_left, basic_set_list_dim, isl_union_set_from_basic_set(isl_basic_set_copy(bset)));
            skip_right++;
            continue;
        }
        // 是否是中间的循环迭代空间夹着的不含循环约束的语句实例
        else if (con_deepth < 1)
        {
            printf("\n\n (con_deepth < 1) is true!!\n\n");
            uset_list_left = isl_union_set_list_insert(uset_list_left, basic_set_list_dim - skip_left, isl_union_set_from_basic_set(isl_basic_set_copy(bset)));
            uset_list_right = isl_union_set_list_insert(uset_list_right, basic_set_list_dim - skip_right, isl_union_set_from_basic_set(isl_basic_set_copy(bset)));
            continue;
        }

        printf("\n\n no continue !\n\n");
        // 获取新的约束,目前默认的separate_dim(要拆分的循环的dim)是0
        amp_bset = amp_get_single_statement_constraints(ctx, bset, separate_dim, rate);
        // 此刻说明，第一个循环迭代空间已经拆分了
        if (outside_of_first_cycle_flag)
            outside_of_first_cycle_flag = isl_bool_false;
        if (!amp_bset)
        {
            printf("\n\033[31m@ERROR:\n       There are some errors because the get_single_statement_constraints function returned amp_bset is NULL , Now will return the original domain in filters way!!! \033[0m\n\n");
            goto error;
        }

        // 逐语句将amp_basic_set的左右分支的isl_basic_set合并成isl_union_set_list
        if (amp_bset->left_flag)
        {
            uset_list_left = isl_union_set_list_insert(uset_list_left, basic_set_list_dim - skip_left, isl_union_set_from_basic_set(amp_bset->left));
        }
        else
        {
            skip_left++;
        }
        if (amp_bset->right_flag)
        {
            uset_list_right = isl_union_set_list_insert(uset_list_right, basic_set_list_dim - skip_right, isl_union_set_from_basic_set(amp_bset->right));
        }
        else
        {
            skip_right++;
        }
    }
    // 获取左右分支union_set_list里面包含的union_set的数量
    isl_size left_size = isl_union_set_list_n_union_set(uset_list_left);
    isl_size right_size = isl_union_set_list_n_union_set(uset_list_right);

    // 对isl_union_set_list求并集,获得左右分支的filter(isl_union_set)
    if (left_size)
        uset_left_filter = isl_union_set_list_union(uset_list_left);
    if (right_size)
        uset_right_filter = isl_union_set_list_union(uset_list_right);

    // 左右分支的filter合并成amp_filters
    if (left_size)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_left_filter);
    if (right_size)
        amp_filters->filters = isl_union_set_list_add(amp_filters->filters, uset_right_filter);

    // 判断类型
    if (left_size && right_size)
        amp_filters->type = amp_filters_both;
    else if (left_size && !right_size)
        amp_filters->type = amp_filters_only_left;
    else if (!left_size && right_size)
        amp_filters->type = amp_filters_only_right;
    else
        goto error;

    // 释放内存
    isl_union_set_free(domain);
    isl_basic_set_list_free(basic_set_list);

    return amp_filters;
error:
    isl_union_set_free(domain);
    isl_basic_set_list_free(basic_set_list);
    isl_union_set_list_free(uset_list_left);
    isl_union_set_list_free(uset_list_right);
    isl_union_set_free(uset_left_filter);
    isl_union_set_free(uset_right_filter);
    isl_union_set_list_free(amp_filters->filters);
    amp_basic_set_free(amp_bset);

    return NULL;
}