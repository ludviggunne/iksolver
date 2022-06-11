/* TODO */
/* Implement constraints */
/* Description */

#ifndef IKSOLVER_H
#define IKSOLVER_H

#define IK_ERROR 0
#define IK_OK    1

/**********************************************
 *                 STRUCTS                    *
 **********************************************/

/*
 * Vector type for representing joint positions
 */
typedef struct { float x, y; } ik_vec2;



/*
 * Struct defining joint, including length of segment 
 *  connecting this joint to it's parent. 
 * Attaching less than 'n_children' before any operation 
 *  on the tree is undefined behavior.
 */
typedef struct ik_joint {

    ik_vec2 position;
    float length;

    int n_children;

    struct ik_joint *parent;
    struct ik_joint *children[0];

} ik_joint;



/*
 * Provides dynamically sized buffer of vertex 
 *  positions for rendering tree as line segments.
 */
typedef struct {
    ik_vec2 *data;
    int size;
    int cap;
} ik_vertex_buffer;



/**********************************************
 *                 FUNCTIONS                  *
 **********************************************/

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Sets up global state required for other functions
 */
void ik_init(void);



/*
 * Creates a new joint with capacity to hold 'n_children'
 *  attached children.
 */
ik_joint *ik_new_joint(float length, int n_children);



/* 
 * Deletes branch beginning at root
 */
void ik_delete_branch(ik_joint *root);



/*
 * Attaches 'child' to 'parent'. 
 * Returns IK_ERROR if maximum children is
 *  reached for 'parent'.
 */
int ik_attach_joint(ik_joint *child, ik_joint *parent);



/*
 * Translates tree by setting root position to (x, y)
 */
void ik_translate(ik_joint *root, float x, float y);



/*
 * Solves IK using FABRIK model.
 */
int ik_solve(ik_joint *effected, float target_x, float target_y);



/*
 * Retrieves vertex positions for rendering tree.
 */
void ik_get_render_data(ik_joint *root, ik_vertex_buffer *buffer);



/*
 * Creates a new vertex buffer.
 */
ik_vertex_buffer ik_new_vertex_buffer(void);



/*
 * Frees buffer data, setting buffer to invalid state.
 */
void ik_free_vertex_buffer(ik_vertex_buffer *buffer);



#ifdef __cplusplus
}
#endif

#endif /* IKSOLVER_H */


#ifdef IKSOLVER_IMPLEMENTATION

/**********************************************
 *                 UTILITY                    *
 **********************************************/

#ifndef IK_MALLOC
# include <stdlib.h>
# define IK_MALLOC malloc
# define IK_FREE   free
#endif

#ifndef IK_SQRT
# include <math.h>
# define IK_SQRT sqrtf
#endif

#ifndef IK_MEMCPY
# include <string.h>
# define IK_MEMCPY memcpy
#endif

#ifndef NULL
# define NULL ((void*)0)
#endif

#ifdef IK_DEBUG
# include <stdio.h>
# define LOG(msg, ...) printf(msg "\n", __VA_ARGS__)
# define INDENT() indent++
# define UNINDENT() indent--
static int indent = 0;
#else
# define LOG(msg, ...)
# define INDENT()
# define UNINDENT()
#endif


/**********************************************
 *             INTERNAL FUNCTIONS             *
 **********************************************/

/*
 * Calculates length of vector (x, y)
 */
static inline float length(float x, float y)
{
    return IK_SQRT(x * x + y * y);
}





/*
 * Stack for pushing ik_joint pointers to during back reach, so that 
 *  we know which path to take when reaching forward. 
 */

#define IK_STACK_SIZE 1024 /* TODO: allow user to define this (?) */
typedef ik_joint **ik_stack_ptr;
struct {
    ik_stack_ptr end;
    ik_joint *data[IK_STACK_SIZE];
} ik_stack;

static inline void ik_stack_push(ik_joint *joint)
{
    /* TODO: bounds-checking (?) */
    LOG("Pushing joint %p", joint);
    *(ik_stack.end++) = joint;
}

static inline ik_joint *ik_stack_pop(void)
{
    if(ik_stack.end == ik_stack.data)
        return NULL;

    ik_joint *joint = *(--ik_stack.end); 
    LOG("Popping joint %p", joint);
    return joint;
}

static inline ik_joint *ik_stack_top(void)
{
    if(ik_stack.end == ik_stack.data)
        return NULL;

    return *(ik_stack.end - 1);
}



/*
 * Struct to represent entries of rotation matrix.
 *
 *  C: cos(abs(angle))
 *  S: sin(abs(angle))
 *  P: sgn(angle)
 *
 */
struct ik_matrix { float C, S, P; };



/* 
 * Aligns branch according to precalculated values.
 *
 *  root: root joint of branch to be aligned
 *  pivot: position to rotate around
 *  mat: entries of rotation matrix
 *  offset: offset to translate
 * 
 */
static void ik_align_branch_precalc(ik_joint *root, ik_vec2 pivot, struct ik_matrix mat, ik_vec2 offset)
{
    LOG("Aligning branch %p", root);

    /* Translate */
    root->position.x += offset.x;
    root->position.y += offset.y;

    /* Rotate */
    root->position.x -= pivot.x;
    root->position.y -= pivot.y;

    float tmp_x = root->position.x;

    /* P is applied using identies cos(-a) = cos(a), sin(-a) = -sin(a) */
    root->position.x = (mat.C        ) * root->position.x + (-mat.P * mat.S) * root->position.y;
    root->position.y = (mat.P * mat.S) * tmp_x            + ( mat.C        ) * root->position.y;

    root->position.x += pivot.x;
    root->position.y += pivot.y;

    /* Recurse */
    /* TODO: Use stack instead of recursion? */
    for(int i = 0; i < root->n_children; i++)
        ik_align_branch_precalc(root->children[i], pivot, mat, offset);
}



/* 
 * This function is used for translating/rotating branches that are not
 *  affected by back or forward reach.
 * 
 *  root: root joint of branch to be aligned
 *  from: previous position of root's parent, relative to it's parent
 *  to: new position of root's parent, relative to it's parent
 * 
 */
static void ik_align_branch(ik_joint *root, ik_vec2 from, ik_vec2 to)
{
    /* Find rotation matrix entries */
    struct ik_matrix mat;

    /* C = cos(angle) */
    float C_denom = length(from.x, from.y) * length(to.x, to.y);
    mat.C = (from.x * to.x + from.y * to.y) / C_denom;
    
    /* S = sin(angle) */
    float S2 = 1 - mat.C * mat.C;
    mat.S = IK_SQRT(S2 > 0.0f ? S2 : 0.0f);

    /* P = sgn(angle) */
    mat.P = from.x * to.y - from.y * to.x > 0.0f ? 1.0f : -1.0f;


    /* Find offset */
    ik_vec2 offset;
    offset.x = to.x - from.x;
    offset.y = to.y - from.y;

    /* Use parents position as pivot, as rotation takes place after translation */
    LOG("Aligning branch %p, C = %f, S = %f, P = %f, offset = (%f, %f)",
            root, mat.C, mat.S, mat.P, offset.x, offset.y);
    ik_align_branch_precalc(root, root->parent->position, mat, offset);
}



/*
 * Translates branch by offset.
 * Used to align branch who's parent is the tree root.
 */
static void ik_align_branch_only_translate(ik_joint *root, float offset_x, float offset_y)
{
    LOG("Translating branch %p", root);
    root->position.x += offset_x;
    root->position.y += offset_y;

    for(int i = 0; i < root->n_children; i++)
        ik_align_branch_only_translate(root->children[i], offset_x, offset_y);
}



/*
 * Moves joint within distance of target.
 */
static inline void ik_move_within_dist(ik_joint *joint, float distance, float target_x, float target_y)
{
    LOG("Moving joint %p within distance %f of (%f, %f)", joint, distance, target_x, target_y);
    float dx = joint->position.x - target_x;
    float dy = joint->position.y - target_y;
    float norm_denom = length(dx, dy);

    if(norm_denom == 0.0f)
    {
        joint->position.x = target_x;
        joint->position.y = target_y;
    } else {
        joint->position.x = target_x + distance * dx / norm_denom;
        joint->position.y = target_y + distance * dy / norm_denom;
    }
}



/*
 * This functions iterates backwards, and must be followed immediatly
 *  by ik_reach_forward, as joint pointers are push to the stack.
 * Distance must be zero when function is called from outside itself.
 * 
 * Since this function traverses tree to the root, the three last arguments provide
 *  access to the root joint, and its original position, which is used for forward reach.
 * 
 */
static void ik_reach_back(ik_joint *effected, float target_x, float target_y, float distance,
                          ik_joint **root, float *root_org_x, float *root_org_y)
{
    LOG("Reach back from joint %p, distance %f", effected, distance);
    float org_x = effected->position.x;
    float org_y = effected->position.y;

    /* If reached root, save original position */
    if(!effected->parent)
    {
        LOG("%s", "Is root");
        *root = effected;
        *root_org_x = org_x;
        *root_org_y = org_y;
    }


    /* Move effected towards target */
    ik_move_within_dist(effected, distance, target_x, target_y);

    /* If effected joints have more than one child, the other */
    /*  children's branches must be aligned according to new  */
    /*  orientation.                                          */
    if(effected->n_children > 1) {

        LOG("Has %d children", effected->n_children);
        ik_joint *path_child = ik_stack_top();

        if(effected->parent)
        {
            LOG("%s", "Aligning child branches");
            ik_vec2 from, to;

            from.x = org_x - effected->parent->position.x;
            from.y = org_y - effected->parent->position.y;

            to.x = effected->position.x - effected->parent->position.x;
            to.y = effected->position.y - effected->parent->position.y;
            
            for(int i = 0; i < effected->n_children; i++)
            {
                ik_joint *child = effected->children[i];

                if(child != path_child)
                    ik_align_branch(child, from, to);
            }
        } else {
            /* Root of whole tree -> no parent to define orientation */
            /*  -> only translate                                    */
            for(int i = 0; i < effected->n_children; i++)
            {
                LOG("%s", "Translating child branches");
                ik_joint *child = effected->children[i];

                if(child != path_child)
                    ik_align_branch_only_translate(
                        child, 
                        effected->position.x - org_x, 
                        effected->position.y - org_y);
            }
        }
    }
    


    /* Terminate if we've reached root */
    if(!effected->parent)
        return;


    /* If 'effected' has siblings, push 'effected' to stack     */
    /*  so that we know which path to take during forward reach */
    if(effected->parent->n_children > 1) {
        LOG("%s", "Push joint to path");
        ik_stack_push(effected);
    }

    /* Recurse */
    /* TODO: Use stack ? */
    ik_reach_back(
        effected->parent, 
        effected->position.x,
        effected->position.y,
        effected->length,
        root, root_org_x, root_org_y);
}



static void ik_reach_forward(ik_joint *root, float distance, float target_x, float target_y)
{
    LOG("Reach forward from joint %p, distance %f", root, distance);
    float org_x = root->position.x;
    float org_y = root->position.y;


    /* Move root towards target */
    ik_move_within_dist(root, distance, target_x, target_y);


    /* Terminate if we've reached leaf */
    if(root->n_children == 0)
        return;


    /* Reach forward through this child */
    ik_joint *path_child;

    /* If 'root' has more than 1 child, all but path_child       */
    /*  must be aligned according to new orientation, provided   */
    /*  that 'root' has a parent to define this orientation, ie. */
    /*  'root' is not the root of the whole tree.                */
    if(root->n_children > 1) {

        LOG("Has %d children", root->n_children);
        path_child = ik_stack_pop();

        if(root->parent)
        {   
            LOG("%s", "Aligning child branches");
            ik_vec2 from, to;

            from.x = org_x - root->parent->position.x;
            from.y = org_y - root->parent->position.y;

            to.x = root->position.x - root->parent->position.x;
            to.y = root->position.y - root->parent->position.y;
            
            for(int i = 0; i < root->n_children; i++)
            {
                ik_joint *child = root->children[i];

                if(child != path_child)
                    ik_align_branch(child, from, to);
            }
        } else {
            /* Root of whole tree -> no parent to define orientation */
            /*  -> only translate                                    */
            LOG("%s", "Translating child branches");
            for(int i = 0; i < root->n_children; i++)
            {
                ik_joint *child = root->children[i];

                if(child != path_child)
                    ik_align_branch_only_translate(
                        child, 
                        root->position.x - org_x, 
                        root->position.y - org_y);
            }
        }
    } else {
        path_child = root->children[0];
    }
    

    /* Recurse */
    /* TODO: Use stack ? */
    ik_reach_forward(
        path_child,
        path_child->length,
        root->position.x,
        root->position.y
    );
}



/*
 * Pushes vector to buffer.
 */
static void ik_vertex_buffer_push(ik_vertex_buffer *buffer, ik_vec2 value)
{
    if(buffer->size == buffer->cap)
    {
        int new_cap = buffer->cap * 2;
        ik_vec2 *new_data = IK_MALLOC(sizeof(ik_vec2) * new_cap);

        memcpy(new_data, buffer->data, sizeof(ik_vec2) * buffer->size);
        free(buffer->data);

        buffer->data = new_data;
        buffer->cap = new_cap;
    }

    buffer->data[buffer->size++] = value;
}


/*
 * Resets buffer data
 */
static void ik_vertex_buffer_reset(ik_vertex_buffer *buffer) 
{
    buffer->size = 0;
}


/*
 * Gets vertex data without reseting buffer, used recursively
 *  from ik_get_render_data.
 */
static void ik_get_render_data_no_reset(ik_joint *root, ik_vertex_buffer *buffer)
{
    for(int i = 0; i < root->n_children; i++)
    {
        ik_joint *child = root->children[i];
        ik_vertex_buffer_push(buffer, root->position);
        ik_vertex_buffer_push(buffer, child->position);
        ik_get_render_data_no_reset(child, buffer);
    }
}


/*
 * Translates branch by vector (dx, dy)
 */
void ik_translate_relative(ik_joint *root, float dx, float dy)
{
    root->position.x += dx;
    root->position.y += dy;
    for(int i = 0; i < root->n_children; i++)
        ik_translate_relative(root->children[i], dx, dy);
}



/**********************************************
 *            INTERFACE FUNCTIONS             *
 **********************************************/

void ik_init(void)
{
    ik_stack.end = ik_stack.data;
}

ik_joint *ik_new_joint(float length, int n_children)
{
    ik_joint *joint = IK_MALLOC(sizeof(ik_joint) + sizeof(ik_joint*) * n_children);

    joint->position.x = 0.0f;
    joint->position.y = 0.0f;
    joint->length = length;
    joint->parent = NULL;
    joint->n_children = n_children;
    for(int i = 0; i < n_children; i++)
        joint->children[i] = NULL;
    
    return joint;
}


void ik_delete_branch(ik_joint *root)
{
    for(int i = 0; i < root->n_children; i++)
        if(root->children[i])
            ik_delete_branch(root->children[i]);

    IK_FREE(root);
}


int ik_attach_joint(ik_joint *child, ik_joint *parent)
{
    for(int i = 0; i < parent->n_children; i++)
        if(parent->children[i] == NULL) {
            parent->children[i] = child;
            child->parent = parent;
            return IK_OK;
        }
    return IK_ERROR;
}


void ik_translate(ik_joint *root, float x, float y)
{
    float dx = x - root->position.x;
    float dy = y - root->position.y;
    ik_translate_relative(root, dx, dy);
}


int ik_solve(ik_joint *effected, float target_x, float target_y)
{
    LOG("%s", "\n *** SOLVE BEGIN ***\n");
    ik_joint *root;
    float root_org_x, root_org_y;

    ik_reach_back(
        effected, 
        target_x, 
        target_y, 
        0.0f, 
        &root, 
        &root_org_x, 
        &root_org_y
    );

    ik_reach_forward(
        root, 
        0.0f, 
        root_org_x, 
        root_org_y
    );

    LOG("%s", "\n *** SOLVE END ***\n");
    return 1;
}


void ik_get_render_data(ik_joint *root, ik_vertex_buffer *buffer)
{
    ik_vertex_buffer_reset(buffer);
    ik_get_render_data_no_reset(root, buffer);
}


ik_vertex_buffer ik_new_vertex_buffer(void)
{
    ik_vertex_buffer buffer;
    buffer.cap = 10;
    buffer.size = 0;
    buffer.data = IK_MALLOC(sizeof(ik_vec2) * buffer.cap);

    return buffer;
}


void ik_free_vertex_buffer(ik_vertex_buffer *buffer)
{
    IK_FREE(buffer->data);
    buffer->cap = 0;
    buffer->size = 0;
}

#endif /* IKSOLVER_IMPLEMENTATION */