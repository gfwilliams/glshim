#include "gl.h"

#ifndef GL_STATE_H
#define GL_STATE_H

#include "eval.h"
#include "texture.h"

typedef struct {
    GLboolean line_stipple,
              blend,
              color_array,
              normal_array,
              tex_coord_array,
              texgen_s,
              texgen_t,
              texture_2d,
              vertex_array;
} EnableState;


typedef struct {
    GLenum S;
    GLenum T;
    GLfloat Tv[4];
    GLfloat Sv[4];
} TexGenState;

typedef struct {
    GLuint unpack_row_length,
           unpack_skip_pixels,
           unpack_skip_rows;
    // TODO: do we only need to worry about GL_TEXTURE_2D?
    GLboolean rect_arb;
    GLtexture *bound;
    khash_t(tex) *list;
} TextureState;


typedef struct {
    GLint size;
    GLenum type;
    GLsizei stride;
    const GLvoid *pointer;
} PointerState;

typedef struct {
    PointerState vertex, color, normal, tex_coord;
} PointerStates;


typedef struct {
    RenderList *active;
    GLboolean compiling;
    GLuint base;
    GLuint name;
    GLenum mode;
} DisplayListState;


typedef struct {
    MapState *vertex3,
             *vertex4,
             *index,
             *color4,
             *normal,
             *texture1,
             *texture2,
             *texture3,
             *texture4;
} MapStates;


typedef struct {
    EnableState enable;
    PointerStates pointers;
    DisplayListState list;
    MapStates map1, map2;
    MapState *map_grid;
    TextureState texture;
    TexGenState texgen;
    GLfloat color[4];
} GLstate;

#endif
