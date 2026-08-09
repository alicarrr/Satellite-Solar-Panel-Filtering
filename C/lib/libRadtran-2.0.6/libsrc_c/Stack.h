/************************************************************************
 * $Id$
 *
 * Corefinder - Identify veiled core; authored by Paul Ockenfuss
 *
 * Copyright (c) 2000-2019 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Paul Ockenfuss
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

#ifndef STACK_H
#define STACK_H

#include <stdlib.h>

typedef struct point_t {
  int   x;
  int   y;
  int   z;
  float cost;
} point_t;

typedef struct my_stack_t {
  point_t* points;
  int      capacity;
  int      size;
} my_stack_t;

size_t     dist2 (const point_t p1, const point_t p2);
my_stack_t create_stack (int capacity);
void       destroy_stack (my_stack_t* stack);
void       free_stack (my_stack_t* stack);
void       stack_extend (my_stack_t* stack);
void       stack_swap (my_stack_t* stack, size_t i, size_t j);
size_t     get_child (size_t i, size_t child);
size_t     get_parent (size_t i);
void       sort_down (my_stack_t* stack, size_t pos);
void       sort_up (my_stack_t* stack);
void       stack_push (my_stack_t* stack, point_t point);
point_t    stack_pop (my_stack_t* stack);
int        stack_empty (const my_stack_t* stack);

#endif /*STACK_H*/
