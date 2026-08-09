#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "errors.h"

// we store all input options going through the lex parser
static char* options_stash = NULL;  // singleton memory to store input options

static size_t stash_size = 128 * 132;  // initial or current stash size, will grow dynamically if more space is needed

static int alloc_options_stash (size_t min_size) {
  if (!options_stash) {  // if there is no options stash yet, just allocate
    options_stash = malloc (stash_size * sizeof (char));
    CHKPOINTEROUT (options_stash, "could not allocate memory for options_stash");
    memset (options_stash, '\0', stash_size * sizeof (char));
  } else {  // otherwise we grow
    size_t new_stash_size = stash_size * 2 + min_size;
    options_stash         = realloc (options_stash, new_stash_size * sizeof (char));
    CHKPOINTEROUT (options_stash, "could not grow memory for options_stash");
    stash_size = new_stash_size;
  }
  return 0;
}

int options_stash_append (char* yytext) {
  int status;
  if (!options_stash) {
    status = alloc_options_stash (stash_size);
    CHKERR (status);
  }

  const size_t current_stash_len = strlen (options_stash);
  const size_t yytext_len        = strlen (yytext);
  int          maxlen            = stash_size - current_stash_len - yytext_len - 1;
  if (maxlen < 0) {  // not enough space to append option
    status = alloc_options_stash (stash_size + yytext_len + 1);
    CHKERR (status);
    options_stash_append (yytext);
    return 0;
  }
  strcat (options_stash, yytext);
  strcat (options_stash, "\n");
  //printf("options_append: %s -> ... \n %s\n", yytext, options_stash);
  return 0;
};

// return reference to a list of options as they were parsed through the lex system
int options_stash_get (char** stash) {
  if (!options_stash) {
    int status = alloc_options_stash (stash_size);
    CHKERR (status);
  }
  *stash = options_stash;
  return 0;
}
