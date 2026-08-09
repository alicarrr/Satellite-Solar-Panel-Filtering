#pragma once

int options_stash_get (char** stash);     // return reference to a list of options as they were parsed through the lex system
int options_stash_append (char* yytext);  // add entry to options stash
