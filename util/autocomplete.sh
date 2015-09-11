#!/bin/bash
if [ -z "$_dials_autocomplete_path" ]; then
 _dials_autocomplete_path=$(libtbx.show_build_path)/dials/autocomplete/
fi

# Predeclare associative arrays, so they can be cached between invocations
# Use global declaration if possible (requires bash 4.2+), local otherwise (bash 4.0+)
# If associative arrays not supported (bash <4), do not declare
declare -gA _dials_autocomplete_flags 2>/dev/null || declare -A _dials_autocomplete_flags 2>/dev/null || :
declare -gA _dials_autocomplete_expansion 2>/dev/null || declare -A _dials_autocomplete_expansion 2>/dev/null || :

# Define two helper functions from bash_completion
# Return 1 if $1 appears to contain a redirection operator.  Handles backslash
# quoting (barely).
_redir_op()
{
  case "$1" in
  *\\'[\<\>]'*) return 1;;
  *[\<\>]*)     return 0;;
  *)            return 1;;
  esac
}

# _redir_test tests the current word ($1) and the previous word ($2) for
# redirection operators and does filename completion if either one contains
# a redirection operator
_redir_test()
{
  if _redir_op "$1" ; then
    COMPREPLY=( $( compgen -f "" ) )
    return 0
  elif _redir_op "$2" ; then
    COMPREPLY=( $( compgen -f "$1" ) )
    return 0
  fi
  return 1
}

type compopt &>/dev/null || {
 function compopt ()
 {
  # unsupported on Mac. Quick fix for now. Proper fix: TODO
  :
 }
}

function _dials_autocomplete ()
{
  # This function provides autocomplete functionality to supported dials commands
  COMPREPLY=()

  local cur prev pprev
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  # Do not attempt autocompletion if user is currently redirecting stdio
  _redir_test "$cur" "$prev" && return 0;

  pprev=""
  if [ "${COMP_CWORD}" -gt "2" ] ; then
    pprev="${COMP_WORDS[COMP_CWORD-2]}"
  fi

#  {
#   echo "1: $1"
#   echo "COMP_WORDS: ${COMP_WORDS[*]}"
#   echo "COMP_CWORD: $COMP_CWORD"
#   echo "CUR: $cur"
#   echo "PREV: $prev"
#   echo "PPREV: $pprev"
#  } > completion.log

  # Load and cache pre-computed hints for the requested command
  if [[ ${_dials_autocomplete_cache} != $1 ]]; then
   # Check if the requested command is supported. If not, use default completion.
   if [ ! -f ${_dials_autocomplete_path}$1 ]; then compopt -o default; return 0; fi

   source ${_dials_autocomplete_path}$1
   _dials_autocomplete_cache=$1
  fi

  if [[ ${cur} == "=" ]]; then
   # initial autocompletion of a choice parameter
   declare -p _dials_autocomplete_flags >/dev/null 2>&1 && \
   if [ ${_dials_autocomplete_flags[${prev}]+exists} ]; then
    COMPREPLY=( ${_dials_autocomplete_flags[${prev}]} )
    return 0
   fi
  fi
  if [[ ${prev} == "=" ]]; then
   # autocompletion of a choice parameter with existing text
   declare -p _dials_autocomplete_flags >/dev/null 2>&1 && \
   if [ ${_dials_autocomplete_flags[${pprev}]+exists} ]; then
    COMPREPLY=( $(compgen -W "${_dials_autocomplete_flags[${pprev}]}" -- "${cur}") )
   fi
   return 0
  fi

  # Identify list of relevant phil parameters
  _dials_autocomplete_hints "${cur}"
  COMPREPLY=( $(compgen -W "${_dials_autocomplete_values}" -- "${cur}") \
              $(compgen -f -X "!*.json" -- "${cur}") \
              $(compgen -f -X "!*.pickle" -- "${cur}") \
              $(compgen -f -X "!*.phil" -- "${cur}") )
  unset -f _dials_autocomplete_values
  compopt -o nospace
  if [[ ${#COMPREPLY[@]} == 1 ]]; then
   # If there's only one option, check if it is expandable
   declare -p _dials_autocomplete_expansion >/dev/null 2>&1 && \
   if [ ${_dials_autocomplete_expansion[${COMPREPLY[0]}]+exists} ] ; then
     COMPREPLY=( ${_dials_autocomplete_expansion[${COMPREPLY[0]}]} )
   fi
   # If the only option is not ending in '.' or '=', then append a space
   if [[ ${COMPREPLY[0]} != *. && ${COMPREPLY[0]} != *= ]] ; then
    compopt +o nospace
   fi
  fi
  return 0
}
