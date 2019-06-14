#!/bin/bash

# This script provides the function _dials_autocomplete(), which is used by
# bash shells to allow command line parameter completion for selected DIALS
# commands
#
# The function contains a debugging routine, which is enabled by setting
# the environmental variable DIALS_AUTOCOMPLETE_DEBUG to any non-empty value.
# Information is written to the file completion.log in the current working
# directory. This can then be monitored in another terminal with e.g.
#   watch -n 0.1 cat completion.log

if [ -z "$_dials_autocomplete_path" ]; then
 _dials_autocomplete_path=$(libtbx.show_build_path)/dials/autocomplete/
fi

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

  [ -n "${DIALS_AUTOCOMPLETE_DEBUG:+1}" ] && {
    echo "1: $1"
    echo "COMP_WORDS: ${COMP_WORDS[*]}"
    echo "COMP_CWORD: $COMP_CWORD"
    echo "CUR: $cur"
    echo "PREV: $prev"
    echo "PPREV: $pprev"
  } > completion.log

  # Load and cache pre-computed hints for the requested command
  if [[ ${_dials_autocomplete_cache} != $1 ]]; then
   # Check if the requested command is supported. If not, use default completion.
   if [ ! -s ${_dials_autocomplete_path}$1 ]; then type compopt &>/dev/null && compopt -o default; return 0; fi

   source ${_dials_autocomplete_path}$1
   _dials_autocomplete_cache=$1
  fi

  if [[ ${cur} == "=" ]]; then
   # initial autocompletion of a choice parameter
   _dials_autocomplete_flags "${prev}"
   if [[ ${_dials_autocomplete_values} != "" ]] ; then
    COMPREPLY=( ${_dials_autocomplete_values} )

    [ -n "${DIALS_AUTOCOMPLETE_DEBUG:+1}" ] && {
     echo "Path is cur="
     echo "COMPREPLY: ${COMPREPLY[@]}"
    } >> completion.log

    return 0
   fi
  fi
  if [[ ${cur} == *= ]]; then
   COMPREPLY=()

   [ -n "${DIALS_AUTOCOMPLETE_DEBUG:+1}" ] && {
    echo "Path is cur*="
    echo "COMPREPLY: ${COMPREPLY[@]}"
   } >> completion.log

   return 0
  fi
  if [[ ${prev} == "=" ]]; then
   # autocompletion of a choice parameter with existing text
  _dials_autocomplete_flags "${pprev}"
   if [[ ${_dials_autocomplete_values} != "" ]] ; then
    COMPREPLY=( $(compgen -W "${_dials_autocomplete_values}" -- "${cur}") )
   fi

   [ -n "${DIALS_AUTOCOMPLETE_DEBUG:+1}" ] && {
    echo "Path is prev="
    echo "COMPREPLY: ${COMPREPLY[@]}"
   } >> completion.log

   return 0
  fi

  # Identify list of relevant phil parameters
  _dials_autocomplete_hints "${cur}"
  COMPREPLY=( $(compgen -W "${_dials_autocomplete_values}" -- "${cur}") \
              $(compgen -f -X "!*.json" -- "${cur}") \
              $(compgen -f -X "!*.expt" -- "${cur}") \
              $(compgen -f -X "!*.pickle" -- "${cur}") \
              $(compgen -f -X "!*.refl" -- "${cur}") \
              $(compgen -f -X "!*.mpack" -- "${cur}") \
              $(compgen -f -X "!*.phil" -- "${cur}") \
              $(compgen -d -S"/" -- "${cur}" ) )
  unset -f _dials_autocomplete_values
  type compopt &>/dev/null && compopt -o nospace
  if [[ ${#COMPREPLY[@]} == 1 ]]; then
   # If there's only one option, check if it is expandable
   _dials_autocomplete_expansion "${COMPREPLY[0]}"
   if [[ ${_dials_autocomplete_values} != "" ]] ; then
     COMPREPLY=( ${_dials_autocomplete_values} )
   fi
   # If the only option is not ending in '.', '=' or '/', then append a space
   if [[ ${COMPREPLY[0]} != *. && ${COMPREPLY[0]} != *= && ${COMPREPLY[0]} != */ ]] ; then
    type compopt &>/dev/null && compopt +o nospace || COMPREPLY=( "${COMPREPLY[0]} " )
   fi
  fi

  [ -n "${DIALS_AUTOCOMPLETE_DEBUG:+1}" ] && {
   echo "Path is other"
   echo "COMPREPLY: ${COMPREPLY[@]}"
  } >> completion.log

  return 0
}
