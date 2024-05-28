#!/usr/bin/env bash

# BATS helper functions.
# To use inside a <pkg>/test/test_XXX.bats do:
#
# bats_load_library wct-bats.sh
#
# See test/docs/bats.org for more info.
#
# Note, these generally exit on error.

shopt -s extglob
shopt -s nullglob

# User may set these to control logging.
#
# Emit log message only if greater or equal to the level.
#
# The following are the interpretations of levels:
#
# - debug :: possibly verbose informaiton for debugging tests themselves
#
# - info :: default, for indicating normal events.
# 
# - warn :: for aberrant events but that do not impede nor negate further running.
#
# - error :: for events which may invalidate test results but testing may continue.
#
# - fatal :: events for which the testing can not continue.
WCT_BATS_LOG_LEVEL=${WCT_BATS_LOG_LEVEL:-info}      

# User may determine where logging goes.
#
# - "capture" is default and allows BATS to capture.
# - "terminal" will bypass BATS capturing.
# - anything else is interpreted as a file to which messages are concatenated.
WCT_BATS_LOG_SINK=${WCT_BATS_LOG_SINK:-capture}

# yell <message>
#
# Echo message to terminal, escaping BATS controlled logging.
function yell () {
    echo "$@" 1>&3
}

# index_of <entry> <element0> <element1> ...
#
# Return array index of entry in one or more elements
function index_of () {
    local entry="$1" ; shift
    local index=0
    for one in "$@"
    do
        if [ "$one" = "$entry" ] ; then
            echo $index
            return
        fi
        index=$(( index + 1 ))
    done
}

declare -g -a log_levels=( debug info warn error fatal )
declare -g log_level
log_level=$(index_of "$WCT_BATS_LOG_LEVEL" "${log_levels[@]}")

# log_ilevel <level>
#
# Emit the numerical value of a log level
function log_ilevel () {
    local lvl=$1 ; shift
    # yell "lvl=$lvl log_level=$log_level WCT_BATS_LOG_LEVEL=$WCT_BATS_LOG_LEVEL"
    # yell "log_levels=${log_levels[@]}"
    index_of "$lvl" "${log_levels[@]}"
}

# log_at_level <level> <message>
#
# Log at given level
function log_at_level () {
    local ilvl lvl="$1" ; shift
    ilvl=$(log_ilevel "$lvl")
    if [ "$ilvl" -lt "$log_level" ] ; then
        return
    fi
    local pre L=${lvl:0:1}
    pre="$(date +"%Y-%m-%d %H:%M:%S.%N") [ ${L^^} ]"

    if [ "$WCT_BATS_LOG_SINK" = capture ] ; then
        echo -e "$pre $*" 1>&2  # stderr
    elif [ "$WCT_BATS_LOG_SINK" = terminal ] ; then
        echo -e "$pre $*" 1>&3  # special avoid capture
    else
        echo -e "$pre $*" >> "$WCT_BATS_LOG_SINK"
    fi
}

# debug <message>
#
# Log at "debug" level
function debug () {
    log_at_level debug "$@"
}


# log <message>
#
# Log at "info" level.
#
# Note, this shadows the GNU info program.
function info () {
    log_at_level info "$@"
}

# warn <message>
#
# Log at warn level
function warn () {
    log_at_level warn "$@"
}

# error <message>
#
# Log at error level
function error () {
    log_at_level error "$@"
}

# fatal <message>
#
# Log at fatal level
function fatal () {
    log_at_level fatal "$@"
}

# die <message>
#
# Log at fatal error and exit
function die () {
    yell "$@"
    info "$@"
    exit 1
}

# file_larger_than <filename> <nbytes>
#
# Assert named file is larger than number of bytes.
function file_larger_than () {
    local filename=$1 ; shift
    local minsize=$1; shift
    echo "check file size for $filename"
    [[ -f "$filename" ]]
    local fsize
    fsize="$(stat -c '%s' "$filename")"
    echo "minsize is $minsize. $filename is size $fsize"
    [[ -n "$minsize" ]]
    [[ "$fsize" -gt "$minsize" ]]
}


# check <command line>
#
# Use BATS run but add some standard boilerplate.
#
# Logs to debug the command, logs $output to info and assert zero
# return status code.
function check () {
    info "RUNNING: $*"
    local status
    run "$@"

    if [ -n "$output" ] ; then
        info "OUTPUT:\n$output"
    else
        info "OUTPUT: (none)"
    fi
    [[ "$status" -eq 0 ]]
}

# run_idempotently [options] -- <command line>
#
# Run command if any target file is missing or older than any source file.
#
# Options:
# -s, --source <filename>
#   Name a source file for idempotency checking (multiple okay).
# -t, --target <filename>
#   Name a target file for idempotency checking (multiple okay).
# --
#   Separate options from command line to run.
# <command line>
#   The command line to run.
function run_idempotently () {

    declare -a src
    declare -a tgt
    while [[ $# -gt 0 ]] ; do
        case $1 in
            -s|--source) src+=( "$2" ) ; shift 2;;
            -t|--target) tgt+=( "$2" ) ; shift 2;;
            --) shift; break;;
            *) die "unknown argument: $1"
        esac
    done

    local need_to_run="no"
    if [ ${#src[@]} -eq 0 ] || [ ${#tgt[@]} -eq 0 ] ; then
        # Always run if sink, source or atomic.
        need_to_run="yes"
        debug "running because source (${src[*]}) or sink (${tgt[*]}) are empty"
    fi
    if [ "$need_to_run" = "no" ] ; then
        # Run if missing any targets
        for one in "${tgt[@]}"  
        do
            if [ ! -f "$one" ] ; then
                need_to_run="yes"
                debug "running because of missing target: $one"
            fi
        done
    fi
    if [ "$need_to_run" = "no" ] ; then
        # Run if any source is newer than any target
        local src_newest tgt_oldest
        # shellcheck disable=SC2012
        src_newest="$(ls -t "${src[@]}" | head -1)"
        # shellcheck disable=SC2012
        tgt_oldest="$(ls -t "${tgt[@]}" | tail -1)"

        if [ "$src_newest" -nt "$tgt_oldest" ] ; then
            need_to_run="yes"
            debug "running because newest src \"$src_newest\" newer than oldest tgt \"$tgt_oldest\""
        fi
    fi
    
    if [ "$need_to_run" = "no" ] ; then
        warn "idempotent: $*"
        return
    fi

    check "$@"
}

# dumpenv
#
# Echo shell environment settings.
function dumpenv () {
    env | grep "=" | sort
    # warn $BATS_TEST_FILENAME
    # warn $BATS_TEST_SOURCE    
}

# tojson
#
# A filter reading stdin for lines like:
#
#  key1=value1
#  key2=value2
#
# and will emit JSON text for the equivalent JSON object.
#
# This requires the "jq" program.
function tojson () {
    declare -a args=("-n")
    while IFS='=' read -r key val
    do
        args+=("--arg" "$key" "$val")
    done < <(printf '%s\n' "$@")
    args+=( "\$ARGS.named" )
    jq "${args[@]}"
}

# topdir
#
# Emit the top level source directory.
#
# Requries this file (wct-bats.sh) to be in its place in the source.
function topdir () {
    # yell "BASH_SOURCE=${BASH_SOURCE[*]}"
    # yell "FUNCNAME=${FUNCNAME[*]}"

    local name="${BASH_SOURCE[0]}"
    [[ -n "$name" ]]
    name="$(realpath "$name")"
    [[ -n "$name" ]]
    name="$(dirname "$name")"
    [[ -n "$name" ]]
    name="$(dirname "$name")"
    [[ -n "$name" ]]
    echo "$name"
}    

# blddir
#
# Emit the biuld directory.
function blddir () {
    echo "$(topdir)/build"
}

# srcdir <pkg>
# 
# Emit the sub-package source directory of hte given pkg name.
function srcdir () {
    local maybe path="$1" ; shift
    maybe="$(topdir)/$path"
    if [ -d "$maybe" ] ; then
        echo "$maybe"
    fi
}


# saveout [options] <src> ...
#
# Save a file out of the test area to the build/tests/ area.
#
# options:
# -c, --category <cat>
#   Specify the category in which to save, default is "output".
#
# -t, --target <tgt>
#   Specify an explicit target file name to save.  Only the first
#   filename is saved out and <tgt> must be a relative path.
#
# <src>
#   The name of a file to save out.
#
# Files are saved to build/tests/<category>/<version>/<test-name>/
#
#   <category> is as given or defaults to "output"
#   <version> is the current wire-cell version
#   <test-name> is from BATS_TEST_FILENAME without the extension.
function saveout () {

    declare -a src=()
    local tgt subdir="output"

    while [[ $# -gt 0 ]] ; do
        case $1 in
            -t|--target) tgt="$2"; shift 2;;
            -c|--category) subdir="$2"; shift 2;;
            -*) die "Unknown option $1";;
            *) src+=("$1"); shift;;
        esac
    done
          
    if [ "${#src[@]}" -eq - ] ; then
        die "saveout: no source given"
    fi

    local name base
    name="$(basename "$BATS_TEST_FILENAME" .bats)"
    base="$(blddir)/tests/$subdir/$(version)/${name}"

    # single, directed target
    if [ -n "$tgt" ] ; then
        tgt="${base}/${tgt}"
        mkdir -p "$(dirname "$tgt")"
        cp "${src[0]}" "$tgt"
        return
    fi

    mkdir -p "${base}"
    # multi, implicit target
    for one in "${src[@]}"
    do
        local tgt
        tgt="${base}/$(basename "$one")"
        mkdir -p "$(dirname "$tgt")"
        cp "$one" "$tgt"
    done
}

# skip_if_no_input <path> [...]
#
# Skip a test if the given path does not exist in the input area.
function skip_if_no_input () {
    declare -a paths=("$@")
    if [ "${#paths[@]}" -eq 0 ] ; then
        paths=( "." )
    fi
    for path in "${paths[@]}"
    do
        local input
        input="$(blddir)/tests/input/$path"
        if [ -f "$input" ] || [ -d "$input" ] ; then
            continue;
        fi
        skip "no input test data at $input"
    done
}

# input_file <relpath> ...
#
# Emit full paths given relative paths in input category
function input_file () {
    local input
    input="$(blddir)/tests/input"
    if [ ! -d "$input" ] ; then
        die "no input data"
    fi

    for path in "$@"
    do
        full="$input/$path"
        if [ -f "$full" ] ; then
            echo "$full"
            continue
        fi
        die "no such file or directory: $full"
    done
}

# downloads
# 
# Emit the download cache directory
function downloads () {
    local d
    d="$(topdir)/downloads"
    mkdir -p "$d"
    echo "$d"
}

# wcb [options]
#
# Run the Wire-Cell build tool in the top directory.
function wcb () {
    local t here
    here="$(pwd)"
    t="$(topdir)"
    cd "$t" || die "failed to cd to $t"
    ./wcb "$@"
    cd "$here" || die "failed to cd to $here"
}

# wcb_env_dump
#
# Emit all Waf environment variable settings.
#
# Note, these are not necessarily varibable settings of the shell
# environment.
function wcb_env_dump () {
    # Keep a cache if running in bats.
    if [ -n "${BATS_FILE_TMPDIR}" ] ; then
        local cache="${BATS_FILE_TMPDIR}/wcb_env.txt"
        if [ ! -f "$cache" ] ; then
            wcb dumpenv | grep 'wcb: ' | sed -e 's/^wcb: '// | sort > "$cache"
        fi
        cat "$cache"
        return
    fi

    # otherwise, do not keep cache as there is no good way to
    # invalidate it.  This suffers about 100ms per call.
    wcb dumpenv | grep 'wcb: ' | sed -e 's/^wcb: '// | sort
}    
    

# wcb_env_vars [variable] ...
#
# Emit key=var lines consisting of Waf build environment variables.
#
# With no arguments, emit all Waf variables, else emit only those named.
#
# User may run "./wcb dumpenv" to dump what is available.
function wcb_env_vars () {
    # Keep a cache if running in bats.
    if [ -z "$1" ] ; then
        wcb_env_dump
    fi
    for one in "$@"
    do
        grep "^${one}=" <(wcb_env_dump)
    done
}


# wcb_env_value <varname>
#
# Emit the value for the Waf build environment variable named <varname>.
wcb_env_value () {
    wcb_env_vars "$1" | sed -e 's/[^=]*=//' | sed -e 's/^"//' -e 's/"$//'
}


# wcb_env <varname> [...]
#
# Evaluate a Waf build environment variable to become a shell
# environment variable.
function wcb_env () {
    eval "$(wcb_env_vars "$@")"
}


# usepkg <pkgname> [...]
#
# Set shell environment variables related to one or more named packages.
#
# This adds test program locations into PATH and defines a <pkg>_src
# variable pointing at the packages source directory.
function usepkg () {
    for pkg in "$@"; do
        local t pkg="$1" ; shift
        t="$(topdir)"
        printf -v "${pkg}_src" "%s/%s" "$t" "$pkg"
        export "${pkg}_src"
        PATH="$(blddir)/$pkg:$PATH"
    done
}

# compile_jsonnet <input> <output>
#
# Compile input Jsonnet file to output JSON file.
#
# Execute the Jsonnet compiler under Bats "run".
# Tries to use wcsonnet then falls back to jsonnet.
# The cfg dir is put into the search path and WIRECELL_PATH is ignored.
function compile_jsonnet () {
    local ifile="$1" ; shift
    local ofile="$1" ; shift
    local cfgdir
    cfgdir="$(topdir)/cfg"

    local orig_wcpath="$WIRECELL_PATH"
    WIRECELL_PATH=""

    # see below for wcsonnet wrapper.  Here we locally DIY.
    local wcsonnet
    wcsonnet="$(wcb_env_value WCSONNET)"
    [[ -n "$wcsonnet" ]]
    if [ -z "$wcsonnet" ] ; then
        fatal "Failed to get WCSONNET!  Cache problem?"
    fi
    declare -a cmd=( "$wcsonnet" -o "$ofile" -P "$cfgdir" "$@" "$ifile" )

    check "${cmd[@]}"
    [[ -s "$ofile" ]]
    WIRECELL_PATH="$orig_wcpath"
}

# wire-cell [options]
# 
# Execute the "wire-cell" CLI with standard handling.
function wire-cell () {
    cli="$(wcb_env_value WIRE_CELL)"
    [[ -n "$cli" ]]
    [[ -x "$cli" ]]
    check "$cli" "$@"
}

# wcsonnet [options]
# 
# Execute the "wcsonnet" CLI with standard handling.
function wcsonnet () {
    cli="$(wcb_env_value WCSONNET)"
    [[ -n "$cli" ]]
    [[ -x "$cli" ]]
    check "$cli" "$@"
}

# wcwires [options]
# 
# Execute the "wcwires" CLI with standard handling.
function wcwires () {
    cli="$(wcb_env_value WCwires)"
    [[ -n "$cli" ]]
    [[ -x "$cli" ]]
    check "$cli" "$@"
}

# version
#
# Emit the wire-cell version string.
function version () {
    cli=$(wcb_env_value WIRE_CELL)
    [[ -n "$cli" ]]
    [[ -x "$cli" ]]
    $cli --version
}

# skip_if_missing <program> ...
#
# Check if the running environment has access to the given programs.
#
# One or more programs may be checked.  First one that is not found
# will lead to the test being skipped.
function skip_if_missing ()  {
    for one in "$@"
    do
        local path
        path="$(which "$one")"
        if [ -n "$path" ] ; then
            continue
        fi
        skip "No program: \"$one\""
        return
    done
}
               

# wcpy <pkg> <command> [options]
#
# Execute the wirecell-<pkg> program with command and optoins.
#
# If the command is not found, the test will be skipped.
function wcpy () {
    local pkg="$1" ; shift
    local wcp
    wcp="$(which "wirecell-$pkg")"
    if [ -z "$wcp" ] ; then
        warn "No wirecell-$pkg found, install wire-cell-python.  Will 'skip' the current test."
        skip "No wirecell-$pkg found"
        return
    fi
    check "$wcp" "$@"
}


# dotify_graph <input> <output>
#
# Output a GraphViz "dot" graph representation of the input.
#
# Input may be Jsonnet or JSON.
#
# Output may be PNG, PDF, DOT
function dotify_graph () {
    local ifile="$1" ; shift
    [[ -n "$ifile" ]]
    [[ -f "$ifile" ]]
    local ofile="$1" ; shift
    [[ -n "$ofile" ]]
    local cfgdir
    cfgdir="$(topdir)/cfg"
    [[ -d "$cfgdir" ]] 

    declare -a cmd
    cmd=( "$(wcb_env_value WCPGRAPH)" )
    if [ "${#cmd[@]}" -eq 0 ] ; then
        die "failed to find wirecell-pgraph from Waf environment"
    fi
    [[ -x "${cmd[0]}" ]]    
    
    cmd+=( dotify -J "$cfgdir" "$@" "$ifile" "$ofile" )

    check "${cmd[@]}"
}



# divine_context
#
# Emit the Bats context.  Will emit: "test", "file" or "run".
function divine_context () {
    if [ -n "$BATS_TEST_NAME" ] ; then
        echo "test"
    elif [ -n "$BATS_TEST_FILENAME" ] ; then
        echo "file"
    elif [ -n "$BATS_RUN_TMPDIR" ] ; then
        echo "run"
    else
        die "apparently not running under bats"
    fi
}


# tmpdir
#
# Emit the Bats temporary directory in the context.
#
# With no arguments, this divines the correct one for the current Bats context.
#
# With a single argument it will emit the temporary directory for that context.
#
# Note, to not purge the tmp call:
#
#   bats --no-tempdir-cleanup
#
# Or, set WCT_BATS_TMPDIR.
function tmpdir () {
    if [ -n "$WCT_BATS_TMPDIR" ] ; then
        echo "$WCT_BATS_TMPDIR"
        return
    fi

    local loc="${1:-test}"
    local ret=""
    case $loc in
        run) ret="$BATS_RUN_TMPDIR";;
        file) ret="$BATS_FILE_TMPDIR";;
        *) ret="$BATS_TEST_TMPDIR";;
    esac
    if [ -z "$ret" ] ; then
        ret="$(mktemp -d /tmp/wct-bats.XXXXX)"
        warn "Not running under bats, will not remove tempdir: $ret"
    fi
    echo "$ret"
}

# cd_tmp [context]
#
# Change to the context temporary directory.
#
# If a context is given, that temporary directory will be used
# otherwise it will be divined.
function cd_tmp () {
    ctx="$1"
    if [ -z "$ctx" ] ; then
        ctx="$(divine_context)"
    fi
    local t
    t="$(tmpdir "$ctx")"
    mkdir -p "$t"
    cd "$t" || die "could not cd to $t"
}

# relative_path <path>
#
# Resolve a file path relative to the path of the BATS test file.
function relative_path () {
    local want="$1"; shift
    [[ -n "$BATS_TEST_FILENAME" ]]
    realpath "$(dirname "$BATS_TEST_FILENAME")/$want"
}

# resolve_pathlist <pathlist> ...
#
# Resolve multiple pathlists to a single list, removing any that do
# not exist.  Consequitive duplicates are removed.
#
# usage:  assume X does not exist
# paths=( $(resolve_pathlist A:X A B:C) )
# A
# B
# C
function resolve_pathlist () {
    IFS=' :' read -r -a parts <<< "$@"
    declare -a exists
    for one in "${parts[@]}"
    do
        if [ -d "$one" ] ; then
            exists+=( "$one" )
        fi
    done
    printf "%s\n" "${exists[@]}" | awk '!seen[$0]++'

}

# resolve_path <path> <pathlist> ...
#
# Resolve a file name in in or more path lists.
# A path list may be a :-separated list
#
# usage:
# local path="$(resolve_path emacs $PATH $MY_OTHER_PATH /sbin:/bin:/etc)"
# /usr/bin/emacs
function resolve_path () {
    local want="$1"; shift

    # already absolute
    if [[ "$want" =~ ^/.* ]] ; then
        echo "$want"
    fi

    for maybe in $(resolve_pathlist "$@") ; do
        if [ -f "${maybe}/$want" ] ; then
            echo "${maybe}/$want"
            return
        fi
    done
}

# find_paths <path> <pathlist>
#
# Emit all matching paths.
#
# usage:
# resolve_all emacs $PATH
# /usr/bin/emacs
# /usr/local/bin/emacso
function find_paths () {
    local want="$1" ; shift

    # already absolute, only one answer
    if [[ "$want" =~ ^/.* ]] ; then
        echo "$want"
    fi
    
    for maybe in $(resolve_pathlist "$@") ; do
        if [ -f "${maybe}/$want" ] ; then
            echo "${maybe}/$want" 
        fi
    done
}


# config_path <path>
#
# Emit an absolute path of a relative path found in config areas.
function config_path () {
    local maybe path="$1" ; shift
    maybe="$(topdir)/cfg/$path"
    if [ -f "$maybe" ] ; then
        echo "$maybe"
        return
    fi
    if [ -n "$WIRECELL_PATH" ] ; then
        resolve_path "$path" "$WIRECELL_PATH"
    fi
}

# resolve_file <file>
#
# Emit the absolute path for a given path.
#
# If absolute, emit as given.  If not found, exit.
#
# The following path lists will be checked for the given relative path
# on a first one wins basis.
# 
# - current working directory
# - relative to directory holding current Bats test file
# - any path lists given on function command line
# - $(topdir)/test/data/
# - $(topdir)/
# - $(topdir)/build/tests/
# - $WIRECELL_PATH
#
# See resolve_path to resolve paths in a policy-free manner.
#
# See relative_path to resolve paths relative to bats test file.
#
# See category_path to resovle paths against specific data repo category.
function resolve_file () {
    [[ -n "$BATS_TEST_FILENAME" ]]
    local want="$1"; shift

    if [ -f "$want" ] || [ -d "$want" ] ; then
        realpath "$want"
        return
    fi

    declare -a paths=(".")
    paths=("$(dirname "${BATS_TEST_FILENAME}")")
    paths+=("$@")
    local t
    t="$(topdir)"
    paths+=("$t/test/data" "$t" "$t/build/tests" "$WIRECELL_PATH")

    for path in $(resolve_pathlist "${paths[@]}")
    do
        local maybe="$path/$want"
        if [ -f "$maybe" ] || [ -d "$maybe" ] ; then
            echo "$maybe"
            return
        fi
    done
    die "No such file or directory: $want"
}



# category_version_paths [options]
#
# Emit all versions paths found for category.
#
# options:
#
# -c, --category <cat>
#   Set the category.  Default is "history".
#
# -d, --dirty
#   If given, will match locally modified versions.
#
# Results are emitted in lexical order of version string.
function find_category_version_paths () {
    local category="history"
    # match "git describe --tags" 
    local matcher='+([0-9])\.+([0-9])\.+([x0-9])'
    local dirty=""
    while [[ $# -gt 0 ]] ; do
        case $1 in
            -c|--category) category="$2"; shift 2;;
            -d|--dirty) dirty="${matcher}-*"; shift;;
            *) break;;          # leave others args in place
        esac
    done

    local catdir
    catdir="$(blddir)/tests/$category"

    declare -a vers
    declare -A pathmap

    for one in "$catdir"/$matcher
    do
        if [ -d "$one" ] ; then
            local ver
            ver="$(basename "$one")"
            vers+=( "$ver" )
            pathmap["$ver"]="$one"
        fi
    done
    if [ -n "$dirty" ] ; then
        for one in "$catdir"/$dirty
        do
            if [ -d "$one" ] ; then
                local ver
                ver="$(basename "$one")"
                vers+=( "$ver" )
                pathmap["$ver"]="$one"
            fi
        done
    fi
    for ver in $(printf '%s\n' "${vers[@]}" | sort -u)
    do
        echo "${pathmap["$ver"]}"
    done
}

# find_category_version <path> ...
#
# Same as category_version_paths but emit only version string.
function find_category_versions () {
    for one in $(category_version_paths "$@")
    do
        basename "$one"
    done
}

# category_path [options] <path> ...
#
# Resolve and emit path to a category path.
#
# options
#
# -c, --category <cat>
#   Set the category, default is "history".
#
# -v, --version <version>
#   Set the version, default uses current version.
#
# <path>
#   A relative path to resolve.
function category_path () {
    local ver cat="history"
    ver="$(version)"
    declare -a paths
    while [[ $# -gt 0 ]] ; do
        case $1 in
            -v|--version) ver="$2"; shift 2;;
            -c|--category) cat="$2"; shift 2;;
            -*) die "unknown argument: $1";;
            *) paths+=( "$1" ); shift;;
        esac
    done
            
    local catdir
    catdir="$(blddir)/tests/$cat/$ver"

    for path in "${paths[@]}"
    do
        if [ -f "$catdir/$path" ] || [ -d "$catdir/$path" ] ; then
            echo "$catdir/$path"
        else
            die "No such file or directory: $catdir/$path"
        fi
    done
}
    

# historical_versions
#
# Emit the versions for which there are known historical files.
#
# Versions are emitted delimited by newlines.
function historical_versions () {
    # shellcheck disable=SC2046
    printf '%s\n' $(wcb_env_value TEST_DATA_VERSIONS)
}


# historical_files [options] <path> ...
#
# Emit cross product of paths across versions.  Paths are
# newline-delimited.
#
# options:
#
# -v, --version <version>
#
#   Set the versions to consider.  Multiple -v options are allowed.
#   If this option is given the default list of versions will not be
#   considered.  The default versions are declared by the current
#   build system.
#
# -l, --last <number>
#
#   Emit the last N versions of the path.
#
# -c, --current
#
#   Include the current software version, even if not a decalred
#   release.
#
function historical_files () {
    local last current
    declare -a paths versions
    while [[ $# -gt 0 ]] ; do
        case $1 in
            -l|--last) last="$2"; shift 2;;
            -v|--version) versions+=( "$2" ); shift 2;;
            -c|--current) current="$(version)"; shift;;
            -*) die "unknown option $1" ;;
            *) paths+=( "$1" ); shift;;
        esac
    done

    local verlines
    if [ "${#versions[@]}" -eq 0 ] ; then
        verlines=$(historical_versions)
    else
        verlines=$(printf '%s\n' "${versions[@]}")
    fi
    if [ -n "$current" ] ; then
        verlines=$(printf '%s\n%s' "$current" "$verlines")
    fi
    verlines=$(echo "$verlines" | sort -u)

    if [ -n "$last" ] ; then
        verlines=$(echo "$verlines" | tail -n "$last")
    fi

    # shellcheck disable=SC2068
    for ver in ${verlines[@]}
    do
        yell "historical_files: ver=|$ver| paths: ${paths[*]}"
        category_path -c history -v "$ver" "${paths[@]}"
    done
}

# skip_if_no_category [options] <category>
#
# Skip test if category is missing.
#
# options:
# -v, --version <version>
#   Set the version of the category, default to current version.
#
# <category>
#   Set the category.  Defaults to "history".
function skip_if_no_category () {
    local ver
    ver="$(version)"
    local cat="history"
    while [[ $# -gt 0 ]] ; do
        case $1 in
            -v|--version) ver="$2"; shift 2;;
            *) category="$1"; shift;;
        esac
    done
    cdir="$(blddir)/tests/$cat/$ver"
    if [ -d "$cdir" ] ; then
        skip "no category: $cat/$ver"
    fi
}

    

# download_file <url> [<tgt>]
#
# Download a file from a URL.
#
# Emit (echo) a path to the target file.  This will be a location in a cache.
# 
# <url>
#   The URL to download.
#
# <tgt>
#   An optional relative path to the file name the URL content should be saved.
function download_file () {
    local url="$1"; shift
    local tgt
    tgt="${1:-$(basename "$url")}"

    local path
    path="$(downloads)/$tgt"
    if [ -s "$path" ] ; then
        echo "$path"
        return
    fi

    wget -O "$path" "$url" 1>&2 || return
    echo "$path"
}


# mv_if_diff <src> <dst>
#
# Move <src> to <dst> if they differ.
#
# This can be useful to initate a chain of run_idempotently steps.
mv_if_diff () {
    local src="$1" ; shift
    local dst="$1" ; shift
    
    if [ ! -f "$dst" ] ; then
        echo "dst does not exist $dst moving $src" 1>&3
        mv "$src" "$dst"
        return
    fi
    if [ -n "$( diff "$src" "$dst" )" ] ; then
        echo "differ: $dst and $src" 1>&3
        mv "$src" "$dst"
        return
    fi
    rm -f "$src"
}



function help_one () {
    local func="$1" ; shift
    awk "/^# $func\>/{flag=1}/^function $func /{flag=0}flag" "$0" | sed -e 's/^ //g' -e 's/^#//g'
}

function help_all () {
    grep ^function "$0" | sed -n 's/^function \(\S*\) .*/\1/p'| sort | while read -r func
    do
        help_one "$func" | head -1
    done
}

if [[ "$(basename -- "$0")" == "wct-bats.sh" ]]; then
    if [ -z "$1" ] ; then
        help_all
    else
        for one in "$@" ; do
            help_one "$one"
        done
    fi
fi
