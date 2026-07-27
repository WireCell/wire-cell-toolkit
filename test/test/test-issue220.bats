#!/usr/bin/env bats

# We eradicate the trace tag "threshold" and put its trace summary on
# "wiener".

bats_load_library wct-bats.sh

# bats file_tags=issue:220

@test "no wiener_threshold_tag in cfg" {

    local cfgdir=$(srcdir cfg)

    # Assure the deprecated parameter is not used in any active configuration.
    # We must ignore commented-out occurrences, so strip Jsonnet comments
    # (// and # line comments plus single-line /* */ block comments) from each
    # file before matching.  Otherwise a mention in a comment yields a false
    # positive.
    local found=""
    local f
    for f in $(find "$cfgdir" -name '*.jsonnet' -print); do
        local hit
        hit=$(sed -E 's://.*$::; s:#.*$::; s:/\*[^*]*\*/::g' "$f" \
                  | grep -nw 'wiener_threshold_tag' || true)
        if [[ -n "$hit" ]]; then
            found+="$f: $hit"$'\n'
        fi
    done
    echo -e "found with wiener_threshold_tag:\n$found"
    [[ -z "$found" ]]

}

