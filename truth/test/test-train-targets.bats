#!/usr/bin/env bats

# pr_ml doc 25 sec 7.4: the training-time label transform (truth/python/train_targets.py) on its synthetic event --
# frame shift on both drift sides, Michel / delta em-daughter distance cut (near kept, far to LE / shower), delta fold,
# d_traj mask (ghosts untouched), a root without points, larg4 folded ids (orig_tid path).

bats_load_library wct-bats.sh

@test "truth train_targets selftest" {
    run python3 "$BATS_TEST_DIRNAME/../python/train_targets.py" --selftest
    echo "$output"
    [ "$status" -eq 0 ]
    [[ "$output" != *FAIL* ]]
}
