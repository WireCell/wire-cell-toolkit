#!/usr/bin/env python3
"""WCT policy guard — a PreToolUse hook.

Deterministically catches edits that would violate firm WCT policy and routes
them to the human (permissionDecision "ask") with a reminder, instead of letting
the agent proceed on its own judgement.

Wire it up in .claude/settings.json:

  "hooks": {
    "PreToolUse": [
      { "matcher": "Edit|Write|MultiEdit",
        "hooks": [ { "type": "command",
                     "command": "python3 .claude/hooks/policy-guard.py" } ] }
    ]
  }

This is a high-signal starting point, not exhaustive. Tune the patterns to taste;
prose policy in CLAUDE.md covers everything the hook does not.
"""
import json
import os
import re
import sys

data = json.load(sys.stdin)
ti = data.get("tool_input", {}) or {}
path = ti.get("file_path", "") or ""
base = os.path.basename(path)

# Gather the text this edit would introduce (Edit/Write/MultiEdit shapes).
added = " ".join(
    str(ti.get(k, ""))
    for k in ("content", "new_string")
) + " ".join(
    str(e.get("new_string", "")) for e in ti.get("edits", []) or []
)

reasons = []

# 1. New external dependency in a build file (team-gated decision).
is_build_file = base in ("CMakeLists.txt", "wscript_build") or base.endswith(".cmake")
if is_build_file and re.search(r"\bfind_package\s*\(", added):
    reasons.append(
        "This build-file edit adds a find_package(...) — a new external dependency. "
        "Per WCT policy, new dependencies are a team decision, not an agent action."
    )

# 2. Any dependency-shaped change under util/ (util takes no new deps, ever).
norm = path.replace("\\", "/")
if ("/util/" in norm or norm.startswith("util/")) and is_build_file \
        and re.search(r"find_package|target_link_libraries", added):
    reasons.append(
        "This touches util/ build dependencies. util/ takes NO new dependencies, ever."
    )

if reasons:
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "ask",
            "permissionDecisionReason": " ".join(reasons)
            + " Confirm only if the WCT team has approved this.",
        }
    }))

sys.exit(0)
