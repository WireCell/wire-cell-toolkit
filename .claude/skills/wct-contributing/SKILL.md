---
name: wct-contributing
description: How to prepare GitHub Issues and Pull Requests for the wirecell/wire-cell-toolkit repo, including the human-preamble rule and the attribution wrapper for any LLM-generated text. Use whenever drafting Issue or PR content.
---

# Contributing Issues & PRs to wire-cell-toolkit

## The preamble is the human's — do not write it

Every Issue and PR opens with a **preamble the human author writes in their own words**:

- **Issue:** describe the problem being solved.
- **PR:** reference the Issue it addresses, and describe how the solution works.

**Do not author this preamble.** If asked to, decline and remind the user that agreed
project practice requires the human to write it themselves. This is a firm policy (see
`CLAUDE.md`), not a stylistic preference.

## What you may contribute

You may provide a **summary** — but only:

1. **Appended** to the human's preamble, never replacing it, and
2. Wrapped so it is unmistakably attributed to a model, using the `<details>` element
   with your model name and version as the `<summary>`:

```
<details><summary>Model Name and Version</summary>

LLM generated content here.

</details>
```

The same wrapper applies to any LLM-generated text posted as an Issue/PR **comment**.

## Quick flow

1. Confirm the human has written (or will write) their own preamble.
2. Offer a summary only if wanted; place it appended and inside the wrapper above.
3. Never present model-generated prose as if authored by the human.
