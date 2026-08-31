# Welcome to the Wire-Cell Toolkit.

This file provides guidance for humans to contribute to the Wire-Cell Toolkit
(WCT) project.

The use of LLM coding agents for such contributions is allowed following the
guidelines in this file.  LLM agents will see CLAUDE.md for further guidance.
Humans and LLMs alike may read README.org for general entry points to other
guidance and documentation.

## Reporting issues

Bugs and feature requests or new ideas can and should be provided using GitHub Issues

  https://github.com/WireCell/wire-cell-toolkit/issues
  
An Issue MUST include a summary of the problem that is written by a human.  Any
Issue that lacks sufficient human context and understanding may be summarily
closed.

An Issue MAY contain LLM generated information following the human description.
If included, it should be wrapped in `<details>` and `<summary>` as in this
example:

```
<details><summary>Model Name and Version</summary>

LLM generated content here.

</details>
```

The WCT repo provides instructions and "skills" to direct popular LLM coding
agents to comply.

## Contributing code

Bug fixes and new features are welcome and should be offered via a GitHub Pull Request.

  https://github.com/WireCell/wire-cell-toolkit/pulls
  
Major contributions from the WCT core developers are also expected to be done
via a PR.  Quick bug fixes or limited ancillary work may be directly pushed to
the "master" branch by core developers on a case-by-case basis.

Like Issues, the initial PR comment must be written by a human and describe WHAT
problem the PR addresses, citing relevant Issues and describing, and HOW the PR
addresses the problem.  

And, like Issues, a PR comment MAY contain LLM generated information in the same
manner as an Issue.

Contributing code indicates that the contributor agrees the code may be
distributed according to WCT's license (see LICENSE file).

