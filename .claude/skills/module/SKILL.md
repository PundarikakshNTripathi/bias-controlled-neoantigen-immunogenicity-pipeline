---
name: module
description: Build one module of the pipeline end to end, test-first, against its contract. Takes a module tag such as M1 or M4.
argument-hint: [M1|M2|M3|M4|M5|M6|M7]
arguments: [module_tag]
allowed-tools: Read Write Edit Grep Glob Bash
---

## Build module $module_tag

Work test-first. The contract is already written; your job is to make it true.

1. **Read the contract.** Open the module's source file and read its docstring
   in full. It lists the numbered functional requirements this module must
   satisfy. Read the matching section of `docs/SRS.md` if present.

2. **Write the tests first.** Create `tests/test_<module>.py` covering every
   requirement in the contract. Each test names the requirement it pins in a
   comment. Tests that assert a degenerate input raises, rather than silently
   returning zero, are required — that property is load-bearing here.

3. **Run them and watch them fail.** Paste the failure output. A test that has
   never failed has not been shown to test anything.

4. **Implement** until they pass. Smallest change that satisfies the contract.

5. **Verify.** Run `python scripts/verify.py` and paste the real output.

6. **Report honestly.** State what is implemented, what is stubbed, and what
   you were unable to do. If a requirement could not be met, say which and why.
   A module that is 80% done and described as 80% done is fine; one described
   as finished is not.

7. **Commit** with a conventional message naming the requirements covered.

Delegate to `data-steward` for M1 and splits, to `ml-engineer` for M3 and M4,
and ask `repro-verifier` to confirm before you call it done.
