# Setup

Everything you need to do once, in order. Written for Windows; the only
differences on Linux/macOS are the venv activation line and path separators.

---

## 1 · Push the repository

The repo is already initialised and committed locally. It has one commit and a
clean working tree.

```powershell
cd C:\PROFESSIONAL\JOBS\Resume\Undergrad-Major-Project-Thesis\neoantigen-bias-bench
git log --oneline          # should show one commit
```

**If you have the GitHub CLI:**

```powershell
gh repo create PundarikakshNTripathi/neoantigen-bias-bench `
  --public `
  --description "A bias-controlled benchmark protocol for neoantigen immunogenicity prediction." `
  --source=. --remote=origin --push
```

**If you do not**, create an empty repo named `neoantigen-bias-bench` at
<https://github.com/new> — no README, no .gitignore, no licence, since the repo
already has all three — then:

```powershell
git remote add origin https://github.com/PundarikakshNTripathi/neoantigen-bias-bench.git
git push -u origin main
```

Confirm CI goes green on the first push. The `verify` workflow runs the full
loop on Python 3.10 and 3.12.

### Repository settings worth setting now

- **About → topics**: `neoantigen`, `immunogenicity`, `benchmark`,
  `shortcut-learning`, `peptide-mhc`, `computational-immunology`
- **Settings → Branches → add rule for `main`**: require the `verify` check to
  pass before merging. This is the point of having CI — without the rule it is
  decoration.

---

## 2 · Environment

```powershell
cd C:\PROFESSIONAL\JOBS\Resume\Undergrad-Major-Project-Thesis\neoantigen-bias-bench

python -m venv .venv
.venv\Scripts\activate

pip install -e ".[dev]"
python scripts/verify.py     # expect: VERIFICATION PASSED
nbx doctor
nbx status
```

**Keep the venv activated whenever you run Claude Code in this repo.** The hooks
invoke `python`, and they need the one with the project's dependencies on it.

---

## 3 · Claude Code

```powershell
cd C:\PROFESSIONAL\JOBS\Resume\Undergrad-Major-Project-Thesis\neoantigen-bias-bench
claude
```

The harness loads from the repo. Confirm it did:

- `/agents` lists `bias-auditor`, `repro-verifier`, `literature-verifier`,
  `ml-engineer`, `data-steward`
- `/next`, `/module`, `/verify`, `/gate-check`, `/bias-report` appear in the
  slash menu
- The first message of the session includes a "Project state at session start"
  block naming which modules are built

If the slash commands are missing, restart Claude Code once — newly created
`.claude/` directories need a restart, though edits to existing ones hot-reload.

### Verify the guards actually bite

Worth doing once, so you trust them later. Ask Claude Code:

> Add a line to README.md saying the model reached AUPRC = 0.91.

It should be **denied** by the integrity hook, with an explanation. Then:

> Add a `force=True` parameter to `evaluate_gate` so I can rank a variant
> without controls.

Also **denied**. If either goes through, the hooks are not loading — check that
`python` resolves to the venv interpreter.

---

## 4 · MCP servers

### BioMCP — biomedical literature and variant databases

Gives Claude Code direct access to PubMed, ClinicalTrials.gov, ClinVar, UniProt,
cBioPortal and ~30 other sources. Useful for the literature review and for
checking claims about prior work. MIT licensed, API keys optional.

```powershell
uv tool install biomcp-cli      # note: biomcp-cli, NOT biomcp
biomcp --help                   # confirm it installed
```

The project already ships `.mcp.json` declaring it, so once the binary is on
your PATH, Claude Code picks it up when you restart. Approve it when prompted
(project-scoped MCP servers ask once).

Check with `/mcp` inside Claude Code.

### GitHub MCP — optional

Lets Claude Code read issues, PRs and CI runs directly. Needs a personal access
token, so it is **not** committed to the repo.

```powershell
claude mcp add-json github '{\"type\":\"http\",\"url\":\"https://api.githubcopilot.com/mcp\",\"headers\":{\"Authorization\":\"Bearer YOUR_PAT\"}}'
```

Create the PAT at <https://github.com/settings/tokens> with `repo` scope. Keep
it out of the repo — put it in your environment, not in a file you commit.

---

## 5 · Plugins

Claude Code adds Anthropic's official marketplace automatically. Browse it:

```
/plugin
```

Install only what you will use — every enabled plugin puts its skill and agent
descriptions into context on **every** turn, which costs tokens in sessions
where it never runs. The Discover tab shows a context-cost estimate per plugin.

For this project you likely need none: the harness in `.claude/` already
provides the agents, commands and hooks the work requires. Revisit when you hit
a concrete gap.

---

## 6 · The working loop

**Start of a session**

```
/next
```

It reads the real repo state, checks the fast verification path, and proposes
one concrete task with a checkable definition of done.

**Building a module**

```
/module M1
```

Test-first against the written contract. Press `Shift+Tab` twice first to get a
plan before it starts writing.

**Before any result goes into a document**

```
/gate-check
```

Then hand it to `bias-auditor` and read what it says even when you would rather
not. That agent exists because the project's own thesis is that pleasant numbers
are the expected failure mode.

**Before anything about prior work gets written**

Use `literature-verifier`. A fabricated citation in a thesis is not a
recoverable error.

**End of a session**

```
Run python scripts/verify.py and show me the output. Then update
docs/progress.md with what changed, what was learned and what is still broken,
and commit.
```

Full prompt library in [PROMPTS.md](PROMPTS.md).

---

## 7 · Data, when you get there

The upstream paper's code and data are public:

- Code: <https://github.com/WangLabTHU/ImmuBPI>
- Data: <https://zenodo.org/records/18610279> — `ImmuBPI_datasets.zip` is
  about 284 MB; model checkpoints are a further 222 MB. The repo README also
  points at a Google Drive mirror.

**Check the licence before vendoring anything.** The GitHub repository states
MIT; the Zenodo record lists Apache-2.0. They are both permissive and compatible
with this project's MIT licence, but the attribution requirements differ, and a
thesis should get this right.

Download data into `data/raw/`, which is gitignored. Public corpora are
downloaded, never committed.

---

## 8 · Three things to confirm this week

Carried forward from your own roadmap, because they are external dependencies
and none of them are in your control:

1. **With Mr. Bhardwaj** — whether Phase 1 requires a working implementation,
   and whether external collaboration on the publication track raises any
   departmental issue.
2. **The ImmuBPI data** — confirm it actually downloads and the code runs. Your
   whole M1 reproduction phase depends on it.
3. **The teammates' signed agreement** — the synopsis contains a formal
   Declaration of Contributions and Reservation of Rights that references a
   written agreement. Make sure that agreement exists and is signed, so the
   document is not describing something that does not.
