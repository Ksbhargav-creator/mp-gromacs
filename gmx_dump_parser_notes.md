# `parse_gmx_dump.py` — design log, difficulties, and rationale

How the `gmx dump` parser evolved from a one-line LJ regex to a format-agnostic
functype tokenizer, and the problems that forced each change. Written so the parser
is maintainable and so the method is reproducible in the paper.

## 1. Why parse `gmx dump` at all

The static operand path reads a `.tpr` with MDAnalysis (pure Python, no GROMACS).
That gives per-atom **charges** and **masses** directly — enough for the
`nonbonded_coulomb` and `integration` kernels. But we confirmed by inspection that
MDAnalysis exposes only *connectivity and type tuples* for bonds/angles, and its
`bond.value`/`bond.length` are **geometry methods** (computed from coordinates), not
the stored force-field constants. It does **not** expose:

- the Lennard-Jones nonbonded matrix (c6/c12),
- bonded force constants and reference values (k, b₀, θ₀, dihedral params),
- constraint lengths.

All of those *are* in the `.tpr`, and the stock tool `gmx dump -s system.tpr` prints
them as text in the `ffparams` section. So `gmx dump` is the authoritative static
source for the LJ, bonded, and constraint kernels — hence a text parser.

## 2. The evolution (four versions)

**v1 — LJ only, line regex.** A regex matching any line containing `c6= <float>` and
`c12= <float>`. Validated against real MDAnalysisTests `.tpr` files (via a bundled
gmx-style dump). Worked, shipped, produced the `nonbonded_lj` records.

**v2 — bonded/constraints, *hardcoded* regexes (the mistake).** I extended it with
single-line regexes keyed on the functype **name** and the field names I *remembered*:
`BONDS … b0A= … cbA=`, `ANGLES … thetaA= … ctA=`, `CONSTR … dA=`. I validated this
against a **synthetic** dump I wrote in that assumed format — it passed. But I could
not test against a real `gmx dump`, because the sandbox has no `gmx` binary and can't
build one. So v2 shipped essentially *blind*.

**v3 — generic tokenizer + `--diagnose` (the fix).** On the real benchMEM dump, v2
returned **0** bonded/constraint entries. Rather than guess again, I made the parser
**format-agnostic**: `parse_functypes()` tokenizes *every* `functype[N]=NAME, key= val,
key= val, …` line into `NAME -> [{field: value}, …]`, regardless of which name or
fields appear. A `--diagnose` mode then prints every family and its field keys — turning
a blind guessing loop into a single round-trip: the user runs `--diagnose` and pastes
the actual structure.

**v4 — locked to the real fields.** The `--diagnose` output of benchMEM revealed the
real structure (below). I mapped the parser to those exact fields, keeping
variant-tolerant key lists so it still works on other versions/force fields.

## 3. What `--diagnose` actually revealed (benchMEM)

```
ANGLES   121   fields: ['ctA','ctB','thA','thB']
CONSTR    61   fields: ['dA','dB']
IDIHS     11   fields: ['cxA','cxB','xiA','xiB']
LJ14     154   fields: ['c12A','c12B','c6A','c6B']
LJ_SR    961   fields: ['c12','c6']
PDIHS     43   fields: ['cpA','cpB','mult','phiA','phiB']
RBDIHS     1   fields: []
SETTLE     1   fields: ['dhh','doh']
```

Three surprises, each of which had broken v2:

1. **Angle reference is `thA`, not `thetaA`.** A pure naming mismatch — v2 looked for a
   field that doesn't exist in this build, so `angle_theta0` came back empty even though
   the data was right there.
2. **The dihedrals are `PDIHS`/`IDIHS`, which v2 didn't parse at all.** v2 only knew
   BONDS/ANGLES/CONSTR, so the entire dihedral operand set was silently dropped.
3. **There is no `BONDS` family.** This system uses `constraints=all-bonds`, so every
   bond was converted to a `CONSTR`. `bond_k`/`bond_b0` are *correctly* empty — the bond
   stiffness lives in the constraint kernel, not the bonded kernel. This is a genuine
   characterization result, not a parser failure.

## 4. The core difficulties

- **No `gmx` in the dev sandbox.** The parser had to be written against a *remembered*
  format and validated only on synthetic input. This is the root cause of the v2 miss:
  synthetic validation proves the *logic*, never the *format assumption*.
- **`gmx dump` format is version- and force-field-dependent.** Field names differ (`thA`
  vs `thetaA`), functype names differ (CHARMM angles are `UREY_BRADLEY`; RB vs proper
  dihedrals), and which families even appear depends on run settings (`constraints=…`).
  Any hardcoded parser is brittle by construction.
- **A `0` count is ambiguous.** Empty output can mean (a) wrong field name, (b) a family
  we don't parse, or (c) a genuine absence. These need *different* fixes, and you can't
  tell which without seeing the raw structure — which is exactly why `--diagnose` exists.
- **Multi-line / array formats.** `RBDIHS` showed `fields: []` — its six Ryckaert–
  Bellemans coefficients are printed in a form our single-line `key=val` scan doesn't
  catch (likely an array or continuation line). Flagged, not yet handled.
- **Zero-value filtering.** Dihedral phases and improper reference angles are often
  exactly `0.0` (or `180.0`); the magnitude reporter drops zeros, so a field can look
  "empty" in the console even though the values were collected. A reporting nuance, not a
  data loss.

## 5. Design decisions (and why)

- **Generic tokenizer over hardcoded regex.** `parse_functypes()` never assumes the set
  of names/fields, so a new force field can't silently break it — worst case it lands in
  a family we haven't mapped, which `--diagnose` makes visible.
- **Variant-tolerant field lists.** Each quantity tries several plausible keys
  (`thetaA`, `thA`, `theta0`, `tA`), so common cross-version/force-field naming
  differences are absorbed without code changes.
- **`--diagnose` as a first-class tool.** It converts an unbounded debugging loop into
  one deterministic round-trip. This is the single most important design choice.
- **Per-distinct-type sampling.** We record one value per *functype* (distinct
  interaction type), matching how the LJ `nbfp` matrix is stored — a modeling choice:
  it characterizes the parameter set, not the per-instance usage frequency. If a
  usage-weighted distribution is wanted later, weight each type by its interaction count.

## 6. Known limitations / TODO

- **`LJ14`** (1–4 scaled LJ) and **`RBDIHS`** are present but not wired in.
- **Multi-line functype formats** (e.g. RB coefficients) aren't stitched.
- **Per-type vs. per-instance** weighting is a documented modeling choice, not measured.
- **CHARMM:** `UREY_BRADLEY` is already handled; **CMAP** (the φ/ψ correction-map grid)
  is a distinct `structured_grids` kernel that will need its own extraction when a
  CHARMM36m system is characterized.

## 7. Outcome

Five of GROMACS's six computational kernels are now obtainable statically from
`.tpr` + `gmx dump`, with no instrumented build: `nonbonded_coulomb`, `nonbonded_lj`,
`bonded`, `lincs_constraint`, `integration`. Only `pme_reciprocal` (runtime, reciprocal-
space) requires the instrumented build. Verify a new dump's format at any time with
`python3 scripts/parse_gmx_dump.py --diagnose system.dump.txt`.
