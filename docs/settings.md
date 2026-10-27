---
title: Settings
nav_order: 7
---

# Settings (`config/config.yaml`)

You don't need to change anything for a normal run. Open the file with
`nano config/config.yaml`. Each setting has a comment explaining it.

| Setting | Default | When to change it |
|---|---|---|
| `samples` | `config/samples.tsv` | To keep several sample sheets |
| `output_dir` | `results` | To keep several runs side by side, e.g. `results_run2` |
| `comparison.mode` | `absent` | Set to `sibling` when the control is a wild-type-looking **sibling** (or a pool of siblings) of a recessive mutant. See below. |
| `comparison.min_alt_reads_mutant` | 2 | Raise to 3–5 when coverage is high (> 20×), to cut false positives |
| `comparison.min_depth_control` | 4 | Raise to 8–10 at high coverage for stricter "absent from control" calls |
| `comparison.mutant_rule` | `any` | Set to `all` with several mutants that should share the same mutation |
| `comparison.min_alt_fraction_mutant` | 0.9 | Sibling mode: share of mutant reads that must carry the variant. Lower to 0.8 for a mutant pool with a few mis-sorted fish |
| `comparison.max_alt_fraction_control` | 0.6 | Sibling mode: highest share of control reads allowed to carry it |
| `comparison.demote_flagged` | `true` | Set to `false` to keep known variants and flagged amino-acid changes at their normal priority |
| `known_variants.enabled` | `true` | Set to `false` to skip the known-variant check (e.g. no internet on the server) |
| `protein_check.enabled` | `true` | Set to `false` to skip the amino-acid check |
| `protein_check.species` | human, rat, chicken, frog, medaka, fugu | Other Ensembl species names to compare with, e.g. `lepisosteus_oculatus` (spotted gar). (Ensembl 113 has no zebrafish–mouse orthologs, so rat stands in for rodents.) |
| `variant_calling.min_qual` | 30 | Lower to 20 for more sensitive (but noisier) calls |
| `variant_calling.min_depth` | 3 | Minimum reads at a site in the mutant |
| `homozygosity_runs.enabled` | `true` | Set to `false` to skip the runs-of-homozygosity step (`bcftools roh`) |
| `homozygosity_runs.min_reads`, `min_reads_each_allele` | `3`, `2` | Genotype calls with fewer reads, or "heterozygous" calls resting on a single read of one allele, are ignored before the run search. At ~3× coverage these single error-prone reads otherwise break every run |
| `homozygosity_runs.az_to_hw` | `1e-10` | Chance of a run ending between two sites; lower gives longer runs |
| `structural_variants.enabled` | `true` | Set to `false` to skip the slowest step |
| `resources.threads_per_job` | 16 | CPU threads for alignment |
| `reference.*` | GRCz11, Ensembl 113 | Only to use a different genome build |

{: .warning }
After changing a setting, run `./run_wgs.sh` again. Only the affected steps are redone.

## Choose a comparison mode

The right mode depends on what your **control** fish are.

| Your control is… | Use | What counts as a candidate |
|---|---|---|
| An unrelated wild-type fish or line (e.g. the parental strain) | `mode: absent` (default) | In the mutant, and **no reads at all** in the control |
| Wild-type-looking **siblings** of a recessive mutant, one fish or a pool | `mode: sibling` | In (nearly) **all** mutant reads; the control may carry it in **some** reads |

**Why siblings need their own mode.** For a recessive mutation, two thirds of the
wild-type-looking siblings carry one copy. One carrier sibling shows the mutation in
about half of its reads, and a sibling pool in about a third. In `absent` mode such a
variant is called *shared* and thrown away, so the real mutation can be missed.
`sibling` mode keeps it, and labels it **mutant-enriched**.

In return, `sibling` mode demands that the mutant carries the variant on **both copies**:
by default at least 90% of the mutant reads (`min_alt_fraction_mutant`). It drops
variants carried by a well-covered control in more than 60% of its reads
(`max_alt_fraction_control`), because no wild-type-looking sibling can be homozygous for
the mutation.

{: .note }
Sibling mode assumes a **recessive** mutation. For a dominant mutation or a new
mutation in one fish, keep `mode: absent`.
