---
title: Tutorial slides
nav_order: 1.5
---

# Tutorial slides

A 67-slide lab tutorial on the whole process: the ideas behind the analysis, running
the pipeline, what each step does, reading the results, and a worked example on our
jaw mutant. Every slide has speaker notes, so anyone can present it.

[Open the slides]({{ '/slides/' | relative_url }}){: .btn .btn-primary .fs-5 .mb-4 }

<iframe src="{{ '/slides/' | relative_url }}" title="Tutorial slides preview"
  style="width:100%;aspect-ratio:16/9;border:1px solid #e4e3df;border-radius:8px"></iframe>

## Presenting

| Key | Does |
|---|---|
| → or Space, ← | Next / previous slide (or click, or swipe on a phone) |
| N | Show or hide the speaker notes |
| F | Full screen |
| Home / End | First / last slide |

Each slide has its own address (for example `slides/#20`), so you can link to a step.
To make a PDF or print handouts, open the slides and use your browser's **Print**
(choose *Save as PDF*): every slide becomes one page.

## The parts

| Part | Slides | Covers |
|---|---|---|
| 1. The ideas | 3–7 | Mutant vs control, key words, sibling genetics, the linked region |
| 2. Running it | 8–17 | Setup, copying files, the sample sheet, comparison mode, settings, running |
| 3. Inside the pipeline | 18–40 | Each of the eight steps, the checks, the math behind each test, and the priority 1–5 ranking |
| 4. Reading the results | 41–45 | Files, the report, the Excel workbook, judging a variant |
| 5. Worked example | 46–58 | The jaw mutant: data, peak, runs of homozygosity, homozygosity mapping, the two-copy block at 24–26 Mb, one fish or a pool, false leads, remaining leads, next steps |
| 6. Appendix: the math | 59–71 | Error models, genotype likelihoods, PL, Bayes, the bcftools caller, depth, homozygosity HMM |

## The math behind each test

Nine slides, labelled *The math* in orange, sit next to the step they explain and
work through each calculation with a real variant from our data:

| Slide | Explains |
|---|---|
| [22]({{ '/slides/#22' | relative_url }}) | Quality scores: what Q20, Q30 and QUAL mean as error chances |
| [24]({{ '/slides/#24' | relative_url }}) | Calling a genotype from read counts (homozygous, heterozygous or reference) |
| [30]({{ '/slides/#30' | relative_url }}) | Why sibling mode uses 90% and 60%, and how often each rule gets it wrong |
| [33]({{ '/slides/#33' | relative_url }}) | Protein alignment scores (BLOSUM62) and the species-counting rules |
| [36]({{ '/slides/#36' | relative_url }}) | The Fisher test behind the Manhattan plot's −log10 p |
| [37]({{ '/slides/#37' | relative_url }}) | Fisher's exact test step by step: dealing the reads out at random |
| [38]({{ '/slides/#38' | relative_url }}) | Reading the Manhattan plot: bands, the depth ceiling, and why chromosome 7 towers |
| [39]({{ '/slides/#39' | relative_url }}) | Window statistics: the peak, homozygosity and copy number |
| [40]({{ '/slides/#40' | relative_url }}) | Homozygosity mapping: the share of mutant reads from another version, and how the best interval is chosen |

Part 6 is an optional appendix with the full derivations, again checked against our
own VCF (the tmppe variant):

| Slide | Explains |
|---|---|
| [60]({{ '/slides/#60' | relative_url }}) | The big picture: likelihood × prior = posterior, then decision rules |
| [61]({{ '/slides/#61' | relative_url }}) | Five Phred numbers (base quality, MAPQ, PL, GQ, QUAL) and what each means |
| [62]({{ '/slides/#62' | relative_url }}) | Three error sources: base quality, mapping quality, base alignment quality (BAQ) |
| [63]({{ '/slides/#63' | relative_url }}) | Genotype likelihoods from single reads |
| [64]({{ '/slides/#64' | relative_url }}) | Why an allele-fraction cut-off is not the calling model; correlated errors |
| [65]({{ '/slides/#65' | relative_url }}) | PL values, checked against our VCF |
| [66]({{ '/slides/#66' | relative_url }}) | Bayes' rule, Hardy–Weinberg priors and genotype quality (GQ) |
| [67]({{ '/slides/#67' | relative_url }}) | The bcftools multiallelic caller, its prior and site QUAL |
| [68]({{ '/slides/#68' | relative_url }}) | Depth: DP vs AD, and Poisson coverage at 3× |
| [69]({{ '/slides/#69' | relative_url }}) | Duplicates, indel normalisation and filters: rules, not probabilities |
| [70]({{ '/slides/#70' | relative_url }}) | Runs of homozygosity as a hidden Markov model, and how the pipeline runs `bcftools roh` |
| [71]({{ '/slides/#71' | relative_url }}) | Six common ways to misread a VCF |

## Updating the slides

The slides live in `tutorial/slides/`, one HTML file per slide, in the order listed
in `tutorial/deck.json`. Formulas are written in LaTeX between `\(` and `\)` (inline)
or `\[` and `\]` (on their own line) and typeset by KaTeX when the page opens.
After editing the slides, rebuild this page and commit:

```bash
python3 tutorial/build_slides.py
```
