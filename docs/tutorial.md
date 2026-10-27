---
title: Tutorial slides
nav_order: 1.5
---

# Tutorial slides

A 51-slide lab tutorial on the whole process: the ideas behind the analysis, running
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

## The five parts

| Part | Slides | Covers |
|---|---|---|
| 1. The ideas | 3–7 | Mutant vs control, key words, sibling genetics, the linked region |
| 2. Running it | 8–17 | Setup, copying files, the sample sheet, comparison mode, settings, running |
| 3. Inside the pipeline | 18–37 | Each of the eight steps, the checks, the math behind each test, and the priority 1–5 ranking |
| 4. Reading the results | 38–42 | Files, the report, the Excel workbook, judging a variant |
| 5. Worked example | 43–51 | The jaw mutant: data, peak, false leads, remaining leads, next steps |

## The math behind each test

Six slides, labelled *The math* in orange, sit next to the step they explain and
work through each calculation with a real variant from our data:

| Slide | Explains |
|---|---|
| [22]({{ '/slides/#22' | relative_url }}) | Quality scores: what Q20, Q30 and QUAL mean as error chances |
| [24]({{ '/slides/#24' | relative_url }}) | Calling a genotype from read counts (homozygous, heterozygous or reference) |
| [30]({{ '/slides/#30' | relative_url }}) | Why sibling mode uses 90% and 60%, and how often each rule gets it wrong |
| [33]({{ '/slides/#33' | relative_url }}) | Protein alignment scores (BLOSUM62) and the species-counting rules |
| [36]({{ '/slides/#36' | relative_url }}) | The Fisher test behind the Manhattan plot's −log10 p |
| [37]({{ '/slides/#37' | relative_url }}) | Window statistics: the peak, homozygosity and copy number |

## Updating the slides

The slides live in `tutorial/slides/`, one HTML file per slide, in the order listed
in `tutorial/deck.json`. After editing them, rebuild this page and commit:

```bash
python3 tutorial/build_slides.py
```
