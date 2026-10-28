# Jaw candidate search, chr7 (completed 2026-10-19)

## Inputs
- Homozygosity peak: chr7 27–33 Mb (GRCz11); broader enrichment ~14–33 Mb. See results/plots/peak.json.
- `chr7_20-35Mb_candidate_genes.tsv`: 110 genes, 234 HIGH/MODERATE variants.
  - Selection: mutant shows only variant reads (at least 2); control shows the variant in at most 60% of its reads.
  - This keeps "shared" variants, allowing for the control being a sibling or a sibling pool.
- `chr7_protein_changing_known_flag.tsv`: the same 234 variants, flagged against Ensembl 113 known zebrafish variation (dbSNP). 118 known, 116 novel.
- `chr7_novel_protein_changing.tsv`: the 116 novel ones.

## Methods
- ZFIN phenotype and wild-type expression files, searched for jaw, cartilage, arch and neural-crest terms. Covered all 110 candidate genes, plus all 464 genes in chr7 20–35 Mb.
- MGI mouse phenotypes for 85 orthologs.
- Literature review by three research passes (PubMed, ZFIN, MGI, OMIM/Monarch).
- Ortholog alignments (Ensembl REST) to test conservation of each residue.
- VEP/SIFT for selected variants.
- Remote query of Ensembl known variation for chr7 20–35 Mb.
- Read-depth checks of exons and of structural-variant calls.
- Unfiltered (QUAL<30) calls checked in the coding exons of the jaw genes.

## Conclusions

**No variant combines strong jaw biology, a damaging effect and novelty.**

### Ruled out or deprioritised
- **aldh1a2** (P393A, H56P, S72A): GRCz11 carries the rare residues H56/S72/P393. Both published wild-type cDNAs (UniProt Q90XS8, A0A8M1PD57) have P56/A72/A393, as the mutant does. These are reference differences, not lesions, verified in this session. No exon loss.
- **slc39a13** A23P has the best jaw biology: the mouse KO has a small mandible and maxilla, and humans with mutations have SCD-EDS. The residue is conserved Ala (human, mouse, medaka, fugu) and SIFT scores it 0.02. But it is a known polymorphism, rs507894324, so it is unlikely to be an ENU lesion.
- **lgr4** I745V/G904A, **adam10a** M194V and **itga11a**: the mutant allele matches the vertebrate or fish consensus, so these are benign.
- **ext2** (dackel): the mutant is not homozygous there (1 homozygous variant), so it lies outside the linked block.
- HIGH-impact calls in unnamed genes are annotation artifacts, a transposon (ENSDARG00000111781), or non-canonical isoforms. The ENSDARG00000073844 splice variant misses the canonical transcript.
- The 10 kb "homozygous deletion" of si:dkey-61m24.1 (pseudogene) is a repeat artifact: MAPQ0 depth is about 10× in both samples. Other DELLY deletions in the peak are intronic or intergenic.

### Remaining weak leads (novel variants)
1. **HERC1**: human biallelic LoF causes MDFPMR (facial dysmorphism, macrocephaly), and the mouse KO has bone defects. Novel variants:
   - +2 splice-donor insertion at 29,203,141 (2 reads; low evidence);
   - Q3550H (7/7 reads in the mutant, 0/5 in the control; high confidence);
   - Q2467L.
2. **plekha7b**: human PLEKHA7 is implicated in cleft lip/palate. Novel A133T (6 reads) and P94S.

### Known jaw genes in the region with no coding lesion found
All exons are covered in the mutant at about 4–5×. Low-QUAL raw calls were checked and none is homozygous in the mutant.
- **cog4** and **thap11** (core, 29.0 Mb): strongly homozygous block.
  - thap11 has a novel mutant-specific 5'UTR variant, c.-217A>T (chr7:29,045,070; 4/4 mutant, 0/4 control).
  - cog4/thap11 also share several novel 3'UTR variants.
- rgp1, b3gat3, cyp26b1, mmp14a (24–25 Mb flank): no novel coding variants.

## Suggested next steps
1. Confirm linkage first. Genotype 3–4 chr7 SNPs across about 20 individual mutants and 20 siblings, e.g. HERC1 Q3550H (29.21 Mb) and markers at 27 and 33 Mb. This narrows the interval.
2. Within the narrowed interval:
   - Sanger-sequence the coding exons and UTRs of cog4, thap11 and HERC1 in mutants.
   - Or test by crispant/complementation: thap11 and cog4 have known jaw phenotypes.
3. Phenotype comparison: Alcian blue at 4–5 dpf vs the published cog4, thap11 and aldh1a2 phenotypes.
4. Deeper paired-end sequencing (≥15×) of a mutant pool vs a sibling pool would resolve the interval far better.

## Pipeline improvements suggested by this analysis (implemented 2026-10-19)
- **Known variation.** Ensembl 113 zebrafish variants (dbSNP/EVA) are normalised once
  (`reference/known_variants/known.vcf.gz`). Their rsIDs are stamped onto our calls by exact
  position and allele (rule `annotate_known`). New columns: `known_id`, `novel`.
  - On chr7 20–35 Mb this reproduces the manual REST-based flags: 118 known, 116 novel, 234/234 agree.
- **Sibling/pool mode.** `comparison: mode: sibling`.
  - The mutant must carry the variant in ≥ 0.9 of its reads (`min_alt_fraction_mutant`).
  - Each well-covered control may carry it in at most 0.6 of its reads (`max_alt_fraction_control`).
  - New categories: `mutant_enriched` (control carries some reads) and `mutant_heterozygous`.
  - New column `af_difference` in all modes.
- **Allele check of missense variants** (`workflow/scripts/protein_check.py`).
  - The zebrafish protein is aligned to the NCBI RefSeq proteins of the same gene and to
    Compara orthologs (human, rat, chicken, *X. tropicalis*, medaka, fugu).
  - Outcomes: `matches_zebrafish_refseq`, `matches_vertebrate_consensus`, `alt_seen_in_orthologs`,
    `conserved`, `variable`, `no_orthologs` or `not_checked`.
  - Known variants and the first two outcomes are moved to the new **priority 5**, with the reason
    in the `warnings` column (`demote_flagged: false` turns this off).
  - Validation on chr7: aldh1a2 P393A → `matches_zebrafish_refseq` (NP_571925 has A);
    lgr4 I745V → `matches_vertebrate_consensus`; slc39a13 A23P → known (rs507894324).
    HERC1 Q3550H stays priority 2 (novel; `variable`).

## Re-run in sibling mode (2026-10-19; control 0038 confirmed as a sibling)
- The peak is now chr7 22–33 Mb, up to 2,582 hom variants/Mb vs a genome mean of 125.
  Inside it: 1 priority-1 variant and 43 priority-2 variants.
- The two novel missense variants that hit residues `conserved` in all aligned vertebrates were
  in the earlier novel list but not highlighted:
  - **tmppe** R370W (29,157,506; 5/5 mutant, 0/6 control; R in human, rat, chicken);
  - **TJP1** P1240L (31,126,090; 3/3 mutant, 2/8 control; P in human, rat, chicken, frog).
- Priority 1: si:dkey-1h4.4 stop-lost (30,929,805; 4/4 mutant, 1/4 control; unnamed gene).
- HERC1 Q3550H stays priority 2. The position is not conserved: zebrafish has Q, while the
  other vertebrates all have E.
- aldh1a2 P393A, mmp14a N338D and b3gat3 S94A are all priority 5 (known variants; aldh1a2 and
  mmp14a also match zebrafish RefSeq).

## Homozygosity mapping (2026-10-27)
New pipeline step (`homozygosity_mapping`, results/homozygosity/mapping_windows.tsv,
plots/homozygosity_mapping.png). At SNVs where the mutant has ≥ 2 alt reads (depth ≥ 4) and the
sibling ≥ 1 ref read (depth ≥ 3), it measures the share of mutant reads showing the reference,
i.e. coming from another version of the chromosome, in 250 kb windows.

- **Genome-wide minimum: chr7 28.75–29.5 Mb, 1.5–1.9%** (error level). Genome typical 45%;
  next-best chromosome 16 at 12%.
- chr7 profile: 21.75–24.0 Mb 13–19%; **24.0–26.45 Mb 35–44% (two copies)**; 26.5–28.75 Mb 4–6%;
  **28.75–29.5 Mb 1.5–1.9%**; 29.5–37.5 Mb mostly 9–21%.
- **The 24.0–26.45 Mb block is real two-copy heterozygosity, not an artefact.** It has sharp edges
  (~23.97 and 26.46 Mb at 20 kb resolution) and normal depth (4.3 vs 4.0–4.1). Allele fractions are
  centred on 0.5. There are no N gaps and no large SV at the edges. Genes keep the same order on
  chr7 in GRCz12tu and GRCz13ab (NCBI Datasets). So the mutation can't be there; btr09 L360V
  (24.52 Mb) drops out, as do rgp1/b3gat3/mmp14a in that flank.
- **The mutant sample is probably a pool.** The in-between values (e.g. 10–15% at 30–33 Mb) can't
  come from one diploid fish. The minority reads are real: at 30–33.2 Mb they are 8.5% where the
  sibling carries ref, vs 0.8% where it carries alt. Allele fractions at those sites are skewed
  (mostly 1 ref read in 6–8), not centred on 0.5. A pool of a few mutants with crossovers near
  the mutation explains this. The 24–26.45 block then means all pooled mutants carry two different
  mutant-chromosome lineages there (e.g. the two carrier parents' chromosomes differ). This
  corrects the earlier guess, from the 50% allele balance alone, that the sample was one fish.
- **Best interval 28.75–29.5 Mb** contains cog4 and thap11 (known jaw genes; thap11 5'UTR
  c.-217A>T at 29,045,070), tmppe R370W, HERC1 (Q3550H, D1201E, splice +2 insertion),
  SLC38A8 C209R, acd T288A, ogfod1 V36I, trpm1a S358C and ENSDARG00000103445 (E252K, S254R).
  Other genes: dok4, dnaaf1, nrn1lb, gfod2, vrk3, fbxl8, tradd, gnao1b, mtbl, exosc6, aph1b,
  dapk2a, sh2d7, rab8b, gtf2a2, foxb1b.
- TJP1 P1240L (31.13 Mb, 12%) and si:dkey-1h4.4 stop-lost (30.93 Mb, 11%) sit in the linked region
  but outside the best interval: less likely if the sample is a pool.
- Caveats: one sample per group at ~4×. The pool inference rests on read fractions, and the number
  of fish is not recorded (ask). A misclassified wild-type embryo in a mutant pool would also raise
  the share everywhere except at the mutation itself.

### Revised next steps
1. Genotype ~20 single mutants and ~20 siblings at chr7 ~23, 25, 28, 29.2 (HERC1 Q3550H) and 31 Mb.
   All mutants should be homozygous at 29.2. Mutants should be heterozygous at 25 if the two-lineage
   explanation is right.
2. Sanger the coding exons and UTRs of cog4, thap11, tmppe and HERC1 in mutants.
3. Find out how many fish went into libraries 0038 and 0039.
