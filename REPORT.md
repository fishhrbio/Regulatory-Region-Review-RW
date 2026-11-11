# Model and scientific report

**Waite Chromatin Atlas — Robert Waite, M.S. candidate in Biotechnology, USF**

## Analytical model

The inputs are two published segmentation functions, `g(x)` and `k(x)`, defined
only on each track's covered hg19 chr22 coordinates. At every common-covered
position they assign an integer state from 1 to 15. The primary eligible set is

`E = {x : g(x) and k(x) are defined, g(x) ∉ {14,15}, k(x) ∉ {14,15}}`.

For state-class set `C`, define `A_C = {x ∈ E : g(x) ∈ C}` and
`B_C = {x ∈ E : k(x) ∈ C}`. We report `J_C = |A_C ∩ B_C| / |A_C ∪ B_C|`,
along with the three exclusive union components and directional retention
`|A_C ∩ B_C| / |A_C|` and `|A_C ∩ B_C| / |B_C|`. Empty denominators are undefined.
The descriptive effect is `J_promoter − J_enhancer`; no significance threshold
or fitted parameter is involved.

An O(n+m) two-pointer intersection enumerates nonoverlapping atomic intervals
whose two state assignments remain constant. Counting their exact integer
lengths is equivalent to counting every base without constructing a genome-sized
per-base table. A 15×15 paired-state count matrix yields state-coverage marginals.
Coverage fractions sum to one separately for each cell line. All definitions
live in `config.json`; low-signal state 13 is not treated as missing.

## Input audit and inclusion flow

| Stage | GM12878 | K562 |
|---|---:|---:|
| Whole-file BED records | 571,339 | 622,257 |
| Records on other chromosomes (not analyzed) | 559,842 | 609,411 |
| chr22 records | 11,497 | 12,846 |
| chr22 covered bp | 34,852,028 | 34,852,028 |
| chr22 internal gaps (bp) | 342,506 | 342,506 |
| chr22 internal gap count | 10 | 10 |
| chr22 duplicate/overlapping records | 0 / 0 | 0 / 0 |
| Malformed/missing records found | 0 | 0 |

Selected coordinates span [16,050,000, 51,244,534). Uncovered terminal regions
are outside this span, not counted as low signal. Within-track overlap and
duplicate checks are chr22-specific; structural parsing covers full files.
There are 21,870 common atoms covering 34,852,028 bp. Removing states 14/15 in
either track excludes 3,384,396 bp (9.71% of common coverage). The remaining
18,523 atoms cover 31,467,632 bp. Gaps have no state and are excluded entirely.

## Main result and sensitivity

| Territory | Enhancer definition | Promoter J | Enhancer J | Difference (pp) |
|---|---|---:|---:|---:|
| Exclude states 14/15 | Strong, 4+5 | 0.473966 | 0.073200 | 40.08 |
| Exclude states 14/15 | Broad, 4–7 | 0.473966 | 0.143112 | 33.09 |
| Retain states 14/15 | Strong, 4+5 | 0.437027 | 0.068056 | 36.90 |
| Retain states 14/15 | Broad, 4–7 | 0.437027 | 0.132916 | 30.41 |

Within the primary territory, 57.8% of GM12878 active-promoter bases and 72.4%
of K562 active-promoter bases are also active promoter in the other line. For
strong enhancers, the corresponding proportions are 11.4% and 16.9%. These
directional fractions differ from Jaccard because their denominators differ.
Exactly 50.95% of eligible bases have the identical numbered state in the two
lines; the paired-state map is exploratory and is not a temporal transition map.

## Biological interpretation

Across 31,467,632 common eligible bases of hg19 chromosome 22, active-promoter
annotations share 265,800 of 560,800 union bases (Jaccard 0.474), whereas
strong-enhancer annotations share 123,400 of 1,685,802 union bases (0.073).
Thus the enhancer class is less shared in this two-line comparison. The
40.08-percentage-point gap remains positive when repetitive/CNV territory is
retained and when the enhancer class includes weak/poised states.

A biological interpretation is that candidate enhancer usage is more dependent
on cellular context than active-promoter annotation in these lines. This is
consistent with the established chromatin-state framework, rather than evidence
of a new mechanism. A competing explanation is that differences in histone
signals, sample quality, transformation or model assignments contribute to the
observed discordance. GM12878 and K562 differ in lineage and disease context;
they do not provide a controlled differentiation experiment.

The analysis is confined to one chromosome, one annotation release, hard state
assignments and a pair of cell lines. It cannot establish enhancer activity,
target genes, temporal switching or remodeling rates. Bases are correlated, and
there are no independent biological replicates for inferential testing. A
proposed next experiment would compare selected shared and line-specific
candidate enhancers in matched functional assays in both contexts. No such
experiment was performed here.

## Software validation

Analytical tests include half-open coordinate examples, zero-length overlap at
touching boundaries, invalid-input rejection, state-class merging, empty unions,
mask symmetry, conservation, and 200 randomized pairs compared with an independent
per-base oracle. Real-file checks compare exact class counts against per-base
reference maps in 100 deterministic 10-kb windows. These are software correctness
checks, not biological replication or predictive validation.

Figure captions are in `README.md` and the notebook. The 15 state labels follow
[UCSC documentation](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeBroadHmm).
Source provenance, paper citations and data terms are in `data/SOURCES.md`.
