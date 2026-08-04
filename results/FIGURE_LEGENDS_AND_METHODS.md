# Figure legends and Methods

Analyses use ANDV Clade III tip-label alignments under `06_v5_label_aln` (plus R600 where indicated). Pairwise p-distances ignore sites with a gap or `N` in either sequence.

---

## Figure legends

**Supplementary Fig. S3. Distance and topological analysis hardly detect evidence of reassortment among Clade III genomes.**  
**(a–c)** Pairwise p-distances between Clade III viruses (including R600) for (a) S vs M, (b) S vs L and (c) M vs L. The dashed line is *y* = λ*x*, where λ is the ML mutation-rate ratio between the indicated segments estimated from non-R600 pairs. Points outlined in black deviate from this expectation at Benjamini–Hochberg *q* < 0.05. Most significant outliers involve R600 (green); among non-R600 pairs, significant deviations mainly involve p1236 and/or p1059 and are examined by sliding-window analysis in (d,e).  
**(d,e)** Sliding-window pairwise p-distances (100-bp window, 10-bp stride) for (d) p1236 vs p1239 and (e) p1059 vs NRC-4/18. Dashed lines indicate overall segment p-distances. In (d), L-segment distance is elevated toward the 3′ half relative to the 5′ half. In (e), the S segment of p1059 is nearly identical to NRC-4/18 (overall p-distance = 0.0016), whereas M and L remain at typical within-clade distances (~0.01).

---

## Results / Discussion (draft)

Overall, pairwise distances among Clade III segments are largely concordant, arguing against frequent reassortment within this clade (Supplementary Fig. S3). Most significant outliers involve R600 and likely reflect its atypical evolutionary history rather than reassortment among human isolates. The remaining signals involve p1236 and p1059. For p1236, elevated L distance is concentrated in the second half of the segment, consistent with local rate heterogeneity rather than a segment swap. For p1059, S is nearly identical to NRC-4/18 while M and L are not, which could indicate recent shared S ancestry but lacks the broader segmental discordance expected for classical reassortment. Together, these patterns provide little evidence that reassortment has shaped Clade III genome relationships.

---

## Methods (for Supplementary Fig. S3)

### Pairwise p-distances
Unless otherwise noted, pairwise nucleotide p-distance between two aligned sequences was calculated as the number of differing sites divided by the number of comparable sites, where a site was considered comparable only if both sequences carried an unambiguous nucleotide (A, C, G or T). Sites containing a gap or ambiguous base (`N`) in either sequence were excluded.

### Detection of segment-discordant distances (Supplementary Fig. S3a–c)
To test whether pairwise genetic distances were concordant across genomic segments, all pairwise p-distances were computed among Clade III viruses for the S, M and L segments, including pairs involving R600. For each isolate pair and each segment comparison (*X*, *Y*) ∈ {(S, M), (S, L), (M, L)}, the number of observed nucleotide differences (*Y<sub>X</sub>*, *Y<sub>Y</sub>*) and the number of comparable sites (*E<sub>X</sub>*, *E<sub>Y</sub>*) were recorded.

Substitution counts were modelled as independent Poisson random variables sharing an isolate-pair rate *r* and a constant segment rate ratio λ:

*Y<sub>X</sub>* ~ Poisson(*E<sub>X</sub>* *r*),  *Y<sub>Y</sub>* ~ Poisson(*E<sub>Y</sub>* λ *r*).

Conditioning on the total number of differences *T* = *Y<sub>X</sub>* + *Y<sub>Y</sub>* eliminates *r* and yields

*Y<sub>Y</sub>* | *T* ~ Binomial(*T*, π),  π = λ *E<sub>Y</sub>* / (*E<sub>X</sub>* + λ *E<sub>Y</sub>*).

The rate ratio λ was estimated by conditional binomial maximum likelihood from non-R600 pairs with expected *Y<sub>Y</sub>* count *T*π > 20, and the fitted relationship was displayed on the p-distance scatter as the line *y* = λ*x*. For each pair meeting the same *T*π > 20 filter, departure from the null expectation under the fitted λ was assessed with a two-sided normal approximation,

*Z* = (*Y<sub>Y</sub>* − *T*π) / √[*T* π (1 − π)],

with *P* = 2 Φ(−|*Z*|). Within each segment comparison, *P*-values were corrected for multiple testing by the Benjamini–Hochberg procedure, and pairs with *q* < 0.05 were considered significant.

### Sliding-window pairwise distances (Supplementary Fig. S3d–e)
To visualise local variation in genetic distance along each segment, pairwise p-distances were calculated in a sliding window of 100 consecutive comparable sites with a stride of 10 comparable sites for selected isolate pairs (p1236 vs p1239; p1059 vs NRC-4/18). Window p-distance was defined as the number of differences within the window divided by 100. Overall segment p-distances were computed from all comparable sites as described above and are indicated as dashed horizontal lines in each panel.