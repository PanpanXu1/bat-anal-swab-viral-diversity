# 06 RdRp Sequence Diversity and Host-Virus Network

This block corresponds to the RdRp-focused manuscript section and related supplementary figures. It contains code-supported summaries for viral-order abundance, RdRp amino-acid identity and host-virus network centrality.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_viral-order-butterfly-plot/` | Summarizes viral-order abundance patterns used in the RdRp sequence-diversity context. |
| 2 | `02_rdrp-amino-acid-identity-barplot/` | Summarizes RdRp amino-acid identity categories by viral order. |
| 3 | `03_rdrp-contig-host-virus-network-centrality/` | Generates the retained host-virus network centrality leaderboard from a prepared centrality table. |

Exact final panel placement is intentionally not hard-coded here because figure assembly may change.

## Method Rationale

This block contains scripts for the tabular-to-figure components that can be reproduced directly from released summary inputs. The repository keeps the reproducible plotting steps for viral-order abundance, amino-acid identity categories and centrality summaries, with each plotted value traceable to a released input table.
