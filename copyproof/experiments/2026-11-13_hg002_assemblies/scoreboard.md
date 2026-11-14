| Assembly | What it is | Size (Gb), sequences | Discordant with the truth (Mb) | Assessable discordances found | High-confidence issues confirmed | Low-confidence issues confirmed | QV by truth k-mers / by read k-mers |
|---|---|---|---|---|---|---|---|
| `verkko_2023` | Verkko HiFi + ONT, trio (2023); names carry no haplotype | 6.08, 3,258 | 149.6 | 3,802 of 4,937 (77%) | 3,380 of 3,390 (100%) | 429 of 551 (78%) | 44.3 / 42.03 |
| `hprc_year1` | HPRC year-1 assembly (hifiasm trio, 2021 data) | 6.02, 1,055 | 59.7 | 290 of 681 (43%) | 254 of 260 (98%) | 50 of 88 (57%) | 52.4 / 49.72 |
| `hprc_curated` | HPRC curated HG002 (accession names; one group) | 6.02, 869 | 59.5 | 273 of 667 (41%) | 233 of 242 (96%) | 60 of 179 (34%) | 52.6 / 49.82 |
| `hifiasm_2021` | hifiasm 0.12, HiFi only (2021) | 5.97, 1,030 | 53.4 | 350 of 838 (42%) | 225 of 234 (96%) | 85 of 135 (63%) | 50.8 / 48.43 |
| `verkko_asm_thic` | Verkko, trio + Hi-C | 5.97, 90 | 23.6 | 89 of 158 (56%) | 77 of 81 (95%) | 5 of 52 (10%) | 55.8 / 50.07 |
| `lc24_nopolish` | ONT-only (LC24), unpolished | 6.06, 891 | 20.6 | 108 of 130 (83%) | 80 of 82 (98%) | 21 of 69 (30%) | 47.0 / 45.23 |
| `hifiasm_ontonly_2025` | hifiasm, ONT only (2025) | 5.98, 163 | 18.8 | 76 of 159 (48%) | 71 of 80 (89%) | 4 of 57 (7%) | 45.7 / 43.89 |
| `hprc_year2_v2_thic` | HPRC year-2 v2, trio + Hi-C | 5.98, 209 | 17.9 | 73 of 155 (47%) | 69 of 70 (99%) | 10 of 45 (22%) | 52.7 / 48.79 |
| `hprc_year2` | HPRC year-2 pipeline, hifiasm trio, polished | 5.98, 302 | 16.5 | 151 of 233 (65%) | 114 of 119 (96%) | 31 of 74 (42%) | 51.4 / 49.21 |
| `v2_trio` | Verkko v2, trio | 6.00, 678 | 15.7 | 71 of 178 (40%) | 69 of 74 (93%) | 7 of 57 (12%) | 50.6 / 49.54 |
| `hg002v0.7` | Q100 draft v0.7 | 6.00, 48 | 2.2 | 2 of 12 (17%) | 4 of 7 (57%) | 2 of 38 (5%) | 56.9 / 51.18 |
| `hg002v0.9` | Q100 draft v0.9 | 6.00, 48 | 2.1 | 0 of 5 (0%) | 3 of 6 (50%) | 0 of 32 (0%) | 58.7 / 51.78 |
| `hg002v1.0.1` | Q100 v1.0.1 | 6.00, 48 | 0.2 | 0 of 2 (0%) | 0 of 3 (0%) | 0 of 43 (0%) | 59.3 / 51.95 |

Summed over the 13 assemblies, by size of the segment:

| Direction | Size | Assessable discordances | Found | Issues | Confirmed | High-confidence issues | Confirmed |
|---|---|---|---|---|---|---|---|
| The truth holds more (collapse) | under 5 kb | 438 | 41 (9%) | 26 | 19 (73%) | 0 |  |
|  | 5 to 20 kb | 537 | 133 (25%) | 155 | 141 (91%) | 140 | 129 (92%) |
|  | 20 to 100 kb | 431 | 215 (50%) | 229 | 214 (93%) | 217 | 207 (95%) |
|  | over 100 kb | 158 | 85 (54%) | 151 | 136 (90%) | 146 | 133 (91%) |
| The truth holds fewer (surplus copy) | under 5 kb | 367 | 65 (18%) | 702 | 26 (4%) | 4 | 4 (100%) |
|  | 5 to 20 kb | 3,515 | 2,778 (79%) | 2,555 | 2,527 (99%) | 2,135 | 2,123 (99%) |
|  | 20 to 100 kb | 2,213 | 1,605 (73%) | 1,856 | 1,844 (99%) | 1,643 | 1,633 (99%) |
|  | over 100 kb | 496 | 363 (73%) | 394 | 376 (95%) | 363 | 350 (96%) |

Sums the record quotes:

- High-confidence copy-number issues confirmed: 4,579 of 4,648 (99%); without `verkko_2023`, 1,199 of 1,258 (95%).
- Per assembly, among the 10 with 20 such issues or more: 89 to 100% confirmed.
- The 3 with fewer than 20 (`hg002v0.7`, `hg002v0.9`, `hg002v1.0.1`): 7 of 16 (44%) confirmed.
- High-confidence issues not confirmed, per assembly: 1 to 10 (69 in all), while the confirmed ones range from 0 to 3,380.
- Found, of assessable discordances of 5 kb and more: surplus copies 73 to 79%, collapses 25 to 54% (the three size classes).
- Under 5 kb: found 9% (collapses) and 18% (surplus copies); 26 of 702 (4%) surplus-copy issues confirmed.
- Segments that fit a depth factor or nothing regular (not reported): 1,310 of 35,372 (4%) confirmed; without `verkko_2023`, 448 of 31,937 (1%).
- Ambiguous segments without support (not reported): 179 of 1,465 (12%) confirmed.
- Sequence discordant with the truth: 0.2 to 149.6 Mb.
- QV by truth k-mers minus QV by read k-mers: 1.1 to 7.3; highest QV by read k-mers 51.95.
- Unread assembly k-mers that are wrong by the truth: 9 to 90%; k-mers wrong by the truth that the reads find: 32 to 92%.
