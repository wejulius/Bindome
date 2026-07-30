# Target Preprocessing

Pipeline from a UniProt organism ID to a filtered `sweep_index.tsv` for `../SweepScripts`.
Each step is a `sbatch NN_run_*.sh` wrapper; modify the conda environment name and pass flags to the underlying `.py`/`.sh` to override defaults.

1. **`01_run_fetch_uniprot_info.sh`** — fetch reviewed UniProt entries + annotations. → `Data/UniProt/uniprot_reviewed_annotated.csv`
2. **`02_run_fetch_afdb_pdbs_and_pae.sh`** — download AlphaFold PDBs + PAE matrices. → `Data/AFDB_v6/{PDB,PAE}/`
3. **`03_run_create_domains_from_pae.sh`** — call domains from PAE (via `pae_to_domains/`). → `Data/AFDB_v6/CSV_Domains/*.domains.csv`
4. **`04_run_segment_and_compute_domain_stats.sh`** — segment per-domain PDBs + compute stats (pLDDT, DSSP, Rg, contacts, membrane, AF-Bind). Needs a working DSSP binary. → `Data/AFDB_v6/PDB_Domains/`, `Data/AFDB_v6/Statistics/domain_stats.csv`
   - **optional `04b_generate_idp_crops.py --design-root <path>`** — crop disordered (IDP) regions for BindCraft. → `Data/AFDB_v6/{PDB_IDP_Crops,Settings_IDP}/`, `idp_crops_stats.csv`, `idp_bindcraft_jobs.txt`
5. **`05_run_filter_domains_and_create_sweep_index.sh`** — filter by quality thresholds, write sweep index (add `--idp` for crops). → `Data/Filtered/filtered_domain_stats.csv`, `sweep_index_batch{1..4}.tsv`

Step 5's output feeds `../SweepScripts/01_create_sweep_dirs.sh`.
