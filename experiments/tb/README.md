# TB Experiments

This folder contains *M. tuberculosis*-specific workflows and settings.

## Files

- `run_tb_gam.py`: run GAM with TB defaults
- `run_tb_lmm.py`: run LMM with TB defaults
- `run_tb_fine_analysis.py`: TB-specific fine-grained SNP analysis
- `fine_analysis_model_tb.py`: TB-specific fine analysis model
- `configs/gam_tb.yaml`: TB GAM column mapping and run settings
- `configs/lmm_tb.yaml`: TB LMM column mapping, covariates, and run settings

## Usage

```bash
python experiments/tb/run_tb_gam.py
python experiments/tb/run_tb_lmm.py
python experiments/tb/run_tb_fine_analysis.py
```
