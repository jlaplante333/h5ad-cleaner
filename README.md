# h5ad-cleaner

Hackathon for Longevity X AI at the Frontier Building!! By Jonathan and team. This is a pipeline to clean, annotate, and document aging biology `.h5ad` datasets for upload to Hugging Face.

h5ad cleaner to clean annotate and process genetic aging data in h5ad files from cellxgene.cziscience created for the Longevity X AI hackathon June 14th at the Frontier Building!!!

---
pretty_name: "Tabula Muris Senis - Bladder"
task_categories:
  - single-cell transcriptomics
license: CC-BY-4.0
---

### Dataset Summary

This is SmartSeq2 single-cell data of mouse bladder from Tabula Muris Senis.

### Usage

```python
from datasets import load_dataset
ds = load_dataset("longevity-db/tabula")
