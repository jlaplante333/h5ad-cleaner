import scanpy as sc
import datasets
import os

class TabulaMurisSenis(datasets.GeneratorBasedBuilder):
    def _info(self):
        return datasets.DatasetInfo(
            description="Cleaned SmartSeq2 single-cell datasets annotated by Claude or OpenAI.",
            features=datasets.Features({
                "cell_id": datasets.Value("string"),
                "source_model": datasets.Value("string"),
                "gene_expression": datasets.Sequence(datasets.Value("float32")),
            }),
            supervised_keys=None,
        )

    def _split_generators(self, dl_manager):
        cleaned_dir = "data/cleaned"  # <-- THIS LINE UPDATED
        files = [os.path.join(cleaned_dir, f) for f in os.listdir(cleaned_dir) if f.endswith(".h5ad")]
        return [
            datasets.SplitGenerator(name="train", gen_kwargs={"filepaths": files})
        ]

    def _generate_examples(self, filepaths):
        idx = 0
        for path in filepaths:
            filename = os.path.basename(path).lower()
            if filename.startswith("claude_"):
                source_model = "claude"
            elif filename.startswith("openai_"):
                source_model = "openai"
            else:
                source_model = "unknown"

            adata = sc.read_h5ad(path)
            for i, row in enumerate(adata.X):
                yield idx, {
                    "cell_id": adata.obs_names[i],
                    "source_model": source_model,
                    "gene_expression": row.tolist()
                }
                idx += 1