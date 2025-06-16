import os
import openai
import numpy as np
import re
from ontology import map_cell_ontology
from metadata import add_metadata
import anthropic
import logging

# Set up logging
logging.basicConfig(filename='llm_cleaning.log', level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

def clean_adata(adata, method="code", **kwargs):
    import scanpy as sc
    if method == "code":
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return adata
    elif method == "openai":
        print("[INFO] Calling OpenAI for cleaning...")
        logging.info("Calling OpenAI for cleaning...")
        return clean_with_openai(adata, **kwargs)
    elif method == "claude":
        print("[INFO] Calling Claude for cleaning...")
        logging.info("Calling Claude for cleaning...")
        return clean_with_claude(adata, **kwargs)
    elif method == "gemini":
        return clean_with_gemini(adata, **kwargs)
    elif method == "anthropic":
        return clean_with_anthropic(adata, **kwargs)
    else:
        raise ValueError(f"Unknown cleaning method: {method}")


def clean_with_openai(adata, **kwargs):
    import scanpy as sc
    # Extract a sample of 100 cells (rows)
    sample_size = min(100, adata.shape[0])
    sample_indices = np.random.choice(adata.shape[0], sample_size, replace=False)
    sample = adata[sample_indices]
    # Convert sample to a DataFrame (obs + a few gene columns)
    sample_df = sample.obs.copy()
    # Add a few gene expression columns (up to 5 genes for brevity)
    gene_cols = sample.var_names[:5]
    for gene in gene_cols:
        sample_df[gene] = sample[:, gene].X.toarray().flatten() if hasattr(sample[:, gene].X, 'toarray') else sample[:, gene].X.flatten()
    # Convert to JSON
    sample_json = sample_df.to_json(orient="split")
    # Prepare prompt
    prompt = (
        "You are an expert in genetic data analysis. "
        "Given the following sample of genetic data from an aging study (in JSON format), "
        "suggest cleaning steps and annotation strategies for the full dataset. "
        "Sample data: " + sample_json
    )
    # Load OpenAI API key
    openai_api_key = os.getenv("OPENAI_API_KEY")
    if not openai_api_key:
        logging.error("OPENAI_API_KEY not set in environment.")
        raise RuntimeError("OPENAI_API_KEY not set in environment.")
    openai.api_key = openai_api_key
    try:
        # Call OpenAI (using GPT-3.5-turbo)
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "system", "content": "You are a helpful assistant for genetic data cleaning."},
                      {"role": "user", "content": prompt}],
            max_tokens=500
        )
        suggestions = response.choices[0].message.content
        print("\n--- OpenAI Cleaning/Annotation Suggestions ---\n", suggestions, "\n--- End Suggestions ---\n")
        logging.info(f"OpenAI suggestions: {suggestions}")
        # Parse and apply common cleaning actions
        adata = apply_cleaning_suggestions(adata, suggestions)
        print("[INFO] Finished applying OpenAI suggestions.")
        logging.info("Finished applying OpenAI suggestions.")
        return adata
    except Exception as e:
        print(f"[ERROR] OpenAI call failed: {e}")
        logging.error(f"OpenAI call failed: {e}")
        return adata


def apply_cleaning_suggestions(adata, suggestions):
    import scanpy as sc
    # Lowercase for easier matching
    text = suggestions.lower()
    # 1. Filter cells by min_genes
    match = re.search(r"remove cells with fewer than (\d+) genes", text)
    if match:
        min_genes = int(match.group(1))
        print(f"Applying: filter_cells(min_genes={min_genes})")
        sc.pp.filter_cells(adata, min_genes=min_genes)
    # 2. Filter genes by min_cells
    match = re.search(r"remove genes expressed in fewer than (\d+) cells", text)
    if match:
        min_cells = int(match.group(1))
        print(f"Applying: filter_genes(min_cells={min_cells})")
        sc.pp.filter_genes(adata, min_cells=min_cells)
    # 3. Normalize total
    if "normalize gene expression" in text or "normalize total" in text:
        print("Applying: normalize_total(target_sum=1e4)")
        sc.pp.normalize_total(adata, target_sum=1e4)
    # 4. Log1p transform
    if "log-transform" in text or "log1p" in text or "log transform" in text:
        print("Applying: log1p()")
        sc.pp.log1p(adata)
    # Warn for unrecognized suggestions
    recognized = any(["remove cells with fewer than" in text,
                      "remove genes expressed in fewer than" in text,
                      "normalize gene expression" in text,
                      "normalize total" in text,
                      "log-transform" in text,
                      "log1p" in text,
                      "log transform" in text])
    if not recognized:
        print("Warning: No recognized cleaning actions found in OpenAI suggestions.")
    return adata


def clean_with_gemini(adata, **kwargs):
    # TODO: Implement Gemini-based cleaning and annotation
    print("🔗 Cleaning with Gemini (not yet implemented)")
    return adata

def clean_with_anthropic(adata, **kwargs):
    # TODO: Implement Anthropic-based cleaning and annotation
    print("🔗 Cleaning with Anthropic (not yet implemented)")
    return adata

def clean_with_claude(adata, **kwargs):
    import scanpy as sc
    # Extract a sample of 100 cells (rows)
    sample_size = min(100, adata.shape[0])
    sample_indices = np.random.choice(adata.shape[0], sample_size, replace=False)
    sample = adata[sample_indices]
    sample_df = sample.obs.copy()
    gene_cols = sample.var_names[:5]
    for gene in gene_cols:
        sample_df[gene] = sample[:, gene].X.toarray().flatten() if hasattr(sample[:, gene].X, 'toarray') else sample[:, gene].X.flatten()
    sample_json = sample_df.to_json(orient="split")
    prompt = (
        "You are an expert in genetic data analysis. "
        "Given the following sample of genetic data from an aging study (in JSON format), "
        "suggest cleaning steps and annotation strategies for the full dataset. "
        "Sample data: " + sample_json
    )
    anthropic_api_key = os.getenv("ANTHROPIC_API_KEY")
    if not anthropic_api_key:
        logging.error("ANTHROPIC_API_KEY not set in environment.")
        raise RuntimeError("ANTHROPIC_API_KEY not set in environment.")
    client = anthropic.Anthropic(api_key=anthropic_api_key)
    try:
        response = client.messages.create(
            model="claude-3-opus-20240229",
            max_tokens=500,
            messages=[{"role": "user", "content": prompt}]
        )
        suggestions = response.content[0].text if hasattr(response.content[0], 'text') else str(response)
        print("\n--- Claude Cleaning/Annotation Suggestions ---\n", suggestions, "\n--- End Suggestions ---\n")
        logging.info(f"Claude suggestions: {suggestions}")
        adata = apply_cleaning_suggestions(adata, suggestions)
        print("[INFO] Finished applying Claude suggestions.")
        logging.info("Finished applying Claude suggestions.")
        return adata
    except Exception as e:
        print(f"[ERROR] Claude call failed: {e}")
        logging.error(f"Claude call failed: {e}")
        return adata

def annotate_adata(adata, method="code", name=None, **kwargs):
    if method == "code":
        adata = map_cell_ontology(adata)
        adata = add_metadata(adata, name=name, source_url="https://cellxgene.cziscience.com", description="Cleaned dataset from hackathon")
        return adata
    elif method == "openai":
        return annotate_with_openai(adata, name=name, **kwargs)
    else:
        raise ValueError(f"Unknown annotation method: {method}")


def annotate_with_openai(adata, name=None, **kwargs):
    # TODO: Implement OpenAI-based annotation
    print("🔗 Annotating with OpenAI (not yet implemented)")
    return adata
