def add_metadata(adata, name, source_url, description):
    adata.uns['dataset_name'] = name
    adata.uns['source'] = source_url
    adata.uns['description'] = description
    return adata



# scripts/metadata.py

def add_metadata(adata, name="unknown", source_url=None, description=None):
    print("🧾 Adding metadata...")
    adata.uns["dataset_name"] = name
    adata.uns["source_url"] = source_url
    adata.uns["description"] = description
    return adata
