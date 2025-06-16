ontology_map = {
    "T cell": "CL:0000084",
    "B cell": "CL:0000236",
    # ...
}


# scripts/ontology.py

def map_cell_ontology(adata):
    # Placeholder: This should map cell types to ontology IDs
    print("📚 Mapping cell ontology (placeholder)")
    if "cell_type" in adata.obs:
        adata.obs["cell_type_ontology"] = adata.obs["cell_type"].apply(lambda x: f"CL:{hash(x) % 10000:04d}")
    else:
        adata.obs["cell_type_ontology"] = "CL:0000000"  # fallback
    return adata
