# Dataset provenance and licenses

## Public Visium validation example

- Dataset: 10x Genomics Human Lymph Node, Space Ranger 1.0.0.
- Landing page: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-0-0
- Files: filtered feature-barcode HDF5 matrix, spatial archive, and metrics CSV.
- License: Creative Commons Attribution 4.0 International (CC BY 4.0), as stated on the dataset page.
- Download date: 2026-07-20.
- Expected input: 4,039 tissue-associated spots. The validated SpatialAtlas output contains 4,037 spots after the default QC thresholds.

## Visium HD validation example

- Dataset: female mouse adrenal gland, Space Ranger output, 8 µm binned matrix.
- Associated publication: Blatkiewicz M, Hryhorowicz S, Szyszka M, et al. *Single-cell spatial transcriptomics reveals sex-dependent gene expression and intercellular signalling in mouse adrenal cortex*. doi:10.1038/s42003-026-10697-9.
- Included input components: filtered HDF5 matrix, Parquet tissue positions, scale-factor JSON, and high-resolution tissue image.
- Input profiles matched to spatial coordinates: 250,669.
- Validated output after default QC: 143,112 profiles and 15,109 genes.
- Redistribution and reuse of this example must follow the provenance and consent terms of the underlying adrenal-gland study. The package does not assert a new license for these biological data.

## Classic Visium adrenal implementation benchmark

- Dataset: four mouse adrenal sections from the authors' classic Visium CytAssist study, represented by 992 tissue-covered spots and 19,465 genes in the retained application object.
- Associated publication: Blatkiewicz M, et al. *Molecular landscape of the mouse adrenal gland and adjacent adipose tissue by spatial transcriptomics*. doi:10.5603/fhc.108988.
- Repository record: GEO `GSE283302`.
- Included source components: the authors' retained Seurat/Shiny application and normalized Seurat RDS, plus a deterministic export script, matched H5AD and extracted histology.
- The matched SpatialAtlas H5AD retains the normalized Seurat `data` layer and is not normalized a second time.
- The comparison also reuses the 143,112-bin, 8 µm Visium HD adrenal output from `examples/02_visium_hd_8um`, associated with doi:10.1038/s42003-026-10697-9. A matched sparse RDS is built directly from the H5AD arrays for a purpose-built R/Shiny reference application; expression values, coordinates and histology are unchanged.
- Reuse remains subject to the originating studies, GEO record, ethics approvals, and institutional data-use terms. The software license does not assign a new license to the biological material.
- File identities and SHA-256 digests are recorded in `examples/05_adrenal_platform_benchmark/provenance.json` and `example-data.sha256`.

Cryptographic checksums are recorded in `checksums.sha256`.
