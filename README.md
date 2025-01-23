# NCBI_AR_Scripts
Scripts for downloading/editing NCBI data related to AR genes

Saute_Remover.py: Removed the Saute (guided) portion of an assembly downloaded from NCBI, as it sometimes has redundant genes.

NCBI_AMR_DB_Maker_Exe.py: Makes an AR database with the resistance types, family, and alleles listed as the gene IDs from the AMRFinder database. Requires a the AMRFinder database (i.e., AMR_CDS.fa) as the ReferenceGeneCatalog.txt as inputs.
