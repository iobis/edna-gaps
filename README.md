# edna-gaps

This is an R package for regional eDNA reference database gap analysis based on OBIS species distributions.

## Reference databases

Reference database creation is documented at <https://github.com/iobis/eDNA_trial_data>, and reference databases were downloaded to `reference_databases` from the LifeWatch server using:

```
rsync -avrm --include='*/' --include='*pga_tax.tsv' --include='*pga_taxa.tsv' --include='*pga_taxon.tsv' --exclude='*' ubuntu@lfw-ds001-i035.lifewatch.dev:/home/ubuntu/data/databases/ ./reference_databases
rm -r ./reference_databases/silva*
```

## Regional species lists

Regional species lists are created using `reference_species.R` and saved to `reference_databases/reference_species.csv.gz`.

## Gap analysis

Species lists by H3 cell of resolution 3 whith columns indicating presence in reference databases are written to `species_lists/species_lists.csv.gz`.
