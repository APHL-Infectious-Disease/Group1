#!/bin/bash 
    
    # Run the SRA query and collect all IDs

    esearch -db sra -query "WGS[Strategy] AND USA AND Wastewater AND Metagenome" | efetch -format native > sra.xml


    # Build xtract parameters
    {
        echo "-pattern EXPERIMENT_PACKAGE"
        echo "-first PRIMARY_ID"
        echo "-element LIBRARY_STRATEGY"
        echo "-element LIBRARY_SELECTION"
        echo "-element TAXON_ID"
        echo "-element SCIENTIFIC_NAME"
        echo "-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"
        echo "-if TAG -equals collection_date -element VALUE"
        echo "-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"
        echo "-if TAG -equals geo_loc_name -element VALUE"
        echo "-block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"
        echo "-if TAG -equals lat_lon -element VALUE"
    } > xtract_dynamic.cmd

    # Build header row
    {
        printf "experiment_id\tlibrary_strategy\tlibrary_selection\ttaxon_id\tscientific_name\tcollection_date\tlocation\tlat_lon"
        printf "\n"
    } > sra_meta.tsv

    xtract -input sra.xml $(cat xtract_dynamic.cmd) | grep "USA" >> sra_meta.tsv

