# Heme auxotroph
A brief description of several commands used for the preparation of our manuscript ["Heme auxotrophy in abundant aquatic microbial lineages"](https://www.biorxiv.org/content/10.1101/2021.01.11.426183v1.full).

## Calculation of the heme biosynthetic pathway

1. KofamScan against the KofamKOALA database using a custom-built hal file (refer to this site; https://github.com/takaram/kofam_scan).
 
```$ ./exec_annotation -o *.KOfam.txt *.faa -p /profiles/HemeBiosynthesis.hal -k /ko_list -f mapper```

2. Combine KofamScan results into a single txt file.

```$ cat *.KOfam.txt > All.KOfam.txt```

3. Calculation of the heme biosynthetic pathway completeness.
 
```$ python3 ./KEGG_decoder_Heme.py -i All.KOfam.txt -o All.KOfam.cal.txt -v static```

The python script "KEGG_decoder_Heme.py" is a modified version of "KEGG_decoder.py" available at https://github.com/bjtully/BioData/blob/master/KEGGDecoder/.

In short, we removed the definition for all pathways from the original script and then inserted the definition for modules and combinations of heme biosynthetic pathway. In addition, we modified lines 293, 294, and 296 of our script to resolve a problem caused by underscores in genome names (accession numbers for the GTDB genomes) as follows (see https://github.com/bjtully/BioData/issues/45).  
info[0].split("\_")[0] --> info[0].rsplit("\_",1)[0]

For detailed instruction on how to use this modified script, refer to the above github repository.
