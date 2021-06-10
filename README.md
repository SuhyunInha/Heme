# Heme auxotroph
A brief description of the lines used in the manuscript ["Heme auxotrophy in abundant aquatic microbial lineages"](https://www.biorxiv.org/content/10.1101/2021.01.11.426183v1.full).

## Calculation of the heme biosynthetic pathway

1. KofamScan against the KofamKOALA database (refer to this site; https://github.com/takaram/kofam_scan).
 
```$ ./exec_annotation -o *.KOfam.txt *.faa -p /profiles/HemeBiosynthesis.hal -k /ko_list -f mapper```

2. Combine KofamScan results into a single txt file.

```$ cat *.KOfam.txt > All.KOfam.txt```

3. Calculate the completeness of the heme biosynthetic pathway.
 
```$ python3 ./KEGG_decoder_Heme.py -i All.KOfam.txt -o All.KOfam.cal.txt -v static```

The python script "KEGG_decoder_Heme.py" is a modifed version of "KEGG_decoder.py" available at https://github.com/bjtully/BioData/blob/master/KEGGDecoder/. 
For detailed instruction on how to use this modified script, refer to the above github page.
Aside from the parts related to heme biosynthetic pathway, we modified lines 293, 294, and 296 to avoid errors caused by underscores in locus names as follows; info[0].split("\_")[0] --> info[0].rsplit("\_",1)[0]
(see https://github.com/bjtully/BioData/issues/45)
