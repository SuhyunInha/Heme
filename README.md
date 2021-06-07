# Heme

## Calculation of the heme biosynthetic pathway

1. KofamScan against the KofamKOALA database (refer to this site; https://github.com/takaram/kofam_scan).
2. 
```$ ./exec_annotation -o *.KOfam.txt *.faa -p /profiles/HemeBiosynthesis.hal -k /ko_list -f mapper```

2. Combine KofamScan results into one txt file.
3. 
```cat *.KOfam.txt > All.KOfam.txt```

3. Calculate the completeness of the heme biosynthetic pathway.
4. 
```$ python3 ./KEGG_decoder_Heme.py -i All.KOfam.txt -o All.KOfam.cal.txt -v static```

The script "KEGG_decoder_Heme.py" is a modifed version of "KEGG_decoder.py" available at https://github.com/bjtully/BioData/blob/master/KEGGDecoder/. 
For detailed instruction on how to use this modified script, refer to the above github page.
Aside from the parts related to heme biosynthetic pathway, we modified the line 293, 294, and 296 of the script to avoid the errors caused by underscore in locus names; info[0].split("\_")[0] --> info[0].rsplit("\_",1)[0]
(refer to this; https://github.com/bjtully/BioData/issues/45)
