# antiSMASH Collator Readme
---
## Install

Download the zip containing the antiSMASH Collator scripts and unpack
```
wget antismash_collator.tar.gz
tar -xvzf antismash_collator.tar.gz
```

Pre-requisites are:
- Singularity >/ 3.5.2
- Singularity image with antiSMASH >/ 5.1 (only for rebuilding output folder)


## Running

To inititate a run, run the 'antismash_collator.py' script as below
```
python antismash_collator.py /path/to/antismash/output/folders job_name /path/to/antismash/singularity/image
```

Calling the `-r, --rebuild` argument will run the utility script `antismash_rebuild.sh` rebuild the antiSMASH output folder using a supplied antiSMASH Singularity image.

## Help & Suggestions

For help, view the help text for 'antismash_collator.py'
```
antismash_collator.py -h
```

The JSON parse functions were written by Jeremy Owen.
antiSMASH Collator is maintained by Kelly Styles.
email **kelly.styles@vuw.ac.nz** for help or suggestions.
