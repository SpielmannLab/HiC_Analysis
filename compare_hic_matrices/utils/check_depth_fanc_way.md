# These are some helper scripts to check several properties of a fanc file

```bash
module load singularity
singularity shell $WORK/singularity/varunkas-fanc-latest.img

cp *.hic $SCRATCH/ && cd $SCRATCH
```

```python
import os
import fanc
import numpy as np

# These totals are the values used by fanc to downsample the data. It is their proxy for number of valid pairs
# When norm = True, it should just be the number of pixels in the matrix
# When norm = False, it should be a proxy for number of valid pairs
def get_valid_pairs(fanc, norm):
   total = int(sum(e.weight for e in fanc.edges(lazy=True, norm=norm)))
   print(f"Valid Pairs: {total}")

# In a normalized matrix, the rowsums should be 1 or close to it.
def get_row_sums(fanc, chrom):
   mat = fanc.matrix((chrom, chrom))
   rowsum = np.sum(mat, axis = 0)
   print(f"Rowsums: {rowsum}")

for file in [x for x in os.listdir() if ".hic" in x ]:
   fanc_data = fanc.load(file)
   print(f"file: {file}")
   get_valid_pairs(fanc_data, norm = False)
   get_row_sums(fanc_data, chrom = "1")

```
