# ensemblrestpy package
```python
from ensembl import Ensembl

ensembl = Ensembl()
print(ensembl.variant_recoder("NP_001361433.1:p.Asp512Asn", vcf_string=True, var_synonyms=True, fields="id", format="xml"))
```