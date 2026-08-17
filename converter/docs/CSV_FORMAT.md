# CSV format reference

The maintained, detailed format and validation specification is in
[`DATASET_EXPORTER.md`](DATASET_EXPORTER.md).

The three original core schemas are preserved exactly:

```text
bus.csv:  bus_id,phases,base_kv,type
ybus.csv: row,col,g,b
load.csv: bus_id,phase,scenario,p_pu,q_pu
```

The Python exporter adds eight compact supporting CSV files plus
`validation.json` and verifies `YV=I` after reading the serialized files back
from disk.

`ybus.csv`, `voltage.csv`, and `load.csv` retain up to 15 decimal places.
Coordinates and equipment tables retain up to eight decimal places.
