# GenePlot
A lightweight Python tool for rendering non-overlapping gene tracks.  
  
This tool generates publication-ready gene visualizations from standard BED12 annotations. Built on Matplotlib, it enables scalable plotting of multiple genes and genomic regions, facilitating seamless integration with other high-throughput genomic data tracks.  

## Installation
Clone this repository and ensure `GenePlot.py` is in your working directory, or install via pip:
```bash
pip install .  
```

## Quick Start  
Dependencies:
* numpy>=1.23.0
* matplotlib>=3.7.0
    
**1. Load and process BED data**  
The GeneInfo class handles file validation and transcript processing. Supports finding gene annotations by gene names and/or by chromosome region
```python
import numpy as np
from GenePlot import GeneInfo
gene_info_gen = GeneInfo(bed_path='annotations.bed', collapse = True)

# Search using a list of gene names
genes = ['GENE1','GENE2','GENE3']
plotting_data = gene_info_gen.get_gene_info(gene_list = genes)

# Or, searching using genomic region
plotting_data = gene_info_gen.get_gene_info(region='Chr1:150000000-150500000')
```
A combination of gene list and region is allowed.

**2. Plot the result**  
The GenePlot class handles all the aesthetics and track logic.
```python
import matplotlib.pyplot as plt
from GenePlot import GenePlot
plotter = GenePlot()
fig, ax = plt.subplots(1, 1, figsize=[20, 5])
ax = plotter.plot_gene_list(plotting_data, ax)
plt.show()
```
<img width="2000" height="500" alt="example" src="https://github.com/user-attachments/assets/8ea0b397-45fa-4b3c-b117-9d24d05564c9" />

Adding a region of interest shading using region_of_interest()
```python
fig, ax = plt.subplots(1, 1, figsize=[20, 5])
ax = plotter.plot_gene_list(plotting_data, ax)
ax = plotter.region_of_interest(ax, 150200000, 150300000, color='green', alpha = 0.3)
plt.show()
```
<img width="2000" height="500" alt="example_roi" src="https://github.com/user-attachments/assets/14843b52-aef1-43d9-9f4e-e3658337c5db" />


## Advanced Options
### GeneInfo  
* get_gene_info(gene_list = gene_list, region = region_string)   
  
| Parameter | Default | Description |
|---|---|---|
| bed_path | None | The file path to the BED12 annotation. |
| collapse | True | Whether multiple records of the same name, usually different transcripts of the same gene, are collapsed to a single entry. |  
  
These options can be set up while initiating GeneInfo (as default) or directly provided to the function: GeneInfo.get_gene_info(bed_path = path, collapse = False)
Input coordinates follow BED12 (0-based, half-open) standards.  
When collapse=True, the tool creates a consensus gene model containing the union of all exons for a given gene name. All internal UTRs are removed for cleaner visual representation.

### GenePlot
* plot_gene_list(plotting_data, ax)
  
| Parameter | Default | Description |
|---|---|---|
| padding | 0.2 | X-axis padding relative to the total region length. |
| track_offset | 1 | Vertical distance between tracks. Use to separate genes clustered together. |
| track_min_gap | 0.05 | Minimum horizontal gap (as fraction of x-axis, e.g., 0.05 = 5%) required to share a track. Smaller value will result in fewer tracks, with higher chance of labels overlapping. |
| tick_density | 50 | Number of strand markers displayed across the visible window. |
| color | blue | Default color of gene plots. |
| global_scaling | None | The width scaling of gene models. Float value to override auto scaling for plots with multiple tracks. |
| font_size | 12 | Font sizes for all text labels. |
| range_arrow | True | Toggle the chromosomal range indicator. |
| arrow_pos | 'top' | Placement of the range indicator: 'top' or 'bottom'. |

* region_of_interest(ax, start, end)

| Parameter | Default | Description |
|---|---|---|
| color | 'orange' | Shading color for region of interest. |
| alpha | 1 | Opacity of shading, 0 is completely transparent. | 
