from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict


class GenePlot:
    default_vals = {'padding' : 0.2,
                    'track_offset' : 1,
                    'color' : 'blue',
                    'tick_density' : 50,
                    'track_min_gap' : 0.05,
                    'global_scaling' : None,
                    'font_size' : 12,
                    'range_arrow' : True,
                    'arrow_pos' : 'top'}

    def __init__(self, **kwargs):
        for key, val in self.default_vals.items():
            setattr(self, key, val)
        self.set_config(**kwargs)

    def get_config(self):
        """Returns a dictionary of all current default plotting parameters."""
        # Filters out any internal methods or private attributes starting with '_'
        return {k: v for k, v in vars(self).items() if not k.startswith('_')}

    def set_config(self, **kwargs):
        """
        Updates default parameters.
        Example: plotter.set_config(color='red', track_offset=1.5)
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                # Optional: Raise warning or error for non-existent attributes
                print(f"Warning: '{key}' is not a valid plotting parameter.")
                
    def _assign_gene_tracks(self, genes_dict, min_gap = None):
        """
        Assigns genes to tracks to avoid overlaps.
        Input: {idx: [start, end], ...}
        Output: {idx: track_number, ...}
        """
        min_gap = min_gap if min_gap is not None else self.track_min_gap
        # 1. Sort gene indices by their start positions
        gene_coords = {k:[v['gene_start'],v['gene_end']] for k,v in genes_dict.items()}
        sorted_indices = sorted(gene_coords.keys(), key=lambda x: gene_coords[x][0])
        
        # 2. Track the 'end' position of the last gene placed in each track
        # track_ends[track_id] = furthest_right_coordinate
        track_ends = [] 
        assignments = {}
    
        for idx in sorted_indices:
            start, end = gene_coords[idx]
            assigned = False
            
            # 3. Try to place the gene in an existing track
            for track_id, last_end in enumerate(track_ends):
                # If the gene starts after the last gene in this track ends
                if start > last_end: 
                    assignments[idx] = track_id
                    track_ends[track_id] = end + min_gap
                    assigned = True
                    break
            
            # 4. If no track is available, create a new one
            if not assigned:
                assignments[idx] = len(track_ends)
                track_ends.append(end + min_gap)
        return assignments
    
    def plot_gene_list(self, gene_list, ax, padding=None, track_offset = None, 
                       color=None, tick_density = None, track_min_gap = None, range_arrow = None, arrow_pos = None):
        """
        Renders a collection of gene models onto a single Matplotlib axes, 
        automatically calculating non-overlapping tracks and scaling line widths.
        """
        padding = padding if padding is not None else self.padding
        track_offset = track_offset if track_offset is not None else self.track_offset
        color = color if color is not None else self.color
        tick_density = tick_density if tick_density is not None else self.tick_density
        track_min_gap = track_min_gap if track_min_gap is not None else self.track_min_gap
        range_arrow = range_arrow if range_arrow is not None else self.range_arrow
        arrow_pos = arrow_pos if arrow_pos is not None else self.arrow_pos
        # Verify all genes belong to the same chromosome to ensure a valid X-axis scale
        if len(set([g['chrom'] for g in gene_list])) > 1:
            raise ValueError('Genes from multiple chromosomes provided. Only one chromosome allowed.')
        
        # Convert list to a dictionary to maintain unique indices during track assignment
        genes_dict = dict(enumerate(gene_list))
        
        # Determine the global genomic range covered by all genes in the list
        left_bound, right_bound = min([g['gene_start'] for g in gene_list]), max([g['gene_end'] for g in gene_list])
        total_length = right_bound - left_bound
        
        # Apply padding to the range (max(..., 0) prevents negative genomic coordinates)
        left_bound = max([left_bound - padding * total_length, 0])
        right_bound = right_bound + padding * total_length
        boundary_ticks = [int(left_bound), int(right_bound)]
        ax.set_xlim([left_bound, right_bound])
        ax.set_xticks([])
        
        # Determine vertical track placement for each gene to avoid visual overlaps
        min_gap = (right_bound - left_bound) * track_min_gap
        track_assignments = self._assign_gene_tracks(genes_dict, min_gap=min_gap)
        total_tracks = len(set(track_assignments.values()))
        
        # Set vertical limits based on the number of tracks and the specified offset
        ax.set_ylim([-1, 1 + (total_tracks - 1) * track_offset])
        
        # Calculate a scaling factor for line thicknesses using exponential decay.
        # This prevents lines from being too thick when many tracks are present.
        if self.global_scaling:
            scaling = self.global_scaling
        else:
            scaling = 0.5 + 0.5 * np.exp(-0.2 * (total_tracks - 1))
        
        # Iterate through genes and draw each one at its assigned vertical position
        for i,gene in genes_dict.items():
            self.plot_gene(gene, ax, y_pos = track_assignments[i] * track_offset, color = color, 
                           scaling = scaling, inherit_xlim=True, tick_density=tick_density)
        
        # Hide Y-axis ticks as they do not represent quantitative data in this track view
        ax.set_yticks([]) 

        if range_arrow:
            ax = self._range_arrow(ax,chrom = gene_list[0]['chrom'], xmin = left_bound, xmax = right_bound, position = arrow_pos)
        
        return ax
        
    def plot_gene(self, plotting_dict, ax, color=None, y_pos = 0, scaling = 1.0, 
                  inherit_xlim = False, tick_density = None):
        """
        Renders a gene model onto a provided Matplotlib Axes object.
        """
        color = color if color is not None else self.color
        tick_density = tick_density if tick_density is not None else self.tick_density
        # Generate strand markers
        if inherit_xlim:
            left_bound, right_bound = ax.get_xlim()
            arrow_gap = (right_bound - left_bound) / tick_density
        else:
            gene_size = plotting_dict['gene_end'] - plotting_dict['gene_start']
            arrow_gap = gene_size / 15
        direction_markers_x = list(np.arange(plotting_dict['gene_start'], plotting_dict['gene_end'], arrow_gap))
        if direction_markers_x:
            direction_markers_x.pop(0)
        direction_markers_y = [y_pos for x in direction_markers_x]
        
        direction_markers_type = '4' if plotting_dict['strand'] == '+' else '3' if plotting_dict['strand'] == '-' else None
    
        # Text label placement
        gene_name = plotting_dict['gene_name']
        midpoint = (plotting_dict['gene_start'] + plotting_dict['gene_end']) / 2
        ax.text(midpoint, y_pos + 0.2, gene_name, ha='center', va='bottom', style='italic', fontsize=self.font_size)
    
        # 1. Draw Intron line
        ax.hlines(y=y_pos, xmin=plotting_dict['gene_start'], xmax=plotting_dict['gene_end'], 
                  ls='-', lw=1, color=color)
        
        # 2. Draw UTR segments
        ax.hlines(y=[y_pos] * len(plotting_dict['utr_starts']), xmin=plotting_dict['utr_starts'], 
                  xmax=plotting_dict['utr_ends'], ls='-', lw=8 * scaling, color=color)
        
        # 3. Draw CDS segments
        ax.hlines(y=[y_pos] * len(plotting_dict['cds_starts']), xmin=plotting_dict['cds_starts'], 
                  xmax=plotting_dict['cds_ends'], ls='-', lw=15 * scaling, color=color)
        
        # 4. Overlay strand markers
        if direction_markers_type:
            ax.plot(direction_markers_x, direction_markers_y, marker=direction_markers_type, 
                    ms=8, linestyle='None', color=color)
        return ax

    def _range_arrow(self, ax, chrom, xmin, xmax, position='top'):
        y = 1.05 if position == 'top' else -0.05
        
        ax.annotate(
            '', 
            xy=(xmin, y), 
            xytext=(xmax, y),
            xycoords=ax.get_xaxis_transform(), 
            textcoords=ax.get_xaxis_transform(),
            annotation_clip=False, 
            arrowprops={
                'arrowstyle': '<->',
                'mutation_scale': 20,
                'color': 'black',
                'lw': 1,
                'shrinkA': 0,
                'shrinkB': 0
            }
        )
        label_text = f'{chrom}:{int(xmin)}-{int(xmax)}'
        mid_x = (xmin + xmax) / 2
        ax.text(
            mid_x, 
            y, 
            label_text,
            transform=ax.get_xaxis_transform(), # Use the same blended transform
            ha='center',                         # Horizontal alignment: center the text on the mid_x point
            va='center',                         # Vertical alignment: center the text on the y point
            bbox=dict(                           # Bounding box for background
                facecolor='white',               # White background
                edgecolor='none',                # No border around the text
                pad=1.0                          # Add a little padding around the text
            ),
            clip_on=False,                        # Allow text to appear outside the main plot area
            fontsize = self.font_size
        )
        return ax

    
    def region_of_interest(self, ax, start, end, color = 'orange', alpha = 1):
        target_start, target_end = int(start), int(end)
        ax.axvspan(target_start,target_end,facecolor = color, alpha = alpha, zorder=0)
        return ax

class GeneInfo:
    def __init__(self, bed_path=None, collapse=True):
        self.bed_path = None
        self.collapse = collapse
        if bed_path:
            self.set_bed_path(bed_path)
    
    def set_bed_path(self, bed_path):
        """Public method to update the BED path with validation."""
        self._verify_bed_file(bed_path)
        self.bed_path = bed_path

    def _verify_bed_file(self, bed_path):
        """Internal helper to validate file existence and basic BED12 format."""
        path = Path(bed_path)
        
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}") 
            
        with open(path, 'r') as f:
            data_found = False
            # Check only first 10 data lines for performance
            count = 0
            for line in f:
                if line.startswith(('#', 'track', 'browser')) or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) != 12:
                    raise ValueError(
                        f"Format Error: BED12 requires 12 columns. "
                        f"Found {len(fields)} in line: {line[:50]}..."
                    )
                data_found = True
                count += 1
                if count >= 10: # Sample size is enough for verification
                    break
        if not data_found:
            raise ValueError(f"BED file {path} appears to be empty.")

    @staticmethod
    def _normalize_chrom(chrom_str):
        s = str(chrom_str)
        if s.lower().startswith('chr'):
            return s[3:]
        return s
    
    @staticmethod
    def _process_gene_info(gene_info):
        # This function translates standard UCSC BED12-style input into
        # a coordinate dictionary optimized for plotting UTR and CDS exons separately.
    
        # Calculate absolute genomic coordinates for all exons (blocks)
        # blockstarts in BED format are relative to gene_start
        exon_starts = np.array(gene_info['blockstarts']) + gene_info['start']
        exon_ends = exon_starts + np.array(gene_info['blocksizes'])
        
        # Define the absolute start and end coordinates of the Coding Sequence (CDS) region
        cds_s, cds_e = gene_info['thickstart'], gene_info['thickend']
        
        # Initialize lists to store the coordinates of the segmented regions
        utr_starts, utr_ends = [], []
        coding_starts, coding_ends = [], []
    
        # Iterate through each individual exon
        for e_s, e_e in zip(exon_starts, exon_ends):
            # Determine the intersection between the current exon [e_s, e_e] 
            # and the overall CDS range [cds_s, cds_e].
            
            # The coding segment starts at the later of the exon start or CDS start
            c_s = max(e_s, cds_s)
            # The coding segment ends at the earlier of the exon end or CDS end
            c_e = min(e_e, cds_e)
            
            if c_s < c_e:
                # If the calculated intersection is valid (start < end), it's a coding exon
                coding_starts.append(c_s)
                coding_ends.append(c_e)
                
                # Identify UTR segments within this exon (Parts of exon outside CDS range)
                if e_s < c_s: 
                    # Part of the exon before the CDS start (5' UTR or part of 3' UTR depending on strand)
                    utr_starts.append(e_s)
                    utr_ends.append(c_s)
                if e_e > c_e: 
                    # Part of the exon after the CDS end (3' UTR or part of 5' UTR depending on strand)
                    utr_starts.append(c_e)
                    utr_ends.append(e_e)
            else:
                # If there is no overlap (c_s >= c_e), the entire exon is non-coding (UTR)
                utr_starts.append(e_s)
                utr_ends.append(e_e)
    
        # Return a structured dictionary tailored for the downstream plotting function
        return {
            'gene_name': gene_info['gene_name'],
            'chrom':gene_info['chrom'],
            'gene_start': gene_info['start'],
            'gene_end': gene_info['end'],
            'utr_starts': utr_starts,
            'utr_ends': utr_ends,
            'cds_starts': coding_starts,
            'cds_ends': coding_ends,
            'strand': gene_info['strand']
        }
    
    def get_gene_info(self,bed_file = None,gene_list = None, region = None, collapse = None):
        bed_file = bed_file if bed_file is not None else self.bed_path
        collapse = collapse if collapse is not None else self.collapse
        if bed_file is None:
            raise ValueError('No BED file provided. Use GeneInfo.set_bed_path(path) to define path or direclty provide path using path=... argument')
        if gene_list is None and region is None:
            raise ValueError('Missing target information. Either gene list or region must be provided.')
        gene_set = set(gene_list) if gene_list else set()
        if region:
            target_chrom, coords = region.split(':')
            target_chrom = self._normalize_chrom(target_chrom)
            target_start, target_end = coords.split('-')
            target_start, target_end = int(target_start), int(target_end)
        out_dict = {}
        with open(bed_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom, start, end, gene_name = fields[0:4]
                chrom = self._normalize_chrom(chrom)
                start, end = int(start), int(end)
                if region:
                    in_region = (chrom == target_chrom) and (start <= target_end) and (end >= target_start)
                else: 
                    in_region = False
                if gene_name not in gene_set and not in_region:
                    continue
                gene_info = {'gene_name': gene_name,
                             'chrom':fields[0],
                             'start':int(fields[1]), 
                             'end':int(fields[2]), 
                             'thickstart':int(fields[6]),
                             'thickend':int(fields[7]),
                             'blockstarts':[int(x) for x in fields[11].strip(',').split(',')],
                             'blocksizes':[int(x) for x in fields[10].strip(',').split(',')],
                             'strand':fields[5]}
                processed_gene_info = self._process_gene_info(gene_info)
                if gene_name not in out_dict:
                    out_dict.update({gene_name:[processed_gene_info]})
                else:
                    out_dict[gene_name].append(processed_gene_info)
            if collapse:
                print('Collapsing transcripts for each gene. Set collapse = False to keep all transcripts separate.')
                out_list = [self._collapse_transcripts(g) for g in out_dict.values()]
            else:
                out_list = self._flatten_transcripts(out_dict)
        return out_list

    @staticmethod
    def _merge_intervals(starts, ends):
        """Standard algorithm to merge overlapping intervals."""
        if not starts: return [], []
        # Sort by start position
        intervals = sorted(zip(starts, ends))
        merged_s, merged_e = [intervals[0][0]], [intervals[0][1]]
        
        for s, e in intervals[1:]:
            if s <= merged_e[-1]: # Overlap
                merged_e[-1] = max(merged_e[-1], e)
            else: # Gap
                merged_s.append(s)
                merged_e.append(e)
        return merged_s, merged_e
    
    def _collapse_transcripts(self, transcript_list):
        """
        Collapses isoforms into a single track, removing internal UTRs.
        Only UTRs outside the global CDS start and global CDS end are kept.
        """
        if not transcript_list: return None
        
        # 1. Basic Metadata
        chrom = transcript_list[0]['chrom']
        strand = transcript_list[0]['strand']
        gene_name = transcript_list[0]['gene_name']
        g_start = min(t['gene_start'] for t in transcript_list)
        g_end = max(t['gene_end'] for t in transcript_list)
        
        # 2. Merge all CDS segments into a single non-overlapping union
        all_cds_s = [s for t in transcript_list for s in t['cds_starts']]
        all_cds_e = [e for t in transcript_list for e in t['cds_ends']]
        merged_cds_s, merged_cds_e = self._merge_intervals(all_cds_s, all_cds_e)
        
        # 3. Define the Global CDS boundaries (The "thick" limits)
        if merged_cds_s:
            global_cds_min = min(merged_cds_s)
            global_cds_max = max(merged_cds_e)
        else:
            # Handling for non-coding genes
            global_cds_min = global_cds_max = None
    
        # 4. Merge all Exon segments (UTR + CDS) to get the total transcript footprint
        all_ex_s = [s for t in transcript_list for s in (t['utr_starts'] + t['cds_starts'])]
        all_ex_e = [e for t in transcript_list for e in (t['utr_ends'] + t['cds_ends'])]
        merged_total_s, merged_total_e = self._merge_intervals(all_ex_s, all_ex_e)
        
        # 5. Clean UTR Logic:
        # We want to keep only the parts of merged_total that are outside [global_cds_min, global_cds_max]
        clean_utr_s, clean_utr_e = [], []
        
        if global_cds_min is not None:
            for s, e in zip(merged_total_s, merged_total_e):
                # Case 1: Exon is entirely before the global CDS start (5' UTR area)
                if e <= global_cds_min:
                    clean_utr_s.append(s)
                    clean_utr_e.append(e)
                # Case 2: Exon spans the global CDS start
                elif s < global_cds_min:
                    clean_utr_s.append(s)
                    clean_utr_e.append(global_cds_min)
                    
                # Case 3: Exon is entirely after the global CDS end (3' UTR area)
                if s >= global_cds_max:
                    clean_utr_s.append(s)
                    clean_utr_e.append(e)
                # Case 4: Exon spans the global CDS end
                elif e > global_cds_max:
                    clean_utr_s.append(global_cds_max)
                    clean_utr_e.append(e)
        else:
            # If the gene is entirely non-coding, all merged exons are "UTRs"
            clean_utr_s, clean_utr_e = merged_total_s, merged_total_e
    
        return {
            'gene_name': f"{gene_name}",
            'chrom': chrom,
            'gene_start': g_start,
            'gene_end': g_end,
            'utr_starts': clean_utr_s, 
            'utr_ends': clean_utr_e,
            'cds_starts': merged_cds_s,
            'cds_ends': merged_cds_e,
            'strand': strand
        }    

    def _flatten_transcripts(self, gene_transcript_dict):
        all_transcripts = [t for t_list in gene_transcript_dict.values() for t in t_list]
        return all_transcripts