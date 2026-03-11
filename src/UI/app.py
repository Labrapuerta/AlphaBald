import ipywidgets as widgets
from IPython.display import display
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import PDBList
import os
import subprocess # Don't forget this for ClustalW!

class CoverageVisualizer:
    def __init__(self, templates, number_of_templates=2):
        self.data = templates
        self.templates_dir = Path("Templates")
        self.target_path = Path("target") / os.listdir("target")[0]
        self.top_templates = templates.head(number_of_templates)

        self.pdb_files = self.download_templates()
        
        # FIX 1: Use the full PDB_ID (e.g. 9H0M_B) as the key so chains don't overwrite each other!
        self.pdb_sequences = {}
        for _, template in self.top_templates.iterrows():
            full_id = template['PDB_ID']
            pdb_4 = full_id[:4].lower()
            chain = full_id[5]
            # Find the downloaded file
            pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
            if pdb_file.exists():
                seq = self.parse_template_sequence(pdb_file, chain)
                if seq:
                    self.pdb_sequences[full_id] = seq # Key is now "9H0M_B"

        self.fasta_for_alignment = Path("temp") / "templates.fasta"
        self.aligned_sequences = Path("alignments") / "aligned_templates.fasta"
        self.target_fasta = SeqIO.read(self.target_path, "fasta")
        self.target_sequence = self.target_fasta.seq
        self.target_length = len(self.target_fasta.seq)

        self.create_fasta_file()
        self.align_sequences()
        
        self.template_data = {}
        self.parse_msa_for_ui()

        self.sliders = {}
        self.ui_render_data = {} 
        self.bg_colors = ['#ADD8E6', '#90EE90', '#FFB6C1', '#FFE4B5', '#D3D3D3'] 
        self.html_output = widgets.HTML() 

        for i, template in self.top_templates.iterrows():
            letter = chr(65 + i)
            color = self.bg_colors[i % len(self.bg_colors)]
            self.create_sliders(template, letter, color)
        
        # Removed the rogue `display(slider)` loop from here so they don't print twice!

    def create_sliders(self, template, letter, color):
        full_id = template['PDB_ID'] # Grab the full ID
        pdb = full_id[0:4]
        chain = full_id[5]
        
        # FIX 2: Look up the MSA info using the full_id
        msa_info = self.template_data.get(full_id, {})

        start = msa_info.get('qstart', template.get('Query_Start', 1))
        end = msa_info.get('qend', template.get('Query_End', self.target_length))
        aligned_seq = msa_info.get('aligned_sequence', "X" * self.target_length)
        
        slider = widgets.IntRangeSlider(
            value=(start, end), 
            min=1, 
            max=self.target_length,
            description=f'Template: {pdb}_{chain}', 
            style={'description_width': 'initial'}, 
            continuous_update=True,
            layout=widgets.Layout(width='80%')
        )
        self.sliders[letter] = slider

        self.ui_render_data[letter] = {
            'sequence': aligned_seq,
            'color': color,
            'pdb': f"{pdb}_{chain}"
        }

    def _update_map(self, **kwargs):
            """Generates the colored HTML sequence alignment map."""
            html = "<div style='font-family: \"Courier New\", monospace; font-size: 14px; letter-spacing: 2px;'>"
            html += "<h4 style='margin-bottom: 2px; margin-top: 10px;'>Target Sequence</h4>"
            html += f"<div style='padding-bottom: 10px; font-weight: bold;'>{self.target_sequence}</div>"
            html += "<h4 style='margin-bottom: 2px;'>Template Alignments</h4>"
            
            for letter, (start, end) in kwargs.items():
                t_info = self.ui_render_data[letter]
                t_seq = t_info['sequence']
                bg_color = t_info['color']
                
                html += f"<div style='margin-bottom: 5px;'><strong>{letter} ({t_info['pdb']}): </strong>"
                for i in range(self.target_length):
                    # Grab the character from our ClustalW alignment
                    t_char = t_seq[i] if i < len(t_seq) else '_' 
                    
                    # Force ClustalW dashes to be underscores for better visibility
                    if t_char == '-':
                        t_char = '_'

                    q_char = self.target_sequence[i]
                    
                    # Check if we are inside the active slider range
                    if start - 1 <= i < end:
                        if t_char == q_char:
                            # Perfect Match
                            html += f"<span style='background-color: {bg_color};'>{t_char}</span>"
                        else:
                            # Mismatch or Structural Gap (now an underscore!)
                            html += f"<span style='background-color: {bg_color}; color: red; font-weight: bold;'>{t_char}</span>"
                    else:
                        # Outside slider range: Display an underscore instead of the old dot
                        html += f"<span style='color: #a0a0a0;'>_</span>"
                        
                html += "</div>"
            html += "</div>"
            self.html_output.value = html
    
    def download_templates(self):
        if not self.templates_dir.exists():
            self.templates_dir.mkdir()
        top_templates = self.top_templates["PDB_ID"].str[:4].tolist()  
        pdbl = PDBList(verbose=False)
        for pdb_id in top_templates:
            pdbl.retrieve_pdb_file(pdb_id, pdir=self.templates_dir, file_format="pdb", overwrite=True)  
        return {pdb_id: self.templates_dir / f"pdb{pdb_id.lower()}.ent" for pdb_id in top_templates}

    def parse_template_sequence(self, pdb_file, chain_id):
        for record in SeqIO.parse(pdb_file, "pdb-atom"):
            chain_letter = record.id.split(':')[-1]
            if chain_letter == chain_id:
                return str(record.seq)
        return None
    
    def create_fasta_file(self):
        with open(self.fasta_for_alignment, 'w') as fasta_file:
            fasta_file.write(f">{self.target_fasta.id}\n{self.target_sequence}\n")
            # FIX 3: Write the full ID to the fasta file
            for full_id, sequence in self.pdb_sequences.items():
                fasta_file.write(f">{full_id}\n{sequence}\n")

    def align_sequences(self):
        cmd = ["clustalw", "-INFILE=" + str(self.fasta_for_alignment), "-OUTFILE=" + str(self.aligned_sequences), "-OUTORDER=INPUT", "-OUTPUT=FASTA"]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def parse_msa_for_ui(self):
        alignment = list(SeqIO.parse(self.aligned_sequences, "fasta"))
        aligned_target = None
        for record in alignment:
            if record.id == self.target_fasta.id:
                aligned_target = str(record.seq)
                self.aligned_target = aligned_target
                break
                
        if not aligned_target:
            raise ValueError(f"Could not find target ID '{self.target_fasta.id}' in the alignment!")

        templates_data = {}
        for record in alignment:
            if record.id == self.target_fasta.id: 
                continue 
                
            aligned_template = str(record.seq)
            target_pos = 0
            start = None
            end = None
            ui_sequence = "" 
            
            for t_char, p_char in zip(aligned_target, aligned_template):
                if t_char != '-':
                    target_pos += 1
                    ui_sequence += p_char 
                if p_char != '-':
                    if start is None:
                        start = target_pos if target_pos > 0 else 1 
                    end = target_pos 
                    
            templates_data[record.id] = {
                'PDB_ID': record.id,
                'qstart': start,
                'qend': end,
                'aligned_sequence': ui_sequence 
            }
        self.template_data = templates_data

    # --- NEW FEATURE: The Submit Button Backend Logic ---
    def _on_submit(self, b):
        """This function runs when the user clicks the Create PDBs button."""
        print("Backend triggered! Preparing to slice PDBs based on sliders...")
        
        for letter, slider in self.sliders.items():
            start_coord, end_coord = slider.value
            pdb_chain = self.ui_render_data[letter]['pdb']
            
            print(f"Creating model for {pdb_chain} -> Cropping to Target Residues: {start_coord} to {end_coord}")
            
            # TODO: Your Biopython structure cropping and MODELLER logic goes here!

    def show(self):
        # 1. Create the Submission Button
        self.submit_btn = widgets.Button(
            description='Create PDBs',
            button_style='success', # Makes it a nice green color
            icon='cogs'             # Adds a little gear icon
        )
        
        # 2. Link the button to our backend function
        self.submit_btn.on_click(self._on_submit)
        
        # 3. Wrap it in an HBox to force it to the right side
        btn_box = widgets.HBox([self.submit_btn], layout=widgets.Layout(justify_content='flex-end'))

        # Link the sliders to the update function
        out = widgets.interactive_output(self._update_map, self.sliders)
        
        # Display sliders on top, HTML visualization in middle, Button on bottom right!
        ui = widgets.VBox([self.html_output, btn_box] + list(self.sliders.values()))
        display(ui)


import IPython.display as display

class DomainVisualizer:
    def __init__(self, target_data, template_data_list):
        """
        target_data: [target_id, target_sequence, target_domains_list]
        template_data_list: list of [pdb_id, sequence, domains_list]
        """
        self.target = target_data
        self.templates = template_data_list
        
        self.colors = [
            "rgba(255, 99, 132, 0.7)", "rgba(54, 162, 235, 0.7)", 
            "rgba(255, 206, 86, 0.7)", "rgba(75, 192, 192, 0.7)", 
            "rgba(153, 102, 255, 0.7)", "rgba(255, 159, 64, 0.7)",
            "rgba(46, 204, 113, 0.7)", "rgba(231, 76, 60, 0.7)"
        ]

    def _generate_sequence_html(self, seq_id, sequence, domains, is_target=False):
        seq_len = len(sequence)
        
        # Styling adjustments
        bg_color = "#f0f8ff" if is_target else "#ffffff"
        border = "2px solid #0056b3" if is_target else "1px solid #ccc"
        label_color = "#0056b3" if is_target else "#333"
        
        # Vertical space for tracks
        bottom_padding = len(domains) * 26 + 10 
        positions = [str(i) for i in range(1, seq_len + 1)]
        pos_html = "<div style='font-family: \"Courier New\", monospace; font-size: 10px; color: #888; padding: 0 4px;'>"
        for i in range(0, seq_len, max(1, seq_len // 10)):
            pos_html += f"<span style='position: absolute; left: {(i / seq_len) * 100}%; transform: translateX(-50%);'>{i+1}</span>"
        pos_html += "</div>"
        html = f"""
        <div style="margin-bottom: 25px; padding: 10px; background-color: {bg_color}; border: {border}; border-radius: 8px; font-family: sans-serif; overflow-x: auto;">
            
            <div style="font-weight: bold; color: {label_color}; margin-bottom: 10px; font-size: 14px; position: sticky; left: 0;">
                {seq_id} <span style="font-weight: normal; font-size: 12px; color: #666;">(Length: {seq_len} aa)</span>
            </div>
            
            <div style="position: relative; display: inline-block; padding-bottom: {bottom_padding}px;">
                
                <div style="font-family: 'Courier New', Courier, monospace; font-size: 14px; letter-spacing: 0; color: #222; background-color: #e9ecef; border-radius: 4px; padding: 0; margin: 0; white-space: nowrap;">
                    {sequence}
                </div>
        """
        html += pos_html
        
        for i, domain in enumerate(domains):
            start = domain.get('start', 1)
            end = domain.get('end', seq_len)
            name = domain.get('name', 'Unknown Domain')
            
            # Safety check: prevent domain boxes from drawing outside the sequence length
            start = max(1, min(start, seq_len))
            end = max(1, min(end, seq_len))
            
            # MATH: Calculate exact percentages based on sequence length
            left_percent = ((start - 1) / seq_len) * 100
            width_percent = ((end - start + 1) / seq_len) * 100
            
            color_idx = hash(domain.get('id', name)) % len(self.colors)
            color = self.colors[color_idx]
            top_pos = 26 + (i * 24) 
            
            html += f"""
                <div title="{name} ({start}-{end})" 
                     style="position: absolute; left: {left_percent}%; top: {top_pos}px; width: {width_percent}%; height: 20px; 
                            background-color: {color}; border: 1px solid rgba(0,0,0,0.4); border-radius: 4px; 
                            cursor: pointer; box-sizing: border-box; font-size: 11px; font-family: sans-serif; 
                            color: #000; line-height: 18px; padding-left: 4px; padding-right: 4px;
                            white-space: nowrap; overflow: hidden; text-overflow: ellipsis; font-weight: bold;">
                    {name}
                </div>
            """
            
        html += """
            </div>
        </div>
        """
        return html

    def show(self):
        final_html = "<div style='width: 100%;'>"
        
        # Render Target First
        final_html += self._generate_sequence_html(self.target[0], self.target[1], self.target[2], is_target=True)
        final_html += "<hr style='border-top: 2px dashed #ccc; margin: 20px 0;'>"
        
        # Render Templates
        for template in self.templates:
            final_html += self._generate_sequence_html(template[0], template[1], template[2])
            
        final_html += "</div>"
        display.display(display.HTML(final_html))