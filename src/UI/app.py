"""
UI visualization components for the homology modeling pipeline.

Provides interactive widgets for template selection and domain visualization.
"""

import ipywidgets as widgets
from IPython.display import display as ipython_display, HTML as IPythonHTML
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select, DSSP
from Bio.PDB.Polypeptide import protein_letters_3to1
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import PDBList
import subprocess
from typing import List, Dict, Any, Optional, Tuple


class ResidueRangeSelect(Select):
    """Select residues within a specific range for PDB cropping."""

    def __init__(self, chain_id: str, start: int, end: int):
        self.chain_id = chain_id
        self.start = start
        self.end = end

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        res_id = residue.get_id()[1]
        return self.start <= res_id <= self.end


class CoverageVisualizer:
    """
    Interactive visualization for template coverage selection.

    Allows users to select coverage regions using sliders and generates
    cropped PDB files and PIR alignments for MODELLER.

    Parameters
    ----------
    templates : pd.DataFrame
        DataFrame with template hits (must have 'PDB_ID' column)
    number_of_templates : int, optional
        Number of templates to visualize. Default 2
    target_domains : list, optional
        List of domain dictionaries for target sequence
    template_domains : dict, optional
        Dict mapping PDB_ID to list of domain dictionaries
    """

    def __init__(
        self,
        templates: pd.DataFrame,
        number_of_templates: int = 2,
        target_domains: Optional[List[Dict]] = None,
        template_domains: Optional[Dict[str, List[Dict]]] = None
    ):
        self.data = templates
        self.templates_dir = Path("Templates")
        self.modeller_dir = Path("Modeller_Templates")
        self.modeller_dir.mkdir(exist_ok=True)
        self.templates_dir.mkdir(exist_ok=True)

        # Store domain data
        self.target_domains = target_domains or []
        self.template_domains = template_domains or {}

        # Find target file
        target_dir = Path("target")
        if target_dir.exists():
            target_files = [f for f in target_dir.iterdir()
                           if f.suffix.lower() in ['.fa', '.fasta', '.faa']]
            self.target_path = target_files[0] if target_files else target_dir / "target.fa"
        else:
            self.target_path = Path("target/target.fa")

        self.top_templates = templates.head(number_of_templates).copy()

        # Download templates
        print("Downloading template PDB files...")
        self.pdb_files = self._download_templates()

        # Parse sequences
        self.pdb_sequences = {}
        self.chain_residue_map = {}

        for _, template in self.top_templates.iterrows():
            full_id = template['PDB_ID']
            pdb_4 = full_id[:4].lower()
            chain = full_id[5] if len(full_id) > 5 else 'A'
            pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
            if pdb_file.exists():
                seq, res_range = self._parse_sequence_with_range(pdb_file, chain)
                if seq:
                    self.pdb_sequences[full_id] = seq
                    self.chain_residue_map[full_id] = res_range

        # Setup paths
        Path("temp").mkdir(exist_ok=True)
        Path("Alignments").mkdir(exist_ok=True)
        self.fasta_for_alignment = Path("temp") / "coverage_templates.fasta"
        self.aligned_sequences = Path("Alignments") / "coverage_aligned.fasta"

        # Load target
        self.target_fasta = SeqIO.read(str(self.target_path), "fasta")
        self.target_sequence = str(self.target_fasta.seq)
        self.target_length = len(self.target_sequence)

        # Create alignment
        print("Creating sequence alignment...")
        self._create_fasta_file()
        self._align_sequences()

        # Parse MSA for UI
        self.template_data = {}
        self._parse_msa_for_ui()

        # Setup widgets
        self.sliders = {}
        self.checkboxes = {}  # NEW: checkboxes for template selection
        self.ui_render_data = {}
        self.bg_colors = ['#ADD8E6', '#90EE90', '#FFB6C1', '#FFE4B5', '#D3D3D3',
                          '#FFDAB9', '#E6E6FA', '#98FB98']
        self.domain_colors = [
            '#FF6B6B', '#4ECDC4', '#FFE66D', '#95E1D3',
            '#DDA0DD', '#F0E68C', '#87CEEB', '#FFA07A'
        ]
        self.html_output = widgets.HTML()
        self.status_output = widgets.HTML()

        # Create sliders and checkboxes for each template
        for i, (_, template) in enumerate(self.top_templates.iterrows()):
            letter = chr(65 + i)
            color = self.bg_colors[i % len(self.bg_colors)]
            self._create_slider(template, letter, color)
            self._create_checkbox(template, letter)

        print(f"CoverageVisualizer ready with {len(self.pdb_sequences)} templates")

    def _download_templates(self) -> Dict[str, Path]:
        """Download PDB files for templates."""
        pdb_files = {}
        pdbl = PDBList(verbose=False)
        unique_pdbs = self.top_templates["PDB_ID"].str[:4].unique().tolist()

        for pdb_id in unique_pdbs:
            pdb_file = self.templates_dir / f"pdb{pdb_id.lower()}.ent"
            if not pdb_file.exists():
                try:
                    pdbl.retrieve_pdb_file(
                        pdb_id,
                        pdir=str(self.templates_dir),
                        file_format="pdb",
                        overwrite=True
                    )
                except Exception as e:
                    print(f"  Warning: Could not download {pdb_id}: {e}")
            if pdb_file.exists():
                pdb_files[pdb_id.upper()] = pdb_file
        return pdb_files

    def _parse_sequence_with_range(self, pdb_file: Path, chain_id: str) -> Tuple[Optional[str], Tuple[int, int]]:
        """Parse sequence and residue number range from PDB."""
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("temp", str(pdb_file))
            if chain_id not in [c.id for c in structure[0].get_chains()]:
                return None, (0, 0)
            chain = structure[0][chain_id]
            residues = [r for r in chain.get_residues() if r.id[0] == ' ']

            if not residues:
                return None, (0, 0)

            first_res = residues[0].id[1]
            last_res = residues[-1].id[1]
            sequence = "".join([protein_letters_3to1.get(r.get_resname(), 'X') for r in residues])

            return sequence, (first_res, last_res)
        except Exception as e:
            print(f"  Warning parsing {pdb_file}: {e}")
            return None, (0, 0)

    def _create_fasta_file(self):
        """Create combined FASTA for alignment."""
        with open(self.fasta_for_alignment, 'w') as f:
            f.write(f">{self.target_fasta.id}\n{self.target_sequence}\n")
            for full_id, sequence in self.pdb_sequences.items():
                f.write(f">{full_id}\n{sequence}\n")

    def _align_sequences(self):
        """Run ClustalW alignment."""
        cmd = [
            "clustalw",
            f"-INFILE={self.fasta_for_alignment}",
            f"-OUTFILE={self.aligned_sequences}",
            "-OUTORDER=INPUT",
            "-OUTPUT=FASTA"
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"Warning: ClustalW alignment failed: {e}")

    def _parse_msa_for_ui(self):
        """Parse MSA to get coverage info for sliders."""
        if not self.aligned_sequences.exists():
            return

        alignment = list(SeqIO.parse(str(self.aligned_sequences), "fasta"))
        aligned_target = None

        for record in alignment:
            if record.id == self.target_fasta.id:
                aligned_target = str(record.seq)
                self.aligned_target = aligned_target
                break

        if not aligned_target:
            print(f"Warning: Target '{self.target_fasta.id}' not found in alignment")
            return

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

            self.template_data[record.id] = {
                'PDB_ID': record.id,
                'qstart': start or 1,
                'qend': end or self.target_length,
                'aligned_sequence': ui_sequence
            }

    def _create_slider(self, template, letter: str, color: str):
        """Create a slider for a template."""
        full_id = template['PDB_ID']
        pdb = full_id[:4]
        chain = full_id[5] if len(full_id) > 5 else 'A'

        msa_info = self.template_data.get(full_id, {})
        start = msa_info.get('qstart', template.get('Query_Start', 1)) or 1
        end = msa_info.get('qend', template.get('Query_End', self.target_length)) or self.target_length
        aligned_seq = msa_info.get('aligned_sequence', '-' * self.target_length)

        start = max(1, min(start, self.target_length))
        end = max(start, min(end, self.target_length))

        slider = widgets.IntRangeSlider(
            value=(start, end),
            min=1,
            max=self.target_length,
            description=f'{pdb}_{chain}:',
            style={'description_width': '100px'},
            continuous_update=False,
            layout=widgets.Layout(width='90%')
        )

        # Observe changes
        slider.observe(self._on_slider_change, names='value')

        self.sliders[letter] = slider
        self.ui_render_data[letter] = {
            'sequence': aligned_seq,
            'color': color,
            'pdb': f"{pdb}_{chain}",
            'full_id': full_id
        }

    def _create_checkbox(self, template, letter: str):
        """Create a checkbox for template selection in PIR file."""
        full_id = template['PDB_ID']
        pdb = full_id[:4]
        chain = full_id[5] if len(full_id) > 5 else 'A'

        checkbox = widgets.Checkbox(
            value=True,  # Selected by default
            description=f'Include {pdb}_{chain}',
            style={'description_width': 'auto'},
            layout=widgets.Layout(width='200px')
        )
        self.checkboxes[letter] = checkbox

    def _on_slider_change(self, change):
        """Handle slider value changes."""
        self._render_alignment()

    def _render_alignment(self):
        """Render the HTML alignment visualization with domain tracks."""
        html = """
        <div style='font-family: "Courier New", monospace; font-size: 12px;
                    background: #f8f9fa; padding: 15px; border-radius: 8px;
                    overflow-x: auto; max-width: 100%;'>
        <h4 style='margin: 0 0 10px 0; color: #333;'>Target Sequence</h4>
        <div style='background: #fff; padding: 5px; border-radius: 4px;
                    margin-bottom: 5px; word-wrap: break-word;'>
        """

        # Render target with position markers
        for i, char in enumerate(self.target_sequence):
            if i % 10 == 0:
                html += f"<span style='color: #999; font-size: 10px;'>{i+1}</span>"
            html += f"<span style='color: #333;'>{char}</span>"
        html += "</div>"

        # Render target domains track
        if self.target_domains:
            html += self._render_domain_track(self.target_domains, "Target Domains")

        html += "<h4 style='margin: 15px 0 10px 0; color: #333;'>Template Coverage</h4>"

        # Render each template with its domain track
        for letter, slider in self.sliders.items():
            start, end = slider.value
            t_info = self.ui_render_data[letter]
            t_seq = t_info['sequence']
            bg_color = t_info['color']
            full_id = t_info.get('full_id', t_info['pdb'])

            html += f"<div style='margin: 10px 0;'>"
            html += f"<strong style='color: #555;'>{letter} ({t_info['pdb']}):</strong> "

            for i in range(self.target_length):
                t_char = t_seq[i] if i < len(t_seq) else '-'
                if t_char == '-':
                    t_char = '_'

                q_char = self.target_sequence[i]

                if start - 1 <= i < end:
                    if t_char == q_char and t_char != '_':
                        html += f"<span style='background: {bg_color};'>{t_char}</span>"
                    else:
                        html += f"<span style='background: {bg_color}; color: red; font-weight: bold;'>{t_char}</span>"
                else:
                    html += f"<span style='color: #ccc;'>·</span>"

            html += "</div>"

            # Render template domains track if available
            template_doms = self.template_domains.get(full_id, [])
            if template_doms:
                html += self._render_domain_track(template_doms, f"{t_info['pdb']} Domains", indent=True)

        html += "</div>"
        self.html_output.value = html

    def _render_domain_track(self, domains: List[Dict], label: str, indent: bool = False) -> str:
        """Render a domain track as colored segments."""
        if not domains:
            return ""

        indent_style = "margin-left: 20px;" if indent else ""
        html = f"<div style='margin: 3px 0; {indent_style}'>"
        html += f"<span style='color: #666; font-size: 10px;'>{label}: </span>"

        # Create position-to-domain mapping
        domain_at_pos = {}
        for i, dom in enumerate(domains):
            start = dom.get('start', 1)
            end = dom.get('end', self.target_length)
            for pos in range(max(0, start - 1), min(self.target_length, end)):
                if pos not in domain_at_pos:
                    domain_at_pos[pos] = (i, dom)

        # Render the track
        for pos in range(self.target_length):
            if pos in domain_at_pos:
                idx, dom = domain_at_pos[pos]
                color = self.domain_colors[idx % len(self.domain_colors)]
                title = f"{dom.get('name', 'Unknown')} ({dom.get('id', '')})"
                html += f"<span style='background: {color};' title='{title}'>\u2588</span>"
            else:
                html += "<span style='color: #ddd;'>\u2591</span>"

        # Add legend for domains
        html += "<br><span style='font-size: 9px; color: #888;'>"
        seen_domains = {}
        for dom in domains:
            dom_id = dom.get('id', dom.get('name', ''))
            if dom_id not in seen_domains:
                idx = domains.index(dom)
                color = self.domain_colors[idx % len(self.domain_colors)]
                name = dom.get('name', 'Unknown')[:20]
                html += f"<span style='background: {color}; padding: 0 3px; margin: 0 2px;'>{name}</span>"
                seen_domains[dom_id] = True
        html += "</span></div>"
        return html

    def _on_submit(self, b):
        """Create cropped PDBs and combined PIR file for MODELLER."""
        self.status_output.value = "<p style='color: blue;'>Processing templates...</p>"
        created_files = []
        selected_templates = []

        # First, delete previous PIR files in Modeller_Templates
        for old_pir in self.modeller_dir.glob("*.pir"):
            old_pir.unlink()

        for letter, slider in self.sliders.items():
            # Check if this template is selected
            if letter in self.checkboxes and not self.checkboxes[letter].value:
                continue  # Skip unselected templates

            target_start, target_end = slider.value
            t_info = self.ui_render_data[letter]
            pdb_chain = t_info['pdb']
            full_id = t_info.get('full_id', pdb_chain)

            pdb_4 = pdb_chain[:4].lower()
            chain_id = pdb_chain[5] if len(pdb_chain) > 5 else 'A'

            pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
            if not pdb_file.exists():
                continue

            # Crop PDB
            cropped_path = self._crop_pdb(pdb_file, chain_id, target_start, target_end, full_id)
            if cropped_path:
                # Store template info for combined PIR
                template_seq = self._get_pdb_sequence(cropped_path)
                target_subseq = self.target_sequence[target_start - 1:target_end]
                selected_templates.append({
                    'pdb_id': full_id,
                    'pdb_4': pdb_4,
                    'chain_id': chain_id,
                    'template_seq': template_seq,
                    'target_subseq': target_subseq,
                    'cropped_path': cropped_path,
                    'target_start': target_start,
                    'target_end': target_end
                })
                created_files.append((pdb_chain, cropped_path))

        if selected_templates:
            # Create combined PIR file with all selected templates
            combined_pir_path = self._create_combined_pir_file(selected_templates)

            html = "<div style='color: green; padding: 10px; background: #d4edda; border-radius: 5px;'>"
            html += "<strong>Templates created for MODELLER:</strong><ul style='margin: 5px 0;'>"
            for pdb_id, pdb_path in created_files:
                html += f"<li><strong>{pdb_id}</strong>: {pdb_path.name}</li>"
            html += f"</ul><p style='margin: 5px 0;'><strong>Combined PIR file:</strong> <code>{combined_pir_path.name}</code></p>"
            html += f"<p style='margin: 5px 0;'>Output: <code>{self.modeller_dir}</code></p></div>"
            self.status_output.value = html
        else:
            self.status_output.value = "<p style='color: red;'>No templates selected.</p>"

    def _create_combined_pir_file(self, selected_templates: List[Dict]) -> Path:
        """Create a combined PIR alignment file with all selected templates for MODELLER."""
        pir_content = ""

        # Add all template sequences first
        for template in selected_templates:
            pdb_4 = template['pdb_4']
            chain_id = template['chain_id']
            template_seq = template['template_seq']

            pir_content += f">P1;{pdb_4}{chain_id}\n"
            pir_content += f"structureX:{pdb_4}{chain_id}::{chain_id}::{chain_id}::::\n"
            pir_content += f"{template_seq}*\n\n"

        # Add target sequence (use the longest target subsequence)
        if selected_templates:
            # Find the range that covers all selected regions
            min_start = min(t['target_start'] for t in selected_templates)
            max_end = max(t['target_end'] for t in selected_templates)
            target_subseq = self.target_sequence[min_start - 1:max_end]

            pir_content += f">P1;{self.target_fasta.id}\n"
            pir_content += f"sequence:{self.target_fasta.id}:::::::0.00:0.00\n"
            pir_content += f"{target_subseq}*\n"

        # Save combined PIR file
        pir_path = self.modeller_dir / "combined_alignment.pir"
        with open(pir_path, 'w') as f:
            f.write(pir_content)

        return pir_path

    def _crop_pdb(self, pdb_file: Path, chain_id: str, target_start: int, target_end: int, pdb_id: str) -> Optional[Path]:
        """Crop PDB to alignment range."""
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("temp", str(pdb_file))
            chain = structure[0][chain_id]
            residues = [r for r in chain.get_residues() if r.id[0] == ' ']

            if not residues:
                return None

            msa_info = self.template_data.get(pdb_id, {})
            qstart = msa_info.get('qstart', 1) or 1

            first_pdb_res = residues[0].id[1]
            offset = target_start - qstart
            template_start = first_pdb_res + max(0, offset)
            template_end = template_start + (target_end - target_start)

            last_pdb_res = residues[-1].id[1]
            template_start = max(first_pdb_res, min(template_start, last_pdb_res))
            template_end = max(template_start, min(template_end, last_pdb_res))

            output_path = self.modeller_dir / f"{pdb_id.replace('_', '')}_cropped.pdb"
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_path), ResidueRangeSelect(chain_id, template_start, template_end))

            return output_path
        except Exception as e:
            print(f"Error cropping {pdb_id}: {e}")
            return None

    def _create_pir_file(self, pdb_id: str, cropped_pdb: Path, target_start: int, target_end: int) -> Path:
        """Create PIR alignment file for MODELLER."""
        pdb_4 = pdb_id[:4].lower()
        chain_id = pdb_id[5] if len(pdb_id) > 5 else 'A'

        template_seq = self._get_pdb_sequence(cropped_pdb)
        target_subseq = self.target_sequence[target_start - 1:target_end]

        pir_content = f""">P1;{pdb_4}{chain_id}
structureX:{pdb_4}{chain_id}::{chain_id}::{chain_id}::::
{template_seq}*

>P1;{self.target_fasta.id}
sequence:{self.target_fasta.id}:::::::0.00:0.00
{target_subseq}*
"""
        pir_path = self.modeller_dir / f"{pdb_id.replace('_', '')}_alignment.pir"
        with open(pir_path, 'w') as f:
            f.write(pir_content)
        return pir_path

    def _get_pdb_sequence(self, pdb_file: Path) -> str:
        """Extract sequence from PDB file."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", str(pdb_file))
        sequence = ""
        for model in structure:
            for chain in model:
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        sequence += protein_letters_3to1.get(residue.get_resname(), 'X')
                break
            break
        return sequence

    def show(self):
        """Display the interactive visualization."""
        # Create submit button
        submit_btn = widgets.Button(
            description='Create MODELLER Templates',
            button_style='success',
            icon='check',
            layout=widgets.Layout(width='250px')
        )
        submit_btn.on_click(self._on_submit)

        # Initial render
        self._render_alignment()

        # Create checkbox row for template selection
        checkbox_widgets = list(self.checkboxes.values())
        checkbox_row = widgets.HBox(
            checkbox_widgets,
            layout=widgets.Layout(justify_content='flex-start', margin='10px 0')
        )

        # Create slider-checkbox pairs
        slider_rows = []
        for letter in self.sliders:
            slider = self.sliders[letter]
            slider_rows.append(slider)

        # Layout
        btn_box = widgets.HBox([submit_btn], layout=widgets.Layout(justify_content='flex-end'))

        ui = widgets.VBox([
            self.html_output,
            widgets.HTML("<hr style='margin: 10px 0;'>"),
            widgets.HTML("<strong>Select templates for PIR file:</strong>"),
            checkbox_row,
            widgets.HTML("<strong>Adjust coverage regions:</strong>"),
            *slider_rows,
            btn_box,
            self.status_output
        ], layout=widgets.Layout(padding='10px'))

        ipython_display(ui)


def run_dssp(pdb_file: Path) -> Dict[int, Dict[str, Any]]:
    """
    Run DSSP to get secondary structure assignments.

    Parameters
    ----------
    pdb_file : Path
        Path to PDB file

    Returns
    -------
    dict
        Mapping of residue number to SS info (ss, phi, psi, rsa)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("temp", str(pdb_file))
    model = structure[0]

    try:
        dssp = DSSP(model, str(pdb_file), dssp='dssp')
        ss_data = {}
        for key in dssp:
            res_id = key[1][1]
            ss_data[res_id] = {
                'ss': dssp[key][2],
                'rsa': dssp[key][3],
                'phi': dssp[key][4],
                'psi': dssp[key][5]
            }
        return ss_data
    except Exception as e:
        print(f"DSSP failed: {e}")
        return {}


# Secondary structure color mapping
SS_COLORS = {
    'H': '#FF6B6B',  # Alpha helix - red
    'G': '#FF8E8E',  # 3-10 helix - light red
    'I': '#FFB0B0',  # Pi helix - lighter red
    'E': '#4ECDC4',  # Beta sheet - cyan
    'B': '#7FDBDB',  # Beta bridge - light cyan
    'T': '#FFE66D',  # Turn - yellow
    'S': '#95E1D3',  # Bend - green
    '-': '#E8E8E8',  # Coil - gray
    ' ': '#E8E8E8',  # Coil - gray
}


class DomainVisualizer:
    """
    Visualize protein domains and secondary structure.

    Parameters
    ----------
    target_data : list
        [target_id, target_sequence, target_domains_list]
    template_data_list : list
        List of [pdb_id, sequence, domains_list]
    show_secondary_structure : bool
        Whether to show DSSP secondary structure
    templates_dir : str
        Directory containing PDB files
    """

    def __init__(
        self,
        target_data: List,
        template_data_list: List[List],
        show_secondary_structure: bool = True,
        templates_dir: str = "Templates"
    ):
        self.target = target_data
        self.templates = template_data_list
        self.show_ss = show_secondary_structure
        self.templates_dir = Path(templates_dir)

        self.domain_colors = [
            "rgba(255, 99, 132, 0.8)", "rgba(54, 162, 235, 0.8)",
            "rgba(255, 206, 86, 0.8)", "rgba(75, 192, 192, 0.8)",
            "rgba(153, 102, 255, 0.8)", "rgba(255, 159, 64, 0.8)",
            "rgba(46, 204, 113, 0.8)", "rgba(231, 76, 60, 0.8)"
        ]

        # Get secondary structure for templates
        self.ss_data = {}
        if self.show_ss:
            for template in self.templates:
                pdb_id = template[0]
                pdb_4 = pdb_id[:4].lower()
                pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
                if pdb_file.exists():
                    self.ss_data[pdb_id] = run_dssp(pdb_file)

    def _generate_ss_track(self, seq_len: int, ss_data: Dict[int, Dict]) -> str:
        """Generate secondary structure track HTML."""
        if not ss_data:
            return ""

        html = "<div style='height: 12px; display: flex; margin-top: 5px; border-radius: 3px; overflow: hidden;'>"

        for i in range(1, seq_len + 1):
            ss = ss_data.get(i, {}).get('ss', '-')
            color = SS_COLORS.get(ss, SS_COLORS['-'])
            title = f"Residue {i}: {ss}"
            html += f"<div title='{title}' style='flex: 1; background: {color};'></div>"

        html += "</div>"
        return html

    def _generate_sequence_html(
        self,
        seq_id: str,
        sequence: str,
        domains: List[Dict],
        is_target: bool = False,
        ss_data: Optional[Dict] = None
    ) -> str:
        """Generate HTML for a single sequence with domains and SS."""
        seq_len = len(sequence)

        bg_color = "#e3f2fd" if is_target else "#ffffff"
        border = "2px solid #1976d2" if is_target else "1px solid #ddd"
        label_color = "#1976d2" if is_target else "#333"
        tag = "TARGET" if is_target else "TEMPLATE"

        # Calculate height needed for domain tracks
        num_tracks = len(domains) + (1 if ss_data else 0)
        bottom_padding = max(40, num_tracks * 28 + 20)

        html = f"""
        <div style="margin-bottom: 20px; padding: 15px; background: {bg_color};
                    border: {border}; border-radius: 10px; font-family: sans-serif;">
            <div style="display: flex; align-items: center; margin-bottom: 10px;">
                <span style="background: {'#1976d2' if is_target else '#666'}; color: white;
                             padding: 2px 8px; border-radius: 4px; font-size: 11px; margin-right: 10px;">
                    {tag}
                </span>
                <span style="font-weight: bold; color: {label_color}; font-size: 14px;">
                    {seq_id}
                </span>
                <span style="color: #666; font-size: 12px; margin-left: 10px;">
                    ({seq_len} aa)
                </span>
            </div>

            <div style="position: relative; padding-bottom: {bottom_padding}px;">
                <div style="font-family: 'Courier New', monospace; font-size: 11px;
                            letter-spacing: 0; color: #333; background: #f5f5f5;
                            border-radius: 4px; padding: 5px; white-space: nowrap;
                            overflow-x: auto;">
                    {sequence[:100]}{'...' if seq_len > 100 else ''}
                </div>
        """

        # Add position markers
        html += """
                <div style="position: relative; height: 15px; margin-top: 3px;">
        """
        for i in range(0, seq_len, max(1, seq_len // 8)):
            left_pct = (i / seq_len) * 100
            html += f"""
                    <span style="position: absolute; left: {left_pct}%; font-size: 9px; color: #999;">
                        {i + 1}
                    </span>
            """
        html += "</div>"

        # Add secondary structure track
        if ss_data:
            html += f"""
                <div style="margin-top: 5px;">
                    <span style="font-size: 10px; color: #666;">Secondary Structure:</span>
                    {self._generate_ss_track(seq_len, ss_data)}
                </div>
            """

        # Add domain tracks
        track_offset = 50 if ss_data else 35
        for i, domain in enumerate(domains):
            start = max(1, domain.get('start', 1))
            end = min(seq_len, domain.get('end', seq_len))
            name = domain.get('name', 'Unknown')
            domain_id = domain.get('id', '')

            left_pct = ((start - 1) / seq_len) * 100
            width_pct = ((end - start + 1) / seq_len) * 100

            color_idx = hash(domain_id or name) % len(self.domain_colors)
            color = self.domain_colors[color_idx]

            top_pos = track_offset + (i * 26)

            html += f"""
                <div title="{name} ({domain_id}): {start}-{end}"
                     style="position: absolute; left: {left_pct}%; top: {top_pos}px;
                            width: {width_pct}%; height: 22px; background: {color};
                            border: 1px solid rgba(0,0,0,0.3); border-radius: 4px;
                            font-size: 11px; color: #000; line-height: 20px;
                            padding: 0 5px; box-sizing: border-box; cursor: pointer;
                            white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
                            font-weight: 600;">
                    {name}
                </div>
            """

        html += """
            </div>
        </div>
        """
        return html

    def show(self):
        """Display the domain visualization."""
        html = "<div style='width: 100%; max-width: 1200px;'>"

        # Legend for secondary structure
        if self.show_ss:
            html += """
            <div style="margin-bottom: 15px; padding: 10px; background: #f8f9fa;
                        border-radius: 8px; font-size: 12px;">
                <strong>Secondary Structure Legend:</strong>
                <span style="margin-left: 10px;">
                    <span style="background: #FF6B6B; padding: 2px 8px; border-radius: 3px;">H: α-helix</span>
                    <span style="background: #4ECDC4; padding: 2px 8px; border-radius: 3px; margin-left: 5px;">E: β-sheet</span>
                    <span style="background: #FFE66D; padding: 2px 8px; border-radius: 3px; margin-left: 5px;">T: Turn</span>
                    <span style="background: #E8E8E8; padding: 2px 8px; border-radius: 3px; margin-left: 5px;">-: Coil</span>
                </span>
            </div>
            """

        # Render target
        target_ss = None  # Target usually doesn't have structure yet
        html += self._generate_sequence_html(
            self.target[0], self.target[1], self.target[2],
            is_target=True, ss_data=target_ss
        )

        html += "<hr style='border: none; border-top: 2px dashed #ddd; margin: 25px 0;'>"

        # Render templates
        for template in self.templates:
            pdb_id = template[0]
            ss_data = self.ss_data.get(pdb_id) if self.show_ss else None
            html += self._generate_sequence_html(
                template[0], template[1], template[2],
                is_target=False, ss_data=ss_data
            )

        html += "</div>"
        ipython_display(IPythonHTML(html))
