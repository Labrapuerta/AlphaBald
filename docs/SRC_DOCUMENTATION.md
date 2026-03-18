# Documentacion completa de `src/` (AlphaBald)

Este documento describe todo el contenido funcional del directorio `src/`, su arquitectura, APIs principales y flujo de uso recomendado.

## 1) Vision general

El paquete `src/` implementa un pipeline de modelado por homologia para proteinas con cinco bloques principales:

1. `Setup`: prepara carpetas y bases de datos (BLAST/HMMER).
2. `Homology`: recupera templates y analiza dominios.
3. `pipeline.py`: orquesta el flujo extremo a extremo.
4. `modeller`: genera scripts para construccion/refinamiento/evaluacion con MODELLER.
5. `Analysis` + `UI`: evaluacion estructural y visualizacion interactiva.

## 2) Mapa del paquete

```text
src/
├── __init__.py
├── pipeline.py
├── Setup/
│   ├── __init__.py
│   ├── setup.py
│   └── preprocessing.py
├── Homology/
│   ├── __init__.py
│   ├── retrieve.py
│   ├── domains.py
│   └── superimpose.py
├── UI/
│   └── app.py
├── modeller/
│   ├── __init__.py
│   └── scripts.py
└── Analysis/
    ├── __init__.py
    ├── assessment.py
    └── visualization.py
```

## 3) API raiz del paquete

Archivo: `src/__init__.py`

Expone una interfaz de alto nivel:

- `HomologyPipeline`
- `run_pipeline`
- `Setup`
- `TargetPreprocessor`
- `TemplateRetriever`
- `TemplateProcessor`
- `SuperimpositionVisualizer`

Version declarada:

- `__version__ = "1.0.0"`

## 4) Orquestador principal

### `src/pipeline.py`

### Clase `HomologyPipeline`

Orquesta el flujo completo:

1. Setup de entorno y bases de datos.
2. Carga/preprocesado de secuencia target.
3. Busqueda de templates por PSI-BLAST y opcionalmente HMMER.
4. Analisis de templates/dominios.
5. Exposicion de resultados y utilidades de visualizacion.

Parametros clave:

- `target_path` (default `target/target.fa`)
- `e_value_threshold` (default `1e-5`)
- `num_iterations` (PSI-BLAST, default `5`)
- `num_templates` (default `5`)
- `skip_setup`
- `run_hmmer`
- `verbose`

Atributos clave:

- `target`: `TargetPreprocessor`
- `retriever`: `TemplateRetriever`
- `processor`: `TemplateProcessor` o `None`
- `homologs`: `pd.DataFrame`
- `templates`: top templates

Metodos utiles:

- `plot_homologs(top_n=10)`
- `visualize(num_templates=None)`
- `get_summary()`
- `export_templates(output_file="templates_summary.csv")`

### Funcion `run_pipeline(...)`

Wrapper de conveniencia para construir y ejecutar el pipeline completo con menos codigo.

## 5) Setup y preprocesado

## `src/Setup/__init__.py`

Re-exporta:

- `Setup`
- `TargetPreprocessor`

## `src/Setup/setup.py`

### Clase `Setup`

Responsable de crear directorios y preparar bases de datos:

- SwissProt BLAST
- PDB sequence BLAST (`pdbaa`)
- Pfam HMM (descarga, descompresion e indexado con `hmmpress`)

Metodos principales:

- `create_directories(clear_temp=True)`
- `setup_swissprot_database()`
- `setup_pdb_database()`
- `setup_pfam_database()`
- `create_pdb_fasta()` (convierte DB BLAST a FASTA para HMMER)

Notas:

- Usa utilidades externas (`update_blastdb.pl`, `blastdbcmd`, `hmmpress`).
- Puede limpiar carpetas temporales (`Alignments`, `Templates`, `temp`).

## `src/Setup/preprocessing.py`

### Clase `TargetPreprocessor`

Carga y valida la entrada objetivo desde FASTA o PDB.

Capacidades:

- Auto-deteccion de archivo target (`.fa/.fasta/.faa/.pdb/.ent`).
- Extraccion de secuencia desde FASTA o desde residuos ATOM de PDB.
- Lectura de metadata basica (incluyendo resolucion si aplica).
- Export a FASTA.

Metodos:

- `_auto_detect_and_load()`
- `_load_from_fasta(...)`
- `_load_from_pdb(...)`
- `to_dict()`
- `write_fasta(output_path=None)`

## 6) Busqueda de templates y dominios

## `src/Homology/__init__.py`

Re-exporta:

- `TemplateRetriever`
- `TemplateProcessor`
- `get_ca_atoms`, `get_ca_mapping`, `superimpose_structures`, `write_structure`
- `SuperimpositionVisualizer`

## `src/Homology/retrieve.py`

### Clase `TemplateRetriever`

Implementa la estrategia de recuperacion de templates:

1. PSI-BLAST contra SwissProt para construir PSSM.
2. PSI-BLAST de PSSM contra PDB (`pdbaa`).
3. Construccion de MSA de homologos top (ClustalW).
4. PSI-BLAST usando MSA contra PDB.
5. Opcional: `jackhmmer` y `hmmbuild/hmmsearch`.
6. Ranking y deduplicacion de hits.

Salidas importantes:

- `top_swissprot_hits`
- `top_psi_hits`
- `top_msa_hits`
- `homologs` final combinado y ordenado

Funciones/metodos destacados:

- `search_pfam_domains(sequence_file=None)`
- `hmmsearch_with_profile(database_fasta=None)`
- `plot_homologs(top_n=10)`
- `get_top_templates(n=5)`

Parsing interno:

- `_parse_pdb_blast_output(...)`
- `_parse_hmmscan_output(...)`
- `_parse_hmmsearch_output(...)`
- `_parse_pdb_title(...)`

## `src/Homology/domains.py`

### Clase `TemplateProcessor`

Procesa templates seleccionados para anotacion de dominios y preparación de datos para visualizacion.

Responsabilidades:

- Descargar estructuras PDB (top templates).
- Extraer secuencias por cadena.
- Construir FASTA combinado target+templates.
- Consultar dominios InterPro/PDBe en templates.
- Ejecutar `hmmscan` (Pfam) para dominios del target y de secuencias arbitrarias.

Metodos clave:

- `align_templates()`
- `run_hmmscan(sequence, output_prefix="query")`
- `search_sequence_domains(sequence, name="query")`
- `get_target_data()`
- `get_template_data_list()`

## `src/Homology/superimpose.py`

Utilidades para superposicion estructural basada en alineamiento.

Funciones:

- `get_ca_atoms(structure, chain_id)`
- `parse_alignment(alignment_file)`
- `get_aligned_positions(target_aligned, template_aligned)`
- `get_ca_mapping(pdb_id, alignment_file, templates_dir, target_id=None)`
- `superimpose_structures(fixed_structure, moving_structure, fixed_atoms, moving_atoms)`
- `superimpose_template_to_target(...)`
- `write_structure(structure, output_path)`
- `get_match_mismatch_positions(target_aligned, template_aligned)`

### Clase `SuperimpositionVisualizer`

- Superpone un conjunto de templates respecto a uno de referencia.
- Guarda PDBs superpuestos.
- Visualiza en `py3Dmol` coloreando match/mismatch/gaps.

Metodos:

- `superimpose_and_save(template_ids, reference_template=None)`
- `visualize_py3dmol(template_ids, width=800, height=600)`

## 7) UI interactiva

## `src/UI/app.py`

### Clase `ResidueRangeSelect`

Selector de residuos para recorte de estructuras PDB en rangos especificos.

### Clase `CoverageVisualizer`

Componente interactivo (ipywidgets) para:

- Visualizar cobertura de templates contra target.
- Ajustar rangos por slider.
- Seleccionar templates por checkbox.
- Recortar PDBs a region de interes.
- Generar `combined_alignment.pir` para MODELLER.

Funciones auxiliares integradas:

- `_crop_pdb(...)`
- `_create_combined_pir_file(...)`
- `_get_pdb_sequence(...)`
- render HTML de alineamientos y dominios.

### Funcion `run_dssp(pdb_file)`

Ejecuta DSSP y devuelve informacion por residuo (`ss`, `phi`, `psi`, `rsa`).

### Clase `DomainVisualizer`

Visualiza target y templates con:

- pistas de dominios,
- opcion de estructura secundaria (DSSP),
- leyendas y representacion HTML para notebook.

## 8) MODELLER

## `src/modeller/__init__.py`

Re-exporta:

- `generate_single_template_script`
- `generate_multi_template_script`
- `generate_loop_refinement_script`
- `generate_evaluation_script`
- `ModellerRunner`

## `src/modeller/scripts.py`

Generadores de scripts Python para MODELLER:

- `generate_single_template_script(...)`
- `generate_multi_template_script(...)`
- `generate_loop_refinement_script(...)`
- `generate_evaluation_script(...)`
- `generate_alignment_checking_script(...)`

### Clase `ModellerRunner`

- Detecta ejecutable MODELLER disponible.
- Ejecuta scripts (`run_script(...)`).
- Permite guardar scripts (`save_script(...)`).

Funcion de alto nivel:

- `create_modeller_scripts_for_pipeline(...)`
  - Genera lote de scripts (`build_model.py`, `evaluate_model.py`, `check_alignment.py`).

## 9) Analisis estructural avanzado

## `src/Analysis/__init__.py`

Re-exporta funciones y clases de evaluacion y visualizacion:

- Evaluacion: familia, HMM, cationes, sitio activo, DSSP, reparacion, `ModelAssessor`.
- Visualizacion: sitio activo, regiones problematicas, SS, sitio de cation, comparacion estructural, `StructureVisualizer`.

## `src/Analysis/assessment.py`

Funciones principales:

- `identify_protein_family(sequence=None, pdb_file=None, output_prefix="query")`
- `create_hmm_profile(alignment_file, output_file="profile.hmm", name="profile")`
- `align_hmm_to_sequence(hmm_file, sequence_file, output_file="alignment.aln")`
- `add_cation_to_structure(pdb_file, cation="CA", output_file=None, binding_residues=None)`
- `identify_functional_residues(alignment_file, output_file=None, conservation_threshold=0.9)`
- `analyze_active_site(pdb_file, active_site_residues=None, reference_pdb=None)`
- `validate_model_regions(pdb_file, method="prosa")`
- `run_dssp_analysis(pdb_file, output_file=None)`
- `fix_model_problems(pdb_file, problems, output_file=None, method="modeller_loop")`
- `extract_sequence_from_pdb(pdb_file)`

### Clase `ModelAssessor`

Interfaz unificada para evaluacion de modelos:

- `identify_family()`
- `create_hmm(...)`
- `run_dssp(...)`
- `add_cation(...)`
- `find_problematic_regions(...)`
- `fix_problems(...)`
- `get_summary()`

## `src/Analysis/visualization.py`

Funciones de visualizacion en PyMOL:

- `generate_pymol_script(...)`
- `run_pymol_script(...)`
- `visualize_active_site(...)`
- `visualize_problematic_regions(...)`
- `visualize_secondary_structure(...)`
- `visualize_cation_binding(...)`
- `compare_structures(...)`

### Clase `StructureVisualizer`

Wrapper orientado a flujo:

- `active_site(...)`
- `problematic_regions(...)`
- `secondary_structure(...)`
- `cation_site(...)`

## 10) Dependencias y herramientas externas usadas en `src`

Python:

- `pandas`, `matplotlib`, `requests`
- `BioPython` (`SeqIO`, `PDBParser`, `DSSP`, etc.)
- `ipywidgets`, `IPython.display`
- `py3Dmol` (visualizacion 3D interactiva)

Herramientas CLI:

- BLAST+: `psiblast`, `blastdbcmd`, `update_blastdb.pl`
- HMMER: `jackhmmer`, `hmmbuild`, `hmmsearch`, `hmmscan`, `hmmpress`
- `clustalw`
- `dssp`
- `pymol`
- MODELLER

## 11) Flujo recomendado de uso

```python
from src.pipeline import HomologyPipeline

pipeline = HomologyPipeline(
    target_path="target/target.fa",
    e_value_threshold=1e-5,
    num_iterations=5,
    num_templates=5,
    skip_setup=False,
    run_hmmer=True,
    verbose=True
)

print(pipeline.get_summary())
pipeline.visualize()
```

Para evaluacion posterior del modelo generado:

```python
from src.Analysis import ModelAssessor

assessor = ModelAssessor("Models/model.pdb")
assessor.identify_family()
assessor.run_dssp()
print(assessor.get_summary())
```

## 12) Notas de mantenimiento

- Ignorar `__pycache__/` y `*.pyc` en documentacion funcional.
- Mantener sincronizados nombres de campos en DataFrames (`PDB_ID`, `E-value`, `Identity`, `Query_Coverage`).
- Si se cambian rutas de bases de datos, actualizar `Setup`, `TemplateRetriever` y `TemplateProcessor` de forma consistente.
- Si falla alguna herramienta externa, verificar primero que el entorno `AlphaBald` este activado.
