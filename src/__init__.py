"""
AlphaBald - Homology Modeling Pipeline

A streamlined pipeline for protein homology modeling using PSI-BLAST,
HMMER, and MODELLER.

Quick Start
-----------
>>> from src import HomologyPipeline
>>> pipeline = HomologyPipeline("target/target.fa")
>>> pipeline.homologs  # View ranked templates
>>> pipeline.visualize()  # Launch interactive UI
"""

from src.pipeline import HomologyPipeline, run_pipeline
from src.Setup.setup import Setup
from src.Setup.preprocessing import TargetPreprocessor
from src.Homology.retrieve import TemplateRetriever
from src.Homology.domains import TemplateProcessor
from src.Homology.superimpose import SuperimpositionVisualizer

__all__ = [
    "HomologyPipeline",
    "run_pipeline",
    "Setup",
    "TargetPreprocessor",
    "TemplateRetriever",
    "TemplateProcessor",
    "SuperimpositionVisualizer",
]

__version__ = "1.0.0"
