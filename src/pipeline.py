"""
Unified Homology Modeling Pipeline.

This module provides a single entry point for the complete homology modeling workflow:
1. Setup (databases, directories)
2. Target preprocessing
3. Template search
4. Template analysis
5. Visualization

Example Usage
-------------
>>> from src.pipeline import HomologyPipeline
>>> pipeline = HomologyPipeline("target/target.fa")
>>> print(pipeline.homologs.head(10))
>>> pipeline.visualize()
"""

from pathlib import Path
from typing import Optional
import pandas as pd

from src.Setup.setup import Setup
from src.Setup.preprocessing import TargetPreprocessor
from src.Homology.retrieve import TemplateRetriever
from src.Homology.domains import TemplateProcessor


class HomologyPipeline:
    """
    Complete homology modeling pipeline.

    This class orchestrates the full workflow from target sequence to
    ranked template list with domain annotations.

    Parameters
    ----------
    target_path : str or Path, optional
        Path to target FASTA file. Default "target/target.fa"
    e_value_threshold : float, optional
        E-value cutoff for template search. Default 1e-5
    num_iterations : int, optional
        PSI-BLAST iterations. Default 5
    num_templates : int, optional
        Number of templates to analyze in detail. Default 5
    skip_setup : bool, optional
        Skip database setup (if already done). Default False
    run_hmmer : bool, optional
        Include HMMER searches in pipeline. Default True
    verbose : bool, optional
        Print progress messages. Default True

    Attributes
    ----------
    target : TargetPreprocessor
        Preprocessed target information
    retriever : TemplateRetriever
        Template search results
    processor : TemplateProcessor
        Template domain analysis
    homologs : pd.DataFrame
        Ranked template candidates
    templates : pd.DataFrame
        Top templates with domain info

    Example
    -------
    >>> pipeline = HomologyPipeline("target/myprotein.fa")
    >>> # Access ranked templates
    >>> pipeline.homologs.head(10)
    >>> # Get domain information
    >>> pipeline.processor.domain_data
    >>> # Visualize with UI
    >>> pipeline.visualize()
    """

    def __init__(
        self,
        target_path: str = "target/target.fa",
        e_value_threshold: float = 1e-5,
        num_iterations: int = 5,
        num_templates: int = 5,
        skip_setup: bool = False,
        run_hmmer: bool = True,
        verbose: bool = True
    ):
        self.target_path = Path(target_path)
        self.e_value_threshold = e_value_threshold
        self.num_iterations = num_iterations
        self.num_templates = num_templates
        self.verbose = verbose

        # Step 1: Setup
        if verbose:
            print("=" * 60)
            print("HOMOLOGY MODELING PIPELINE")
            print("=" * 60)

        if not skip_setup:
            if verbose:
                print("\n[1/4] Setting up directories and databases...")
            self.setup = Setup(skip_databases=False, clear_temp=True)
        else:
            if verbose:
                print("\n[1/4] Skipping setup (using existing databases)...")
            self.setup = Setup(skip_databases=True, clear_temp=True)

        # Step 2: Preprocess target
        if verbose:
            print("\n[2/4] Loading target sequence...")
        self.target = TargetPreprocessor(target_dir=str(self.target_path.parent))
        if verbose:
            print(f"  Target: {self.target.name} ({len(self.target.sequence)} aa)")

        # Step 3: Search for templates
        if verbose:
            print("\n[3/4] Searching for templates...")
        self.retriever = TemplateRetriever(
            target_path=str(self.target_path),
            e_value_threshold=e_value_threshold,
            num_iterations=num_iterations,
            num_hits=20,  # Use more hits for MSA
            run_hmmer=run_hmmer,
            verbose=verbose
        )

        # Step 4: Analyze top templates
        if verbose:
            print("\n[4/4] Analyzing top templates...")
        if len(self.retriever.homologs) > 0:
            self.processor = TemplateProcessor(
                homologs=self.retriever.homologs,
                num_templates=num_templates,
                target_path=str(self.target_path),
                verbose=verbose
            )
        else:
            self.processor = None
            if verbose:
                print("  Warning: No templates found!")

        # Expose key results
        self.homologs = self.retriever.homologs
        self.templates = self.retriever.get_top_templates(num_templates)

        if verbose:
            print("\n" + "=" * 60)
            print("PIPELINE COMPLETE")
            print(f"  Found {len(self.homologs)} potential templates")
            if len(self.homologs) > 0:
                best = self.homologs.iloc[0]
                print(f"  Best hit: {best['PDB_ID']} (E={best['E-value']:.2e}, {best['Identity']:.1f}% identity)")
            print("=" * 60)

    def plot_homologs(self, top_n: int = 10):
        """
        Plot distribution of homolog protein types and organisms.

        Parameters
        ----------
        top_n : int
            Number of top entries to show
        """
        self.retriever.plot_homologs(top_n)

    def visualize(self, num_templates: Optional[int] = None):
        """
        Launch interactive visualization widgets.

        Parameters
        ----------
        num_templates : int, optional
            Number of templates to visualize. Uses pipeline default if None.
        """
        try:
            from src.UI.app import CoverageVisualizer, DomainVisualizer
        except ImportError:
            print("UI module not available. Install ipywidgets for visualization.")
            return

        n = num_templates or self.num_templates

        if len(self.homologs) == 0:
            print("No templates available for visualization.")
            return

        # Prepare domain data for CoverageVisualizer
        target_domains = []
        template_domains = {}

        if self.processor:
            target_domains = self.processor.target_domains or []
            # Convert domain_data list to dict keyed by PDB_ID
            for template_data in self.processor.domain_data:
                pdb_id = template_data[0]
                domains = template_data[2]
                template_domains[pdb_id] = domains

        print("\n--- Coverage Visualization ---")
        cv = CoverageVisualizer(
            self.homologs,
            number_of_templates=n,
            target_domains=target_domains,
            template_domains=template_domains
        )
        cv.show()

        if self.processor and self.processor.domain_data:
            print("\n--- Domain Visualization ---")
            target_data = self.processor.get_target_data()
            template_data = self.processor.get_template_data_list()
            dv = DomainVisualizer(target_data, template_data)
            dv.show()

    def get_summary(self) -> dict:
        """
        Get pipeline summary as dictionary.

        Returns
        -------
        dict
            Summary with target info, template counts, and best hit
        """
        summary = {
            "target_name": self.target.name,
            "target_length": len(self.target.sequence),
            "total_templates": len(self.homologs),
            "search_evalue": self.e_value_threshold,
        }

        if len(self.homologs) > 0:
            best = self.homologs.iloc[0]
            summary["best_hit"] = {
                "pdb_id": best["PDB_ID"],
                "evalue": best["E-value"],
                "identity": best["Identity"],
                "coverage": best["Query_Coverage"],
                "protein_name": best.get("Protein Name", "Unknown")
            }

        return summary

    def export_templates(self, output_file: str = "templates_summary.csv"):
        """
        Export template list to CSV file.

        Parameters
        ----------
        output_file : str
            Output file path
        """
        self.homologs.to_csv(output_file, index=False)
        print(f"Templates exported to {output_file}")

    def __repr__(self) -> str:
        return (
            f"HomologyPipeline(\n"
            f"  target='{self.target.name}' ({len(self.target.sequence)} aa),\n"
            f"  templates={len(self.homologs)},\n"
            f"  e_value<{self.e_value_threshold}\n"
            f")"
        )


def run_pipeline(
    target_path: str = "target/target.fa",
    e_value: float = 1e-5,
    num_templates: int = 5
) -> HomologyPipeline:
    """
    Convenience function to run the complete pipeline.

    Parameters
    ----------
    target_path : str
        Path to target FASTA file
    e_value : float
        E-value threshold
    num_templates : int
        Number of templates to analyze

    Returns
    -------
    HomologyPipeline
        Completed pipeline instance
    """
    return HomologyPipeline(
        target_path=target_path,
        e_value_threshold=e_value,
        num_templates=num_templates
    )
