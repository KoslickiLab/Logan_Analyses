# Makefile for Component Metadata Analysis
# Makes it easy to run common analysis tasks

.PHONY: all quick full stats explore queries clean help

# Default target - show help
help:
	@echo "Component Metadata Analysis - Available targets:"
	@echo ""
	@echo "  make all       - Run complete analysis pipeline"
	@echo "  make quick     - Quick insights summary only (fast)"
	@echo "  make full      - Comprehensive metadata analysis"
	@echo "  make stats     - Statistical enrichment testing"
	@echo "  make explore   - Outlier detection and exploration"
	@echo "  make queries   - Example queries"
	@echo "  make clean     - Remove all generated output files"
	@echo "  make help      - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make quick     # Fast overview"
	@echo "  make all       # Complete analysis"

# Run everything
all:
	python run_all_analyses.py

# Quick overview
quick:
	python quick_insights.py

# Individual analysis steps
full:
	python analyze_component_metadata.py

stats:
	python statistical_enrichment_analysis.py

explore:
	python explore_components_interactive.py

queries:
	python example_queries.py

# Clean up generated files
clean:
	@echo "Removing generated files..."
	rm -f component_metadata_summary.csv
	rm -f merged_components_metadata.csv
	rm -f statistical_enrichment_results.csv
	rm -f pairwise_field_interactions.csv
	rm -f geographic_outliers.csv
	rm -f bioproject_mixtures.csv
	rm -f temporal_patterns.csv
	rm -f *.png
	@echo "Clean complete!"

# Check that input files exist
check:
	@echo "Checking for required input files..."
	@test -f components.json || (echo "ERROR: components.json not found" && exit 1)
	@test -f accessions_mbases_geq_10.txt || (echo "ERROR: accessions_mbases_geq_10.txt not found" && exit 1)
	@test -f metagenome_metadata_with_geo.duckdb || (echo "ERROR: metagenome_metadata_with_geo.duckdb not found" && exit 1)
	@echo "All required files found!"
