rm -r reports/data/alpha
uv run src/simulation.py
zip -r reports/data/alpha.zip reports/data/alpha ; rm -r reports/.observablehq