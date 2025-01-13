output_dir="alpha"

uv run src/simulation.py --output-dir reports/data/$output_dir
uv run src/metrics.py --input-file reports/data/$output_dir/lattice.json

cd reports/data/$output_dir
zip -r data.zip *
cd ../../../
rm -r reports/.observablehq