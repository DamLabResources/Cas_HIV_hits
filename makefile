test:
	pytest tests/*.py
    
prepare:
	dvc pull data/GRCh38.fna.gz.dvc