DIR_RB_READS := ../../nonn-lab/rachel-test-crispr/

py:
	find grna_extraction bin -name '*.py' -type f

run_script:
	cd grna_extraction; python extraction_attempt.py

reads:
	find $(DIR_RB_READS)/reads

hello:
	ls -l $(DIR_RB_READS)

clean:
	rm -f grna_extraction/.ipynb_checkpoints/extraction_attempt-checkpoint.py
