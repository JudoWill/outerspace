DIR_RB_READS := ../../nonn-lab/rachel-test-crispr/

py:
	find grna_extraction -name '*.py' -type f

run_script:
	cd grna_extraction; python extraction_attempt.py

reads:
	find $(DIR_RB_READS)/reads

hello:
	ls -l $(DIR_RB_READS)

