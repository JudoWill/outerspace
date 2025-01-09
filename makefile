DIR_RB_READS := ../../nonn-lab/rachel-test-crispr/

run_script:
	cd grna_extraction; python extraction_attempt.py

reads:
	find $(DIR_RB_READS)/reads

hello:
	ls -l $(DIR_RB_READS)

