DIR_RB_READS := /data/share/nonn-lab/rachel-test-crispr/reads
P1 := $(DIR_RB_READS)/409-4_S1_L001_R1_001.fastq.gz
P2 := $(DIR_RB_READS)/409-4_S1_L001_R2_001.fastq.gz

run:
	bin/main.py $(P1) $(P2)

files:
	@ls $(P1)
	@ls $(P2)

py:
	find original tests grna_extraction bin -name '*.py' -type f | grep -v checkpoint

md:
	find . -type f -name '*.md'

vim:
	vim -p grna_extraction/extraction_attempt.py bin/main.py

run_script:
	cd grna_extraction; python extraction_attempt.py

reads:
	find $(DIR_RB_READS)/reads

hello:
	ls -l $(DIR_RB_READS)

clean:
	rm -f grna_extraction/.ipynb_checkpoints/extraction_attempt-checkpoint.py
	rm -f testing.cfg
