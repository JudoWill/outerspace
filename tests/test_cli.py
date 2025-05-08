#!/usr/bin/env python3
"""testing cli"""
from sys import argv
from os.path import exists
from logging import debug
from logging import DEBUG
from logging import basicConfig
# basicConfig (level = DEBUG)
from pytest import raises
from grna_extraction.cli import Cli
from grna_extraction.cli import main
from tests.pkgtest.utils import get_filename

# cdef test_cli_empty():
#    """testing no arguments"""
## TODO add new test where there is no arguments, system exit
#    with raises(SystemExit) as excinfo:
#        cli=Cli(args=[])
#        args = cli.args
#        print(args)
#    ## testing for error code 1 "systemexit 1" = missing read files
#    ## to see if there are NO args -- we want an error here
#    ## args with no FASTQ files should result in an error
#    assert excinfo.value.code == 1
#    print("TEST PASSED")


def test_cfg():
    """testing all cases of configuration - was it given? does it exist?"""
    filenamecfg = get_filename("tests/configs/grnaquery.cfg")
    args = [
        filenamecfg,
        '-1', get_filename("tests/data/409-4_S1_L002_R1_001.fastq.gz"),
        '-2', get_filename("tests/data/409-4_S1_L002_R2_001.fastq.gz"),
        # TODO: consider putting output to tests in a working directory that will be deleted
        '-o', get_filename("409-4_S1_L002_R1_R2_output.csv")
    ]
    cli = Cli(args)
    cli.run()
    # filename exists
    assert exists(filenamecfg)
    # did we get a config out of it
    assert cli.top is not None
    # did the filename get transferred into the object
    assert filenamecfg == cli.top.cfg.filename
    
def test_main():
    """testing main"""
    print(f'ARGV: {argv}')
    with raises(SystemExit) as excinfo:
        main(args=[])
        assert excinfo.value.code == 0 


if __name__ == '__main__':
    #test_cli()
    #test_cfg()  # TODO: Uncomment this after test
    test_main()
