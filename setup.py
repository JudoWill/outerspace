"""For building and installing outerspace"""

__copyright__ = 'Copyright (C) 2025-present, DV Klopfenstein, PhD. All rights reserved'
__author__ = "DV Klopfenstein, PhD"


from os.path import join
from os.path import abspath
from os.path import dirname
from setuptools import setup

PACKAGES = [
    'grna_extraction'
]

PACKAGE_DIRS = {p:join(*p.split('.')) for p in PACKAGES}

def get_long_description():
    """Return the contents of the README.md as a string"""
    dir_cur = abspath(dirname(__file__))
    with open(join(dir_cur, 'README.md'), 'rb') as ifstrm:
        return ifstrm.read().decode("UTF-8")

# These are executables on left and module name:function name
CONSOLE_SCRIPTS = [
    'findseq=grna_extraction.cli:main'
]

REQUIRES = [
    'tomlkit',
    'tqdm',
    'Bio',
    'regex',
    'umi_tools'
]

#KEYWORDS = [
    # TODO: fill with key words that describe this code and it's applications
#]

setup(
    # The name of the project on PyPi
    name='outerspace',
    # https://peps.python.org/pep-0440/
    version='0.0a1',
    author='SB, DVK, REB, WD',
    author_email='wnd22@drexel.edu',
    packages=PACKAGES,
    package_dir=PACKAGE_DIRS,
    entry_points={
        'console_scripts': CONSOLE_SCRIPTS,
    },
    # https://pypi.org/classifiers/
    classifiers=[
        'Programming Language :: Python',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    url='http://github.com/dvklopfenstein/timetracking',
    description="A code to extract sequences from fastq files "
                "when the surrounding sequence is known",
    long_description=get_long_description(),
    long_description_content_type='text/markdown',
    install_requires=REQUIRES,
    # Needed for assignment expressions & the walrus operator
    python_requires='>=3.8',
#    keywords=KEYWORDS,
# TOD0: uncomment out keywords when the keywords are filled in
)

# Copyright (C) 2025-present, DV Klopfenstein, PhD. All rights reserved
