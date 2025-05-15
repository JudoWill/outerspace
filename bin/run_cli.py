#!/usr/bin/env python3
"""running cli"""
from outerspace.cli import get_args

def main():
    """running cli"""
    args = get_args()
    print(args)

if __name__ == '__main__':
    main()
