#!/usr/bin/env python

from distutils.core import setup

setup(name="seqstats",
      version="0.0.1",
      author="Edmund Lau",
      author_email="edmundlth95@gmail.com",
      packages=["seqstats"],
      entry_points={
          "console_scripts":["seqstats = seqstats.seqstats:main"]
          },
      url="",
      licence="LICENCE",
      install_requires=["Biopython"],
      )
