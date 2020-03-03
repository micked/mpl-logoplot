#!/usr/bin/env python

import io
import os
import re

from setuptools import setup


classifiers =[
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: BSD",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Visualization",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
]


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_long_description():
    return _read('README.rst')


setup(
    name='mpl_logoplot',
    author='Michael Schantz Klausen',
    author_email='sch@ntz.nu',
    version='0.1',
    license='BSD',
    description='Sequence logo plotting for matplotlib',
    long_description=get_long_description(),
    keywords=['matplotlib', 'logoplot', 'visualization'],
    url='https://github.com/themicked/mpl-logoplot',
    packages=['mpl_logoplot'],
    zip_safe=False,
    classifiers=classifiers,
    install_requires=[
        'numpy',
        'matplotlib'
    ],
)
