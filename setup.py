#!/usr/bin/env python

import io
import os
import re

from setuptools import setup


classifiers =[
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: MIT",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Visualization",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
]


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read('mpl_logoplot.py'),
        re.MULTILINE).group(1)
    return version


def get_long_description():
    return _read('README.rst')


setup(
    name='mpl_logoplot',
    author='Michael Schantz Klausen',
    author_email='sch@ntz.nu',
    version=get_version(),
    license='BSD',
    description='Sequence logo plotting for matplotlib',
    long_description=get_long_description(),
    keywords=['matplotlib', 'logoplot', 'visualization'],
    url='https://github.com/themicked/mpl-logoplot',
    py_modules=['mpl_logoplot'],
    zip_safe=False,
    classifiers=classifiers,
    install_requires=[
        'numpy',
        'matplotlib'
    ],
)
