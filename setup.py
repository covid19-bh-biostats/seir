#!/usr/bin/env python3

import seir as my_pkg
from setuptools import setup, find_packages, Command
import os
import subprocess
import re
import sys


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class MypyCommand(Command):
    """ Class for running the MYPY static type checker on the package """

    description = 'Runs the MYPY static type checker on ' + my_pkg.__name__
    user_options = []

    def initialize_options(self):
        ...

    def finalize_options(self):
        ...

    def run(self):
        import mypy.api
        res = mypy.api.run(['-p', my_pkg.__name__, '--ignore-missing-imports'])
        print(res[0])
        return res[1]

def get_requirements():
    with open('requirements.txt', 'r') as f:
        return f.readlines()

setup(
    name=my_pkg.__name__,
    author=my_pkg.__author__,
    author_email=my_pkg.__author_email__,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: Other/Proprietary License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7'
    ],
    cmdclass={
        'mypy': MypyCommand,
    },
    data_files=[],
    description="Modeling of epidemics using the SEIR model",
    entry_points = {
        'console_scripts': [
            'seir-model = seir.scripts.cli:main'
        ]
    },
    install_requires=get_requirements(),
    license=my_pkg.__license__,
    long_description=read('README.rst'),
    packages=find_packages(),
    python_requires='>=3.7',
    test_suite='nose2.collector.collector',
    tests_require=['nose2', 'mypy'],
    url='https://https://github.com/covid19-bh-biostats/seir',
    version=my_pkg.__version__,
    zip_safe=True)
