# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import os.path

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    readme = f.read()

with open(os.path.join(here, 'LICENSE'), encoding='utf-8') as f:
    my_license = f.read()

setup(
    name='dia_analyze',
    version='0.1',  
    packages=find_packages(exclude=('tests', 'docs')),
    url='https://github.com/seerbio/dia_analyze.git',
    license=Creative Commons Attribution-ShareAlike 4.0 International Public License,
    author='Seer, Inc',
    long_description=readme,
    description='',
    entry_points={'console_scripts': ['dia_compare=dia_analyze.scripts.compare:main',
                                      'dia_convert=dia_analyze.scripts.convert:main'],
                  },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,

    install_requires=[
        'crozes_base>=5.0.3',

    ]
)
