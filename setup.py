# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

# !/usr/bin/env python
import re
import ast
from sys import argv, exit
from json import dumps
from setuptools import setup, find_packages, convert_path


requirements = ["OpenEye-orionplatform==3.1.0",
                "OpenEye-ensemble2img==0.1.1",
                "OpenEye-floereport==0.1.7",
                "OpenEye_oetrajanalysis==0.7.1",
                "openeye-oemdtoolbox==1.0.7",
                "Openeye-toolkits==2020.1.1"]

_version_re = re.compile(r'__version__\s+=\s+(.*)')
version_file = convert_path("MDOrion/__init__.py")
with open(version_file, 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(f.read().decode(
        'utf-8')).group(1)))

if argv[-1] == "--requires":
    print(dumps(requirements))
    exit()

setup(
    name="OpenEye-MD-Floes",
    version=version,
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    include_package_data=True,
    author="Gaetano Calabro, Christopher Bayly",
    author_email="gcalabro@eyesopen.com",
    description='Orion cubes to perform MD and MD analysis',
    install_requires=requirements,
    license='Other/Proprietary License',
    url='https://github.com/oess/openmm_orion',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ],
    entry_points='''
        [console_scripts]
        mdocli=MDOcli.command_line:main
    ''',
    zip_safe=False
)
