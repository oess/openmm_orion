#!/usr/bin/env python
from __future__ import absolute_import, print_function
import re, ast, os, yaml, json
from os.path import relpath, join
from pip.req import parse_requirements
from setuptools import setup, find_packages
from collections import OrderedDict
import yaml

def from_env_yml(filename):
    with open(filename, 'r') as fp:
        yamlstr = fp.read()
        data = yaml.load(yamlstr)
        yml = { 'name': None, 'channels': None, 'dependencies': None }
        for k, v in yml.items():
            yml[k] = data[k]
        return yml

env_yml = from_env_yml('.environment.yml')
man_json = {
    "name": "OpenMMCubes",
    "version": "0.1.6",
    "python_version": "3.5",
    "conda_channels": env_yml['channels'],
    "conda_dependencies": env_yml['dependencies'],
    "requirements": "requirements.txt"
    }

with open('manifest.json', 'w') as outjson:
    json.dump(man_json, outjson, indent=4, ensure_ascii=False)

def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]

try:
    install_reqs = get_reqs(parse_requirements(man_json['requirements']))
except TypeError:
    from pip.download import PipSession
    install_reqs = get_reqs(
        parse_requirements(man_json['requirements'], session=PipSession())
    )

def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('OpenMMCubes/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version



setup(
    name=man_json['name'],
    version=man_json['version'],
    packages=find_packages(include=['examples'], exclude=['tests*']),
    include_package_data=True,
    package_data={ 'examples' : find_package_data('examples/data', 'examples')},
    author="Christopher Bayly, Nathan M. Lim, John Chodera",
    author_email="bayly@eyesopen.com",
    description='Prepare complex for MD with OpenMM',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ]
)
