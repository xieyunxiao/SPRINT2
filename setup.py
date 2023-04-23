"""

SPRINT2

"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path


setup(
    name='sprint2',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'argparse',
    ],
    entry_points={
        'console_scripts': [
            'MaskAGref=MaskAGref:main',
            'getDsRNA=getDsRNA:run',
            'getRES=getRES:run',
        ],
    },
    author='Feng Zhang, Yunxiao Xie',
    author_email='20210700107@fudan.edu.cn',
    description='SPRINT-2.0: An enhanced tool to identify RNA editing sites',
    url='https://github.com/xieyunxiao/SPRINT2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)

