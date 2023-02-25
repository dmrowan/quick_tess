from setuptools import setup, find_packages
from distutils.util import convert_path

setup(name = "quick_tess",
    version = 0.1,
    description = "Color-Magnitude Diagram Analysis tools for binary stars",
    author = "Dom Rowan",
    author_email = "",
    url = "https://github.com/dmrowan/quick_tess",
    packages = find_packages(include=['quick_tess', 'quick_tess.*']),
    package_data = {'quick_tess':['data/*']},
    include_package_data = True,
    classifiers=[
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["astropy", "cmasher", "matplotlib", "numpy", "pandas", "scipy"]
)
