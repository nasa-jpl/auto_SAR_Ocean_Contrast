"""
Installation version listing.  Adapted from Rydiqule.
"""

import platform
import OilClassification
import inspect
import os
import sys
import psutil
from pathlib import Path
from importlib.metadata import version


def about():
    """About box describing auto_SAR_Ocean_Contrast and its core dependencies.

    Prints human readable strings of information about the system.

    Examples
    --------
    >>> import OilClassification as oc
    >>> oc.about()
    <BLANKLINE>
      auto_SAR_Ocean_Contrast
    ===========================
    <BLANKLINE>
    Version:              1.3
    Installation Path:    ~/src/auto_SAR_Ocean_Contrast/OilClassification
    <BLANKLINE>
           Dependencies
    ===========================
    <BLANKLINE>
    gdal:                 3.6.2
    jupyter:              1.0.0
    jupyterlab:           4.2.4
    matplotlib:           3.9.1.post1
    netCDF4:              1.7.1.post1
    numpy:                1.26.4
    PyYAML:               6.0.2
    scikit-image:         0.24.0
    scipy:                1.14.0
    spectral:             0.23.1
    Python:               3.12.4
    Python Install Path:  ~/miniconda3/envs/newenv/bin
    Platform Info:        Linux (x86_64)
    CPU Count:            6
    Total System Memory:  16 GB
    """

    home = Path.home()
    install_path = inspect.getsourcefile(OilClassification)
    assert install_path is not None
    this_install_path = Path(install_path).parent
    try:
        this_path = '~' / this_install_path.relative_to(home)
    except ValueError:
        this_path = this_install_path

    python_install_path = Path(sys.executable).parent
    try:
        py_path = '~' / python_install_path.relative_to(home)
    except ValueError:
        py_path = python_install_path

    header = """
      auto_SAR_Ocean_Contrast
    ===========================
    """
    print(header)
    print(f'Version:              {OilClassification.io.get_most_recent_git_tag():s}')

    print(f'Installation Path:    {this_path}')
    dep_header = """
           Dependencies
    ===========================
    """
    print(dep_header)
    print(f'gdal:                 {version("gdal"):s}')
    print(f'jupyter:              {version("jupyter"):s}')
    print(f'jupyterlab:           {version("jupyterlab"):s}')
    print(f'matplotlib:           {version("matplotlib"):s}')
    print(f'netCDF4:              {version("netCDF4"):s}')
    print(f'numpy:                {version("numpy"):s}')
    print(f'PyYAML:               {version("PyYAML"):s}')
    print(f'scikit-image:         {version("scikit-image"):s}')
    print(f'scipy:                {version("scipy"):s}')
    print(f'spectral:             {version("spectral"):s}')
    print(f'Python:               {platform.python_version():s}')
    print(f'Python Install Path:  {py_path}')
    print(f'Platform Info:        {platform.system():s} ({platform.machine():s})')
    print(f'CPU Count:            {os.cpu_count()}')
    print(f'Total System Memory:  {psutil.virtual_memory().total/1024**3:.0f} GB')


if __name__ == "__main__":
    about()
