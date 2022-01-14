from ._version import get_versions
__version__ = "NotMe"
del __version__
from ._parse_version import *
del get_versions
