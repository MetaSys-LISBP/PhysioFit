import importlib.metadata
import logging

__version__ = importlib.metadata.version("physiofit")
logger = logging.getLogger("physiofit")
logger.setLevel(logging.DEBUG)