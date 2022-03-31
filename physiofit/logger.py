"""Logger module containing the PhysioFit logger setup"""

import logging

# Setup base logger

mod_logger = logging.getLogger("PhysioFit_Logger")
mod_logger.setLevel(logging.DEBUG)


def initialize_fitter_logger(verbose):
    logger = logging.getLogger(f"PhysioFit_Logger.base.fitter.PhysioFitter")
    # fh = logging.FileHandler(f"{self.run_name}.log")
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    # For debugging purposes
    if verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    if not logger.hasHandlers():
        logger.addHandler(handler)

    return logger
