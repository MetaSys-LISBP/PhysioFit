"""Logger module containing the PhysioFit logger setup"""

import logging

# Setup base logger

mod_logger = logging.getLogger("PhysioFit_Logger")
mod_logger.setLevel(logging.DEBUG)


def initialize_fitter_logger(verbose):
    logger = logging.getLogger(f"PhysioFit_Logger.base.fitter.PhysioFitter")
    logger.setLevel(logging.DEBUG)

    return logger
