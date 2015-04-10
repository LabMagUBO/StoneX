import logging
from StoneX.constants import *


def set_colors():
    """
        Define the colors in the console.
    """
    logging.addLevelName( logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
    logging.addLevelName( logging.DEBUG, "\033[1;34m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
    logging.addLevelName( logging.WARNING, "\033[1;93m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
    logging.addLevelName( logging.ERROR, "\033[1;91m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

def init_log(name, console_level='debug', file_level='info', log_file=main_file, mode='a'):
    """
        Initialize the logger.

        Mode : append by default. In order to clean the log, need to set the fist logger to 'w' mode.

        Usage: initialize a new logger at each level.
            logger = init_log(__name__)
        then
            logger.info() / logger.warn() / ...
    """
    # Create the logger
    logger = logging.getLogger(name)

    # Deactivate default console output
    #logger.propagate = False
    if not logger.handlers:

        # Set the logger's level, to write everything
        logger.setLevel(logging.DEBUG)

        # Formatter
        file_formatter = logging.Formatter("[%(asctime)s] %(name)s - %(module)s.%(funcName)s -:- %(message)s")
        console_formatter = logging.Formatter("%(levelname)s\t%(name)s - %(module)s.%(funcName)s \t-> %(message)s")

        # Console handler, with its own level
        console_handler = logging.StreamHandler()
        console_handler.setLevel(set_level(console_level))
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        # Create a file handler (log file), with write mode, and applying the formatter.
        # Append the handler to the logger
        file_handler = logging.FileHandler(log_file + '.log', mode=mode)
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(set_level(file_level))
        logger.addHandler(file_handler)

        #Setting colors
        set_colors()

    return logger

def set_level(level):
    if level == 'critical':
        return logging.CRITICAL
    elif level == 'error':
        return logging.ERROR
    elif level == 'warning':
        return logging.WARNING
    elif level == 'info':
        return logging.INFO
    elif level == 'debug':
        return logging.DEBUG
    else:
        #not critical
        logger.warn("WARNING: Logging level is not set : {}".format(__name__))
        return logging.DEBUG
