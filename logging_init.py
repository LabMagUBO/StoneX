import logging
import __main__ as main
print(main.__file__)

def clean_logs():
    print("cleaning")


def set_colors():
    """
        Define the colors in the console
    """
    logging.addLevelName( logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
    logging.addLevelName( logging.DEBUG, "\033[1;34m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
    logging.addLevelName( logging.WARNING, "\033[1;93m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
    logging.addLevelName( logging.ERROR, "\033[1;91m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

def init_log(name, console_level='debug', file_level='info', log_file='logging.log'):
    """
        Initialize the logger.
    """
    # Create the logger
    logger = logging.getLogger(name)

    # Deactivate default console output
    #logger.propagate = False

    # Set the logger's level, to write everything
    logger.setLevel(logging.DEBUG)

    # Formatter
    file_formatter = logging.Formatter("[%(asctime)s] %(name)s - %(module)s.%(funcName)s -:- %(message)s")
    console_formatter = logging.Formatter("%(levelname)s  \t\t %(name)s - %(module)s.%(funcName)s | %(message)s")

    # Console handler, with its own level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # Create a file handler (log file), with write mode, and applying the formatter.
    # Append the handler to the logger
    file_handler = logging.FileHandler(log_file, mode='a')
    file_handler.setFormatter(file_formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)



    #Setting colors
    set_colors()

    return logger

def init_log_bis(name, level='debug'):
    """
        Initialize the logger.
    """
    #Create the logger
    logger = logging.getLogger(name)

    # create a file handler
    handler = logging.FileHandler('logging.log', mode='w')

    # create a logging format

    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    #handler.setFormatter(formatter)

    # add the handlers to the logger

    logger.addHandler(handler)

    #Setting the level
    if level == 'critical':
        logging.basicConfig(level=logging.CRITICAL)
    elif level == 'error':
        logging.basicConfig(level=logging.ERROR)
    elif level == 'warning':
        logging.basicConfig(level=logging.WARNING)
    elif level == 'info':
        logging.basicConfig(level=logging.INFO)
    elif level == 'debug':
        logging.basicConfig(level=logging.DEBUG)
    else:
        #not critical
        logger.warn("WARNING: Logging level is not set : {}".format(__name__))

    #Setting colors
    set_colors()

    return logger
