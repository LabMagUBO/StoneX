import logging
import logging.config
import time


def init_log(log_conf):
    logging.addLevelName( logging.INFO, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.INFO))
    logging.addLevelName( logging.DEBUG, "\033[1;34m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
    logging.addLevelName( logging.WARNING, "\033[1;93m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
    logging.addLevelName( logging.ERROR, "\033[1;91m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

    logging.config.fileConfig(log_conf, disable_existing_loggers=False)

def launch_log():
    init_log('logging.ini')
    logger = logging.getLogger(__name__)
    logger.info("Modules Stoner chargé")


def launch_time(record_time):
    """
        Record the starting time.
    """

    if record_time:
        return time.time()

def end_time(start_time):
    """
        Print the time
    """
    logger = logging.getLogger(__name__)

    timeLapse = time.time()-start_time
    minutes = int(timeLapse / 60)
    seconds = timeLapse % 60
    logger.info("Temps écoulé : %s min %s s", minutes, seconds)


## Function to launched at start
start_time = launch_time(True)
launch_log()
