import time
from StoneX.logging_init import *

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
    logger = init_log(__name__)
    timeLapse = time.time()-start_time
    minutes = int(timeLapse / 60)
    seconds = timeLapse % 60
    logger.info("Temps écoulé : %s min %s s", minutes, seconds)


## Function to launched at start
start_time = launch_time(True)
