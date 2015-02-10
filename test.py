from logging_init import *
from module import *

logger = init_log(__name__, console_level='debug', file_level='info')
logger.info('test')
logger.info('Start reading database')
# read database here

records = {'john': 55, 'tom': 66}
logger.debug('Records: %s', records)
logger.info('Updating records ...')
# update records here

amel = Personne('Sediri')

#amel.presente()

#myFunction()

logger.info('Finish updating records')
