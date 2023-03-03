import logging
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(filename='log_file.log', format=FORMAT)
logger.setLevel(logging.DEBUG)