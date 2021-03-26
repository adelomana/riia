import logging, datetime, socket, time
#import data_load_module, coexpression_module

# def coexpression(expression_data):
#
#     """
#     This function builds coexpression clusters
#     """
#
#     coexpression_clusters = coexpression_module.main(expression_data)
#
#     return coexpression_clusters
#
# def data_load(expression_path):
#
# 	expression_data = data_load_module.main(expression_path)
#
# 	return expression_data

###
### MAIN
###

hostname = socket.gethostname()
now = time.strftime('%Y-%m-%d %H:%M:%S')
working_format = '{} {} %(levelname)s | %(message)s'.format(hostname, now)
logging.basicConfig(format=working_format)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

logging.info('riia version 0.0.0')
