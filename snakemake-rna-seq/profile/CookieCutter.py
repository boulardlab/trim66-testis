#
# Based on lsf CookieCutter.py
#
import os
import json
from datetime import date

d = os.path.dirname(__file__)
with open(os.path.join(d, "settings.json")) as fh:
    settings = json.load(fh)


class CookieCutter:

    SBATCH_DEFAULTS = settings['SBATCH_DEFAULTS']
    CLUSTER_NAME = settings['CLUSTER_NAME']
    CLUSTER_CONFIG = settings['CLUSTER_CONFIG']
    ADVANCED_ARGUMENT_CONVERSION = settings['ADVANCED_ARGUMENT_CONVERSION']
    LOG_FOLDER = settings['LOG_FOLDER']

    @staticmethod
    def get_cluster_option() -> str:
        cluster = CookieCutter.CLUSTER_NAME
        if cluster != "":
            return f"--cluster={cluster}"
        return ""

    @staticmethod
    def get_advanced_argument_conversion() -> bool:
        val = {"yes": True, "no": False}[
            CookieCutter.ADVANCED_ARGUMENT_CONVERSION
        ]
        return val

    @staticmethod
    def get_formatted_log_folder() -> str:
        """ Get log folder from config and append date"""
        log_folder = CookieCutter.LOG_FOLDER
        if log_folder != "":
            today = date.today().isoformat()
            log_folder = os.path.join(log_folder, today)
        return log_folder
