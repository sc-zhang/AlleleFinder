import time


class Message:
    def __init__(self):
        pass

    @staticmethod
    def info(info):
        print("\033[32m%s\033[0m %s" % (time.strftime('[%H:%M:%S]', time.localtime(time.time())), info))

    @staticmethod
    def warn(info):
        print("\033[33m%s\033[0m %s" % (time.strftime('[%H:%M:%S]', time.localtime(time.time())), info))

    @staticmethod
    def error(info):
        print("\033[31m%s\033[0m %s" % (time.strftime('[%H:%M:%S]', time.localtime(time.time())), info))
