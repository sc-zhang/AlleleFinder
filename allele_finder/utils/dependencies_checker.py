from subprocess import Popen, PIPE


class DependChecker:
    def __init__(self):
        pass

    @staticmethod
    def check(cmd):
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding='utf-8')
        res, _ = p.communicate()

        if res:
            return True
        return False
