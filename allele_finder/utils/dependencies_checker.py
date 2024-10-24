from subprocess import Popen, PIPE
from allele_finder.utils.message import Message


class DependChecker:
    def __init__(self):
        self.__res = ""
        self.__err = ""

    def check(self, cmd):
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding='utf-8')
        self.__res, self.__err = p.communicate()

        if self.__res:
            return True
        return False

    def exit_program_with_errors(self):
        Message.error("\tStdout:")
        for _ in self.__res.strip().split('\n'):
            Message.error("\t%s" % _)
        Message.error("\tStderr:")
        for _ in self.__err.strip().split('\n'):
            Message.error("\t%s" % _)
        exit(-1)
