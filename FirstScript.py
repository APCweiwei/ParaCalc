#staticmethod example

import time

class TimeTest(object):

    def __init__(self, hour, minute, second):
        self.hour = hour
        self.minute = minute
        self.second = second

    @staticmethod
    def showTime(cls):
        return time.strftime("%H:%M:%S", time.localtime())

print(TimeTest.showTime())
t = TimeTest(2, 10, 10)
nowTime = t.showTime()
print(nowTime)
