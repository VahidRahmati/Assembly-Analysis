class vahidSucks(object):
    def __init__(self, vahid_sucks=True):
        self.vahid_sucks = vahid_sucks

    def print_absolute_truth(self,):
        print('Does Vahid suck?', self.vahid_sucks)

class howMuchDoesHeSuck(vahidSucks):
    def __init__(self, vahid_sucks=True, suckiness=100):
        # This next line calls the __init__ of the parent class, so vahid_sucks
        # is saved to self too.
        super(howMuchDoesHeSuck, self).__init__(vahid_sucks)
        self.suckiness = suckiness

    def print_suckiness(self,):
        self.print_absolute_truth()
        print('How much does he suck?', self.suckiness)


v = vahidSucks(vahid_sucks=True)
vv = howMuchDoesHeSuck(True, 1000)