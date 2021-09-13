class ArgParser:
    def __init__(self, args: dict):
        self.args = args
        self.val = args

    def parse(self, argv):
        for arg in argv[1:]:
            if self.args[arg]:
                self.val[arg] = True
        for key in self.val:
            if type(self.val[key]) == str:
                self.val[key] = False
        return self.val
