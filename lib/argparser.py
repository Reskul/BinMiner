import sys


class ArgParser:
    def __init__(self, args: list):
        self.args = args
        self.val = {}

    def parse(self, argv):
        for arg in argv[1:]:
            if arg in self.args:
                self.val[arg] = True
        for key in self.args:
            if key not in self.val:
                self.val[key] = False
        return self.val


if __name__ == "__main__":
    expected_args = ['-d', '-t']
    ap = ArgParser(expected_args)

    given_args = ap.parse(sys.argv)
    print(given_args)
