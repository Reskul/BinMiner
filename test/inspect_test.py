import inspect as i


class Test:
    def __init__(self):
        self.x = 42

    @staticmethod
    def lol():
        print(f"{i.currentframe().f_code.co_name}::#############")
        print(i.stack())
        cf = i.currentframe()
        print(type(cf))
        print(f"CurrentFrame:{cf}\nLineNo:{cf.f_lineno}\nCode:{cf.f_code.co_name}")


    class Nested:
        def __init__(self):
            print(i.stack())
            cf = i.currentframe()
            print(f"{cf}\n{cf.f_lineno}\n{cf.f_code.co_name}")


if __name__ == '__main__':
    print(i.stack())
    cf = i.currentframe()
    print(type(cf))
    print(f"CurrentFrame:{cf}\nLineNo:{cf.f_lineno}\nCode:{cf.f_code}")
    print("#############")
    # t = Test()
    Test.lol()
    print("#############")
    tn = Test.Nested()
