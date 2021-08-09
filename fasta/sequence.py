##
# This Class is derived from str to use the sequences as if they were just strings, but also have additional information
class Sequence(str):
    def __new__(cls, seq, unique_id, header=None):
        obj = str.__new__(cls, seq)
        obj.unique_id = unique_id
        if header is not None:
            obj.header = header
        return obj

    def add_header(self, header):
        self.header = header