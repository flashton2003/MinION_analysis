__author__ = 'flashton'

__version__ = '0.1.0'

class BlastRes:
    def __init__(self):
        self.read_name = str
        self.hits = []
        self.read_len = int

class BlastHit:
    def __init__(self):
        self.sbjct = int
        self.score = int
        self.match_len = int
        self.match_pos = int
        self.match_gap = int
        self.query_start = int
        self.query_stop = int
        self.sbjct_start = int
        self.sbjct_stop = int
        self.query_coordinates = []
        self.sbjct_coordinates = []
        self.orientation = str

    def print_res(self, read_name, read_len):
        self.query_start = min(self.query_coordinates)
        self.query_stop = max(self.query_coordinates)
        self.sbjct_start = min(self.sbjct_coordinates)
        self.sbjct_stop = max(self.sbjct_coordinates)
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (read_name, read_len, self.sbjct, self.score, self.match_len, self.match_pos, self.match_gap, self.query_start, self.query_stop, self.sbjct_start, self.sbjct_stop)


def print_res_dict(res_dict):
    """

    :param res_dict:
    """
    #outhandle = open('/Users/flashton/Dropbox/H58_from_iMac/H58/blast_txt', 'w')
    print 'query    read len	subject	orientation	score	match len	match pos	match gap	q start	q stop	s start	s stop'
    for read in res_dict:
    #print each, res_dict[each]
    #print res_dict[each].hits
        for every in res_dict[read].hits:
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (res_dict[read].read_name, res_dict[read].read_len,
                                                                    every.sbjct, every.orientation, every.score, every.match_len, every.match_pos, every.match_gap, every.query_start, every.query_stop, every.sbjct_start, every.sbjct_stop))
    #outhandle.close()