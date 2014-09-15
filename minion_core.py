__author__ = 'flashton'

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


class MinionRead:
    def __init__(self):
        pass

def print_res_dict(res_dict):
    """

    :param res_dict:
    """
    print 'query\tquery_len	subject	orientation	score	match len	match pos	match gap	q start	q stop	s start	s stop'
    for read in res_dict:
    #print each, res_dict[each]
    #print res_dict[each].hits
        for every in res_dict[read].hits:
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (res_dict[read].read_name, res_dict[read].read_len, every.sbjct, every.orientation, every.score, every.match_len, every.match_pos, every.match_gap, every.query_start, every.query_stop, every.sbjct_start, every.sbjct_stop)


def find_best_hits(res_dict):

    for read in res_dict:

        hit_dict = {}
        #print res_dict[read].read_name

        for hit in res_dict[read].hits:
            #print (hit.query_start, hit.query_stop)
            hit_dict[hit] = range(hit.query_start, hit.query_stop)


        #print hit_dict

        overlap_dict = {}

        for hit1 in hit_dict:
            overlap_dict[hit1] = {}
            for hit2 in hit_dict:
                if hit1 != hit2:
                    if hit2 in overlap_dict:
                        if hit1 not in overlap_dict[hit2]:
                            overlap_dict[hit1][hit2] = len(set(hit_dict[hit1]).intersection(hit_dict[hit2]))

        #print overlap_dict
        to_delete = []
        if len(overlap_dict) > 1:

            for hit1 in overlap_dict:
                for hit2 in overlap_dict[hit1]:
                    if hit1 != hit2:
                        #print hit1, hit1.match_len, hit2, hit2.match_len, overlap_dict[hit1][hit2]



                        hits = [hit1, hit2]

                        for h in hits:
                        #    #print h.match_len, overlap_dict[hit1][hit2] * 0.75
                            #print h.match_len
                            if h.match_len * 0.75 < overlap_dict[hit1][hit2]:

                                if hit1.score < hit2.score:
                                    to_delete.append(hit1)

                                if hit2.score < hit1.score:
                                    to_delete.append(hit2)

                                if hit1.score == hit2.score:
                                    to_delete = []
                                    to_delete.append(hit2)



        #print res_dict[read].hits
        #print 'to delete', to_delete
        for h in to_delete:
            for i, hit in enumerate(res_dict[read].hits):
                if hit == h:
                    #print len(res_dict[read].hits)
                    del res_dict[read].hits[i]
                    #print 'left', res_dict[read].read_name, res_dict[read].hits

    return res_dict
