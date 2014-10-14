__author__ = 'flashton'





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


