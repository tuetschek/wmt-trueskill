#!/usr/bin/env python
#encoding: utf-8

__author__ = "Keisuke Sakaguchi"
__version__ = "0.1"

#Input: JUDGEMENTS.csv which must contain one language-pair judgements.
#Output: *_mu_sigma.json: Mu and Sigma for each system
#        *.count: number of judgements among systems (for generating a heatmap) if -n is set to 2 and -e.

import sys
import os
import argparse
import random
import json
import numpy as np
import math
import scripts.next_comparison
from itertools import combinations
from collections import defaultdict
from csv import DictReader
from trueskill import *
from math import factorial

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('prefix', help='output ID (e.g. fr-en0)')
arg_parser.add_argument('input_file', nargs='?', help='input file (default to stdin)', type=str)
arg_parser.add_argument('-n', '--freeN', type=int,
        help='Free-for-All N (2-5)', required=True)
arg_parser.add_argument('-d', '--dp', type=int,
        help='Number of judgments to use (0 == all)', required=True)
arg_parser.add_argument('-p', '--dp_pct', type=float, default=1.0,
        help='Percentage of judgments to use (default: 1.0)')
arg_parser.add_argument('-s', '--num_systems', type=int, default=5,
        help='Number of systems in one ranking in CSV file (default=5)')
arg_parser.add_argument('-e', '--heat', default=False, action="store_true",
        help='Produce a file for generating a heatmap (default=False)')
args = arg_parser.parse_args()

#######################################
### Global Variables and Parameters ###
param_sigma = 0.5
param_tau = 0.
draw_rate = 0.25

# You can set arbitrary number(s) for record (dp is the number assigned by -d).
#num_record = [int(args.dp*0.9), args.dp]
num_record = [ args.dp ]
#e.g. num_record = [args.dp*0.125, args.dp*0.25, args.dp*0.5, args.dp*0.9, args.dp]
#e.g. num_record = [400, 800, 1600, 3200, 6400, 11520, 12800]

# When -n is set to 2, you can set beginning and ending between (0 and 1) for counting the number of comparisons among systems.
# This is used for generating a heatmap.
# e.g. "count_begin=0.4 and count_end=0.6" records the number of comparisons from 40% to 60% of total comparisons.
count_begin = 0.8
count_end = 1.0
if count_begin > count_end:
    raise
#######################################

comparison_d = defaultdict(list)


def parse_csv(fh=sys.stdin):
    """Parse the WMT-formatted CSV file and return system names and rank(1-5)
    for each sentence."""
    all_systems = []
    sent_sys_rank = defaultdict(list)
    for i,row in enumerate(DictReader(fh)):
        sentID = int(row.get('segmentId'))
        systems = []
        ranks = []
        for num in range(1, args.num_systems+1):
            if row.get('system%dId' % num) in all_systems:
                pass
            else:
                all_systems.append(row.get('system%dId' % num))
            systems.append(row.get('system%dId' % num))
            ranks.append(int(row.get('system%drank' % num)))
        if -1 in ranks:
            pass
        else:
            sent_sys_rank[sentID].append({'systems': systems, 'ranks': ranks})
    return all_systems, sent_sys_rank


def fill_comparisons(sent_sys_rank):
    """Convert all multi-way ranking data to freeN-wise comparisons (using all subset
    combinations for all sentences). Store this in the global comparison_d variable."""
    sentIDs = sent_sys_rank.keys()
    for sid in sentIDs:
        for rand_sid in sent_sys_rank[sid]:
            system_list = list(combinations(rand_sid['systems'], args.freeN))
            rank_list = list(combinations(rand_sid['ranks'], args.freeN))
            for system_tuple, rank_tuple in zip(system_list, rank_list):
                comparison_d["_".join(tuple(sorted(set(system_tuple))))].append((system_tuple, rank_tuple))


def sort_by_mu(sys_rate):
    sortlist = []
    for k, v in sys_rate.items():
        mu = v.mu
        sortlist.append((mu, k))
    sortlist.sort(reverse=True)
    return sortlist


def get_counts(s_name, c_dict, n_play):
    c_list = np.zeros((len(s_name), len(s_name)))
    total = sum(c_dict.values())
    for i, s_a in enumerate(s_name):
        for j, s_b in enumerate(s_name):
            c_list[i][j] = (c_dict[s_a + '_' + s_b] / float(sum(c_dict.values()))) *2
    return c_list.tolist()


def comb(n, k):
    """Combination number"""
    return factorial(n) / factorial(k) / factorial(n - k)


def estimate_by_number():
    #Format of rating by one judgement:
    #  [[r1], [r2], [r3], [r4], [r5]] = rate([[r1], [r2], [r3], [r4], [r5]], ranks=[1,2,3,3,5])

    for num_iter_org in num_record:
        # determine how many comparisons to use
        if num_iter_org == 0:
            ### by # of pairwise judgements
            num_rankings = 0
            for key in comparison_d.keys():
                num_rankings += len(comparison_d[key])
            data_points = num_rankings / comb(args.freeN, 2) + 1
        else:
            data_points = num_iter_org  # by # of matches
        num_iter = int(args.dp_pct * data_points)
        print >> sys.stderr, "Sampling %d / %d pairwise judgments" % (num_iter, data_points)
        # initialize TrueSkill environment
        param_beta = param_sigma * (num_iter/40.0)
        env = TrueSkill(mu=0.0, sigma=param_sigma, beta=param_beta, tau=param_tau, draw_probability=draw_rate)
        env.make_as_global()
        system_rating = {}
        counter_dict = defaultdict(int)  # keep track of # of comparisons A vs B in here
        for s in all_systems:
            system_rating[s] = Rating()

        # run TrueSkill
        for num_play in xrange(num_iter):
            systems_to_compare = scripts.next_comparison.get(system_rating, args.freeN)
            obs = random.choice(comparison_d["_".join(sorted(systems_to_compare))])  # observation -- randomly choose a sentence (systems, rank)
            systems_name_compared = obs[0]  # actually the same as systems_to_compare, but in different order
            partial_rank = obs[1]  # ranks in the same order as sytems_name_compared

            if args.freeN == 2:  # keep track of # of comparisons for heat map
                if (num_play >= (num_iter * count_begin)) and (num_play <= (num_iter * count_end)):
                    sys_a = obs[0][0]
                    sys_b = obs[0][1]
                    counter_dict[sys_a + '_' + sys_b] += 1
                    counter_dict[sys_b + '_' + sys_a] += 1

            ratings = []  # gather current TS rankings for the systems to compare
            for s in systems_name_compared:
                ratings.append([system_rating[s]])
            updated_ratings = rate(ratings, ranks=partial_rank)  # update rankings using TS
            for s, r in zip(systems_name_compared, updated_ratings):
                system_rating[s] = r[0]

        # write the resulting mus & sigmas
        f = open(args.prefix + '_mu_sigma.json', 'w')
        t = {name: (rating.mu, rating.sigma ** 2) for name, rating in system_rating.iteritems()}
        t['data_points'] = [data_points, args.dp_pct]
        json.dump(t, f)
        f.close()

        # write heatmaps
        if (args.freeN == 2) and (num_iter_org == num_record[-1]) and args.heat:
            f = open(args.prefix + '-' + str(count_begin) + '-' + str(count_end) + '_count.json', 'w')
            sys_names = zip(*sort_by_mu(system_rating))[1]
            counts = get_counts(sys_names, counter_dict, num_play)
            outf = {}
            outf['sysname'] = sys_names
            outf['counts'] = counts
            json.dump(outf, f)
            f.close()

if __name__ == '__main__':
    all_systems, sent_sys_rank = parse_csv(open(args.input_file) if args.input_file else sys.stdin)
    fill_comparisons(sent_sys_rank)
    estimate_by_number()

