#!/usr/bin/env python
# -"- encoding: utf-8 -"-

"""
# Input: JUDGEMENTS.csv which must contain one language-pair judgements.
# Output: *_mu_sigma.json: Mu and Sigma for each system
#         *.count: number of judgements among systems (for generating a heatmap) if -n is set to 2 and -e.
"""

import sys
import argparse
import random
import json
import numpy as np
import scripts.next_comparison
from itertools import combinations
from collections import defaultdict
from csv import DictReader
from scored_trueskill import ScoreBasedTrueSkill


__author__ = "Keisuke Sakaguchi, Ondrej Dusek"
__version__ = "0.1"


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('prefix', help='output ID (e.g. fr-en0)')
arg_parser.add_argument('input_file', nargs='?', help='input file (default to stdin)', type=str)
arg_parser.add_argument('-d', '--dp', type=int,
                        help='Number of judgments to use (0 == all)', required=True)
arg_parser.add_argument('-p', '--dp_pct', type=float, default=1.0,
                        help='Percentage of judgments to use (default: 1.0)')
arg_parser.add_argument('-s', '--num_systems', type=int, default=5,
                        help='Number of systems in one ranking in CSV file (default=5)')
args = arg_parser.parse_args()

#######################################
# Global Variables and Parameters     #
param_sigma = 0.5
param_tau = 0.
draw_rate = 0.25

# You can set arbitrary number(s) for record (dp is the number assigned by -d).
# num_record = [int(args.dp*0.9), args.dp]
num_record = [args.dp]
# e.g. num_record = [args.dp*0.125, args.dp*0.25, args.dp*0.5, args.dp*0.9, args.dp]
# e.g. num_record = [400, 800, 1600, 3200, 6400, 11520, 12800]

# When -n is set to 2, you can set beginning and ending between (0 and 1) for counting the number of comparisons among systems.
# This is used for generating a heatmap.
# e.g. "count_begin=0.4 and count_end=0.6" records the number of comparisons from 40% to 60% of total comparisons.
count_begin = 0.8
count_end = 1.0
if count_begin > count_end:
    raise
#######################################


def parse_csv(fh=sys.stdin):
    """Parse the WMT-formatted CSV file and return system names and rank(1-5)
    for each sentence."""
    all_systems = []
    sent_sys_rank = defaultdict(list)
    for i, row in enumerate(DictReader(fh)):
        sentID = int(row.get('segmentId'))
        systems = []
        ranks = []
        for num in range(1, args.num_systems + 1):
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
    """Convert all multi-way ranking data to pairwise comparisons (using all subset
    combinations for all sentences). Store this in the global comparison_d variable."""
    comparison_d = defaultdict(list)
    for sid in sent_sys_rank.keys():
        for rand_sid in sent_sys_rank[sid]:
            system_list = list(combinations(rand_sid['systems'], 2))
            rank_list = list(combinations(rand_sid['ranks'], 2))
            for system_tuple, rank_tuple in zip(system_list, rank_list):
                comparison_d["_".join(tuple(sorted(set(system_tuple))))].append((system_tuple, rank_tuple))
    return comparison_d


def compute_gamma(sent_sys_rank):
    """Compute overall score variance."""
    all_ranks = []
    for sent_ratings in sent_sys_rank.itervalues():
        for rating in sent_ratings:
            all_ranks.extend(rating['ranks'])
    return np.var(all_ranks)


def estimate_by_number(all_systems, comparison_d, variance):

    for num_iter_org in num_record:
        # determine how many comparisons to use
        if num_iter_org == 0:
            num_rankings = sum(len(val) for val in comparison_d.itervalues())
            data_points = num_rankings + 1
        else:
            data_points = num_iter_org  # by # of matches
        num_iter = int(args.dp_pct * data_points)
        print >> sys.stderr, "Sampling %d / %d pairwise judgments" % (num_iter, data_points)
        # initialize TrueSkill environment
        param_beta = param_sigma * (num_iter / 40.0)
        sts = ScoreBasedTrueSkill(mu=0.0, sigma=param_sigma, beta=param_beta, tau=param_tau, gamma=gamma, draw_probability=draw_rate)
        system_rating = {}
        for s in all_systems:
            system_rating[s] = sts.create_rating()

        # run TrueSkill
        for num_play in xrange(num_iter):
            systems_to_compare = scripts.next_comparison.get(system_rating, 2)
            obs = random.choice(comparison_d["_".join(sorted(systems_to_compare))])  # observation -- randomly choose a sentence (systems, rank)
            sys_a, sys_b = obs[0]  # actually the same as systems_to_compare, but in different order
            score_a, score_b = obs[1]  # ranks in the same order as sytems_name_compared

            # update rankings using score-based TS
            sts.rate_score_based_1vs1(system_rating[sys_a], system_rating[sys_b], score_a - score_b)

        # write the resulting mus & sigmas
        f = open(args.prefix + '_mu_sigma.json', 'w')
        t = {name: (rating.mu, rating.sigma ** 2) for name, rating in system_rating.iteritems()}
        t['data_points'] = [data_points, args.dp_pct]
        t['gamma'] = variance
        json.dump(t, f)
        f.close()


if __name__ == '__main__':
    all_systems, sent_sys_rank = parse_csv(open(args.input_file) if args.input_file else sys.stdin)
    comparison_d = fill_comparisons(sent_sys_rank)
    gamma = compute_gamma(sent_sys_rank)
    estimate_by_number(all_systems, comparison_d, gamma)
