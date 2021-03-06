#!/usr/bin/env python
#encoding: utf-8

__author__ = ""
__version__ = ""
__copyright__ = ""
__license__ = ""
__descripstion__ = ""
__usage__ = ""

import sys
import os
import math
import random_sample
import json

def get(ratings, n):
    """Get the next group of systems to be rated.
    Select the 1st system as the one with the currently largest sigma, then
    sample the rest according to a distribution of draw probabilities, where
    p_draw_ab = 1 / exp(abs(mu_a - mu_b)).

    @param ratings: dict of system name -> TrueSkill Rating
    @param n: number of systems to draw out for comparison
    @return: a tuple of system names to be compared
    """
    systems_compared = []
    # select the 1st system -- the one with the largest sigma
    sys_a = sorted(ratings.items(), key=lambda r: r[1].sigma, reverse=True)[0][0]
    systems_compared.append(sys_a)
    # establish p_draw against sys_a for all others
    sys_chance = [[], []]
    for k, v in ratings.items():
        if k != sys_a:
            sys_chance[0].append(k)
            sys_chance[1].append(1./math.exp(abs(ratings[sys_a].mu - ratings[k].mu)))
    # sample n-1 other systems according to p_draw
    while len(systems_compared) != n:
        sys_b = random_sample.choose(sys_chance[0], sys_chance[1])
        systems_compared.append(sys_b)
        systems_compared = list(set(systems_compared))
    # return the result
    return tuple(set(systems_compared))

if __name__ == '__main__':
    f = open(sys.argv[1], 'r')
    sample = json.load(f)
    f.close()
    print get(sample, int(sys.argv[2]))

    # unit test
    #sample = {"newstest2013.de-en.umd.2922": [0.04438640281747725, 0.24863199120996113], "newstest2013.de-en.uedin-syntax.2605": [0.2910463746313596, 0.2486317088666395], "newstest2013.de-en.MES:.2916": [0.19995113646670695, 0.2486319554372286], "newstest2013.de-en.Shef-wproa.2761": [-0.5510601821890464, 0.24863128613105073], "newstest2013.de-en.LIMSI-Ncode-SOUL-primary.2591": [-0.1654550673971147, 0.24863127737879956], "newstest2013.de-en.online-B": [0.39098637483956117, 0.24863160528166164], "newstest2013.de-en.desrt.2704": [-0.6721125315154654, 0.2486320743852226], "newstest2013.de-en.CNGL_DCU.2703": [-0.2554593630399819, 0.24863165949540617], "newstest2013.de-en.uedin-wmt13.2636": [0.015619564035093384, 0.24863126504934419], "newstest2013.de-en.online-A": [0.45771054302365144, 0.24863130455560303], "newstest2013.de-en.MES-Szeged-reorder-split-primary.2682": [0.16552429527837645, 0.24863211796423232], "newstest2013.de-en.cu-zeman.2720": [-0.19890084464280997, 0.24863075248532357], "newstest2013.de-en.RWTH-Jane-primary.2615": [0.24888930698114511, 0.24863155212763818], "newstest2013.de-en.JHU.2887": [-0.4021356113199055, 0.24863094991568652], "newstest2013.de-en.QUAERO_primary.2601": [0.19104152412885064, 0.24863187337869821], "newstest2013.de-en.KIT_primary.2653": [0.21994473742860132, 0.24863187401140185], "newstest2013.de-en.TUBITAK.2613": [0.02002381307448493, 0.24863111203905974]}
    #print get(sample, 3)

