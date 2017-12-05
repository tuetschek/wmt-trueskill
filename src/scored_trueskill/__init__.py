#!/usr/bin/env python

from __future__ import unicode_literals
from trueskill import TrueSkill
from trueskill import MU, SIGMA, BETA, TAU, DRAW_PROBABILITY


# variance of the game outcome (should be possible to compute from data)
GAMMA = 0.1


class ScoreBasedTrueSkill(TrueSkill):
    """A TrueSkill environment to compute pairwise matches including score/ratings."""

    def __init__(self, mu=MU, sigma=SIGMA, beta=BETA, tau=TAU, gamma=GAMMA,
                 draw_probability=DRAW_PROBABILITY, backend=None):
        super(ScoreBasedTrueSkill, self).__init__(mu, sigma, beta, tau, draw_probability, backend)
        self.gamma = gamma

    def rate_score_based_1vs1(self, pl1, pl2, score_diff):
        """Closed-form update for score-based TrueSkill (1-vs-1 play only).

        @param pl1: player 1 rating
        @param pl2: player 2 rating
        @param score_diff: player 1 score - player 2 score (positive if player 1 won)
        @todo: global tau parameter is ignored so far (~ 0).
        """
        # Gaussian parameters: pi  = 1 / (sigma ^ 2), tau = mu / (sigma ^ 2)
        pi1_new = pl1.pi + (1 / (2 * self.beta ** 2 + self.gamma ** 2 + pl2.sigma ** 2))
        pi2_new = pl2.pi + (1 / (2 * self.beta ** 2 + self.gamma ** 2 + pl1.sigma ** 2))

        tau1_new = pl1.tau + ((pl2.mu + score_diff) / (2 * self.beta ** 2 + self.gamma ** 2 + pl2.sigma ** 2))
        tau2_new = pl2.tau + ((pl1.mu - score_diff) / (2 * self.beta ** 2 + self.gamma ** 2 + pl1.sigma ** 2))

        (pl1.pi, pl1.tau) = (pi1_new, tau1_new)
        (pl2.pi, pl2.tau) = (pi2_new, tau2_new)
