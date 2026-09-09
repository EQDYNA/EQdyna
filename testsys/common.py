"""Shared helpers for the testsys tiers (see PROJECT_RULES.md rule 3).

Kept deliberately tiny: printing a result is not a gate, so every tier
funnels its outcome through here into a single non-zero-on-failure exit,
per rule 3 ("the script ran" and "the script passed" are different
questions).
"""
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def banner(msg):
    print(f'\n==== testsys: {msg} ====')
