#!/usr/bin/env python
import sys
import os
import argparse
import allele_backbone as ab
import allele_gmap as ag
import allele_blast as abl


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument('-