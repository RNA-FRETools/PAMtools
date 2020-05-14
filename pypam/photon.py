#!/usr/bin/env python3

import phconvert as phc

def ptu_meta(filename):
    d, meta = phc.loader.nsalex_pq(filename)
    return meta['tags']['File_Comment']['data']

