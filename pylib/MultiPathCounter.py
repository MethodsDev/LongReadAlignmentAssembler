#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
import logging
import MultiPath

logger = logging.getLogger(__name__)

class MultiPathCountPair:

    def __init__(self, multipath):
        self._multipath = multipath
        self.reset_count()

    def increment(self, increment=1):
        self._count += increment

    def get_multipath_and_count(self):
        return(self._multipath, self._count)

    def include_read_type(self, read_type):
        self._multipath.include_read_type(read_type)

    def include_read_name(self, read_name):
        self._multipath.include_read_name(read_name)

    def reset_count(self):
        self._count = len(self._multipath.get_read_names())
    
        
    def __repr__(self):
        ret_text = "{}\t{}".format(str(self._multipath), self._count)
        return ret_text
    
    
class MultiPathCounter:
    
    
    def __init__(self):

        self._multipath_counter = dict()
        
        return


    def add(self, multipath_obj):

        assert len(multipath_obj.get_read_names()) > 0
        
        assert type(multipath_obj) == MultiPath.MultiPath
        multipath_key = str(multipath_obj.get_simple_path())
        
        if multipath_key in self._multipath_counter:
            orig_mp_count_pair = self._multipath_counter[multipath_key]
            orig_mp_count_pair.include_read_name(multipath_obj.get_read_names())
            orig_mp_count_pair.include_read_type(multipath_obj.get_read_types())
            orig_mp_count_pair.reset_count()
            

        else:
            self._multipath_counter[multipath_key] = MultiPathCountPair(multipath_obj)


    
    def get_all_MultiPathCountPairs(self):
        return(self._multipath_counter.values())

    
    def __repr__(self):

        ret_text = "\t"

        for multipath_count_pair in self._multipath_counter.values():
            ret_text += str(multipath_count_pair) + "\n"

        return ret_text
            

    
