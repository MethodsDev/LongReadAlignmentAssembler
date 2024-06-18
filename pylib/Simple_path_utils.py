#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER
from Splice_graph import Splice_graph
from GenomeFeature import *
from Util_funcs import coordpairs_overlap
import logging

logger = logging.getLogger(__name__)

# Namespace: Simple_path_utils
# includes basic functions for evaluating relationships between simple paths in the graph.



## Scenarios:
#
## Scenario 1
## path_A  ==============
## path_B        ===============  (idx_B == 0)
## or
#
## Scenario 2
## path_A        ================= (idx_A == 0)
## path_B  ==============
## or
#
## Scenarion 3
## path_A  ======================= (idx_B == 0)
## path_B         =======
## or
#
#  Scenario 4
## path_A         =======
## path_B  ======================= (idx_A == 0)
## or
#
# Scenario 5
## path_A       ==========     or   ============   or =====        (either are idx 0)
## path_B       ==========          ====              ===========



def path_A_contains_path_B(simple_path_A, simple_path_B):

    if len(simple_path_B) > len(simple_path_A):
        return False

    if simple_path_B[0] not in simple_path_A:
        return False
    
    idx_A = simple_path_A.index(simple_path_B[0])
    if idx_A < 0:
        return False

    if idx_A + len(simple_path_B) > len(simple_path_A):
        return False

    idx_B = 0
    while idx_B < len(simple_path_B):
        idx_B += 1
        idx_A += 1
        if (idx_B < len(simple_path_B) and
            idx_A < len(simple_path_A) and
            simple_path_A[idx_A] != simple_path_B[idx_B]):
            return False

    return True


def are_overlapping_and_compatible_NO_gaps_in_overlap(simple_path_A, simple_path_B):
    
    ## find first non-spacer match between two paths.  Ensure remaining parts of paths are identical
    for idx_A in range(0, len(simple_path_A)):
        A_node_id = simple_path_A[idx_A]
        if A_node_id == SPACER:
            continue
        if A_node_id in simple_path_B:
            idx_B = simple_path_B.index(A_node_id)

            ## one of the indexes needs to start at zero or there'll be some unmatched upstream nodes.
            if (idx_B != 0 and idx_A != 0):
                return False

            # ensure remainder of paths overlap, no gaps allowed.
            idx_A += 1
            idx_B += 1
            while (idx_A < len(simple_path_A) and idx_B < len(simple_path_B)):
                if simple_path_A[idx_A] == SPACER or simple_path_A[idx_A] != simple_path_B[idx_B]:
                    return False
                idx_A += 1
                idx_B += 1
                
            # if made it here, the remainder of the paths are gap-free and identical
            return True

    return False # no matching index between paths



def merge_simple_paths(simple_path_A, simple_path_B):
    
    if not are_overlapping_and_compatible_NO_gaps_in_overlap(simple_path_A, simple_path_B):
        raise RuntimeException("cannot merge paths that are not compatible in overlapping region")


    ## find first non-spacer match between two paths, then merge.
    for idx_A in range(0, len(simple_path_A)):
        A_node_id = simple_path_A[idx_A]
        if A_node_id == SPACER:
            continue
        if A_node_id in simple_path_B:
            idx_B = simple_path_B.index(A_node_id)

            merged_path = None
            if idx_A == 0: # scenarios 2,4,5
                if idx_B > 0: # scenario 2 or 4
                    # begin path with B prefix
                    merged_path = simple_path_B
                    # if A extends past B, need to include that.
                    #  path A:        0 1 2 3 4 5 6
                    #  path B:    0 1 2 3 4 5 6
                    extension_idx_A = len(simple_path_B) - idx_B
                    if extension_idx_A < len(simple_path_A):
                        merged_path.extend(simple_path_A[extension_idx_A: ])
                    return merged_path
                else: #scenario 5, return longer path
                    if len(simple_path_A) >= len(simple_path_B):
                        return simple_path_A
                    else:
                        return simple_path_B
            
            else:
                # idx_A != 0, so idx_B must be == 0
                assert(idx_B == 0)
                # scenarios 1,3
                # begin path with A prefix
                merged_path = simple_path_A
                # if A extends past B, need to include that.
                #  path A:   0 1 2 3 4 5 6
                #  path B:         0 1 2 3 4 5 6
                extension_idx_B = len(simple_path_A) - idx_A
                if extension_idx_B < len(simple_path_B):
                    merged_path.extend(simple_path_B[extension_idx_B: ])
                return merged_path
            
    raise RuntimeException("Error, could not merge simple paths {} and {} ... bug... ".format(simple_path_A, simple_path_B))



def count_exons_in_simple_path(simple_path):
    num_exons = 0

    for node_id in simple_path:
        if re.match("E:", node_id):
            num_exons += 1

    return num_exons


def merge_adjacent_segments(segment_list):

    assert len(segment_list) > 0

    segment_list = sorted(segment_list, key=lambda x: x[0])
    
    ret_segments = list()
    ret_segments.append(list(segment_list.pop(0)))

    while len(segment_list) > 0:
        next_seg = list(segment_list.pop(0))

        prev_seg = ret_segments[-1]
        if next_seg[0] == prev_seg[1] + 1:
            # merge them
            assert(next_seg[1] > prev_seg[1])
            prev_seg[1] = next_seg[1]
        elif prev_seg[1] == next_seg[1] and prev_seg[0] <= next_seg[0]:
            # identical or contained feature - skip it.
            pass
        elif next_seg[0] > prev_seg[1]:
            ret_segments.append(next_seg)
        else:
            raise RuntimeError("Error, not sure hwo to merge adjacent segments: {} and {}".format(prev_seg, next_seg))

    
    return ret_segments




def _convert_path_to_nodes_with_coords_list(sg:Splice_graph, simple_path:list, ret_node_objs=False) -> list:

    simple_path = trim_terminal_spacers(simple_path.copy())
    
    node_ids_with_coords_list = list()
    found_spacer = False
    for i, node_id in enumerate(simple_path):
        if node_id != SPACER:
            node_obj = sg.get_node_obj_via_id(node_id)
            lend, rend = node_obj.get_coords()
            if ret_node_objs:
                node_id_with_coords = [node_obj, lend, rend]
            else:
                node_id_with_coords = [node_id, lend, rend]
                
            node_ids_with_coords_list.append(node_id_with_coords)
        else:
            found_spacer = True
            node_ids_with_coords_list.append([SPACER, -1, -1])


    if found_spacer:
        # adjust spacer coordinates to neighboring bounds of nodes.
        #print(node_ids_with_coords_list)
        for i, node_info_list in enumerate(node_ids_with_coords_list):
            if node_info_list[0] == SPACER:
                node_info_list[1] = node_ids_with_coords_list[i-1][2] + 1
                node_info_list[2] = node_ids_with_coords_list[i+1][1] -1

    return node_ids_with_coords_list



def fraction_read_overlap(sg:Splice_graph, read_simple_path:list, transcript_simple_path:list) -> float:

    nodes_in_common = list()
    for node in read_simple_path:
        if node in transcript_simple_path:
            nodes_in_common.append(node)

    read_simple_path_exons = simple_path_to_exon_objs(sg, read_simple_path)
    read_align_len = 0
    for exon in read_simple_path_exons:
        lend, rend = exon.get_coords()
        exon_len = rend - lend + 1
        read_align_len += exon_len
    
    
    exons_in_common = simple_path_to_exon_objs(sg, nodes_in_common)
    common_region_len = 0
    for exon in exons_in_common:
        lend, rend = exon.get_coords()
        exon_len = rend - lend + 1
        common_region_len += exon_len

    fraction_read_aligned_to_transcript = common_region_len / read_align_len

    return fraction_read_aligned_to_transcript




def simple_path_to_exon_objs(sg:Splice_graph, simple_path:list) -> list:

    exons = list()
    for node_id in simple_path:
        if node_id != SPACER:
            obj = sg.get_node_obj_via_id(node_id)
            if type(obj) == Exon:
                exons.append(obj)
                
    return exons

def simple_path_to_intron_objs(sg:Splice_graph, simple_path:list) -> list:

    introns = list()
    for node_id in simple_path:
        if node_id != SPACER:
            obj = sg.get_node_obj_via_id(node_id)
            if type(obj) == Intron:
                introns.append(obj)
                
    return introns




def _split_spacers_with_coords(simple_path_w_coords:list) -> list:

    adj_list = list()

    for node_coordset in simple_path_w_coords:
        node_id, lend, rend = node_coordset
        if node_id == SPACER:
            adj_list.append([SPACER, lend, lend])
            adj_list.append([SPACER, rend, rend])
        else:
            adj_list.append(node_coordset)

    return adj_list




def merge_simple_paths_containing_spacers(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> list:

    """
    Remove redundancies and adjust SPACER coordinates
    """

    
    A_list = _convert_path_to_nodes_with_coords_list(sg, simple_path_A)
    A_list = _split_spacers_with_coords(A_list)
    
    B_list = _convert_path_to_nodes_with_coords_list(sg, simple_path_B)
    B_list = _split_spacers_with_coords(B_list)

    
    merged_list = A_list + B_list
    merged_list = sorted(merged_list, key=lambda x: (x[1], x[2]) )

    adj_merged_list = list()
    adj_merged_list.append(merged_list.pop(0))
    
    for entry in merged_list:
        prev_entry = adj_merged_list[-1]
        prev_node_id, prev_lend, prev_rend = prev_entry

        curr_node_id, curr_lend, curr_rend = entry

        if prev_node_id == curr_node_id:
            if curr_node_id == SPACER:
                if curr_rend > prev_rend:
                    prev_entry[2] = curr_rend
            else:
                assert(prev_lend == curr_lend and prev_rend == curr_rend)

        else:
            adj_merged_list.append(entry)


    # need simple path to return
    simple_path_ret = list()
    for entry in adj_merged_list:
        simple_path_ret.append(entry[0])

            
    return simple_path_ret



def simple_paths_overlap(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool:

    # only checks bounding coordinates for overlaps
    
    simple_path_A = trim_terminal_spacers(simple_path_A.copy())
    simple_path_B = trim_terminal_spacers(simple_path_B.copy())

    path_A_lend = sg.get_node_obj_via_id(simple_path_A[0]).get_coords()[0]
    path_A_rend = sg.get_node_obj_via_id(simple_path_A[-1]).get_coords()[1]

    path_B_lend = sg.get_node_obj_via_id(simple_path_B[0]).get_coords()[0]
    path_B_rend = sg.get_node_obj_via_id(simple_path_B[-1]).get_coords()[1]

    return path_A_lend < path_B_rend and path_A_rend > path_B_lend


def trim_terminal_spacers(simple_path:list) -> list:


    # trim any terminal spacer
    while len(simple_path) > 0  and simple_path[-1] == SPACER:
        simple_path.pop()
        
    # trim any initial spacer
    while len(simple_path) > 0 and simple_path[0] == SPACER:
        simple_path.pop(0)
    
    return simple_path



def simple_path_A_within_bounds_of_simple_path_B(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool:

    # only checks bounding coordinates for overlaps
    
    path_A_lend = sg.get_node_obj_via_id(simple_path_A[0]).get_coords()[0]
    path_A_rend = sg.get_node_obj_via_id(simple_path_A[-1]).get_coords()[1]

    path_B_lend = sg.get_node_obj_via_id(simple_path_B[0]).get_coords()[0]
    path_B_rend = sg.get_node_obj_via_id(simple_path_B[-1]).get_coords()[1]

    return path_A_lend >= path_B_lend and path_A_rend <= path_B_rend
    


def simple_path_A_contains_and_compatible_with_simple_path_B_spacefree_region_path_A(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool: 

    if simple_path_A == simple_path_B:
        return True
    elif path_A_contains_path_B(remove_spacers_from_path(simple_path_A), remove_spacers_from_path(simple_path_B)):
        return True
    else:
        return simple_path_A_within_bounds_of_simple_path_B(sg, simple_path_B, simple_path_A) and simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, simple_path_A, simple_path_B)



def simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool: 

    if simple_path_A == simple_path_B:
        return True
    elif path_A_contains_path_B(remove_spacers_from_path(simple_path_A), remove_spacers_from_path(simple_path_B)):
        return True
    else:
        return simple_path_A_within_bounds_of_simple_path_B(sg, simple_path_B, simple_path_A) and simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, simple_path_A, simple_path_B)


    

def simple_paths_overlap_and_compatible_spacefree_region_path_A(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool:

    ## perform compatibility analysis, require no spaces in path A region of overlap, but ok if in path B

    if not simple_paths_overlap(sg, simple_path_A, simple_path_B):
        return False

    # begin in spacer mode because they are unlikely to be aligned.
    my_path_spacer_mode = True
    other_path_spacer_mode = True

    my_path = simple_path_A.copy()
    other_path = simple_path_B.copy()

    alignment_started = False
    
    while(len(my_path) > 0 and len(other_path) > 0):
        my_node_id = my_path[0]
        if my_node_id == SPACER:
            my_path.pop(0)
            my_path_spacer_mode = True
            if alignment_started:
                return False
            continue

        other_node_id = other_path[0]
        if other_node_id == SPACER:
            other_path.pop(0)
            other_path_spacer_mode = True
            continue

        if my_node_id == other_node_id:
            # great! advance and continue
            my_path.pop(0)
            other_path.pop(0)
            my_path_spacer_mode = False
            other_path_spacer_mode = False
            alignment_started = True
            continue

        my_node_obj = sg.get_node_obj_via_id(my_node_id)
        other_node_obj = sg.get_node_obj_via_id(other_node_id)

        if coordpairs_overlap(my_node_obj.get_coords(), other_node_obj.get_coords()):
            # uh oh, they overlap but they're not the same.
            return False

        else:
            if not (my_path_spacer_mode or other_path_spacer_mode):
                return False

            # advance the side that's before the other
            if other_path_spacer_mode and my_node_obj.get_coords()[0] < other_node_obj.get_coords()[0]:
                my_path.pop(0)
            elif my_path_spacer_mode and my_node_obj.get_coords()[0] > other_node_obj.get_coords()[0]:
                other_path.pop(0)
            else:
                return False

    if not alignment_started:
        return False
            
    return True # no conflicts detected
    


def simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg:Splice_graph, simple_path_A:list, simple_path_B:list) -> bool:

    ## perform compatibility analysis allowing spacers to appear in both paths and within overlapping regions of paths
    
    if not simple_paths_overlap(sg, simple_path_A, simple_path_B):
        return False


    # begin in spacer mode because they are unlikely to be aligned.
    my_path_spacer_mode = True
    other_path_spacer_mode = True

    my_path = simple_path_A.copy()
    other_path = simple_path_B.copy()

    alignment_started = False
    
    while(len(my_path) > 0 and len(other_path) > 0):
        my_node_id = my_path[0]
        if my_node_id == SPACER:
            my_path.pop(0)
            my_path_spacer_mode = True
            continue

        other_node_id = other_path[0]
        if other_node_id == SPACER:
            other_path.pop(0)
            other_path_spacer_mode = True
            continue

        if my_node_id == other_node_id:
            # great! advance and continue
            my_path.pop(0)
            other_path.pop(0)
            my_path_spacer_mode = False
            other_path_spacer_mode = False
            alignment_started = True
            continue

        my_node_obj = sg.get_node_obj_via_id(my_node_id)
        other_node_obj = sg.get_node_obj_via_id(other_node_id)

        if coordpairs_overlap(my_node_obj.get_coords(), other_node_obj.get_coords()):
            # uh oh, they overlap but they're not the same.
            return False

        else:
            if not (my_path_spacer_mode or other_path_spacer_mode):
                return False

            # advance the side that's before the other
            if other_path_spacer_mode and my_node_obj.get_coords()[0] < other_node_obj.get_coords()[0]:
                my_path.pop(0)
            elif my_path_spacer_mode and my_node_obj.get_coords()[0] > other_node_obj.get_coords()[0]:
                other_path.pop(0)
            else:
                return False

    if not alignment_started:
        return False
            
    return True # no conflicts detected


def try_fill_spacers_via_splicegraph(sg, simple_path):



    # first, try to insert missing introns
    new_path = try_fill_spacers_via_splicegraph_insert_missing_introns(sg, simple_path)

    # remove spacers between connected nodes.
    new_path = remove_spacers_between_connected_nodes(sg, new_path)

    return new_path



def try_fill_spacers_via_splicegraph_insert_missing_introns(sg, simple_path):

    new_path = simple_path.copy()

    for i in range(len(simple_path)):
        if simple_path[i] == SPACER:
            if i == 0 or i == len(simple_path)-1:
                # cant be terminal
                continue

            prev_node = simple_path[i-1]
            next_node = simple_path[i+1]

            if prev_node == SPACER or next_node == SPACER:
                continue

            prev_node_obj = sg.get_node_obj_via_id(prev_node)
            next_node_obj = sg.get_node_obj_via_id(next_node)

            # require exons on both sides of spacer
            if not ( type(prev_node_obj) == Exon and type(next_node_obj) == Exon):
                continue

            intron_lend = prev_node_obj._rend + 1
            intron_rend = next_node_obj._lend - 1

            intron_obj = sg.get_intron_node_obj(intron_lend, intron_rend)
            if intron_obj:
                intron_id = intron_obj.get_id()
                new_path[i] = intron_id # replace spacer


    return new_path
                


def remove_spacers_between_connected_nodes(sg, simple_path):

    simple_path = simple_path.copy()
    
    deleted_spacer = "DELETED_SPACER"
    
    for i in range(1, len(simple_path)-1):
        if simple_path[i] == SPACER:
            prev_node_id = simple_path[i-1]
            next_node_id = simple_path[i+1]
            prev_node_obj = sg.get_node_obj_via_id(prev_node_id)
            next_node_obj = sg.get_node_obj_via_id(next_node_id)

            if next_node_obj in sg.get_successors(prev_node_obj):
                # connected. Remove SPACER
                simple_path[i] = deleted_spacer
                logger.debug("-deleting spacer between connected nodes {} and {} in path {}".format(prev_node_id, next_node_id, simple_path))

    new_path = simple_path
    if deleted_spacer in simple_path:
        new_path = [x for x in simple_path if x != deleted_spacer]

    return new_path



def split_path_at_spacers(simple_path):

    if SPACER not in simple_path:
        return [simple_path]

    else:
        simple_path = trim_terminal_spacers(simple_path)
        paths = [ [] ]
        
        for node_id in simple_path:
            if node_id == SPACER:
                paths.append([])
            else:
                paths[-1].append(node_id)


    return paths

def remove_spacers_from_path(simple_path):

    new_path = list()
    for node_id in simple_path:
        if node_id != SPACER:
            new_path.append(node_id)

    return new_path



def refine_TSS_simple_path(splice_graph, simple_path):

    #print("Simple path: {}".format(simple_path))
    contig_strand = splice_graph.get_contig_strand()

    nodes_with_coords_list = _convert_path_to_nodes_with_coords_list(splice_graph, simple_path)

    nodes_with_coords_list = sorted(nodes_with_coords_list, key=lambda x: (x[1], x[2]))
    
    
    TSS_indices = list()
    for i, node_w_coords in enumerate(nodes_with_coords_list):
        if re.match("TSS:", node_w_coords[0]):
            TSS_indices.append(i)

    if len(TSS_indices) == 0:
        return(simple_path)

    logger.debug("Found TSS in path {} at indices {}".format(nodes_with_coords_list, TSS_indices))

    # trim off region beyond candidate TSS if within allowed distance
    if contig_strand == '+' and TSS_indices[0] != 0:
        TSS_index = TSS_indices[0]
        lend_coord =  nodes_with_coords_list[0][1]
        rend_coord =  nodes_with_coords_list[0][2]
        TSS_coord = nodes_with_coords_list[TSS_index][1]

        # special case where first exon is a single base segment and aligns to TSS
        if lend_coord == rend_coord and lend_coord == TSS_coord and TSS_index == 1:
            # swap order
            (nodes_with_coords_list[0], nodes_with_coords_list[1]) = (nodes_with_coords_list[1], nodes_with_coords_list[0])
        
        elif TSS_coord - lend_coord <= LRAA_Globals.config['max_dist_between_alt_TSS_sites']:
            nodes_with_coords_list = nodes_with_coords_list[TSS_index:]

    elif contig_strand == '-' and TSS_indices[-1] != len(nodes_with_coords_list)-1:
        TSS_index = TSS_indices[-1]
        lend_coord = nodes_with_coords_list[-1][1]
        rend_coord = nodes_with_coords_list[-1][2]
        TSS_coord = nodes_with_coords_list[TSS_index][2]

        # special case where last exon is a single base segment and overlaps TSS
        if lend_coord == rend_coord and rend_coord == TSS_coord and TSS_index == len(nodes_with_coords_list)-2:
            # swap
            (nodes_with_coords_list[-2], nodes_with_coords_list[-1]) = (nodes_with_coords_list[-1], nodes_with_coords_list[-2])
        
        
        elif rend_coord - TSS_coord <= LRAA_Globals.config['max_dist_between_alt_TSS_sites']:
            nodes_with_coords_list = nodes_with_coords_list[0:TSS_index+1]
        
            
    # remove intervening TSS annotations
    if contig_strand == '+':
        idx_low = 1
        idx_high = len(nodes_with_coords_list) -1
    else: # (-)strand
        idx_low = 0
        idx_high = len(nodes_with_coords_list) -2

    for i in range(idx_high, idx_low-1, -1):
        if re.match("TSS:", nodes_with_coords_list[i][0]):
            nodes_with_coords_list.pop(i)

    # regenerate the simple path
    simple_path = [ x[0] for x in nodes_with_coords_list]

    logger.debug("TRIMMED to TSS {} {}".format(contig_strand, simple_path))
    
    return simple_path



def refine_PolyA_simple_path(splice_graph, simple_path):

    #print("\nSimple path input: {}".format(simple_path))
    contig_strand = splice_graph.get_contig_strand()

    nodes_with_coords_list = _convert_path_to_nodes_with_coords_list(splice_graph, simple_path)

    nodes_with_coords_list = sorted(nodes_with_coords_list, key=lambda x: (x[2], x[1]))
        
    polyA_indices = list()
    for i, node_w_coords in enumerate(nodes_with_coords_list):
        if re.match("POLYA:", node_w_coords[0]):
            polyA_indices.append(i)

    if len(polyA_indices) == 0:
        return(simple_path)

    logger.debug("Found POLYA in path {} at indices {}".format(nodes_with_coords_list, polyA_indices))

    # trim off region beyond candidate polyA if within allowed distance

    if contig_strand == '+' and polyA_indices[-1] != len(nodes_with_coords_list)-1:
        polyA_index = polyA_indices[-1]
        lend_coord = nodes_with_coords_list[-1][1]
        rend_coord = nodes_with_coords_list[-1][2]
        polyA_coord = nodes_with_coords_list[polyA_index][2]

        # special case where last exon is a single base segment and overlaps polyA
        if lend_coord == rend_coord and rend_coord == polyA_coord and polyA_index == len(nodes_with_coords_list)-2:
            # swap
            (nodes_with_coords_list[-2], nodes_with_coords_list[-1]) = (nodes_with_coords_list[-1], nodes_with_coords_list[-2])


        elif rend_coord - polyA_coord <= LRAA_Globals.config['max_dist_between_alt_polyA_sites']:
            nodes_with_coords_list = nodes_with_coords_list[0:polyA_index+1]

            
    elif contig_strand == '-' and polyA_indices[0] != 0:
        polyA_index = polyA_indices[0]
        lend_coord =  nodes_with_coords_list[0][1]
        rend_coord = nodes_with_coords_list[0][2]
        polyA_coord = nodes_with_coords_list[polyA_index][1]

        # special case where first exon is a single base segment and aligns to polyA
        if lend_coord == rend_coord and lend_coord == polyA_coord and polyA_index == 1:
            # swap order
            (nodes_with_coords_list[0], nodes_with_coords_list[1]) = (nodes_with_coords_list[1], nodes_with_coords_list[0])
        
        elif polyA_coord - lend_coord <= LRAA_Globals.config['max_dist_between_alt_polyA_sites']:
            nodes_with_coords_list = nodes_with_coords_list[polyA_index:]

            
    # remove intervening polyAsite annotations
    if contig_strand == '+':
        idx_low = 0
        idx_high = len(nodes_with_coords_list) -2
    else: # (-)strand
        idx_low = 1
        idx_high = len(nodes_with_coords_list) -1
        

    for i in range(idx_high, idx_low-1, -1):
        if re.match("POLYA:", nodes_with_coords_list[i][0]):
            nodes_with_coords_list.pop(i)

    #print(nodes_with_coords_list)
            
    # regenerate the simple path
    simple_path = [ x[0] for x in nodes_with_coords_list]

    logger.debug("TRIMMED to POLYA {} {}".format(contig_strand, simple_path))

    #print("\nSimple path output: {}".format(simple_path))
    
    return simple_path

    

def trim_TSS_and_PolyA(simple_path, strand):

    assert strand in ('+','-')

    TSS_id = None
    polyA_id = None
    
    if strand == '+':
        if re.match("TSS:", simple_path[0]):
            TSS_id = simple_path[0]
            simple_path = simple_path[1:]
        if re.match("POLYA:", simple_path[-1]):
            polyA_id = simple_path[-1]
            simple_path = simple_path[0:-1]

    else:
        # - strand
        if re.match("POLYA:", simple_path[0]):
            polyA_id = simple_path[0]
            simple_path = simple_path[1:]
        if re.match("TSS:", simple_path[-1]):
            TSS_id = simple_path[-1]
            simple_path = simple_path[0:-1]
        
    return (simple_path, TSS_id, polyA_id)



def add_spacers_between_disconnected_nodes(splice_graph, simple_path):

    if len(simple_path) < 2:
        return simple_path # nothing to do here.
    
    new_path = [ simple_path[0] ]

    spacer_added = False
    
    for i in range(1, len(simple_path)):
        prev_node_id = new_path[-1]
        node_id = simple_path[i]

        if prev_node_id != SPACER and node_id != SPACER:
            # verify connection
            prev_node = splice_graph.get_node_obj_via_id(prev_node_id)
            node = splice_graph.get_node_obj_via_id(node_id)
            if node in splice_graph.get_successors(prev_node):
                # connection verified. No spacer needed
                new_path.append(node_id)
            else:
                # must add spacer between them:
                new_path.extend([SPACER, node_id])
                spacer_added = True
                
        else:
            # one must be a spacer already, so no connection to verify here.
            new_path.append(node_id)

    if spacer_added:
        logger.debug("-spacer added between unconnected nodes. Prev path:\n{}\nnew path\n{}".format(simple_path, new_path))
    
    return new_path

    
###############
# unit tests ##
###############


def test_are_overlapping_and_compatible_NO_gaps_in_overlap():
    
    # Tests that should return True

    path_a = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    path_b =             ["n2", "n3", "n4", "n5", "n6", "n7", "n8"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    path_b =             ["n2", "n3", "n4", "n5"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)

    
    path_a =              ["n2", "n3", "n4", "n5"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)


    path_a =             ["n2", "n3", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)

    

    ################################
    # Tests that should return False

    path_a =             ["n2", "n10", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)
    
        
    path_a =             ["n2", "X10", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)

    
    path_a = ["n2", "n10"]
    path_b =               ["n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)
    


def test_merge_simple_paths():
    
    ###################
    ## Test merging paths

    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2", "n3", "n4", "n5"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3", "n4", "n5"])



    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2", "n3"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3"])


    path_a =              ["n3", "n4", "n5"]
    path_b = ["n1", "n2", "n3"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3", "n4", "n5"])


def test_path_A_contains_path_B():
    
    ####################
    ## Test containments
    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n1", "n2", "n3"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n1", "n2", "n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)

    path_a = ["n1", "n2", "n3"]
    path_b = ["n0", "n1", "n2", "n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)



def  __get_dummy_splice_graph():
    
    sg = Splice_graph()

    Exon.reset_counter()
    
    
    #   E1:100-200   E2:300-400      E3:500-600          E4:700-800    E5:900-1000
    #    [-----]     [--------]      [---------]         [--------]    [---------]
    #             
    
    e1 = Exon("contig", 100, 200, '+', 1)
    e1_ID = e1.get_id()
    sg._node_id_to_node[ e1_ID ] = e1
    
    
    e2 = Exon("contig", 300, 400, '+', 1)
    e2_ID = e2.get_id()
    sg._node_id_to_node[ e2_ID ] = e2

    e3 = Exon("contig", 500, 600, '+', 1)
    e3_ID = e3.get_id()
    sg._node_id_to_node[ e3_ID ] = e3

    e4 = Exon("contig", 700, 800, '+',  1)
    e4_ID = e4.get_id()
    sg._node_id_to_node[ e4_ID ] = e4

    e5 = Exon("contig", 900, 1000, '+', 1)
    e5_ID = e5.get_id()
    sg._node_id_to_node[ e5_ID ] = e5

    #print(str(sg._node_id_to_node))
        
    return sg

def test_overlapping_n_compatible_spacers_aware_both_paths():

    sg = __get_dummy_splice_graph()

    p1 = ['E:1', 'E:2', 'E:3']
    p2 = ['E:2', 'E:3', 'E:4'] 
    
    # test compatible paths - no spacers
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, p1, p2)  == True)
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, p2, p1) == True)

    # test incompatible paths - no spacers
    p3 = ['E:1', 'E:2', 'E:4']
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, p1, p3) == False)
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, p3, p1) == False)
    
    # test compatible paths with spacers
    sp1 = ['E:1', SPACER, 'E:4'] 
    sp2 = ['E:1', SPACER, 'E:3', 'E:4'] 
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, sp1, sp2) == True)

    sp3 = ['E:1', SPACER, 'E:3']
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, sp1, sp3) == True)

    # test incompatible paths with spacers
    p4 = ['E:2', 'E:3', 'E:5']
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, sp2, p4) == False)
    
    # test multiple spacers
    sp4 = ['E:1', SPACER, 'E:2', 'E:3', 'E:4', 'E:5']
    sp5 = ['E:1', 'E:2', 'E:3', 'E:4', SPACER, 'E:5']
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, sp4, sp5) == True)
    
    # test incompatible multipel spacers
    p5 = ['E:1', 'E:2', 'E:5']
    assert(simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, p5, sp4) == False)



def test_simple_paths_overlap_and_compatible_spacefree_region_path_A():
    
    sg = __get_dummy_splice_graph()

    p1 = ['E:1', 'E:2', 'E:3']
    p2 = ['E:2', 'E:3', 'E:4'] 
    
    # test compatible paths - no spacers
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p1, p2)  == True)
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p2, p1) == True)

    # test incompatible paths - no spacers
    p3 = ['E:1', 'E:2', 'E:4']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p1, p3) == False)
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p3, p1) == False)
    
    # test compatible paths with spacers
    sp1 = ['E:1', SPACER, 'E:4'] 
    sp2 = ['E:1', SPACER, 'E:3', 'E:4'] 
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp1, sp2) == False) ## 

    sp3 = ['E:1', SPACER, 'E:3']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp1, sp3) == False) ##

    # test incompatible paths with spacers
    p4 = ['E:2', 'E:3', 'E:5']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp2, p4) == False)
    
    # test multiple spacers
    sp4 = ['E:1', SPACER, 'E:2', 'E:3', 'E:4', 'E:5']
    sp5 = ['E:1', 'E:2', 'E:3', 'E:4', SPACER, 'E:5']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp4, sp5) == False) ##
    
    # test incompatible multiple spacers
    p5 = ['E:1', 'E:2', 'E:5']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p5, sp4) == False)


    # test additional specialized cases
    
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp5, p1) == True)

    sp6 = ['E:2', 'E:3', 'E:4', SPACER, 'E:5']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp6, p1) == True)

    p6 = ['E:3', 'E:4', 'E:5']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, p6, sp6) == True)
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp6, p6) == False)

    p7 = ['E:2']
    assert(simple_paths_overlap_and_compatible_spacefree_region_path_A(sg, sp3, p7) == False)
    
    
def test_merge_simple_paths_containing_spacers():

    sg = Splice_graph()
    Exon.reset_counter()
    
    e1 = Exon("contig", 100, 200, '+', 1)
    e1_ID = e1.get_id()
    sg._node_id_to_node[ e1_ID ] = e1

    e2 = Exon("contig", 300, 400, '+', 1)
    e2_ID = e2.get_id()
    sg._node_id_to_node[ e2_ID ] = e2

    e3 = Exon("contig", 500, 600, '+', 1)
    e3_ID = e3.get_id()
    sg._node_id_to_node[ e3_ID ] = e3

    e4 = Exon("contig", 700, 800, '+', 1)
    e4_ID = e4.get_id()
    sg._node_id_to_node[ e4_ID ] = e4

    e5 = Exon("contig", 900, 1000, '+', 1)
    e5_ID = e5.get_id()
    sg._node_id_to_node[ e5_ID ] = e5
    

    sp1 = [e1_ID, e2_ID, SPACER, e3_ID]
    sp2 = [e1_ID, SPACER, e2_ID, SPACER, e3_ID, SPACER, e5_ID]

    conv_sp1 = _convert_path_to_nodes_with_coords_list(sg, sp1)
    conv_sp2 = _convert_path_to_nodes_with_coords_list(sg, sp2)

    #print(str(conv_sp1))
    #print(str(conv_sp2))

    # mewrge them:
    merged = merge_simple_paths_containing_spacers(sg, sp1, sp2)
    print(str(merged))

    #assert (merged == [['E:1', 100, 200], ['???', 201, 299], ['E:2', 300, 400], ['???', 401, 499], ['E:3', 500, 600], ['???', 601, 899], ['E:5', 900, 1000]] )
    assert (merged == ['E:1', '???', 'E:2', '???', 'E:3', '???', 'E:5'] )



    
    

def test_trim_terminal_spacers():

    simple_path_A = [SPACER, "a", "b", "c"]

    assert(trim_terminal_spacers(simple_path_A) == ["a", "b", "c"] )

    simple_path_B = ["a", "b", "c", SPACER]

    assert(trim_terminal_spacers(simple_path_B) == ["a", "b", "c"] )

    simple_path_C = [SPACER, "a", "b", "c", SPACER]

    assert(trim_terminal_spacers(simple_path_C) == ["a", "b", "c"] )




def test_split_path_at_spacers():

    p = ['a', SPACER, 'b']
    assert(split_path_at_spacers(p) == [ ['a'], ['b'] ])

    p = [SPACER]
    assert(split_path_at_spacers(p) == [[]])

    p = [SPACER, 'c']
    assert(split_path_at_spacers(p) == [['c']])
    
    p = ['d', SPACER]
    assert(split_path_at_spacers(p) == [['d']])

    p = [SPACER, 'd', SPACER, 'e', SPACER, 'f', 'g', SPACER]
    assert(split_path_at_spacers(p) == [['d'], ['e'], ['f', 'g']])

    
