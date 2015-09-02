#! /usr/bin/env python

import argparse
import sys
from collections import defaultdict
from yeti.genomics.roitools import Transcript, SegmentChain
import pandas as pd
import numpy as np
import multiprocessing as mp
import subprocess as sp
import os

parser = argparse.ArgumentParser()
parser.add_argument('tfamstem', help='Transcript family information generated by make_tfams.py. '
                                     'Both TFAMSTEM.txt and TFAMSTEM.bed should exist or an error will result.')
parser.add_argument('orfstore', help='Path to pandas HDF store containing ORFs to regress; generated by find_ORFs.py')
parser.add_argument('cdsstore', help='File to which to output table with CDS information. Formatted as pandas HDF store (preferred extension is .h5; '
                                     'table name is "annot_CDSs").')
parser.add_argument('--inbed', type=argparse.FileType('rU'), default=sys.stdin, help='Transcriptome BED-file. Annotated CDSs are assumed to be bona '
                                                                                     'fide CDSs, unless --ignoreannotations is set. (Default: stdin)')
parser.add_argument('--ignoreannotations', action='store_true', help='If flag is set, CDS annotations in INBED will be ignored. Typically used in '
                                                                     'conjunction with --extracdsbeds')
parser.add_argument('--extracdsbeds', nargs='+', help='Extra bed file(s) containing additional annotated CDSs beyond (or instead of) those in inbed. '
                                                      'If transcript names are repeated across these files, sources of annotated CDSs may become '
                                                      'ambiguous, but it will not result in an error or warning.')
parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes to run. Defaults to 1 but recommended to use more (e.g. 12-16)')
opts = parser.parse_args()

# hash transcripts by ID for easy reference later
bedlinedict = {line.split()[3]: line for line in opts.inbed}
if not bedlinedict:
    raise EOFError('Insufficient input or empty file provided')

tfamtids = defaultdict(list)
with open('%s.txt' % opts.tfamstem, 'rU') as tfamtable:
    for line in tfamtable:
        ls = line.strip().split()
        tfamtids[ls[1]].append(ls[0])

with open('%s.bed' % opts.tfamstem, 'rU') as tfambed:
    tfambedlines = {line.split()[3]: line for line in tfambed}

if not opts.ignoreannotations:
    annot_tfam_lookups = [tfamtids]
    annot_tid_lookups = [bedlinedict]
else:
    annot_tfam_lookups = []
    annot_tid_lookups = []
if opts.extracdsbeds:
    import pybedtools  # to handle identifying which tfams get the extra CDSs - otherwise would need to replicate a lot of intersection functionality
    tfambedtool = pybedtools.BedTool('%s.bed' % opts.tfamstem)
    for cdsbedfile in opts.extracdsbeds:
        with open(cdsbedfile, 'rU') as cdsbed:
            annot_tid_lookups.append({line.split()[3]: line for line in cdsbed})  # as usual, hash bed lines by transcript ID
        annot_tfam_lookups.append(defaultdict(list))
        for line in tfambedtool.intersect(pybedtools.BedTool(cdsbedfile), split=True, s=True, wa=True, wb=True):
            annot_tfam_lookups[-1][line[3]].append(line[15])
# after this has finished, each element of annot_tfam_lookup will be a dictionary mapping tfams to lists of transcript IDs in the annotation bed files
# similarly, each element of annot_tid_lookup will map transcript IDs to BED lines

tfams_with_annots = set(sum([x.keys() for x in annot_tfam_lookups], []))

def find_annot_CDSs(tfam, tfam_ORFs):
    currtfam = SegmentChain.from_bed(tfambedlines[tfam])
    chrom = currtfam.chrom
    strand = currtfam.strand
    tfam_genpos = np.array(currtfam.get_position_list(stranded=True))
    # annot_CDS_dfs = []
    found_CDS_info = []
    unfound_CDS_info = []
    for (annot_tfam_lookup, annot_tid_lookup) in zip(annot_tfam_lookups, annot_tid_lookups):
        if tfam in annot_tfam_lookup:
            # annot_CDS_dfs.append(pd.DataFrame.from_items([('tfam', tfam),
            #                                               ('tid', annot_tfam_lookup[tfam]),
            #                                               ('chrom', chrom),
            #                                               ('gcoord', 0),
            #                                               ('gstop', 0),
            #                                               ('strand', strand),
            #                                               ('AAlen', 0),
            #                                               ('ORF_name', '')]))  # ORF_name left blank if not found
            for (annot_tidx, annot_tid) in enumerate(annot_tfam_lookup[tfam]):
                curr_trans = Transcript.from_bed(annot_tid_lookup[annot_tid])
                if curr_trans.cds_start is not None and curr_trans.cds_end is not None:
                    found = False
                    curr_gcoord = curr_trans.get_genomic_coordinate(curr_trans.cds_start)[1]
                    # annot_CDS_dfs[-1].loc[annot_tidx, 'gcoord'] = curr_gcoord
                    curr_gstop = curr_trans.get_genomic_coordinate(curr_trans.cds_end-1)[1]+(strand == '+')*2-1
                    # annot_CDS_dfs[-1].loc[annot_tidx, 'gstop'] = curr_gstop
                    shared_start = (tfam_ORFs['gcoord'] == curr_gcoord)
                    shared_stop = (tfam_ORFs['gstop'] == curr_gstop)
                    curr_cds_pos = curr_trans.get_cds().get_position_set()
                    # annot_CDS_dfs[-1].loc[annot_tidx, 'AAlen'] = len(curr_cds_pos)/3 - 1
                    shared_len = (tfam_ORFs['tstop']-tfam_ORFs['tcoord'] == len(curr_cds_pos))
                    possible_ORFs = shared_start & shared_stop & shared_len
                    if possible_ORFs.any() and curr_cds_pos.issubset(tfam_genpos):
                        # curr_cds_genidx = np.flatnonzero(np.in1d(tfam_genpos, list(curr_cds_pos), assume_unique=True))
                        for (ORF_name, tid, tcoord, tstop) in tfam_ORFs.loc[possible_ORFs, ['ORF_name', 'tid', 'tcoord', 'tstop']].itertuples(False):
                            curr_ORF_pos = SegmentChain.from_bed(bedlinedict[tid]).get_genomic_coordinate(np.arange(tcoord, tstop))[1]
                            if curr_cds_pos.issubset(curr_ORF_pos):  # this actually tests for equality, as they are guaranteed the same length
                                found = True
                                # annot_CDS_dfs[-1].loc[annot_tidx, 'ORF_name'] = ORF_name
                                found_CDS_info.append((tfam,
                                                       annot_tid,
                                                       chrom,
                                                       curr_gcoord,
                                                       curr_gstop,
                                                       strand,
                                                       len(curr_cds_pos)/3-1,
                                                       ORF_name))
                                break
                    if not found:
                        unfound_CDS_info.append((tfam,
                                                 annot_tid,
                                                 chrom,
                                                 curr_gcoord,
                                                 curr_gstop,
                                                 strand,
                                                 len(curr_cds_pos)/3-1))
    return (pd.DataFrame(found_CDS_info, columns=['tfam', 'tid', 'chrom', 'gcoord', 'gstop', 'strand', 'AAlen', 'ORF_name']),
            pd.DataFrame(unfound_CDS_info, columns=['tfam', 'tid', 'chrom', 'gcoord', 'gstop', 'strand', 'AAlen']))
    # return pd.concat(annot_CDS_dfs, ignore_index=True)


def find_annot_CDSs_by_chrom(chrom_to_do):
    named_ORFs = pd.read_hdf(opts.orfstore, 'all_ORFs', where="chrom == '%s' and tstop > 0" % chrom_to_do,
                             mode='r', columns=['tfam', 'tmap', 'tid', 'tcoord', 'tstop', 'chrom', 'gcoord', 'gstop', 'strand', 'ORF_name']) \
        .drop_duplicates('ORF_name')
    named_ORFs = named_ORFs[named_ORFs['tfam'].isin(tfams_with_annots)]  # don't bother looking if there aren't any annotated CDSs to be found
    if not named_ORFs.empty:
        return tuple(
            [pd.concat(dfs, ignore_index=True) for dfs in zip(*[find_annot_CDSs(tfam, tfam_ORFs) for (tfam, tfam_ORFs) in named_ORFs.groupby('tfam')])])
    else:
        return pd.DataFrame(), pd.DataFrame()
    # return pd.concat([find_annot_CDSs(tfam, tfam_ORFs) for (tfam, tfam_ORFs) in named_ORFs.groupby('tfam')], ignore_index=True)

with pd.get_store(opts.orfstore, mode='r') as orfstore:
    chroms = orfstore.select('all_ORFs/meta/chrom/meta').values  # because saved as categorical, this is a list of all chromosomes
workers = mp.Pool(opts.numproc)
# annot_CDSs = pd.concat(workers.map(find_annot_CDSs_by_chrom, chroms), ignore_index=True)
(found_CDSs, unfound_CDSs) = [pd.concat(dfs, ignore_index=True) for dfs in zip(*workers.map(find_annot_CDSs_by_chrom, chroms))]
workers.close()

for catfield in ['chrom', 'strand']:
    found_CDSs[catfield] = found_CDSs[catfield].astype('category')  # saves disk space and read/write time
    unfound_CDSs[catfield] = unfound_CDSs[catfield].astype('category')  # saves disk space and read/write time

origname = opts.cdsstore+'.tmp'
with pd.get_store(origname, mode='w') as outstore:
    outstore.put('found_CDSs', found_CDSs, format='t', data_columns=True)
    outstore.put('unfound_CDSs', unfound_CDSs, format='t', data_columns=True)
# annot_CDSs.to_hdf(origname, 'annot_CDSs', format='t', data_columns=True)
sp.call(['ptrepack', origname, opts.cdsstore])  # repack for efficiency
os.remove(origname)
