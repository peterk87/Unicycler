#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This script is a stand-alone tool for running iterative Pilon polishing rounds. It is quite
similar to the Pilon polishing stage of Unicycler, and depends on Unicycler being installed. It is
intended mainly for Pilon polishing of long read assemblies for which a Unicycler hybrid assembly
is not appropriate (e.g. very poor Illumina coverage).

Note that it takes a GFA graph as input (so it can tell which sequences are circular).

Usage:
  pilon_polish.py --input unpolished.gfa --output polished -1 reads_1.fastq.gz -2 reads_2.fastq.gz

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

import unicycler.assembly_graph
import unicycler.log
import unicycler.pilon_func
import unicycler.unicycler
from unicycler.blast_func import find_start_gene, CannotFindStart
from unicycler.misc import int_to_str, float_to_str, print_table, gfa_path


def rotate_completed_replicons(graph,
                               start_genes=None,
                               start_gene_id=90.0,
                               start_gene_cov=95.0,
                               out='.',
                               makeblastdb_path='makeblastdb',
                               tblastn_path='tblastn',
                               keep=1):
    # get default start_genes.fasta from script root dir
    if start_genes is None:
        start_genes = default_start_genes_fasta()

    completed_replicons = graph.completed_circular_replicons()
    if len(completed_replicons) > 0:
        unicycler.log.log_section_header('Rotating completed replicons')
        unicycler.log.log_explanation('Any completed circular contigs (i.e. single contigs which have one '
                                      'link connecting end to start) can have their start position changed '
                                      'without altering the sequence. For consistency, Unicycler now '
                                      'searches for a starting gene (dnaA or repA) in each such contig, and '
                                      'if one is found, the contig is rotated to start with that gene on '
                                      'the forward strand.')

        rotation_result_table = [['Segment', 'Length', 'Depth', 'Starting gene', 'Position',
                                  'Strand', 'Identity', 'Coverage']]
        blast_dir = os.path.join(out, 'blast')
        if not os.path.exists(blast_dir):
            os.makedirs(blast_dir)
        completed_replicons = sorted(completed_replicons, reverse=True,
                                     key=lambda x: graph.segments[x].get_length())
        rotation_count = 0
        for completed_replicon in completed_replicons:
            segment = graph.segments[completed_replicon]
            sequence = segment.forward_sequence

            try:
                seg_name = str(segment.number)
            except AttributeError:
                seg_name = segment.full_name

            unicycler.log.log('Segment ' + seg_name + ':', 2)
            rotation_result_row = [seg_name, int_to_str(len(sequence)),
                                   float_to_str(segment.depth, 2) + 'x']
            try:
                blast_hit = find_start_gene(sequence, start_genes, start_gene_id,
                                            start_gene_cov, blast_dir, makeblastdb_path,
                                            tblastn_path)
            except CannotFindStart:
                rotation_result_row += ['none found', '', '', '', '']
            else:
                rotation_result_row += [blast_hit.qseqid, int_to_str(blast_hit.start_pos),
                                        'reverse' if blast_hit.flip else 'forward',
                                        '%.1f' % blast_hit.pident + '%',
                                        '%.1f' % blast_hit.query_cov + '%']
                segment.rotate_sequence(blast_hit.start_pos, blast_hit.flip)
                rotation_count += 1
            rotation_result_table.append(rotation_result_row)

        unicycler.log.log('', 2)
        print_table(rotation_result_table, alignments='RRRLRLRR', indent=0,
                    sub_colour={'none found': 'red'})
        if rotation_count and keep > 0:
            out_dir = Path(out)
            graph.save_to_gfa(str(out_dir / 'rotated.gfa'), newline=True)
        if keep < 3 and os.path.exists(blast_dir):
            shutil.rmtree(blast_dir, ignore_errors=True)


def default_start_genes_fasta():
    out = None
    path = Path(__file__).resolve().absolute().parent.parent
    for x in path.rglob('**/start_genes.fasta'):
        if x:
            out = x
            break
    if out is None:
        raise FileNotFoundError(f'"start_genes.fasta" not found in {path}')
    return out


def main():
    args = get_arguments()
    unicycler.log.logger = unicycler.log.Log(None, args.verbosity)
    check_dependencies(args)
    polish_dir = os.getcwd()

    graph = unicycler.assembly_graph.AssemblyGraph(args.input, None)
    if not args.no_polish:
        args.keep = 1
        insert_size_1st, insert_size_99th = unicycler.pilon_func.get_insert_size_range(graph, args,
                                                                                       polish_dir)
        args.keep = 3

        for i in range(10):
            if i > 0:
                graph.rotate_circular_sequences()

            round_num = i + 1
            change_count = unicycler.pilon_func.polish_with_pilon(graph, args, polish_dir,
                                                                  insert_size_1st, insert_size_99th,
                                                                  round_num, 'all')
            input_filename = str(round_num) + '_polish_input.fasta'
            paired_bam_filename = str(round_num) + '_paired_alignments.bam'
            unpaired_bam_filename = str(round_num) + '_unpaired_alignments.bam'
            for f in [input_filename,
                      input_filename + '.1.bt2', input_filename + '.2.bt2',
                      input_filename + '.3.bt2', input_filename + '.4.bt2',
                      input_filename + '.rev.1.bt2', input_filename + '.rev.2.bt2',
                      paired_bam_filename, paired_bam_filename + '.bai',
                      unpaired_bam_filename, unpaired_bam_filename + '.bai']:
                try:
                    os.remove(os.path.join(polish_dir, f))
                except FileNotFoundError:
                    pass

            if not change_count:
                break
    else:
        unicycler.log.log('Skipping Pilon polishing!')
    if args.do_rotate:
        rotate_completed_replicons(graph)
    graph.save_to_gfa(args.output + '.gfa')
    graph.save_to_fasta(args.output + '.fasta')
    unicycler.log.log('')


def get_arguments():
    parser = argparse.ArgumentParser(description='Pilon polishing tool for long read assemblies')
    parser.add_argument('-i', '--input', required=True,
                        help='Input GFA to be polished')
    parser.add_argument('-o', '--output', required=True,
                        help='Output prefix for GFA and FASTA files')

    parser.add_argument('-1', '--short1', required=False,
                        help='FASTQ file of first short reads in each pair')
    parser.add_argument('-2', '--short2', required=False,
                        help='FASTQ file of second short reads in each pair')
    parser.add_argument('-s', '--unpaired', required=False,
                        help='FASTQ file of unpaired short reads')

    parser.add_argument('-t', '--threads', type=int, required=False, default=16,
                        help='Number of threads used')

    parser.add_argument('--min_polish_size', type=int, default=10000,
                        help='Contigs shorter than this value (bp) will not be polished '
                             'using Pilon')

    parser.add_argument('--bowtie2_path', type=str, default='bowtie2',
                        help='Path to the bowtie2 executable')
    parser.add_argument('--bowtie2_build_path', type=str, default='bowtie2-build',
                        help='Path to the bowtie2_build executable')
    parser.add_argument('--pilon_path', type=str, default='pilon',
                        help='Path to a Pilon executable or the Pilon Java archive file')
    parser.add_argument('--java_path', type=str, default='java',
                        help='Path to the java executable')
    parser.add_argument('--samtools_path', type=str, default='samtools',
                        help='Path to the samtools executable')
    parser.add_argument('--no-polish', action='store_true')
    parser.add_argument('--do-rotate', action='store_true')

    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if (args.short1 and not args.short2) or (args.short2 and not args.short1):
        unicycler.misc.quit_with_error('you must use both --short1 and --short2 or neither')

    if not args.short1 and not args.short2 and not args.unpaired:
        unicycler.misc.quit_with_error('no input reads provided (--short1, --short2, --unpaired')

    # Change some arguments to full paths.
    if args.short1:
        args.short1 = os.path.abspath(args.short1)
    if args.short2:
        args.short2 = os.path.abspath(args.short2)
    if args.unpaired:
        args.unpaired = os.path.abspath(args.unpaired)

    pilon_path, _, _ = unicycler.misc.pilon_path_and_version(args.pilon_path, args.java_path, args)

    args.verbosity = 2
    args.keep = 3

    return args


def check_dependencies(args):
    unicycler.log.log('\nDependencies:')
    if args.verbosity <= 1:
        program_table = [['Program', 'Version', 'Status']]
    else:
        program_table = [['Program', 'Version', 'Status', 'Path']]

    bowtie2_build_path, bowtie2_build_version, bowtie2_build_status = \
        unicycler.misc.bowtie2_build_path_and_version(args.bowtie2_build_path)
    bowtie2_path, bowtie2_version, bowtie2_status = \
        unicycler.misc.bowtie2_path_and_version(args.bowtie2_path)
    samtools_path, samtools_version, samtools_status = \
        unicycler.misc.samtools_path_and_version(args.samtools_path)
    java_path, java_version, java_status = \
        unicycler.misc.java_path_and_version(args.java_path)
    pilon_path, pilon_version, pilon_status = \
        unicycler.misc.pilon_path_and_version(args.pilon_path, args.java_path, args)

    bowtie2_build_row = ['bowtie2-build', bowtie2_build_version, bowtie2_build_status]
    bowtie2_row = ['bowtie2', bowtie2_version, bowtie2_status]
    samtools_row = ['samtools', samtools_version, samtools_status]
    java_row = ['java', java_version, java_status]
    pilon_row = ['pilon', pilon_version, pilon_status]
    if args.verbosity > 1:
        bowtie2_build_row.append(bowtie2_build_path)
        bowtie2_row.append(bowtie2_path)
        samtools_row.append(samtools_path)
        java_row.append(java_path)
        pilon_row.append(pilon_path)
    program_table.append(bowtie2_build_row)
    program_table.append(bowtie2_row)
    program_table.append(samtools_row)
    program_table.append(java_row)
    program_table.append(pilon_row)

    row_colours = {}
    for i, row in enumerate(program_table):
        if 'not used' in row:
            row_colours[i] = 'dim'
        elif 'too old' in row or 'not found' in row or 'bad' in row or 'Python problem' in row:
            row_colours[i] = 'red'

    unicycler.misc.print_table(program_table, alignments='LLLL', row_colour=row_colours,
                               max_col_width=60, sub_colour={'good': 'green'})
    unicycler.log.log('')

    quit_if_dependency_problem(bowtie2_build_status, bowtie2_status, samtools_status, java_status,
                               pilon_status, args)


def quit_if_dependency_problem(bowtie2_build_status, bowtie2_status, samtools_status, java_status,
                               pilon_status, args):
    if all(x == 'good' or x == 'not used'
           for x in [bowtie2_build_status, bowtie2_status, samtools_status, java_status,
                     pilon_status]):
        return

    unicycler.log.log('')
    if bowtie2_build_status == 'not found':
        unicycler.misc.quit_with_error('could not find bowtie2-build - specify its location using '
                                       '--bowtie2_build_path')
    if bowtie2_status == 'not found':
        unicycler.misc.quit_with_error('could not find bowtie2 - specify its location using '
                                       '--bowtie2_path')
    if samtools_status == 'not found':
        unicycler.misc.quit_with_error('could not find samtools - specify its location using '
                                       '--samtools_path')
    if java_status == 'not found':
        unicycler.misc.quit_with_error('could not find java - specify its location using '
                                       '--java_path')
    if java_status == 'bad':
        unicycler.misc.quit_with_error('Java did not run correctly - specify its location '
                                       'using --java_path')
    if pilon_status == 'not found':
        unicycler.misc.quit_with_error('could not find pilon or pilon*.jar - specify its location '
                                       'using --pilon_path')
    if pilon_status == 'bad':
        unicycler.misc.quit_with_error('Pilon was found (' + args.pilon_path + ') but does not '
                                                                               'work - either fix it or specify a different location '
                                                                               'using --pilon_path')

    # Code should never get here!
    unicycler.misc.quit_with_error('Unspecified error with dependencies')


if __name__ == '__main__':
    main()
