#!/usr/bin/env python

import os, sys
import string
import argparse
import re

parser = argparse.ArgumentParser(description="""
creates mut file from snpEFF vcf
""", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--vcf', metavar = '<file.vcf>', required=True,
help="""VCF file.  <REQUIRED>\n\n""")

parser.add_argument('--depth', metavar = '<file.vcf>', default = 10,
help="""Sample minimum DP (default:10)\n\n""")

parser.add_argument('--targets', metavar = '<targets.txt>',
help="""List of targets to create scoring for\n\n""")

parser._optionals.title = "Program Options"

args = parser.parse_args()

VHEADER = re.compile('^##', re.IGNORECASE)
SHEADER = re.compile('^#', re.IGNORECASE)
VAR_GT  = re.compile('(0/1|1/1|1/0)', re.IGNORECASE)
DP = re.compile(':(\d+):', re.IGNORECASE)
GT_CHECK = re.compile('GT:DP', re.IGNORECASE)
OUT_HEADER = "chr\tstart\tend\tsample\ttype\tgene\tallele\ttranscript\tHGVS.c\tcDNA position\tCDS position\tHGVS.p\tProtein position\tLOF\tNMD\n"
SCORE_EFFECT = {'MODIFIER' : 0, 'LOW' : 1, 'MODERATE' : 2, 'HIGH' : 4}


def check_file(file):
    if not os.path.isfile(file):
        sys.stderr.write("could not locate {}\n".format(file))
        raise SystemExit(1)


def target_list(targets, score):
    with open(targets) as topen:
        for line in topen:
            line = line.rstrip()
            score[line] = {}


def make_mut(vcf, d, score):
    samples = []
    score['global'] = {}
    mut_out = open('./mut_out.mut', 'w')
    with open(vcf) as vopen:
        for line in vopen:
            if VHEADER.match(line):
                continue
            line   = line.rstrip()
            fields = line.split("\t")
            if SHEADER.match(line):
                mut_out.write(OUT_HEADER)
                for i in range(9, len(fields)):
                    samples.append(fields[i])
                    for r in score:
                        score[r][fields[i]] = 0
            else:
                chr   = fields[0]
                start = fields[1]
                stop  = int(fields[1]) + len(fields[3]) - 1
                info = { i[0] : i[1] for i in [f.split("=") for f in fields[7].split(";")] };
                #ANN is last field in INFO 
                ann = info['ANN']
                #most impactful effect is listed first
                eff = ann.split(",")
                eff1 = eff[0].split("|");
                type = eff1[1]
                gene = eff1[3]
                # if chr in score, if start >= min && start <= max
                effect = eff1[2]
                escore = SCORE_EFFECT[effect]

                allele = eff1[0];
                transcript = eff1[6];
                hgvs_c = eff1[9]
                hgvs_p = eff1[10]
                cdna_pos = eff1[11]
                cds_pos = eff1[12]
                protein_pos = eff1[13]
                if 'LOF' in info:
                    lof = info['LOF']
                else:
                    lof = ""
                if 'NMD' in info:
                    nmd = info['NMD']
                else:
                    nmd = ""

                for i in range(9, len(fields)):
                    if VAR_GT.search(fields[i]):
                        sample = samples[i - 9]
                        check = 1
                        if DP.search(fields[i]):
                            if DP.search(fields[i]).group(1) >= d:
                                check = 1
                            else:
                                check = 0
                        if check:
                            score['global'][sample] += escore
                            if gene in score:
                                score[gene][sample] += escore
                            field = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                      chr,
                                                                      start,
                                                                      stop,
                                                                      sample,
                                                                      type,
                                                                      gene,
                                                                      allele,
                                                                      transcript,
                                                                      hgvs_c,
                                                                      cdna_pos,
                                                                      cds_pos,
                                                                      hgvs_p,
                                                                      protein_pos,
                                                                      lof,
                                                                      nmd,
                                                                      )
                            mut_out.write(field)
    mut_out.close()


if __name__ == "__main__":
    vcf = args.vcf
    depth = args.depth
    targets = args.targets
    score = {}
    results = {}
    check_file(vcf)
    if targets:
        check_file(targets)
        target_list(targets, score)
    make_mut(vcf, depth, score)
    header = ['sample']
    for r in score:
        if r == 'global':
            continue
        header.append(str(r + "_score"))
        for s in score[r]:
            rscore = str(score[r][s])
            if s in results:
                results[s].append(rscore)
            else:
                results[s] = [s, rscore]

    header.append('global_score')
    for s in score['global']:
        if not s in results:
            rscore = str(score['global'][s])
            results[s] = [s, rscore]
        else:
            rscore = str(score['global'][s])
            results[s].append(rscore)

    file = "attributes.txt"
    r_out = open(file, 'w')
    r_out.write("\t".join(header) + "\n")
    for s in results:
        r_out.write("\t".join(results[s]) + "\n")
    r_out.close()

