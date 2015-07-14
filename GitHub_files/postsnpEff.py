#!python3
# Manipulation of snpEff output

import os.path
import urllib.request  # Required for get_keggid/path()
from bs4 import BeautifulSoup as bs  # Required for get_keggid/path()
m = __import__('basic')


def extract_effs(snpEff_vcf, tag=''):
    effcomm, effs = m.read_file(snpEff_vcf, header=False)
    effs = [l[0:7] + l[7].split(',')[0].split('|') for l in effs if tag in l[7]]
    return effcomm, effs

def select_fa_descr(annots, fagenes):
    ph_fa = []
    phgenes = sorted([x[3] for x in annots])
    phmrnas = sorted([x[6] for x in annots])
    for gen in fagenes:
        if 'gene=' in gen:
            geneid = gen.split('gene=')[1].split(' CDS')[0]
            if geneid in phgenes:
                phgenes.remove(geneid)
                ph_fa.append(gen)
        else:
            geneid=gen.lstrip('>').split(' CDS')[0]
            if geneid in phmrnas:
                phmrnas.remove(geneid)
                ph_fa.append(gen)
    return ph_fa

def select_fasta(annots, fastafile):
    s, seqdict = 0, {}

    with open(fastafile) as fa:
        fagenes=[x for x in fa if x.startswith('>')]
    ph_fa = select_fa_descr(annots, fagenes)
    
    with open(fastafile) as fa:
        for line in fa:
            if line.startswith('>'):
                s = 0
                if line in ph_fa:
                    seq_id = line.lstrip('>')
                    s, seqdict[seq_id] = 1, ''
            elif s == 1:
                seqdict[seq_id] += line
            else:
                continue
    return seqdict

def write_subfasta(snpEff_vcf, fastafile, outPath, outFile, tag=''):
    annlist = [l[7:] for l in extract_effs(snpEff_vcf, tag=tag)[1]]
    fas_dict = select_fasta(annlist, fastafile)
    fas_list = []
    for seq in sorted(fas_dict.keys()):
        fas_list.append('>' + seq)
        fas_list.append(fas_dict[seq])
    m.write_output(fas_list, outPath, outFile)

def get_coords(snplist, snpEff_vcf):
    effs = extract_effs(snpEff_vcf)[1]
    coords = {l[0]: {int(l[1]): l[2]} for l in effs if l[2] in snplist}
    return coords

def write_mapchart_input(coords, outPath, outfile):
    ends = {'LG1': 20954957, 'LG2': 24538926, 'LG3': 32069524, 
    'LG4': 27214541, 'LG5': 28438568, 'LG6': 39347594, 'LG7': 22556666}
    map = []
    for lg in sorted(ends.keys()):
        map.append('group ' + lg + '\n')
        map.append('Start\t0\n')
        for pos in sorted(coords[lg]):
            map.append(coords[lg][pos]+'\t'+str(round(pos/1E6,2))+'\n')
        map.append('End\t' + str(round(ends[lg]/1E6,2)) + '\n')
    m.write_output(map, outPath, outfile)

def mapchart_by_eff(snpEff_vcf, outPath, outfile, tag=''):
    snplist = [l[2] for l in extract_effs(snpEff_vcf, tag=tag)[1]]
    coords = get_coords(snplist, snpEff_vcf)
    write_mapchart_input(coords, outPath, outfile)

def write_calls_by_eff(snpEff_vcf, callsfile, outPath, eff='HIGH', tas_out=False):
    snplist = [l[2] for l in extract_effs(snpEff_vcf, tag=tag)[1]]
    with open(callsfile) as c:
        header = [c.readline()]
        effcalls = [header] + [line for line in c if line.split()[0] in snplist]
    m.write_output(effcalls, outPath, eff+'EffPH.calls.txt')
    if tas_out == True:
        m.write_calls4tassel(eff+'EffPH.calls.txt', outPath)

def get_keggpath(gene_id):
    path = ''
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/www_bget?'
    response = urllib.request.urlopen(URL + FUN + gene_id)
    html = bs(response.read())
    for tr in html.find_all('tr'):
        if 'Pathway' in tr.text:
            p = tr.text
            p2 = [x.split(gene_id.split(':')[0])[0] for x in p.split(
                  '\xa0')[1:] if x != '']
            path = ', '.join(p2)
    return path

def get_keggid(unip_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/www_bconv?uniprot+'
    response = urllib.request.urlopen(URL + FUN + unip_id)
    html = bs(response.read())
    for tr in html.find_all('tr'):
        if 'Entry' in tr.text:
            e = tr.text
            org = e.split('Organism')[1].split()[0]
            ent = e.split('Entry')[1].split()[0]
            return org + ':' + ent

def create_u2map(blast_table):
        u2map = {}
        with open(blast_table) as f:
            h = next(f)
            for line in f:
                keggid = get_keggid(line.split('\t')[3])
                if keggid:
                    u2map[line.split('\t')[3]] = get_keggpath(keggid)
        return u2map

def get_name(name):
    if 'RecName' in name:
        name = name.split('; ')[0].lstrip('RecName: Full=')
    return name

def annotate_genes(snpEff, blast_table, tag='HIGH', paths=False):
    header = ['predicted_gene\tprobeset_id\ttype_effet\tHGVS.c\tHGVS.p' \
             + '\t#_hits\tLowest_E-value\tAccession\tDescription\n']
    if paths == True:
        u2map = create_u2map(blast_table)
        header = [header[0].rstrip() + '\tPathway\n']
    effs = extract_effs(snpEff, tag=tag)[1]
    with open(blast_table) as fa:
        f = fa.readlines()
        data = [x.rstrip().replace('"','').split('\t') for x in f]
        bl_res = {d[0]: d[1:4] + [get_name(d[4])] for d in data}

    for snp in effs:
        inf = [snp[10], snp[2], snp[8], snp[16], snp[17]] + bl_res[snp[13]]
        if paths == True:
            if bl_res[snp[13]][2] in u2map:
                inf.append(u2map[bl_res[snp[13]][2]])
            else:
                inf.append('')
        header.append('\t'.join(inf) + '\n')
    m.write_output(header, os.path.dirname(snpEff), 'annotation.txt')

