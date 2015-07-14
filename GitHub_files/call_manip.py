#!python3
#Handling of calls from Genotyping Console

import os.path
m = __import__('basic')

def PHvar_calls(PHpath, callspath, outPath):
    PolyHigh = m.file2set(PHpath)
    comments, varnames, calls = m.read_file(callspath)

    varTypes = {'probeset_id', 'Hibrid', 'Dover', 'Camarosa', '122_', 'a550466'}
    varList = [var for var in varnames if True in [c in var for c in varTypes]]
    idxs = [varnames.index(x) for x in varnames if x in varList]

    phlist = ['\t'.join([v for v in varList]) + '\n']
    phlist.extend(['\t'.join([c[i] for i in idxs])+'\n' for c in calls if c[0] in PolyHigh])
    m.write_output(phlist, outPath, 'polyhigh.calls.txt')
   
def join_dupl_calls(PHfile, outPath):
    dupl_vars = ['Camarosa', 'Hibrid', '_Dover', '122_48', '122_50']
    comments, varnames, calls = m.read_file(PHfile)

    EDovidx = varnames.index('1F_E_Dover.CEL')
    varnames.remove(varnames[EDovidx])
    varidx = {name: [] for name in dupl_vars}
    for name in dupl_vars:
        vars = [var for var in varnames if name in var]
        for var in vars:
            varidx[name].append(varnames.index(var))
            varnames.remove(var)
    varnames = varnames + ['E_Dover'] + dupl_vars
    
    snpcalls = {}
    nonsegrSNPs = ['probeset_id']
    a = 0
    while a < len(calls):
        line = calls[a]
        EDov = line.pop(EDovidx)
        for name in dupl_vars:
            snpcalls[name] = [line.pop(idx) for idx in varidx[name]]
            dif_calls = sorted(list(set(snpcalls[name])))
            if len(dif_calls) == 1:    #All calls identical
                snpcalls[name] = dif_calls[0]
            elif len(dif_calls) == 2 and '-1' in dif_calls:  # With undefined call
                snpcalls[name] = dif_calls[1]
            elif len(snpcalls[name]) > 2 and len(dif_calls) == 2 and int(
                dif_calls[1])-int(dif_calls[0]) == 1: #AA AB or AB BB if one twice
                snpcalls[name] = max(set(snpcalls[name]), key=snpcalls[name].count)
            else:
                snpcalls[name] = '-1'

        if len(set(line[1:len(line)])) < 3:   # non-polymorphic SNPs
            nonsegrSNPs.append(line[0])
            calls.remove(line)
            continue

        line.extend([EDov] + [snpcalls[name] for name in dupl_vars])
        a += 1
        
    m.change_dir(outPath)
    with open('jmPHcalls.txt', mode='w') as jm:
        with open('clPHcalls.txt', mode='w') as cl:
            jm.writelines('\t'.join([v for v in varnames[0:len(varnames)-2]]) + '\n')
            cl.writelines('\t'.join([v for v in varnames]) + '\n')
            for line in calls:
                jm.writelines('\t'.join([c for c in line[0:len(line)-2]]) + '\n')
                cl.writelines('\t'.join([c for c in line]) + '\n')
    with open('nonsegSNPs.txt', 'w') as ns:
        ns.writelines('\n'.join([snp for snp in nonsegrSNPs]))

def match_mapsnps2ref(annots, map, coreidpos=1):
    fvbs = m.single_dict(annots, store=5)  # d[snp] = Fvb
    lgs = m.single_dict(annots, store=3)   # d[snp] = LG
    mappedsnps = {}  # d[SNP] = [HG, cM]
    unmatched = ['probeset_id\tfvb\tLG\tcore\n']
    for line in map:
        if line[coreidpos-1] not in fvbs:
            continue
        if fvbs[line[coreidpos-1]].strip('Fvb') == line[6].strip("LG").split("-")[0]:
            mappedsnps[line[coreidpos-1]] = [line[6], line[7]]
        elif '/' in line[6]:
            mappedsnps[line[coreidpos-1]] = ['LG'+fvbs[line[coreidpos-1]].strip('Fvb'
                                            )+'-?', line[7]]
        else:
            unmatched.append('\t'.join([line[coreidpos-1], fvbs[line[coreidpos-1]], 
                lgs[line[coreidpos-1]], line[6]])+ '\n')
    return mappedsnps, unmatched  # dict and list

def write_unmatch(annotFile, coremap, outdir, coreidpos=1):
    annots = m.read_file(annotFile)[2]
    map = m.read_file(coremap)[2]
    mappedsnps, unmatched = match_mapsnps2ref(annots, map, coreidpos=coreidpos)
    m.write_output(unmatched, outdir, 'coreLGmismatch.txt')

def prepare_calls4clustering(annotFile, coremap, PH_callsfile, outdir, coreidpos=1):
    annots = m.subfile(PH_callsfile, annotFile)
    map = m.subfile(PH_callsfile, coremap, f2=coreidpos)
    phcomm, phhead, phcalls = m.read_file(PH_callsfile)
    fvbs = m.single_dict(annots, store=5)  # d[snp] = Fvb

    mappedsnps, unmatched = match_mapsnps2ref(annots, map, coreidpos=coreidpos)

    callshead = '\t'.join([h for h in phhead]) + '\n'
    map_calls = {'Fvb' + str(chr): [callshead] for chr in range(1,8)}
    unmap_calls = {'Fvb' + str(chr): [callshead] for chr in range(1,8)}
    for call in phcalls:
        fv = fvbs[call[0]]
        if call[0] in mappedsnps:   # SNP HG Pos calls
            newcall = [call[0]] + mappedsnps[call[0]] + [c for c in call[1:]]
            map_calls[fv].append('\t'.join(newcall) + '\n')
        else:       # SNP calls
            unmap_calls[fv].append('\t'.join([c for c in call]) + '\n')
    
    for x in range(1,8):  # 7 chr
        m.write_output(map_calls['Fvb'+str(x)], outdir, 'LG'+str(x)+'.txt')
        m.write_output(unmap_calls['Fvb'+str(x)], outdir, 'Fvb'+str(x)+'.txt')
    m.write_output(unmatched, outdir, 'coreLGmismatch.txt')

    top = ''.join(['Calls distribution:\nTotal PH:\t', str(len(phcalls)), '\nAnnot PH:\t', 
        str(len(annots)-1), '\nTotal core PH:\t', str(len(map)-1), '\nUsed core PH:\t', 
        str(len(mappedsnps)), '\nUnmatched core:\t', str(len(unmatched)-1), '\n\n'])
    chrs = '\t'.join(['Chr']+[chr for chr in sorted(map_calls.keys())]) + '\n'
    mapped = '\t'.join(['Core']+[str(len(map_calls[fv])-1) for fv in sorted(map_calls.keys())]) + '\n'
    unmap = '\t'.join(['PH']+[str(len(unmap_calls[fv])-1) for fv in sorted(unmap_calls.keys())]) + '\n'
    m.write_output([top, chrs, mapped, unmap], outdir, 'stats.txt')

def rephase_calls(indir, perc='90'):
    inv = {'a': 'b', 'b': 'a', 'h': 'h', '-': '-', 'NORMAL': 'INVERS', 
            'INVERS': 'NORMAL'}
    for i in range(1,8):
        lg = 'LG' + str(i)
        path = os.path.join(indir, lg, perc)
        jm_in = os.path.join(path, 'jm_step1.txt')
        phas_in = os.path.join(path, 'phasing.txt')
        header, data = m.read_file(jm_in)[1:]
        phase = m.single_dict(m.read_file(phas_in)[2])
        for call in data:
            call.remove(call[4])
            if call[0] == 'INVERS':
                snp = call[2].split('_')[1].split('LG')[0]
                call[4:] = [inv[c] for c in call[4:]]
                phase[snp] = inv[phase[snp]]
        header.remove(header[4])
        newdata = ['\t'.join(header[2:]) + '\n'] + \
                  ['\t'.join([c for c in call[2:]]) + '\n' for call in data]
        m.write_output(newdata, path, 'calls4step2.txt')
        m.write_output(['probeset_id\tphase\n'] + ['\t'.join([key, phase[key]]
            ) + '\n' for key in sorted(phase.keys())], path, 'phasing_2.txt')

def totstats(lg, lg_frags):
    cols = list(zip(*lg_frags))
    lg = [[lg, len(lg_frags), max(cols[1]), sum(cols[3]), sum(cols[4])]]
    for hg in sorted(set(cols[0])):
        sub = list(zip(*[x for x in lg_frags if x[0] == hg]))
        lg.append([hg, len(sub[0]), max(sub[1]), sum(sub[3]), sum(sub[4])])
    return ['\t'.join([str(x) for x in l]) + '\n' for l in lg]

def fragstats(frag):
    homs = ['LG' + x[1].split('LG')[1] for x in frag]
    homeol = max(set(homs), key=homs.count)
    others = ''
    if len(set(homs)) > 1:
        hg_set = set(homs)
        hg_set.remove(homeol)
        others = ', '.join(list(hg_set))
    n = len(frag)
    a = len([x for x in frag if x[1].startswith('A_')])
    return [homeol, n, others, n-a, a]

def mapstats(indir, perc='90'):
    stats = ['Group\t#Fragments\tMax_#SNPs\tMapped\tAssoc\n']
    tot_frag = ['Homeolog\tOther_hom\t#_SNPs\tMapped\tAssoc\n']
    for i in range(1,8):
        frag = []
        lg_frags = []
        lg = 'LG' + str(i)
        mappath = os.path.join(indir, lg, perc, 'jm_step2.txt')
        with open(mappath) as f:
            for line in f:
                if len(line.rstrip()) == 0:
                    tot_frag.append('\t'.join(str(x) for x in fragstats(frag[1:]))+'\n')
                    lg_frags.append(fragstats(frag[1:]))
                    frag = []
                    continue
                frag.append(line.rstrip().split('\t'))
        stats.extend(totstats(lg, lg_frags))
    m.write_output(stats, indir, 'jm_stats.txt')
    m.write_output(tot_frag, indir, 'jm_fragments.txt')
        
def nonclus_PH(indir, perc='90', call2jm = False):
    for i in range(1,8):
        lg = 'LG' + str(i)
        path = os.path.join(indir, lg, perc)
        clus_in = os.path.join(path, 'clusterout.txt')
        ph_in = os.path.join(indir, 'Fvb' + str(i) + '.txt')
        clus = m.file2set(clus_in, field=3)
        with open(ph_in) as f:
            uncl = [x for x in f if x.split()[0] not in clus]
        m.write_output(uncl, path, 'unclus_ph.txt')
        if call2jm == True:
            m.GCcalls2Joinmap(os.path.join(path, 'unclus_ph.txt'))
