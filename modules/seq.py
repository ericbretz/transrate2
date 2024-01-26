import pysam
from multiprocessing import Pool

def process_task(args):
    split = args[0]
    bamfile = args[1]
    threads = args[2]
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refs = bam.references[split::threads]
    tmp_dct = {ref: {'seq': {'true':0, 'count':0}, 'stats': {'p_seqtrue': 0}} for ref in refs}

    for ref in refs:
        for fetch in bam.fetch(reference=ref):
            try:
                name_r  = fetch.reference_name
                len_r   = fetch.reference_length
                nm      = fetch.get_tag('NM', with_value_type = True)
                scale   = (len_r - 35) / len_r  # ! Why 35?
                tmp_dct[name_r]['seq']['count']      += 1
                tmp_dct[name_r]['seq']['true']       += (((len_r - nm[0]) / len_r) - scale) * (1 / (1 - scale))
                seqtrue = round(tmp_dct[name_r]['seq']['true'] / tmp_dct[name_r]['seq']['count'], 6)
                tmp_dct[name_r]['stats']['p_seqtrue'] = seqtrue
            except Exception as e:
                continue
    return tmp_dct

def seq(bamfile, bamdct, threads, single):

    with Pool(threads) as p:
        results = p.map(process_task, [[i, bamfile, threads] for i in range(threads)])

    for result in results:
        for key in result:
            if key in bamdct:
                for subkey in result[key]:
                    if isinstance(bamdct[key][subkey], dict):
                        bamdct[key][subkey].update(result[key][subkey])
                    else:
                        bamdct[key][subkey] += result[key][subkey]
            else:
                bamdct[key] = result[key]

    for k in bamdct:
        if bamdct[k]['stats']['p_seqtrue'] == 0:
            bamdct[k]['stats']['p_seqtrue'] = 1
    return bamdct