import sys
import pysam
import re


if len(sys.argv) != 3:
    print(f'Usage: python {sys.argv[0]} input.vcf(.gz) output.vcf(.gz)')
    print('Please index the inputted vcf first')
    sys.exit(1)


def get_new_samples(vcf):
    new_samples = []
    for sample in vcf.header.samples:
        if re.search(r'-(M|P).+$', sample): ## the format of sample names???
            new_sample = '-'.join(sample.split('-')[:-1])
        else:
            new_sample = sample.split('_')[0]
        
        if new_sample not in new_samples:
            new_samples.append(new_sample)
    return new_samples


f_in = pysam.VariantFile(sys.argv[1])
new_header = pysam.VariantHeader()
for i in f_in.header.records:
#    if i.key == 'INFO' and list(i.values())[0] == 'AF':
#        new_header.add_meta('INFO', items=[('ID', 'AF'), ('Number', '.'), ('Type', 'Float'), ('Description', 'Estimated allele frequency in the range (0,1]')])
#    elif i.key == 'INFO' and list(i.values())[0] == 'AT':
#        new_header.add_meta('INFO', items=[('ID', 'AT'), ('Number', '.'), ('Type', 'String'), ('Description', 'Allele Traversal as path in graph')])
#    else:
#        new_header.add_record(i)
    new_header.add_record(i)
new_header.add_samples(get_new_samples(f_in))
f_out = pysam.VariantFile(sys.argv[-1], 'w', header=new_header)


for rec in f_in.fetch():
    new_samples_gt = {sample:[] for sample in new_header.samples}
    for sample, value_dict in rec.samples.items():
        if re.search(r'-(M|P).+$', sample):
            new_sample = '-'.join(sample.split('-')[:-1])
        else:
            new_sample = sample.split('_')[0]
        new_samples_gt[new_sample].append(value_dict['GT'][0])
        
    number = 0
    for sample in new_samples_gt:
        if new_samples_gt[sample] != [None, None]:
            number += 1
    new_info = dict(rec.info)
    
    new_info['AT'] = tuple(new_info['AT'][:len(new_info['AC']) + 1])
    new_info['AN'] = len(['1' for value_dict in rec.samples.values() if value_dict['GT'][0] != None])
    new_info['AF'] = tuple([i / new_info['AN'] for i in new_info['AC']])
    new_info['NS'] = number
    
    new_rec = new_header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=rec.alleles, id=rec.id, qual=rec.qual, filter=rec.filter, info=new_info)
    for sample in new_samples_gt:
        new_rec.samples[sample]['GT'] = tuple(new_samples_gt[sample])
        new_rec.samples[sample].phased = True
    f_out.write(new_rec)

f_in.close()
f_out.close()
