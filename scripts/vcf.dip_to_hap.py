import sys
import pysam


if len(sys.argv) != 3:
    print(f'Usage: python {sys.argv[0]} input.vcf(.gz) output.vcf(.gz)')
    print('Please index the inputted vcf first')
    sys.exit(1)


def get_new_samples(vcf):
    new_samples = []
    rec = next(vcf.fetch())
    for sample, value_dict in rec.samples.items():
        if len(value_dict['GT']) == 1:
            new_samples.append(sample)
        else:
            for i, value in enumerate(value_dict['GT']):
                new_samples.append(f'{sample}.hap{i+1}')
    return new_samples


f_in = pysam.VariantFile(sys.argv[1])
new_header = pysam.VariantHeader()
for i in f_in.header.records:
    new_header.add_record(i)
new_header.add_samples(get_new_samples(f_in))
f_out = pysam.VariantFile(sys.argv[-1], 'w', header=new_header)

for rec in f_in.fetch():
    new_info = dict(rec.info)
    new_info['NS'] = rec.info['AN']
    new_rec = new_header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=rec.alleles, id=rec.id, qual=rec.qual, filter=rec.filter, info=new_info)

    for sample, value_dict in rec.samples.items():
        if len(value_dict['GT']) == 1:
            new_rec.samples[sample]['GT'] = value_dict['GT']
        else:
            for i, value in enumerate(value_dict['GT']):
                new_rec.samples[f'{sample}.hap{i+1}']['GT'] = (value, )
    
    f_out.write(new_rec)

f_in.close()
f_out.close()
