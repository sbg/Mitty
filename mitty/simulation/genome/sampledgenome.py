import pysam
from numpy.random import RandomState


def assign_random_gt(input_vcf, outname, sample_name="HG", default_af=0.01, seed=None):
    vcf_pointer = pysam.VariantFile(filename=input_vcf)
    new_header = vcf_pointer.header.copy()
    if "GT" not in new_header.formats:
        new_header.formats.add("GT", "1", "String", "Consensus Genotype across all datasets with called genotype")
        new_header.samples.add(sample_name)

    default_probs = [1 - default_af * (1 + default_af), default_af, default_af * default_af]
    rng = RandomState(seed)

    with open(outname, 'w') as out_vcf:
        out_vcf.write(str(new_header))
        for rec in vcf_pointer.fetch():
            rec_copy = rec.copy()
            if "GT" not in rec_copy.format.keys():
                if "AF" not in rec_copy.info.keys():
                    gt_probs = default_probs
                else:
                    af = rec_copy.info["AF"]
                    gt_probs = [1 - af * (1 + af), af, af * af]
                c = rng.choice(["0/0", "0/1", "1/1"], p=gt_probs)
                out_vcf.write("\t".join([str(rec_copy)[:-1], "GT", c]) + "\n")
    vcf_pointer.close()
