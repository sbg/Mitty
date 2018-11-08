import pysam
from numpy.random import RandomState


def assign_random_gt(input_vcf, output, sample_name="HG", default_af=0.01, seed=None):
    vcf_pointer = pysam.VariantFile(filename=input_vcf)
    new_header = vcf_pointer.header.copy()
    if "GT" not in new_header.formats:
        new_header.formats.add("GT", "1", "String", "Consensus Genotype across all datasets with called genotype")
        new_header.samples.add(sample_name)
    output.write(str(new_header))

    default_probs = [1 - default_af * (1 + default_af), default_af/2, default_af/2, default_af * default_af]
    rng = RandomState(seed)
    previous_pos = 0
    for rec in vcf_pointer.fetch():
        rec_copy = rec.copy()
        if "GT" not in rec_copy.format.keys():
            if rec_copy.pos == previous_pos:
                c = "0|0"
            else:
                if "AF" not in rec_copy.info.keys():
                    gt_probs = default_probs
                else:
                    af = rec_copy.info["AF"]
                    gt_probs = [1 - af * (1 + af), af/2, af/2, af * af]
                c = rng.choice(["0|0", "0|1", "1|0", "1|1"], p=gt_probs)
            output.write("\t".join([str(rec_copy)[:-1], "GT", c]) + "\n")
        previous_pos = rec_copy.pos

    vcf_pointer.close()
