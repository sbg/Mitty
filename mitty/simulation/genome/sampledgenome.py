import pysam
from numpy.random import RandomState


def assign_random_gt(input_vcf, output, sample_name="HG", default_af=0.01, seed=None):
    vcf_pointer = pysam.VariantFile(filename=input_vcf)
    new_header = vcf_pointer.header.copy()
    if "GT" not in new_header.formats:
        new_header.formats.add("GT", "1", "String", "Consensus Genotype across all datasets with called genotype")
        new_header.samples.add(sample_name)
    output.write(str(new_header))

    gt = ["0|0", "0|1", "1|0", "1|1"]
    default_probs_all_gt = [1 - default_af * (1 + default_af), default_af/2, default_af/2, default_af * default_af]
    default_probs_no_hom_alt = [1 - default_af, default_af]
    rng = RandomState(seed)
    prev_locus = -1, -1
    for rec in vcf_pointer.fetch():
        r_copy = rec.copy()
        if "GT" not in r_copy.format.keys():
            locus = r_copy.pos
            ok_alts = (locus > prev_locus[0]), (locus > prev_locus[1])
            if ok_alts[0] or ok_alts[1]:
                hom_ok = ok_alts[0] and ok_alts[1]
                if "AF" not in r_copy.info.keys():
                    gt_p = default_probs_all_gt if hom_ok else default_probs_no_hom_alt
                else:
                    af = r_copy.info["AF"]
                    gt_p = [1 - af * (1 + af), af/2, af/2, af * af] if hom_ok else [1 - af, af]
                c = rng.choice(gt if hom_ok else (["0|.", "1|."] if ok_alts[0] else [".|0", ".|1"]), p=gt_p)
                prev_locus = (locus + r_copy.rlen - 1) if c[0] == '1' else locus, \
                             (locus + r_copy.rlen - 1) if c[2] == '1' else locus
            else:
                c = ".|."
                prev_locus = locus, locus

            output.write("\t".join([str(r_copy)[:-1], "GT", c]) + "\n")

    vcf_pointer.close()
