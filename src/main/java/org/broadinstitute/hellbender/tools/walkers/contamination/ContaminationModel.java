package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.List;

public class ContaminationModel {
    private boolean isInfiniteContaminant;  // the ContEst model of infinite contaminants
    private double eps;
    private double[] c;

    public ContaminationModel(final boolean isInfiniteContaminant, final double eps, final double[] c) {
        this.isInfiniteContaminant = isInfiniteContaminant;
        this.eps = eps;
        this.c = c;
    }

    /**
     * Array of likelihoods over the four sample genotypes
     * @param site
     * @param maf   the local minor allele fraction
     * @param eps   the base error rate
     * @param c     array of contaminant fractions eg [0.03, 0.4] means two contaminants with 3% and 4% contamination
     * @return a double[4], never null
     */
    private double[] genotypeLikelihoods(final PileupSummary site, final double maf, final double eps, final double[] c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();
        final int numContaminants = c.length;

        // u for "uncontaminated" part eg 1 - sum of contaminant fractions
        final double u = 1 - MathUtils.sum(c);

        // sample genotypes in order hom ref, alt minor, alt major, hom alt
        final double[] samplePriors = sampleGenotypePriors(f);
        final double[] sampleAFs = sampleGenotypeAlleleFractions(maf, eps);

        // contaminant genotypes in order hom ref, het, hom alt
        final double[] contaminantPriors = new double[] {(1-f)*(1-f), 2*f*(1-f), f*f};
        final double[] contaminantAFs = new double[] {eps, 0.5, 1 - eps};

        return new IndexRange(0, 4).mapToDouble(sg -> { // sg for "sample genotype"
            double sumOverContaminants = 0;
            if (isInfiniteContaminant) {
                sumOverContaminants = binom(k,n, u * sampleAFs[sg] + c[0] * f);
            } else {
                for (int[] contaminantGenotypes : new MultipleContaminantIterator(numContaminants)) {
                    final double weightedAF = u * sampleAFs[sg] + new IndexRange(0, numContaminants).sum(i -> c[i] * contaminantAFs[contaminantGenotypes[i]]);
                    final double contaminantPrior = new IndexRange(0, numContaminants).product(i -> contaminantPriors[contaminantGenotypes[i]]);
                    sumOverContaminants += contaminantPrior * binom(k, n, weightedAF);
                }
            }
            return samplePriors[sg] * sumOverContaminants;
        });

    }

    private static double binom(final int k, final int n, final double p) {
        return new BinomialDistribution(null, n, p).probability(k);
    }

    private double[] sampleGenotypeLikelihoods(PileupSummary site, double maf) {
        return genotypeLikelihoods(site, maf, eps, c);
    }

    public double probability(final PileupSummary site, final double maf, final SampleGenotype sampleGenotype) {
        final double[] likelihoods = sampleGenotypeLikelihoods(site, maf);
        return likelihoods[sampleGenotype.ordinal()] / MathUtils.sum(likelihoods);
    }

    public double logLikelihood(final List<PileupSummary> segment, final double maf) {
        return segment.stream().mapToDouble(site -> FastMath.log(MathUtils.sum(sampleGenotypeLikelihoods(site, maf)))).sum();
    }

    public double logLikelihood(final List<List<PileupSummary>> segments, final List<Double> mafs) {
        Utils.validate(segments.size() == mafs.size(), " Must have one MAF per segment");
        return new IndexRange(0, segments.size()).sum(n -> logLikelihood(segments.get(n), mafs.get(n)));
    }

    public enum SampleGenotype {
        HOM_REF, ALT_MINOR, ALT_MAJOR, HOM_ALT;
        public final static int NUMBER_OF_SAMPLE_GENOTYPES = values().length;
    }

    private static double[] sampleGenotypePriors(final double f) {
        final double[] result = new double[SampleGenotype.NUMBER_OF_SAMPLE_GENOTYPES];
        result[SampleGenotype.HOM_REF.ordinal()] = (1 - f) * (1 - f);
        result[SampleGenotype.ALT_MINOR.ordinal()] = f * (1 - f);
        result[SampleGenotype.ALT_MAJOR.ordinal()] = f * (1 - f);
        result[SampleGenotype.HOM_ALT.ordinal()] = f * f;
        return result;
    }

    private static double[] sampleGenotypeAlleleFractions(final double maf, final double errorRate) {
        final double[] result = new double[SampleGenotype.NUMBER_OF_SAMPLE_GENOTYPES];
        result[SampleGenotype.HOM_REF.ordinal()] = errorRate;
        result[SampleGenotype.ALT_MINOR.ordinal()] = maf;
        result[SampleGenotype.ALT_MAJOR.ordinal()] = 1 - maf;
        result[SampleGenotype.HOM_ALT.ordinal()] = 1 - errorRate;
        return result;
    }

    // iterate over all tuples of contaminant genotype ordinals.  For example, there are three contaminant genotypes,
    // hom ref, het, and hom alt, so iterating over two contaminants goes over
    // [0,0]; [0,1]; [0,2]; [1,0]; [1,1]; [1,2]; [2,0]; [2,1]; [2,2]
    private static class MultipleContaminantIterator implements Iterator<int[]>, Iterable<int[]> {
        private final int[] indices;
        private final int numContaminants;
        private int overallIndex = 0;
        private int maxOverallIndex;
        private static final int NUM_CONTAMINANT_GENOTYPES = 3;

        public MultipleContaminantIterator(final int numContaminants) {
            this.numContaminants = numContaminants;
            indices = new int[numContaminants];
            maxOverallIndex = 1;
            for (int n = 0; n < numContaminants; n++) {
                maxOverallIndex *= NUM_CONTAMINANT_GENOTYPES;
            }
            maxOverallIndex--;
        }

        @Override
        public boolean hasNext() { return overallIndex < maxOverallIndex; }

        @Override
        public int[] next() {
            overallIndex++;
            // try to increment the first index -- if it's at 2, zero it and try to increment the second, etc.
            for (int n = 0; n < numContaminants; n++) {
                if (++indices[n] < NUM_CONTAMINANT_GENOTYPES) {
                    return indices;
                } else {
                    indices[n] = 0;
                }
            }

            // if we got here, all the indices were at 2 and next() shouldn't have been called
            throw new GATKException.ShouldNeverReachHereException("Iteration shoul dbe done by now");
        }

        @Override
        public Iterator<int[]> iterator() { return this; }

    }
}
