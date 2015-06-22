package ec.vector;

import ec.*;
import ec.util.*;

/* 
 * FloatVectorSpecies.java
 * 
 * Created: Tue Feb 20 13:26:00 2001
 * By: Sean Luke
 */

/**
 * FloatVectorSpecies is a subclass of VectorSpecies with special
 * constraints for floating-point vectors, namely FloatVectorIndividual and
 * DoubleVectorIndividual.
 *
 * <p>FloatVectorSpecies can specify a number of parameters globally, per-segment, and per-gene.
 * See <a href="VectorSpecies.html">VectorSpecies</a> for information on how to this works.
 *
 * <p>FloatVectorSpecies defines a minimum and maximum gene value.  These values
 * are used during initialization and, depending on whether <tt>mutation-bounded</tt>
 * is true, also during various mutation algorithms to guarantee that the gene value
 * will not exceed these minimum and maximum bounds.
 *
 * <p>
 * FloatVectorSpecies provides support for five ways of mutating a gene.
 * <ul>
 * <li><b>reset</b> Replacing the gene's value with a value uniformly drawn from the gene's
 * range (the default behavior).</li>
 * <li><b>gauss</b>Perturbing the gene's value with gaussian noise; if the gene-by-gene range 
 * is used, than the standard deviation is scaled to reflect each gene's range. 
 * If the gaussian mutation's standard deviation is too large for the range,
 * than there's a large probability the mutated value will land outside range.
 * We will try again a number of times (100) before giving up and using the 
 * previous mutation method.</li>
 * <li><b>polynomial</b> Perturbing the gene's value with noise chosen from a <i>polynomial distribution</i>,
 * similar to the gaussian distribution.  The polynomial distribution was popularized
 * by Kalyanmoy Deb and is found in many of his publications (see http://www.iitk.ac.in/kangal/deb.shtml).
 * The polynomial distribution has two options.  First, there is the <i>index</i>.  This
 * variable defines the shape of the distribution and is in some sense the equivalent of the
 * standard deviation in the gaussian distribution.  The index is an integer.  If it is zero,
 * the polynomial distribution is simply the uniform distribution from [1,-1].  If it is 1, the
 * polynomial distribution is basically a triangular distribution from [1,-1] peaking at 0.  If
 * it is 2, the polynomial distribution follows a squared function, again peaking at 0.  Larger
 * values result in even more peaking and narrowness.  The default values used in nearly all of
 * the NSGA-II and Deb work is 20.  Second, there is whether or not the value is intended for
 * <i>bounded</i> genes.  The default polynomial distribution is used when we assume the gene can
 * take on literally any value, even beyond the min and max values.  For genes which are restricted
 * to be between min and max, there is an alternative version of the polynomial distribution, used by
 * Deb's team but not discussed much in the literature, desiged for that situation.  We assume boundedness
 * by default, and have found it to be somewhat better for NSGA-II and SPEA2 problems.  For a description
 * of this alternative version, see "A Niched-Penalty Approach for Constraint Handling in Genetic Algorithms"
 * by Kalyanmoy Deb and Samir Agrawal.  Deb's default implementation bounds the result to min or max;
 * instead ECJ's implementation of the polynomial distribution retries until it finds a legal value.  This
 * will be just fine for ranges like [0,1], but for smaller ranges you may be waiting a long time.
 * <li><b>integer-reset</b> Replacing the gene's value with a value uniformly drawn from the gene's range
 * but restricted to only integers.
 * <li><b>integer-random-walk</b> Replacing the gene's value by performing a random walk starting at the gene
 * value.  The random walk either adds 1 or subtracts 1 (chosen at random), then does a coin-flip
 * to see whether to continue the random walk.  When the coin-flip finally comes up false, the gene value
 * is set to the current random walk position.
 * </ul>
 *
 * <p>
 * FloatVectorSpecies provides support for two ways of initializing a gene.  The initialization procedure
 * is determined by the choice of mutation procedure as described above.  If the mutation is floating-point
 * (<tt>reset, gauss, polynomial</tt>), then initialization will be done by resetting the gene
 * to uniformly chosen floating-point value between the minimum and maximum legal gene values, inclusive.
 * If the mutation is integer (<tt>integer-reset, integer-random-walk</tt>), then initialization will be done
 * by performing the same kind of reset, but restricting values to integers only.
 * 
 * 
 * <p>
 * <b>Parameters</b><br>
 * <table>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>min-gene</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>min-gene</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>min-gene</tt>.<i>gene-number</i><br>
 <font size=-1>0.0 &lt;= double &lt;= 1.0 </font></td>
 <td valign=top>(probability that a gene will get mutated over default mutation)</td></tr>
 * <font size=-1>double (default=0.0)</font></td>
 * <td valign=top>(the minimum gene value)</td>
 * </tr>
 * 
 <tr><td>&nbsp;
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>max-gene</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>max-gene</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>max-gene</tt>.<i>gene-number</i><br>
 <font size=-1>0.0 &lt;= double &lt;= 1.0 </font></td>
 <td valign=top>(probability that a gene will get mutated over default mutation)</td></tr>
 * <font size=-1>double &gt;= <i>base</i>.min-gene</font></td>
 * <td valign=top>(the maximum gene value)</td>
 * </tr>
 * 
 <tr><td>&nbsp;
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-type</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>mutation-type</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-prob</tt>.<i>gene-number</i><br>
 * <font size=-1><tt>reset</tt>, <tt>gauss</tt>, <tt>polynomial</tt>, <tt>integer-reset</tt>, or <tt>integer-random-walk</tt> (default=<tt>reset</tt>)</font></td>
 * <td valign=top>(the mutation type)</td>
 * </tr>
 * 
 <tr><td>&nbsp;
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-stdev</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>mutation-stdev</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-stdev</tt>.<i>gene-number</i><br>
 * <font size=-1>double &ge; 0</font></td>
 * <td valign=top>(the standard deviation or the gauss perturbation)</td>
 * </tr>
 * 
 * <tr>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>distribution-index</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>distribution-index</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>distribution-index</tt>.<i>gene-number</i><br>
 * <font size=-1>int &ge; 0</font></td>
 * <td valign=top>(the mutation distribution index for the polynomial mutation distribution)</td>
 * </tr>
 * 
 <tr><td>&nbsp;
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>alternative-polynomial-version</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>alternative-polynomial-version</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>alternative-polynomial-version</tt>.<i>gene-number</i><br>
 *  <font size=-1>boolean (default=true)</font></td>
 *  <td valign=top>(whether to use the "bounded" variation of the polynomial mutation or the standard ("unbounded") version)</td>
 * </tr>
 *
 <tr><td>&nbsp;
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>random-walk-probability</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>random-walk-probability</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>random-walk-probability</tt>.<i>gene-number</i><br>
 <font size=-1>0.0 &lt;= double &lt;= 1.0 </font></td>
 *  <td valign=top>(the probability that a random walk will continue.  Random walks go up or down by 1.0 until the coin flip comes up false.)</td>
 * </tr>
 * 
 * <tr>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-bounded</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>segment</tt>.<i>segment-number</i>.<tt>mutation-bounded</tt>&nbsp;&nbsp;&nbsp;<i>or</i><br>
 <tr><td valign=top style="white-space: nowrap"><i>base</i>.<tt>mutation-bounded</tt>.<i>gene-number</i><br>
 *  <font size=-1>boolean (default=true)</font></td>
 *  <td valign=top>(whether mutation is restricted to only being within the min/max gene values.  Does not apply to SimulatedBinaryCrossover (which is always bounded))</td>
 * </tr>
 * 
 <tr><td>&nbsp;
 * <td valign=top><i>base</i>.<tt>out-of-bounds-retries</tt><br>
 *  <font size=-1>int &ge; 0 (default=100)</font></td>
 *  <td valign=top>(number of times the gaussian mutation got the gene out of range 
 *  before we give up and reset the gene's value; 0 means "never give up")</td>
 * </tr>
 *
 * </table>
 * @author Sean Luke, Gabriel Balan, Rafal Kicinger
 * @version 2.0
 */
 
public class FloatVectorSpecies extends VectorSpecies
    {
    public final static String P_MINGENE = "min-gene";
    public final static String P_MAXGENE = "max-gene";
    public final static String P_MUTATIONTYPE = "mutation-type";
    public final static String P_STDEV = "mutation-stdev";
    public final static String P_MUTATION_DISTRIBUTION_INDEX = "mutation-distribution-index";
    public final static String P_POLYNOMIAL_ALTERNATIVE = "alternative-polynomial-version";
    public final static String P_RANDOM_WALK_PROBABILITY = "random-walk-probability";
    public final static String P_OUTOFBOUNDS_RETRIES = "out-of-bounds-retries";
    public final static String P_MUTATION_BOUNDED = "mutation-bounded";

    public final static String V_RESET_MUTATION = "reset";
    public final static String V_GAUSS_MUTATION = "gauss";
    public final static String V_POLYNOMIAL_MUTATION = "polynomial";
    public final static String V_INTEGER_RANDOM_WALK_MUTATION = "integer-random-walk";
    public final static String V_INTEGER_RESET_MUTATION = "integer-reset";
    
    // Set to true when setup() is called.  This is only used for the representation invariant.
    private boolean isSetup;
    
    public static final int DEFAULT_OUT_OF_BOUNDS_RETRIES = 100;
                
    static final double SIMULATED_BINARY_CROSSOVER_EPS = 1.0e-14;   

    @Override
    public DoubleVectorMutator mutator(final int gene) {
        assert(gene >= 0);
        assert(gene < mutators.length);
        return (DoubleVectorMutator) mutators[gene];
    }
    
    public double maxGene(int gene)
        {
        final VectorMutator[] m = mutators;
        if (m.length <= gene)
            gene = m.length - 1;
        return ((DoubleVectorMutator)m[gene]).maxGene();
        }

    public double minGene(int gene)
        {
        final VectorMutator[] m = mutators;
        if (m.length <= gene)
            gene = m.length - 1;
        return ((DoubleVectorMutator)m[gene]).minGene();
        }


    private boolean inNumericalTypeRange(double geneVal)
        {
        if (i_prototype instanceof FloatVectorIndividual)
            return (geneVal <= Float.MAX_VALUE && geneVal >= -Float.MAX_VALUE);
        else if (i_prototype instanceof DoubleVectorIndividual)
            return true; // geneVal is valid for all double
        else
            return false; // dunno what the individual is...
        }


    public void setup(final EvolutionState state, final Parameter base)
        {
        Parameter def = defaultBase();
        
        setupGenome(state, base);
        
        mutators = new VectorMutator[genomeSize];
        // Set the global mutation parameters
        globalMutator = getMutatorForType(null, state, base, def, "");
        
        // The call to super will populate the mutator array for per-gene
        // and per-segment parameters via calls to loadParametersForGene().
        super.setup(state, base);
        }
    
    @Override
    protected void loadParametersForGene(final EvolutionState state, final int index, final Parameter base, final Parameter def, final String postfix)
        {       
        super.loadParametersForGene(state, index, base, def, postfix);
        // The mutator may already have been set up with segment-level parameters.  If so,
        // inherit the segment parameters.  If not, inherit the global parameters.
        final VectorMutator parent = (mutators[index] == null) ? globalMutator : mutators[index];
        mutators[index] = getMutatorForType(parent, state, base, def, postfix);
        }
    
    /** Load global, per-segment, or per-gene parameters for a mutation operator. 
     * @return A DoubleVectorMutator that encodes a mutation operator. */
    protected DoubleVectorMutator getMutatorForType(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix)
    {            
        String mutationType = state.parameters.getStringWithDefault(base.push(P_MUTATIONTYPE).push(postfix), def.push(P_MUTATIONTYPE).push(postfix), (parent == null) ? null : parent.mutationType());
        if (mutationType == null)
        {
            mutationType = V_RESET_MUTATION;
            state.output.warning(String.format("No mutation type given for %s, assuming '%s' mutation", this.getClass().getSimpleName(), V_RESET_MUTATION),
                base.push(P_MUTATIONTYPE).push(postfix), def.push(P_MUTATIONTYPE).push(postfix));
        }
        if (mutationType.equals(V_RESET_MUTATION))
            return new ResetMutator(parent, state, base, def, postfix);
        else if (mutationType.equals(V_POLYNOMIAL_MUTATION))
            return new PolynomialMutator(parent, state, base, def, postfix);
        else if (mutationType.equals(V_GAUSS_MUTATION))
            return new GaussianMutator(parent, state, base, def, postfix);
        else if (mutationType.equals(V_INTEGER_RESET_MUTATION))
        {
            state.output.warnOnce(String.format("Integer Reset Mutation used in %s.  Be advised that during initialization these genes will only be set to integer values.", this.getClass().getSimpleName()));
            return new IntegerResetMutator(parent, state, base, def, postfix);
        }
        else if (mutationType.equals(V_INTEGER_RANDOM_WALK_MUTATION))
        {
            state.output.warnOnce(String.format("Integer Random Walk Mutation used in %s.  Be advised that during initialization these genes will only be set to integer values.", this.getClass().getSimpleName()));
            return new IntegerRandomWalkMutator(parent, state, base, def, postfix);
        }
        else {
            state.output.fatal(String.format("%s given a bad mutation type: %s", this.getClass().getSimpleName(), mutationType), base.push(P_MUTATIONTYPE).push(postfix), def.push(P_MUTATIONTYPE).push(postfix));
            return null;
        }
    }
    
    /** Set the value of a gene from an Individual that is either a 
     * FloatVectorIndividual or a DoubleVectorIndividual.  This hack allows us to
     * write general code that works with both individual types.
     */
    private static void setGeneValue(final Individual ind, final int index, double value) {
        assert(ind != null);
        assert(!Double.isNaN(value));
        if (ind instanceof FloatVectorIndividual)
            ((FloatVectorIndividual)ind).genome[index] = (float) value;
        else
            ((DoubleVectorIndividual)ind).genome[index] = value;
    }
    
    /** Retrieve the value of a gene from an Individual that is either a 
     * FloatVectorIndividual or a DoubleVectorIndividual. This hack allows us to
     * write general code that works with both individual types.
     */
    private static double getGeneValue(final Individual ind, final int index) {
        assert(ind != null);
        if (ind instanceof FloatVectorIndividual)
            return ((FloatVectorIndividual)ind).genome[index];
        else
            return ((DoubleVectorIndividual)ind).genome[index];
    }

    /** Representation invariant.  This is used for asserts and unit tests,
     * to ensure that the class is always in a valid state.
     * 
     * @return Always true.  If this is ever false, it means that a fault has caused the class to be in an invalid state.
     */
    public final boolean repOK() {
        return V_GAUSS_MUTATION != null
                && !V_GAUSS_MUTATION.isEmpty()
                && V_GAUSS_MUTATION.equals(V_GAUSS_MUTATION.toLowerCase())
                && V_INTEGER_RANDOM_WALK_MUTATION != null
                && !V_INTEGER_RANDOM_WALK_MUTATION.isEmpty()
                && V_INTEGER_RESET_MUTATION != null
                && !V_INTEGER_RESET_MUTATION.isEmpty()
                && V_POLYNOMIAL_MUTATION != null
                && !V_POLYNOMIAL_MUTATION.isEmpty()
                && V_RESET_MUTATION != null
                && !V_RESET_MUTATION.isEmpty()
                && DEFAULT_OUT_OF_BOUNDS_RETRIES >= 0
                && !isSetup
                || (mutators != null
                && mutators.length == genomeSize
                && !containsNulls(mutators));
    }
    
    private static boolean containsNulls(final Object[] array) {
        for (final Object o : array)
            if (o == null)
                return true;
        return false;
    }
    
    /** This class stores the parameters and algorithm for a specific mutation operator
     * that operates on FloatVectorIndividuals and DoubleVectorIndividuals. */
    protected abstract static class DoubleVectorMutator extends VectorMutator {
        final protected double minGene;
        final protected double maxGene;
        final protected boolean mutationIsBounded;
        
        public DoubleVectorMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
            
            // If the parent is of this type, inherit its parameters
            final double defaultMin = (parent instanceof DoubleVectorMutator) ? ((DoubleVectorMutator)parent).minGene : 0;
            final double defaultMax = (parent instanceof DoubleVectorMutator) ? ((DoubleVectorMutator)parent).maxGene : defaultMin;
            final boolean defaultBounded = (parent instanceof DoubleVectorMutator) ? ((DoubleVectorMutator)parent).mutationIsBounded : true;
            
            minGene = state.parameters.getDoubleWithDefault(base.push(P_MINGENE).push(postfix), def.push(P_MINGENE).push(postfix), defaultMin);
            maxGene = state.parameters.getDoubleWithDefault(base.push(P_MAXGENE).push(postfix), def.push(P_MAXGENE).push(postfix), defaultMax);
            if (maxGene < minGene)
                state.output.fatal(String.format("%s must have a %s which is <= the %s", FloatVectorSpecies.class.getSimpleName(), P_MINGENE, P_MAXGENE),
                        base.push(P_MAXGENE).push(postfix), def.push(P_MAXGENE).push(postfix));
            if (!state.parameters.exists(base.push(P_MUTATION_BOUNDED).push(postfix), def.push(P_MUTATION_BOUNDED).push(postfix)))
                state.output.warning(String.format("%s is using gaussian, polynomial, or integer random walk mutation as its global mutation type, but '%s' is not defined.  Assuming 'true'", FloatVectorSpecies.class.getSimpleName(), P_MUTATION_BOUNDED));
            mutationIsBounded = state.parameters.getBoolean(base.push(P_MUTATION_BOUNDED).push(postfix), def.push(P_MUTATION_BOUNDED).push(postfix), defaultBounded);
        }
        
        /** @return True if the gene that this operator manipulates can only assume integer values. */
        public abstract boolean isIntegerType();

        public double minGene() { return minGene; }

        public double maxGene() { return maxGene; }

        public boolean mutationIsBounded() { return mutationIsBounded; }
    }
        
    /** Stores the parameters and algorithm for the integer-random-walk mutation method. */
    protected static final class IntegerRandomWalkMutator extends DoubleVectorMutator {
        public static final double MAXIMUM_INTEGER_IN_DOUBLE = 9.007199254740992E15;
        
        final private double randomWalkProbability;
        
        public double randomWalkProbability() { return randomWalkProbability; }
        
        public IntegerRandomWalkMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
            // If the parent is of this type, inherit its parameters
            if (parent instanceof IntegerRandomWalkMutator)
                randomWalkProbability = state.parameters.getDoubleWithDefault(base.push(P_RANDOM_WALK_PROBABILITY).push(postfix), def.push(P_RANDOM_WALK_PROBABILITY).push(postfix), ((IntegerRandomWalkMutator)parent).randomWalkProbability);
            else
                randomWalkProbability = state.parameters.getDouble(base.push(P_RANDOM_WALK_PROBABILITY).push(postfix), def.push(P_RANDOM_WALK_PROBABILITY).push(postfix));
            if (randomWalkProbability < 0 || randomWalkProbability > 1.0)
                state.output.fatal(String.format("'%s' is set to %f for '%s' mutation, but must be a probability between 0 and 1.", P_RANDOM_WALK_PROBABILITY, randomWalkProbability, V_INTEGER_RANDOM_WALK_MUTATION),
                        base.push(P_RANDOM_WALK_PROBABILITY).push(postfix), def.push(P_RANDOM_WALK_PROBABILITY).push(postfix));
            assert(repOK());
        }
        
        @Override
        public void mutate(final EvolutionState state, final Individual individual, final MersenneTwisterFast random, final int index) {
            assert(state != null);
            assert(individual != null);
            assert(random != null);
            assert(index > 0);
            final double min = mutationIsBounded() ? minGene() : -MAXIMUM_INTEGER_IN_DOUBLE;
            final double max = mutationIsBounded() ? maxGene() : MAXIMUM_INTEGER_IN_DOUBLE;
            int g = (int) getGeneValue(individual, index);
            do
            {
                int n = (int)(random.nextBoolean() ? 1 : -1);
                if ((n == 1 && g < max) ||
                    (n == -1 && g > min))
                    g += n;
                else if ((n == -1 && g < max) ||
                    (n == 1 && g > min))
                    g -= n;
                else
                    state.output.fatal(String.format("Illegal state reached in %s.mutate().", IntegerRandomWalkMutator.class.getSimpleName()));
            }
            while (random.nextBoolean(randomWalkProbability));
            setGeneValue(individual, index, g);
        }

        @Override
        public String mutationType() {
            return V_INTEGER_RANDOM_WALK_MUTATION;
        }
        
        public final boolean repOK() {
            return minGene >= 0.0
                   && maxGene >= minGene
                    && !Double.isNaN(randomWalkProbability)
                    && randomWalkProbability >= 0.0
                    && randomWalkProbability <= 1.0;
        }

        @Override
        public boolean isIntegerType() { return true; }
    }
    
    /** Stores the parameters and algorithm for the float-reset mutation method. */
    protected static final class ResetMutator extends DoubleVectorMutator {
        
        public ResetMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
            assert(repOK());
        }
        
        @Override
        public void mutate(final EvolutionState state, final Individual individual, final MersenneTwisterFast random, final int index) {
            assert(state != null);
            assert(individual != null);
            assert(random != null);
            assert(index > 0);
            final double val = (float)(minGene + random.nextFloat(true, true) * (maxGene - minGene));
            setGeneValue(individual, index, val);
        }

        @Override
        public String mutationType() {
            return V_RESET_MUTATION;
        }
        
        public final boolean repOK() {
            return minGene >= 0.0
                   && maxGene >= minGene;
        }

        @Override
        public boolean isIntegerType() { return false; }
    }
    
    /** Stores the parameters and algorithm for the Gaussian mutation methods. */
    protected static final class GaussianMutator extends DoubleVectorMutator {
        
        final private double stdev;
        final private int outOfBoundsRetries;
        
        public double stdev() { return stdev; }
        public int outOfBoundsRetries() { return outOfBoundsRetries; }
        
        public GaussianMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
            // If the parent is of this type, inherit its parameters
            final double defaultStdev = (parent instanceof GaussianMutator) ? ((GaussianMutator)parent).stdev : 0;
            final int defaultOutOfBoundsRetries = (parent instanceof GaussianMutator) ? ((GaussianMutator)parent).outOfBoundsRetries : DEFAULT_OUT_OF_BOUNDS_RETRIES;
            stdev = state.parameters.getDoubleWithDefault(base.push(P_STDEV).push(postfix),def.push(P_STDEV).push(postfix), defaultStdev);
            outOfBoundsRetries = state.parameters.getIntWithDefault(base.push(P_OUTOFBOUNDS_RETRIES).push(postfix), def.push(P_OUTOFBOUNDS_RETRIES).push(postfix), defaultOutOfBoundsRetries);
            if(outOfBoundsRetries<0)
                state.output.fatal("Out of bounds retries must be >= 0.", base.push(P_OUTOFBOUNDS_RETRIES).push(postfix), def.push(P_OUTOFBOUNDS_RETRIES).push(postfix));
            assert(repOK());
        }
        
        @Override
        public void mutate(final EvolutionState state, final Individual individual, final MersenneTwisterFast random, final int index) {
            assert(state != null);
            assert(individual != null);
            assert(random != null);
            assert(index > 0);
            final double oldVal = getGeneValue(individual, index);
            double val;
            int outOfBoundsLeftOverTries = outOfBoundsRetries;
            boolean givingUpAllowed = (outOfBoundsRetries != 0);
            do
            {
                val = random.nextGaussian() * stdev + oldVal;
                outOfBoundsLeftOverTries--;
                if (mutationIsBounded() && (val > maxGene || val < minGene))
                {
                    if (givingUpAllowed && (outOfBoundsLeftOverTries == 0))
                    {
                        val = minGene + random.nextDouble() * (maxGene - minGene);
                        state.output.warnOnce("The limit of 'out-of-range' retries for gaussian or polynomial mutation (" + outOfBoundsRetries + ") was reached.");
                        break;
                    }
                } 
                else break;
            }
            while (true);
            setGeneValue(individual, index, val);
        }

        @Override
        public String mutationType() {
            return V_GAUSS_MUTATION;
        }
        
        public final boolean repOK() {
            return minGene >= 0.0
                   && maxGene >= minGene
                    && !Double.isNaN(stdev)
                    && stdev >= 0.0
                    && outOfBoundsRetries >= 0;
        }

        @Override
        public boolean isIntegerType() { return false; }
    }
    
    /** Stores the parameters and algorithm for the polynomial mutation method. */
    protected static final class PolynomialMutator extends DoubleVectorMutator {
        
        final private int mutationDistributionIndex;
        final private boolean polynomialIsAlternative;
        final private int outOfBoundsRetries;
        
        public double mutationDistributionIndex() { return mutationDistributionIndex; }
        public int outOfBoundsRetries() { return outOfBoundsRetries; }
        public boolean polynomialIsAlternative() { return polynomialIsAlternative; }
        
        public PolynomialMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
            // If the parent is of this type, inherit its parameters
            if (parent instanceof PolynomialMutator) {
                final PolynomialMutator pParent = (PolynomialMutator)parent;
                outOfBoundsRetries = state.parameters.getIntWithDefault(base.push(P_OUTOFBOUNDS_RETRIES).push(postfix), def.push(P_OUTOFBOUNDS_RETRIES).push(postfix), pParent.outOfBoundsRetries);
                mutationDistributionIndex = state.parameters.getIntWithDefault(base.push(P_MUTATION_DISTRIBUTION_INDEX).push(postfix), def.push(P_MUTATION_DISTRIBUTION_INDEX).push(postfix), pParent.mutationDistributionIndex);
                polynomialIsAlternative = state.parameters.getBoolean(base.push(P_POLYNOMIAL_ALTERNATIVE).push(postfix), def.push(P_POLYNOMIAL_ALTERNATIVE).push(postfix), pParent.polynomialIsAlternative);
            }
            else {
                outOfBoundsRetries = state.parameters.getIntWithDefault(base.push(P_OUTOFBOUNDS_RETRIES), def.push(P_OUTOFBOUNDS_RETRIES), DEFAULT_OUT_OF_BOUNDS_RETRIES);
                mutationDistributionIndex = state.parameters.getInt(base.push(P_MUTATION_DISTRIBUTION_INDEX), def.push(P_MUTATION_DISTRIBUTION_INDEX), 0);
                if (!state.parameters.exists(base.push(P_POLYNOMIAL_ALTERNATIVE), def.push(P_POLYNOMIAL_ALTERNATIVE)))
                    state.output.warning(String.format("%s is using polynomial mutation as its global mutation type, but '%s' is not defined.  Assuming 'true'", FloatVectorSpecies.class.getSimpleName(), P_POLYNOMIAL_ALTERNATIVE));
                polynomialIsAlternative = state.parameters.getBoolean(base.push(P_POLYNOMIAL_ALTERNATIVE), def.push(P_POLYNOMIAL_ALTERNATIVE), true);
            }
            if(outOfBoundsRetries<0)
                state.output.fatal("Out of bounds retries must be >= 0.", base.push(P_OUTOFBOUNDS_RETRIES).push(postfix), def.push(P_OUTOFBOUNDS_RETRIES).push(postfix));
            if (mutationDistributionIndex < 0)
                state.output.fatal(String.format("If %s is going to use polynomial mutation as its global mutation type, the global distribution index must be defined and >= 0.", FloatVectorSpecies.class.getSimpleName()),
                    base.push(P_MUTATION_DISTRIBUTION_INDEX).push(postfix), def.push(P_MUTATION_DISTRIBUTION_INDEX).push(postfix));
            
            assert(repOK());
        }
        
        @Override
        public void mutate(final EvolutionState state, final Individual individual, final MersenneTwisterFast random, final int index) {
            assert(state != null);
            assert(individual != null);
            assert(random != null);
            assert(index > 0);
            final double oldVal = getGeneValue(individual, index);

            double rnd, delta1, delta2, mut_pow, deltaq;
            double y, yl, yu, val, xy;
            double y1;

            y1 = y = oldVal;  // ind[index];
            yl = minGene; // min_realvar[index];
            yu = maxGene; // max_realvar[index];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);

            int totalTries = outOfBoundsRetries;
            int tries;
            for(tries = 0; tries < totalTries || totalTries == 0; tries++)  // keep trying until totalTries is reached if it's not zero.  If it's zero, go on forever.
                {
                rnd = random.nextFloat();
                mut_pow = 1.0/(mutationDistributionIndex+1.0);
                if (rnd <= 0.5)
                    {
                    xy = 1.0-delta1;
                    val = 2.0*rnd + (polynomialIsAlternative ? (1.0-2.0*rnd)*(Math.pow(xy,(mutationDistributionIndex+1.0))) : 0.0);
                    deltaq =  Math.pow(val,mut_pow) - 1.0;
                    }
                else
                    {
                    xy = 1.0-delta2;
                    val = 2.0*(1.0-rnd) + (polynomialIsAlternative ? 2.0*(rnd-0.5)*(Math.pow(xy,(mutationDistributionIndex+1.0))) : 0.0);
                    deltaq = 1.0 - (Math.pow(val,mut_pow));
                    }
                y1 = y + deltaq*(yu-yl);
                if (!mutationIsBounded || (y1 >= yl && y1 <= yu)) break;  // yay, found one
                }

            // at this point, if tries is totalTries, we failed
            if (totalTries != 0 && tries == totalTries)
                {
                // just randomize
                y1 = (minGene + random.nextFloat(true, true) * (maxGene - minGene));
                state.output.warnOnce("The limit of 'out-of-range' retries for gaussian or polynomial mutation (" + outOfBoundsRetries + ") was reached.");
                }
            setGeneValue(individual, index, y1);
        }

        @Override
        public String mutationType() {
            return V_POLYNOMIAL_MUTATION;
        }
        
        public final boolean repOK() {
            return minGene >= 0.0
                   && maxGene >= minGene
                    && mutationDistributionIndex >= 0
                    && outOfBoundsRetries >= 0;
        }

        @Override
        public boolean isIntegerType() { return false; }
    }
    
    /** Stores the parameters and algorithm for the integer-reset mutation method. */
    protected static final class IntegerResetMutator extends DoubleVectorMutator {
        
        public IntegerResetMutator(final VectorMutator parent, final EvolutionState state, final Parameter base, final Parameter def, final String postfix) {
            super(parent, state, base, def, postfix);
            assert(state != null);
            assert(base != null);
             assert(repOK());
        }
        
        @Override
        public void mutate(final EvolutionState state, final Individual individual, final MersenneTwisterFast random, final int index) {
            assert(state != null);
            assert(individual != null);
            assert(random != null);
            assert(index > 0);
            final int val = randomValueFromClosedInterval((int) minGene, (int) maxGene, random);
            setGeneValue(individual, index, val);
        }
        
        private int randomValueFromClosedInterval(final int min, final int max, final MersenneTwisterFast random)
        {
            if (max - min < 0) // we had an overflow
            {
                int l = 0;
                do l = random.nextInt();
                while(l < min || l > max);
                return l;
            }
            else return min + random.nextInt(max - min + 1);
        }

        @Override
        public String mutationType() {
            return V_INTEGER_RESET_MUTATION;
        }
        
        public final boolean repOK() {
            return minGene >= 0.0
                   && maxGene >= minGene;
        }

        @Override
        public boolean isIntegerType() { return true; }
    }
    
    }
