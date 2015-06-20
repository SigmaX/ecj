package ec.vector;

import ec.EvolutionState;
import ec.Evolve;
import ec.Population;
import ec.Subpopulation;
import ec.select.FirstSelection;
import ec.simple.SimpleEvolutionState;
import ec.simple.SimpleFitness;
import ec.util.Parameter;
import ec.util.ParameterDatabase;
import ec.vector.breed.VectorMutationPipeline;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Unit tests to verify the mutation-stdev-scale-to-range parameter of FloatVectorSpecies.
 * 
 * @author Eric O. Scott <escott8@gmu.edu>
 */
public class FloatVectorSpeciesTest {
    private final static Parameter BASE = new Parameter("pop.subpop.0.species");
    private final static int NUM_GENES = 3;
    private final static double GLOBAL_STDEV = 0.5;
    private final static double DEFAULT_MIN = 0;
    private final static double DEFAULT_MAX = 100;
    
    private EvolutionState state;
    
    public FloatVectorSpeciesTest() {
    }
    
    @Before
    public void setUp() {
        state = getFreshState();
    }
    
    private static EvolutionState getFreshState() {
        final EvolutionState state = new SimpleEvolutionState();
        // Set up just the parameters needed for the SUT to initialize itself
        state.parameters = getParams();
        
        // We need errors to throw exceptions (rather than exit the program) so we can verify them
        state.output = Evolve.buildOutput();
        state.output.setThrowsErrors(true);
        
        // Set up an empty population
        state.population = new Population();
        state.population.subpops = new Subpopulation[] { new Subpopulation() };
        return state;
    }
    
    private static ParameterDatabase getParams() {
        final ParameterDatabase parameters = new ParameterDatabase();
        parameters.set(BASE, FloatVectorSpecies.class.getSimpleName());
        parameters.set(BASE.push(FloatVectorSpecies.P_INDIVIDUAL), DoubleVectorIndividual.class.getCanonicalName());
        parameters.set(BASE.push(FloatVectorSpecies.P_GENOMESIZE), ""+NUM_GENES);
        parameters.set(BASE.push(FloatVectorSpecies.P_MINGENE), ""+DEFAULT_MIN);
        parameters.set(BASE.push(FloatVectorSpecies.P_MAXGENE), ""+DEFAULT_MAX);
        parameters.set(BASE.push(FloatVectorSpecies.P_MUTATIONPROB), "0.5");
        parameters.set(BASE.push(FloatVectorSpecies.P_MUTATIONTYPE), "gauss");
        parameters.set(BASE.push(FloatVectorSpecies.P_STDEV), ""+GLOBAL_STDEV);
        parameters.set(BASE.push(FloatVectorSpecies.P_PIPE), VectorMutationPipeline.class.getCanonicalName());
        parameters.set(BASE.push(FloatVectorSpecies.P_PIPE).push(VectorMutationPipeline.P_SOURCE).push("0"), FirstSelection.class.getCanonicalName());
        parameters.set(BASE.push(FloatVectorSpecies.P_FITNESS), SimpleFitness.class.getCanonicalName());
        
        return parameters;
    }

    @Test
    public void testGlobalParameters1() {
        System.out.println("setup (global parameters)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        instance.setup(state, BASE);
        
        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(1), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(1), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(2), 0.000001);
        // Mutation bounding defaults to true
        assertTrue(instance.mutationIsBounded(0));
        assertTrue(instance.mutationIsBounded(1));
        assertTrue(instance.mutationIsBounded(2));
    }

    @Test
    public void testGlobalParameters2() {
        System.out.println("setup (global parameters, bounded)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MUTATION_BOUNDED), "false");
        instance.setup(state, BASE);
        
        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(1), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(1), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(2), 0.000001);
        assertFalse(instance.mutationIsBounded(0));
        assertFalse(instance.mutationIsBounded(1));
        assertFalse(instance.mutationIsBounded(2));
    }
    
    @Test
    public void testGeneParameters1() {
        System.out.println("setup (gene parameters, initialization bounds)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MINGENE).push("1"), ""+min1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MAXGENE).push("1"), ""+max1);
        instance.setup(state, BASE);

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(min1, instance.minGene(1), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(max1, instance.maxGene(1), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(2), 0.000001);
        assertTrue(instance.mutationIsBounded(0));
        assertTrue(instance.mutationIsBounded(1));
        assertTrue(instance.mutationIsBounded(2));
    }
    
    @Test
    public void testGeneParameters2() {
        System.out.println("setup (gene parameters, mutation bounds)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MINGENE).push("1"), ""+min1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MAXGENE).push("1"), ""+max1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_MUTATION_BOUNDED).push("1"), "false");
        instance.setup(state, BASE);

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(min1, instance.minGene(1), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(max1, instance.maxGene(1), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(2), 0.000001);
        assertTrue(instance.mutationIsBounded(0));
        assertFalse(instance.mutationIsBounded(1));
        assertTrue(instance.mutationIsBounded(2));
    }
    
    @Test
    public void testSegmentParameters1() {
        System.out.println("setup (segment parameters, initialization bounds)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        state.parameters.set(BASE.push(FloatVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT_TYPE), FloatVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("0").push(FloatVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MAXGENE), ""+max1);
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(min1, instance.minGene(1), 0.000001);
        assertEquals(min1, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(max1, instance.maxGene(1), 0.000001);
        assertEquals(max1, instance.maxGene(2), 0.000001);
        assertTrue(instance.mutationIsBounded(0));
        assertTrue(instance.mutationIsBounded(1));
        assertTrue(instance.mutationIsBounded(2));
    }
    
    @Test
    public void testSegmentParameters2() {
        System.out.println("setup (segment parameters, mutation bounds, mutation type)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        state.parameters.set(BASE.push(FloatVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT_TYPE), FloatVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("0").push(FloatVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        // Set segment-specific initialization and mutation bounds
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MAXGENE), ""+max1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MUTATION_BOUNDED), "false");
        // Set a segment-specific mutation operator
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MUTATIONTYPE), FloatVectorSpecies.V_INTEGER_RANDOM_WALK_MUTATION);
        final double rProb = 0.5;
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_RANDOM_WALK_PROBABILITY), ""+rProb);
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(min1, instance.minGene(1), 0.000001);
        assertEquals(min1, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(max1, instance.maxGene(1), 0.000001);
        assertEquals(max1, instance.maxGene(2), 0.000001);
        assertEquals(Double.NaN, instance.randomWalkProbability(0), 0);
        assertEquals(rProb, instance.randomWalkProbability(1), 0.000001);
        assertEquals(rProb, instance.randomWalkProbability(2), 0.000001);
        assertTrue(instance.mutationIsBounded(0));
        assertFalse(instance.mutationIsBounded(1));
        assertFalse(instance.mutationIsBounded(2));
    }
    
    @Test
    public void testSegmentParameters3() {
        System.out.println("setup (segment parameters, mutation bounds, no mutation type)");
        final FloatVectorSpecies instance = new FloatVectorSpecies();
        state.parameters.set(BASE.push(FloatVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT_TYPE), FloatVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("0").push(FloatVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MAXGENE), ""+max1);
        state.parameters.set(BASE.push(FloatVectorSpecies.P_SEGMENT).push("1").push(FloatVectorSpecies.P_MUTATION_BOUNDED), "false");
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(0), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(1), 0.000001);
        assertEquals(GLOBAL_STDEV, instance.gaussMutationStdev(2), 0.000001);
        assertEquals(DEFAULT_MIN, instance.minGene(0), 0.000001);
        assertEquals(min1, instance.minGene(1), 0.000001);
        assertEquals(min1, instance.minGene(2), 0.000001);
        assertEquals(DEFAULT_MAX, instance.maxGene(0), 0.000001);
        assertEquals(max1, instance.maxGene(1), 0.000001);
        assertEquals(max1, instance.maxGene(2), 0.000001);
        assertTrue(instance.mutationIsBounded(0));
        assertFalse(instance.mutationIsBounded(1));
        assertFalse(instance.mutationIsBounded(2));
    }
}
