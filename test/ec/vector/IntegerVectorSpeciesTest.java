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
 *
 * @author Eric O. Scott
 */
public class IntegerVectorSpeciesTest {
    private final static Parameter BASE = new Parameter("pop.subpop.0.species");
    private final static int NUM_GENES = 3;
    private final static double GLOBAL_WALK_PROB = 0.5;
    private final static double DEFAULT_MIN = 0;
    private final static double DEFAULT_MAX = 100;
    
    private EvolutionState state;
    
    public IntegerVectorSpeciesTest() {
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
        parameters.set(BASE, IntegerVectorSpecies.class.getSimpleName());
        parameters.set(BASE.push(IntegerVectorSpecies.P_INDIVIDUAL), IntegerVectorIndividual.class.getCanonicalName());
        parameters.set(BASE.push(IntegerVectorSpecies.P_GENOMESIZE), ""+NUM_GENES);
        parameters.set(BASE.push(IntegerVectorSpecies.P_MINGENE), ""+DEFAULT_MIN);
        parameters.set(BASE.push(IntegerVectorSpecies.P_MAXGENE), ""+DEFAULT_MAX);
        parameters.set(BASE.push(IntegerVectorSpecies.P_MUTATIONPROB), "0.5");
        parameters.set(BASE.push(IntegerVectorSpecies.P_MUTATIONTYPE), IntegerVectorSpecies.V_RANDOM_WALK_MUTATION);
        parameters.set(BASE.push(IntegerVectorSpecies.P_RANDOM_WALK_PROBABILITY), ""+GLOBAL_WALK_PROB);
        parameters.set(BASE.push(IntegerVectorSpecies.P_PIPE), VectorMutationPipeline.class.getCanonicalName());
        parameters.set(BASE.push(IntegerVectorSpecies.P_PIPE).push(VectorMutationPipeline.P_SOURCE).push("0"), FirstSelection.class.getCanonicalName());
        parameters.set(BASE.push(IntegerVectorSpecies.P_FITNESS), SimpleFitness.class.getCanonicalName());
        
        return parameters;
    }

    @Test
    public void testGlobalParameters1() {
        System.out.println("setup (global parameters)");
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        instance.setup(state, BASE);
        
        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MUTATION_BOUNDED), "false");
        instance.setup(state, BASE);
        
        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MINGENE).push("1"), ""+min1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MAXGENE).push("1"), ""+max1);
        instance.setup(state, BASE);

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MINGENE).push("1"), ""+min1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MAXGENE).push("1"), ""+max1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_MUTATION_BOUNDED).push("1"), "false");
        instance.setup(state, BASE);

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT_TYPE), IntegerVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("0").push(IntegerVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MAXGENE), ""+max1);
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT_TYPE), IntegerVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("0").push(IntegerVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        // Set segment-specific initialization and mutation bounds
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MAXGENE), ""+max1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MUTATION_BOUNDED), "false");
        // Set a segment-specific mutation operator
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MUTATIONTYPE), IntegerVectorSpecies.V_RESET_MUTATION);
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RESET_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RESET_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(Double.NaN, instance.randomWalkProbability(1), 0.000001);
        assertEquals(Double.NaN, instance.randomWalkProbability(2), 0.000001);
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
    
    @Test
    public void testSegmentParameters3() {
        System.out.println("setup (segment parameters, mutation bounds, no mutation type)");
        final IntegerVectorSpecies instance = new IntegerVectorSpecies();
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_NUM_SEGMENTS), "2");
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT_TYPE), IntegerVectorSpecies.P_SEGMENT_START);
        
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("0").push(IntegerVectorSpecies.P_SEGMENT_START), ""+0);
        
        final double min1 = 1.0;
        final double max1 = 32.0;
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_SEGMENT_START), ""+1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MINGENE), ""+min1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MAXGENE), ""+max1);
        state.parameters.set(BASE.push(IntegerVectorSpecies.P_SEGMENT).push("1").push(IntegerVectorSpecies.P_MUTATION_BOUNDED), "false");
        instance.setup(state, BASE);
        

        assertEquals(NUM_GENES, instance.genomeSize);
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(0));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(1));
        assertEquals(IntegerVectorSpecies.C_RANDOM_WALK_MUTATION, instance.mutationType(2));
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(0), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(1), 0.000001);
        assertEquals(GLOBAL_WALK_PROB, instance.randomWalkProbability(2), 0.000001);
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
