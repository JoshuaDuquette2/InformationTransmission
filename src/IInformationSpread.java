import java.util.Collection;

/**
 * @author ericfouh
 */
public interface IInformationSpread
{

    /**
     * 
     */
    public static final Object VISITED  = -5;
    /**
     * 
     */
    public static final int    INFINITY = Integer.MAX_VALUE;


    /**
     * Create a graph representation of the dataset. The first line of the file
     * contains the number of nodes add 1 to the number of nodes in the graph
     * since there is no node with id 0
     * 
     * @param filePath the path of the data
     * @param tau the transmissibility, probability of
     *            infection given contact between a susceptible and infected individual
     * @return the number of entries (nodes) in the dataset (graph)
     */
    int loadGraphFromDataSet(String filePath, double tau);

    /**
     * Return the neighbors ids of a specific node
     * 
     * @param id the id of the page
     * @return the array of neighbor(s)
     */
    int[] getNeighbors(int id);

    /**
     * return the shorthest path between two nodes
     * include the source and destination nodes in your collection
     * @param source      - the id of the origin node
     * @param destination - the id of the destination node
     * @return collection of nodes to follow to go from source to destination
     */
    Collection<Integer> path(int source, int destination);

    /**
     * Compute the average degree of the graph
     */
    double avgDegree();


    /**
     * Compute the basic reproduction number R0
     * R0 = TRANSMISSIBILITY (tau) * average_degree
     * @return the basic reproduction number R0
     */
     double rNumber();


    /**
     * Given a specific node id (seed) this method will return the number of
     * "spreadLevels" necessary to reach a percentage (threshold) of the nodes
     * in the graph
     * 
     * @param seed      - the id of the seed page
     * @param threshold - the percentage of nodes to reach
     * @return the number of spread Levels necessary to reach threshold percent
     *         nodes in the graph or -1 if the seed is not in the graph
     */
    int generations(int seed, double threshold);

    // -- Degree
    /**
     * @param n the node
     * @return the degree of the node
     */
    int degree(int n);
    
    /**
     * @param d the degree
     * @return all the node with degree d
     */
    Collection<Integer> degreeNodes(int d);
    
    /**
     * Given a specific node id (seed) this method will return the number of
     * "generations" necessary to reach a percentage (threshold) of the nodes
     * in the graph when all the nodes with a given degree d are removed
     * 
     * @param seed      - the id of the seed page
     * @param threshold - the percentage of nodes to reach
     * @param d        - the degree of the nodes to be removed
     * @return the number of spread Levels necessary to reach threshold percent
     *         nodes in the graph
     */
    int generationsDegree(int seed, double threshold, int d);

    /**
     * Compute the basic reproduction number R0 when
     * all the nodes with a given degree d are removed
     * R0 = tau * average_degree
     *
     * @param d        - the degree of the nodes to be removed
     * @return the basic reproduction number
     */
    double rNumberDegree(int d);
    
    // -- CLustering Coefficient

    /**
     * nodes with degree 0 or 1 have a cc of 0
     * @param n the node
     * @return the  clustering coefficient of n
     */
    double clustCoeff(int n);

    /**
     * precision: 0.01 (use when comparing CC values)
     * @param low - the lower bound (inclusive) of the cc range
     * @param high - the upper bound (inclusive) of the cc range
     * @return a collection of nodes with a clustering coefficient 
     * within [low, high]
     */
    Collection<Integer> clustCoeffNodes(double low, double high);



    /**
     * precision: 0.01
     * Given a specific node id (seed) this method will return the number of
     * "generations" necessary to reach a percentage (threshold) of the nodes
     * in the graph when all the nodes with a clustering coefficient within the 
     * range [low, high] are removed
     * 
     * @param seed      - the id of the seed page
     * @param threshold - the percentage of nodes to reach
     * @param low - the lower bound (inclusive) of the cc range
     * @param high - the upper bound (inclusive) of the cc range
     * @return the number of spread Levels necessary to reach threshold percent
     *         nodes in the graph
     */
    int generationsCC(int seed, double threshold, double low, double high);


    /**
     * Compute the basic reproduction number R0 when
     * all the nodes with a clustering coefficient within the
     *  range [low, high] are removed
     * R0 = tau * average_degree
     *
     * @param low - the lower bound (inclusive) of the cc range
     * @param high - the upper bound (inclusive) of the cc range
     * @return the basic reproduction number
     */
    double rNumberCC(double low, double high);
    
    // high degree low cc
    /**
     * precision: 0.01
     * @param lowBoundDegree - the lower bound (inclusive) of the degree
     * @param upBoundCC - the upper bound (inclusive) of the cc 
     * @return a collection of nodes with degree >= lowBoundDegree and
     *  clustering coefficient <= upBoundCC
     */
    Collection<Integer> highDegLowCCNodes(int lowBoundDegree, double upBoundCC);
    
    /**
     * Given a specific node id (seed) this method will return the number of
     * "generations" necessary to reach a percentage (threshold) of the nodes
     * in the graph when all the nodes with a clustering coefficient below a 
     * given value and a degree above a given value are removed.
     * 
     * @param seed      - the id of the seed page
     * @param threshold - the percentage of nodes to reach
     * @param lowBoundDegree - the lower bound (inclusive) of the degree
     * @param upBoundCC - the upper bound (inclusive) of the cc
     * @return the number of spread Levels necessary to reach threshold percent
     *         nodes in the graph
     */
    int generationsHighDegLowCC(int seed, double threshold, int lowBoundDegree, double upBoundCC);

    /**
     * Compute the basic reproduction number R0 when
     * all the nodes with a clustering coefficient below a
     * given value and a degree above a given value are removed.
     * R0 = TRANSMISSIBILITY * average_degree
     *
     * @param lowBoundDegree - the lower bound (inclusive) of the degree
     * @param upBoundCC - the upper bound (inclusive) of the cc
     * @return the basic reproduction number
     */
    double rNumberDegCC(int lowBoundDegree, double upBoundCC);
}
