import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.sql.Array;
import java.util.*;

public class InformationSpread implements IInformationSpread {


    private double tau;
    Set<Integer> indices = new HashSet<Integer>();
    Graph graph;

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
    @Override
    public int loadGraphFromDataSet(String filePath, double tau) {
        graph = new GraphL();
        this.tau = tau;
        int nodeCount = 1; //Graph includes one node implicitly
        int nodes = 0;
        int edges = 0;

        try {
            BufferedReader reader = new BufferedReader(new FileReader(filePath));
            String   fileData = reader.readLine();
            String[] fileDataArr = fileData.split(" ");

            nodes = Integer.parseInt(fileDataArr[0]);
            edges = Integer.parseInt(fileDataArr[1]);
            graph.init(nodes + 1); //may need to be checked later down the line


            while(true){
                String data = reader.readLine();

                if(data == null){
                    break;
                }
                String[] dataArr = data.split(" ");
                int nodeStart = Integer.parseInt(dataArr[0]);
                int nodeEnd   = Integer.parseInt(dataArr[1]);
                double weight    = Double.parseDouble(dataArr[2]);
                if(weight >= tau) {
                    int multWeight = (int)(weight * 100);
                    graph.addEdge(nodeStart, nodeEnd, multWeight);
                    graph.addEdge(nodeEnd, nodeStart, multWeight);
                    indices.add(nodeStart);
                    indices.add(nodeEnd);
                    nodeCount++;
                }
            }
        } catch (IOException e){
            e.printStackTrace();
        }

        System.out.println("Created graph with nodes = " + graph.nodeCount() + " and edges = " + graph.edgeCount());

        return nodeCount;
    }
    
    
    /**
     * Return the neighbors ids of a specific node
     * 
     * @param id the id of the page
     * @return the array of neighbor(s)
     */
    @Override
    public int[] getNeighbors(int id) {
        List<Integer> neighbors = new ArrayList<Integer>();
        for(int i = 1; i < graph.nodeCount(); i++){
            if(graph.hasEdge(id, i) || graph.hasEdge(i, id)){
                neighbors.add(i);
            }
        }

        if(neighbors.contains(0)){
            neighbors.clear();
        }

        int[] output = new int[neighbors.size()];
        for(int i = 0; i < neighbors.size(); i++){
            output[i] = neighbors.get(i);
        }

        return output;
    }

    
    /**
     * return the shorthest path between two nodes
     * include the source and destination nodes in your collection
     * @param source      - the id of the origin node
     * @param destination - the id of the destination node
     * @return collection of nodes to follow to go from source to destination
     */
    @Override
    public Collection<Integer> path(int source, int destination) {

        int n = graph.nodeCount();
        for (int i = 1; i < n; i++) {
            graph.setValue(i, false);
        }

        List<Integer> path = new ArrayList<Integer>();
        Integer[] distances = new Integer[n];
        Integer[] predictions = new Integer[n];

        for(int i = 0; i < n; i++){
            predictions[i] = 0;
        }

        for (int i = 1; i < n; i++) {
            distances[i] = Integer.MAX_VALUE;
        }
        distances[source] = 0;

        Queue<Integer> queue = new LinkedList<Integer>();
        graph.setValue(source, true);
        queue.add(source);

        while (queue.size() > 0) {
            Integer currentIndex = queue.poll();

            for (Integer neighborIndex : graph.neighbors(currentIndex)) {
                if (!((boolean)graph.getValue(neighborIndex))) {
                    queue.add(neighborIndex);
                    graph.setValue(neighborIndex, true);
                }
                if (distances[currentIndex] + graph.weight(currentIndex, neighborIndex) <= distances[neighborIndex]) {
                    predictions[neighborIndex] = currentIndex;
                    distances[neighborIndex] = distances[currentIndex] + graph.weight(currentIndex, neighborIndex);
                }
            }
        }

        for (int j = destination; predictions[j] != 0; j = predictions[j]) {
            path.add(j); //Implicit addition of destination
        }

        path.add(source); //Addition of source
        Collections.reverse(path);

        return path;
    }
    
    
    
    /**
     * Compute the average degree of the graph
     */
    @Override
    public double avgDegree() {
        int totalDegree = 0;
        for(Integer id : indices){
            totalDegree += degree(id);
        }
        return (double)totalDegree / (graph.nodeCount() - 1);
    }
    
    
    
    /**
     * Compute the basic reproduction number R0
     * R0 = TRANSMISSIBILITY (tau) * average_degree
     * @return the basic reproduction number R0
     */
    @Override
    public double rNumber() {
        return 1 * tau * avgDegree();
    }
    
    
    
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

    @Override
    public int generations(int seed, double threshold) {

        int n = graph.nodeCount();
        if (seed <= 0 || seed >= n || threshold < 0 || threshold > 1) {
            return -1;
        } else if (threshold == 0) {
            return 0;
        }

        for (int i = 1; i < n; i++) {
            graph.setValue(i, false);
        }

        int level = 0;
        int currentNeighborsCount = 0;
        int currentLevel = 1;
        int nodeCount = 1;


        if (((double)nodeCount / (n - 1)) >= threshold) {
            return level;
        }

        Queue<Integer> queue = new LinkedList<Integer>();
        graph.setValue(seed, true);
        queue.add(seed);

        while (!queue.isEmpty()) {
            int current = queue.poll();
            currentNeighborsCount += getNeighbors(current).length;
            int[] neighbors = getNeighbors(current);

            for (Integer v : neighbors) {
                if ((boolean) graph.getValue(v)) {
                    currentNeighborsCount--;
                } else {
                    queue.add(v);
                    nodeCount++;
                    graph.setValue(v, true);
                }
            }

            currentLevel -= 1;
            if (currentLevel == 0) {
                level++;
                if (((double) nodeCount / (graph.nodeCount() - 1)) >= threshold) {
                    return level;
                }
                currentLevel += currentNeighborsCount;
                currentNeighborsCount = 0;
            }
        }
        return -1;
    }

    /**
     * @param n the node
     * @return the degree of the node
     */
    @Override
    public int degree(int n) {
        return graph.neighbors(n).length;
    }

    /**
     * @param d the degree
     * @return all the node with degree d
     */
    @Override
    public Collection<Integer> degreeNodes(int d) {
        List<Integer> degreeMatch = new ArrayList<Integer>();
        for(Integer id : indices){
            if(degree(id) == d) {
                degreeMatch.add(id);
            }
        }
        return degreeMatch;
    }
    
    
    
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
    @Override
    public int generationsDegree(int seed, double threshold, int d) {
        Collection<Integer> nodesToRemove = degreeNodes(d);

        if(seed <= 0 || seed >= graph.nodeCount() || threshold < 0 || threshold > 1){
            return -1;
        }

        if(nodesToRemove.isEmpty()){
            return -1;
        } else if(nodesToRemove.contains(seed)){
            return 0;
        }
        for(Integer currentIndexToRemove : nodesToRemove) {
            int[] neighbors = graph.neighbors(currentIndexToRemove);
            for(Integer edge : neighbors){
                if(nodesToRemove.contains(edge)) {
                    graph.removeEdge(currentIndexToRemove, edge);
                    graph.removeEdge(edge, currentIndexToRemove);
                }
            }
        }
        return generations(seed, threshold);
    }
    
    /**
     * Compute the basic reproduction number R0 when
     * all the nodes with a given degree d are removed
     * R0 = tau * average_degree
     *
     * @param d        - the degree of the nodes to be removed
     * @return the basic reproduction number
     */
    @Override
    public double rNumberDegree(int d) {
        Collection<Integer> nodesToRemove = degreeNodes(d);
        for(Integer currentIndexToRemove : nodesToRemove) {
            int[] neighbors = graph.neighbors(currentIndexToRemove);
            for(Integer edge : neighbors){
                graph.removeEdge(currentIndexToRemove, edge);
                graph.removeEdge(edge, currentIndexToRemove);
            }
        }
        return rNumber();
    }
    
    
    
    // -- CLustering Coefficient
    /**
     * nodes with degree 0 or 1 have a cc of 0
     * @param n the node
     * @return the  clustering coefficient of n
     */
    @Override
    public double clustCoeff(int n) {
        int degree = degree(n);
        if(n <= 0 || n >= graph.nodeCount()){
            return -1;
        }
        if(degree < 2) {
            return 0;
        }
        int[] neighbors = graph.neighbors(n);
        int connectedness = 0;
        for(int i = 0; i < neighbors.length; i++){
            int connection1 = neighbors[i];
            for(int j = 0; j < neighbors.length; j++) {
                int connection2 = neighbors[j];
                //System.out.println("Checking edge between " + connection1 + " and " + connection2 + " = " + graph.hasEdge(connection1, connection2));
                if (graph.hasEdge(connection1, connection2)) {
                    connectedness++;
                }
            }

        }
        return ((double)connectedness) / (degree * (degree - 1));
    }
    
    
    
    /**
     * precision: 0.01 (use when comparing CC values)
     * @param low - the lower bound (inclusive) of the cc range
     * @param high - the upper bound (inclusive) of the cc range
     * @return a collection of nodes with a clustering coefficient 
     * within [low, high]
     */
    @Override
    public Collection<Integer> clustCoeffNodes(double low, double high) {
        List<Integer> coeffs = new ArrayList<Integer>();
        for(Integer node : indices){
            double clustCoeff = clustCoeff(node);
            if(((int)(clustCoeff * 100)) >= ((int)(low*100)) && ((int)(clustCoeff * 100)) <= ((int)(high * 100))){
                coeffs.add(node);
            }
        }
        return coeffs;
    }
    
    
    
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
    @Override
    public int generationsCC(int seed, double threshold, double low, double high) {
        if(seed <= 0 || seed >= graph.nodeCount() || threshold < 0 || threshold > 1){
            return -1;
        }
        Collection<Integer> nodesToRemove = clustCoeffNodes(low, high);

        if(nodesToRemove.isEmpty()){
            return -1;
        } else if(nodesToRemove.contains(seed)){
            return 0;
        }

        for(Integer current : indices){
            if(nodesToRemove.contains(current)){
                for(Integer neighbor : getNeighbors(current)){
                    graph.removeEdge(current, neighbor);
                    graph.removeEdge(neighbor, current);
                }
            }
        }
        return generations(seed, threshold);
    }
    
    
    
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
    @Override
    public double rNumberCC(double low, double high) {
        Collection<Integer> nodesToRemove = clustCoeffNodes(low, high);
        for(Integer currentIndexToRemove : nodesToRemove) {
            int[] neighbors = graph.neighbors(currentIndexToRemove);
            for(Integer edge : neighbors){
                if(nodesToRemove.contains(edge)) {
                    graph.removeEdge(currentIndexToRemove, edge);
                    graph.removeEdge(edge, currentIndexToRemove);
                }
            }
        }
        return rNumber();
    }
    
    
    
    // high degree low cc
    /**
     * precision: 0.01
     * @param lowBoundDegree - the lower bound (inclusive) of the degree
     * @param upBoundCC - the upper bound (inclusive) of the cc 
     * @return a collection of nodes with degree >= lowBoundDegree and
     *  clustering coefficient <= upBoundCC
     */
    @Override
    public Collection<Integer> highDegLowCCNodes(int lowBoundDegree, double upBoundCC) {
        List<Integer> nodesToRemove = new ArrayList<Integer>();


        for(Integer node : indices){
            double degree = degree(node);
            BigDecimal clusterCoeff = new BigDecimal(clustCoeff(node)).setScale(2, RoundingMode.DOWN);
            if(degree >= lowBoundDegree && clusterCoeff.doubleValue() <= upBoundCC){
                nodesToRemove.add(node);
            }
        }
        return nodesToRemove;
    }
    
    
    
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
    @Override
    public int generationsHighDegLowCC(int seed, double threshold, int lowBoundDegree, double upBoundCC) {
        Collection<Integer> nodesToRemove = highDegLowCCNodes(lowBoundDegree, upBoundCC);

        if(seed <= 0 || seed >= graph.nodeCount() || threshold < 0 || threshold > 1){
            return -1;
        }

        if(nodesToRemove.isEmpty()){
            return -1;
        } else if(nodesToRemove.contains(seed)){
            return 0;
        }
        for(Integer currentIndexToRemove : nodesToRemove) {
            int[] neighbors = graph.neighbors(currentIndexToRemove);
            for(Integer edge : neighbors){
                if(nodesToRemove.contains(edge)) {
                    graph.removeEdge(currentIndexToRemove, edge);
                    graph.removeEdge(edge, currentIndexToRemove);
                }
            }
        }
        return generations(seed, threshold);
    }
    
    
    
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
    @Override
    public double rNumberDegCC(int lowBoundDegree, double upBoundCC) {
        Collection<Integer> nodesToRemove = highDegLowCCNodes(lowBoundDegree, upBoundCC);
        for(Integer currentIndexToRemove : nodesToRemove) {
            int[] neighbors = graph.neighbors(currentIndexToRemove);
            for(Integer edge : neighbors){
                graph.removeEdge(currentIndexToRemove, edge);
                graph.removeEdge(edge, currentIndexToRemove);
            }
        }
        return rNumber();
    }

}