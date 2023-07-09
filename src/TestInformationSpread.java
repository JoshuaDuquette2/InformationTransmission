import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class TestInformationSpread {
    public IInformationSpread createInformationSpread(){
        return new InformationSpread();
    }

    @Test
    public void testloadGraphFromDataSet(){
        IInformationSpread infoSpread = createInformationSpread();
        int nodes = infoSpread.loadGraphFromDataSet("datasets/disconnected.mtx", 0.55);
        //assertEquals(8, nodes);
        System.out.println(infoSpread.avgDegree());
        System.out.println(infoSpread.rNumber());
        System.out.println(infoSpread.path(3,7));
        int[] init = infoSpread.getNeighbors(1);
        ArrayList<Integer> round2 = new ArrayList<Integer>();
        ArrayList<Integer> round3 = new ArrayList<Integer>();
        for(Integer node : init){
            round2.addAll(Arrays.stream(infoSpread.getNeighbors(node)).boxed().toList());
        }
        System.out.println(init.length + round2.size());
        for(Integer node : round2){
            round3.addAll(Arrays.stream(infoSpread.getNeighbors(node)).boxed().toList());
        }
        System.out.println(init.length + round2.size() + round3.size());
        System.out.println(round3.size());
        System.out.println("Level = " + infoSpread.generations(1, 1));
        System.out.println(infoSpread.clustCoeff(3));
        System.out.println(infoSpread.clustCoeffNodes(0.5, 1.0));
        System.out.println(infoSpread.generationsDegree(1, 0.9, 5));
        System.out.println(infoSpread.generationsCC(1, 0.9, 0.75, 1.0));
        System.out.println(infoSpread.degree(1));
        System.out.println(infoSpread.generationsHighDegLowCC(1, 0.9, 25, 0.75));
        System.out.println(infoSpread.rNumberCC(0.5, 1.0));
    }

}
