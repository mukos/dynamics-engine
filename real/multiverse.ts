import {listToMatrix} from "./helper";
import Universe from "./universe";

const SIM_OFFSET = 0;
const SIM_INCREMENT = 1;

export class Multiverse {
    size: number;
    resultsMap = new Map()

    constructor(size: number) {
        this.size = size
    }

    generation(){
        const numberOfElements = this.size*2;
        const matrixSize = numberOfElements * numberOfElements;
        const SIM_COMBINATIONS = Math.pow(2, matrixSize);
        for(let i = SIM_OFFSET; i < SIM_COMBINATIONS; i += SIM_INCREMENT){
            const numString = i.toString(2).padStart(matrixSize, '0').split('').map(sn => parseInt(sn));
            const activationMatrix = listToMatrix(numString, numberOfElements);

            // TODO: Row activation input variations results compare to total row activation results

            // Total activation vector
            let inputVec = []
            for(let k = 0; k < numberOfElements; k++){
                inputVec.push(k);
            }

            const universe = new Universe(this.size, activationMatrix, inputVec);
            universe.simulation();

            const uniState = universe.combinedCycleType();
            const resultValue = this.resultsMap.get(uniState) || 0;
            this.resultsMap.set(uniState, resultValue+1);
        }
    }

    analysis(){
        console.log('RESULTS MAP:', this.resultsMap);

        // TODO: Mathematical analysis of results to find patterns and formulas to predict different simulations
    }
}
