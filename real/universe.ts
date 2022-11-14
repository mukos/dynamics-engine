import Element from "./element";

const SANITY_COUNT = 10;

export default class Universe {
    activeElements: Element[] // COMPUTER
    activationMatrix: number[][] // PROGRAM
    activationSet: Set<number> // INPUT

    history: any[][] = []

    // TODO: Refactor data structures above and make code more methodic and procedural

    constructor(numOfElements: number, activationMatrix: number[][], activationVector: number[]) {
        this.activeElements = []
        for(let i = 1; i <= numOfElements; i++){
            this.activeElements.push(new Element(i))
        }
        if(activationVector.length > this.activeElements.length * 2){
            throw new Error('invalid activation vector length');
        }
        if(activationMatrix.length !== this.activeElements.length * 2){
            throw new Error('invalid activation matrix length');
        }
        this.activationMatrix = activationMatrix;
        // console.log('Initial activation matrix: ', this.activationMatrix);
        this.activationSet = new Set(activationVector.sort());
        // console.log('Initial activation vector: ', this.activationVector);
    }

    getCycleType(elementIndex: number, countType: number){
        if(this.history.length < SANITY_COUNT){
            return '|'
        }
        for (let i = 0; i < this.history.length - 2; i++) {
            if(this.history[i][elementIndex][countType] === this.history[this.history.length-1][elementIndex][countType]){
                return '0'
            }
        }
        return countType === 0 ? '1' : '-1'
    }

    combinedCycleType(){
        const uniState = [];
        for(let i = 0; i < this.activeElements.length; i++){
            const countType = this.getCycleType(i, 0);
            const anticountType = this.getCycleType(i, 1);
            uniState.push([ countType, anticountType ]);
        }
        return uniState.toString();
    }

    simulation(){
        let sanity = SANITY_COUNT;
        while(this.activationSet.size > 0 && sanity-- > 0) {
            const change = this.activation();
            if(!change){
                return;
            }
            this.activationSet = this.evaluation();
        }
    }

    activation(){
        let valueChange = false;
        for (let rowIndex = 0; rowIndex < this.activationMatrix.length; rowIndex++) {
            if (this.activationSet.has(rowIndex)) {
                const row = this.activationMatrix[rowIndex];
                for (let columnIndex = 0; columnIndex < row.length; columnIndex++) {
                    const value = row[columnIndex];
                    if(value){
                        valueChange = true;
                        const activeElementIndex = Math.floor(columnIndex / 2);
                        const element = this.activeElements[activeElementIndex];
                        if (columnIndex % 2 == 0) {
                            element.count += value;
                        } else {
                            element.antiCount -= value;
                        }
                    }
                }
            }
        }
        return valueChange;
    }

    evaluation(){
        const stateVector = [];
        let newActivationSet = new Set<number>();
        for(let i = 0; i < this.activeElements.length; i++){
            const element = this.activeElements[i];
            element.evaluate();
            stateVector.push([element.count, element.antiCount])
            element.log();
            // TODO: Multi activation of matrix from evaluated elements (Vector sum instead of Set)
            if(element.sum > 0){
                newActivationSet.add(i*2);
            } else if(element.sum < 0){
                newActivationSet.add(i*2+1);
            }
        }
        this.history.push(stateVector);
        return newActivationSet;
    }
}
