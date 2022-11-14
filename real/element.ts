export default class Element {
    index: number | undefined;
    count = 0
    antiCount = 0
    sum = 0

    constructor(index: number) {
        this.index = index;
    }

    evaluate() {
        this.sum = this.count + this.antiCount;
        switch(true) {
            case this.sum > 0: {
                this.count = this.sum;
                this.antiCount = 0;
                break;
            }
            case this.sum < 0: {
                this.count = 0;
                this.antiCount = this.sum;
                break;
            }
            default: {
                this.reset()
            }
        }
    }

    reset(){
        this.count = 0;
        this.antiCount = 0;
        this.sum = 0;
    }

    log(){
        // console.log('Count::Anticount', this.count, this.antiCount);
    }
}
