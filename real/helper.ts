export function listToMatrix(list: any[], elementsPerSubArray: number) {
    var matrix = [], i, k;

    for (i = 0, k = -1; i < list.length; i++) {
        if (i % elementsPerSubArray === 0) {
            k++;
            matrix[k] = [];
        }

        // @ts-ignore
        matrix[k].push(list[i]);
    }

    return matrix;
}
