import {Multiverse} from "./multiverse";

console.time('Session');

const multiverse = new Multiverse(2);

multiverse.generation();
multiverse.analysis();

console.timeEnd('Session');

