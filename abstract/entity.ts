// For real world simulation and artificial intelligence

abstract class Entity {}

export class Action extends Entity {
    name: string = ''

    constructor(name: string) {
        super();
        this.name = name;
    }
}

export class Element extends Entity {
    name: string = ''
    count: number = 0

    constructor(name?: string, initialCount?: number) {
        super();
        this.name = name ?? '';
        this.count = (initialCount && initialCount > 0) ? initialCount : 0;
    }

    add(){
        this.count++;
    }

    remove(){
        if(this.count > 0){
            this.count--;
        }
    }
}

abstract class Collection<T extends Entity> extends Element {
    // [Symbol.iterator](){
    //     // todo
    // }
}

export class Set<T extends Entity> extends Collection<T> {
    elements = new global.Set<T>();

    set(element: T){
        this.elements.add(element)
    }

    empty(){
        this.elements.clear()
    }
}

export class Vector<T extends Entity> extends Collection<T> {
    elements: Array<T> = []

    push(element: T){
        this.elements.push(element)
    }

    get(index: number){
        if (index < this.elements.length) {
            return this.elements[index];
        }
        return undefined;
    }
}

export class Relation<T extends Entity, K extends Entity> {}

export class Event extends Relation<Action, Element> {
    element: Element
    action: Action
    count: number = 0

    constructor(action: Action, element: Element) {
        super();
        this.action = action;
        this.element = element
    }
}

export class Cause extends Relation<Event, Event> {
    eventInit: Event
    eventResult: Event

    constructor(init: Event, result: Event) {
        super();
        this.eventInit = init;
        this.eventResult = result;
    }
}

export class System extends Element {
    // NOTE - clock: system with two states, inclusive, define through initial state

    states = new Vector<Set<Element>>()
    events: Set<Event> // current events, all listed with numeric differences

    // WAY 1 (OLD)
    causes: Set<Cause>; // all system event causalities

    // WAY 2 (NEW) - TODO
    // CC, CD, DC, DD - CAUSE matrices

    constructor(causes?: Set<Cause>, events?: Set<Event>) {
        super();
        this.causes = causes || new Set<Cause>()
        this.events = events || new Set<Event>()
    }

    addCause(event1: Event, event2: Event){
        this.causes.set(new Cause(event1, event2));
    }

    run(){
        // this.states.push()
        // for(const cause in this.causes){
        //
        // }
    }
}

// abstract class Start extends Action {}
// abstract class Stop extends Action {}


