#!/usr/bin/env python3
import sys
import random
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--N_founders",
        type=int,
        required=True,
        help="number of founder samples",
    )

    parser.add_argument(
        "--sample_map",
        type=str,
        required=True,
        help="map file to sample founders from"
    )

    parser.add_argument(
        "--props",
        required=True,
        type=str,
        nargs="+",
        help="proportion of each ancestry population (pop=fraction)"
    )

    parser.add_argument(
        "--N_gen",
        type=int,
        default=8,
        help="number of generations",
    )

    parser.add_argument(
        "--growth_rate",
        type=float,
        default=1.0,
        help="growth rate per generation",
    )

    parser.add_argument(
        "--max_samples",
        type=int,
        default=1000,
        help="max samples to simulate",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1234)

    parser.add_argument(
        "--outbase",
        default='mating',
        type=str,
    )

    return parser.parse_args(args=None if argv else ['--help'])


def parse_props(props):
    dic = {}
    for pop in props:
        (name, val) = pop.split("=")
        dic[name] = float(val)
    return(dic)

def sample_map(sample_map, props, n_founders):
    founders_per_pop = {pop: int(n_founders*prop) for pop, prop in props.items()}
    data = {}
    with open(sample_map, 'r') as fh:
        for line in fh:
            fields = line.rstrip().split()
            try:
                data[fields[1]].append(fields[0])
            except KeyError:
                data[fields[1]] = [fields[0]]

    samples_per_pop = {}
    for pop, sample_list in data.items():
        n_founders = founders_per_pop[pop]
        indices = random.sample(range(len(sample_list)), n_founders)
        samples_per_pop[pop] = [sample_list[i] for i in sorted(indices)]
    return(samples_per_pop)

def main(argv):
    args = parse_args(argv)
    props = parse_props(args.props)
    random.seed(args.seed)


    log = open(args.outbase+".log", 'w')
    print("##ARGUMENTS:", file=log)
    for k,v in vars(args).items():
        print(f"{k} {v}", file=log)
    print("##END", file=log)

    samples_per_pop = sample_map(args.sample_map, props, args.N_founders)
    with open(args.outbase+".founders", 'w') as fh:
        for pop, sample_list in samples_per_pop.items():
            for samplename in sample_list:
                print(samplename, pop, file=fh)


    parents_size = args.N_founders
    curr_gen = 0
    with open(args.outbase+".pedigree", 'w') as fh:
        while curr_gen < args.N_gen:
            parents = list(range(parents_size))
            children = []
            children_size = min(int(parents_size * args.growth_rate), args.max_samples)
            for i in range(children_size):
                if i % parents_size == 0:
                    ## shuffle in the beginning and when i >= parent size
                    print("shuffling", curr_gen, i, file=log)
                    random.shuffle(parents)
                p1 = parents[i % parents_size]
                p2 = parents[(i+1) % parents_size]
                if i == 0:
                    print(f"{p1} {p2}", end='', file=fh)
                else:
                    print(f" {p1} {p2}", end='', file=fh)
                children.append(i)
            print("", file=fh)

            parents_size = len(children)
            parents = children ## not needed
            curr_gen += 1

    log.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
